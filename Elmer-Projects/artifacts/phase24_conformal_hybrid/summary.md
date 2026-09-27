# Phase24 完全共形プリズムメッシュ + HYPRE PCG：COMSOL/MUMPS parity 検証

実施日: 2026-09-25。物性・TES法則・回路定数・熱浴条件は変更していない。変えたのはメッシュ生成と線形ソルバー設定だけ。

## 結論

1. **HYPRE と MUMPS は完全に一致した。** 同じ完全共形メッシュで、定常電流の差は 0.0001%。100 µs の過渡波形の差は最大 1.9e-5 µA だった。
2. **以前の「HYPRE の過渡応答が平坦になる」現象は、HYPRE の偽収束が原因だった。** 相対許容値 2e-8 では、回路の電力更新やパルスによる右辺の変化が熱浴 Dirichlet 行の支配する ‖b‖ に埋もれ、`Required iterations 0` のまま解が更新されない。mortar のない対称正定値の系に PCG+BoomerAMG を使い、許容値を 1e-12 にすると解消した。
3. **Phase24 の四面体メッシュは熱コンダクタンス G を数 % 過大評価する。** 厚さ 1 µm の薄膜を四面体で分割しているためで、h を変えると非単調に ±3% ばらつく。z=192 µm 界面だけを細分化する案（next_controlled_fix.md）は 1.0433 → 1.0412 にしかならず、効果はなかった。
4. **完全共形のプリズムメッシュでは、G は単調に収束する。** 連続体の G は旧基準（prod_v2）の約 0.955〜0.96 倍と外挿される。旧基準は TES/膜の mortar によって G が約 4% 上乗せされていた。
5. COMSOL との比較結果は下表のとおり。RMSE、t10、t50 は Gate4 の基準を満たした。最大差だけが 0.080 µA で、基準の 0.05 µA を超えた。原因は基準電流の −2.36%（G が低いこと）と考えられる。

## メッシュ

`generate_hybrid_prism_geometry.py --conformal-tes-stack --conformal-abs` で生成した（レシピは `project.json` の `mesh_hybrid_fullconf_h12` を参照）。

- 積層部の底面に TES の長方形と Stycast の円板を刻み込み、全層をプリズム要素で押し出した。TES/膜と TES/Stycast は同じ OCC 面を共有する。
- 吸収体の底面に円板を刻み込み、`setPeriodic`（恒等変換）で Stycast 上面と同じ表面メッシュにした。そのうえで ElmerGrid `-merge 1e-10` を実行した。
- 共有節点数は TES/膜 2252、TES/Stycast 1716、Stycast/吸収体 1716。**mortar は 0 組**。

## 電力固定での G_eff（MUMPS、旧基準=1）

| メッシュ | G 比 |
|---|---:|
| Phase24 四面体 対照 h=10 | 1.0433 |
| 上記 + z=192 界面 10.5 µm | 1.0412 |
| 四面体 h=7 / h=7（Netgen 最適化なし）/ h=6 | 1.0178 / 1.0161 / 1.0486 |
| 共形プリズム h=17（層 2/2/1） | 0.9785 |
| 共形プリズム h=12（層 2/2/1） | 0.9716 |
| 共形プリズム h=12（層 6/3/3） | 0.9694 |
| 共形プリズム h=8（層 4/2/2） | 0.9645 |

## 定常（初期温度 T0 から独立に計算、stage11 バイナリ）

| バックエンド | 電流 [µA] | MUMPS との差 |
|---|---:|---:|
| MUMPS | 139.676348 | — |
| HYPRE PCG+BoomerAMG 1e-12 | 139.676200 | 0.0001% |
| HYPRE FlexGMRES 2e-8（旧運用値） | 142.961209 | +2.35%（偽収束） |

### 過去の「HYPRE が COMSOL に近かった」定常値の正体

過去の Gate3 の HYPRE 結果を見直すと、COMSOL（143.055 µA）に近いものはすべて同じ経過をたどっていた。HYPRE の反復数が最後に 0 になり、非線形反復 3〜4 回で「非線形残差 0.0」として終了している。

| ケース | HYPRE [µA] | 同じメッシュの MUMPS | HYPRE の反復数の推移 |
|---|---:|---:|---|
| conformal_refine20（2e-8） | 143.522 | 144.711 | 71 → 9 → 2 → 0 |
| s32m2（2e-6） | 143.568 | 218.6（steady_branch_diagnosis） | 3991 → 15 → 1 → 0 |
| s32 / strip（2e-8） | 143.534 | — | 281 → 19 → 3 → 0 |
| hybfull_h12（2e-8） | 142.961 | 139.676 | → 0 |

どのケースも、初期条件の TES 温度は `T_0 = 0.16857 K` で、これは旧 MUMPS 解の値である。完全共形メッシュで初期温度だけを変えて、同じ HYPRE 2e-8 を実行した（`initial_guess_probe/`）。

| 初期温度 | HYPRE の最終 TES 温度 | HYPRE の電流 |
|---:|---:|---:|
| 0.1680 K | 0.168021 K | 163.08 µA |
| 0.16857 K | — | 142.96 µA |
| 0.1690 K | 0.169013 K | 128.49 µA |

許容値 2e-8 の HYPRE では、TES 温度は初期値から約 20 µK しか動かない。電流は初期温度しだいで決まり、MUMPS の真の解（139.68 µA）とは関係がない。過去に HYPRE が COMSOL に近く見えたのは、旧基準（≒ COMSOL）の TES 温度を初期値に使っていたためである。

## 100 µs 過渡（同じ MUMPS 定常解から再始動、0–99.376 µs）

| | 最大差 [µA] | RMSE [µA] | t10 差 [µs] | t50 差 [µs] | 基準電流 |
|---|---:|---:|---:|---:|---:|
| HYPRE（完全共形） | 0.0803 | 0.0258 | −0.23 | +0.98 | −2.36% |
| MUMPS（完全共形） | 0.0803 | 0.0258 | −0.23 | +0.98 | −2.36% |
| 旧基準 MUMPS（prod_v2） | 0.0287 | 0.0159 | −0.46 | +0.01 | +0.51% |

Gate4 の基準は最大差 0.05 µA、RMSE 0.03 µA、時定数の差 1 µs。

## COMSOL 相当の解像度での再現（2026-09-25 追加）

COMSOL の設定を確認した結果は次のとおり。

- 伝熱の離散化は「2 次ラグランジュ」
- 全体の要素サイズは「普通」（最大 600 µm、最小 108 µm）
- メンブレンのスイープは厚さ方向に 7 分割
- ドメイン構成は Elmer と同じ。SC は Stycast、Al 配線は無効、熱源は TES 全体

Elmer でも同じ条件（完全共形、面内サイズ一様、Si1 を 7 層、ElmerGrid `-increase` で 2 次要素化）を作り、G を比べた。

| 面内サイズ | 1 次要素の G 比 | 2 次要素の G 比 |
|---|---:|---:|
| 100 µm | 1.0596 | **0.9900** |
| 50 µm | 1.0327 | 0.9816 |
| 収束値（外挿） | ≈0.955 | ≈0.955 |

COMSOL の定常電流に対応する G 比は約 0.994。2 次要素・100 µm での非線形定常計算（MUMPS）では **142.098 µA（COMSOL 比 −0.67%）** となり、完全共形 h=12・1 次要素の −2.36% から大きく COMSOL に近づいた。

COMSOL の値そのものも、TES の縁の面内解像度が粗いため、G が約 4% 高めに出ていると判断する。旧 MUMPS は TES/膜の mortar によって G が同じ程度上乗せされ、結果として COMSOL と一致していた。

## 75 ms 過渡の高速化（h=8、1次要素、2026-09-26〜27）

実行環境：WSL ディストリ ElmerGPU（G:、Ubuntu 24.04）、CPU 版 HYPRE+Elmer（/opt/elmer-gpu/cpu）、Ryzen 5 5600X。

| 実行 | 時間積分 | 実時間 | 線形解き | ピーク | COMSOL との差（最大 / RMSE） |
|---|---|---|---|---|---|
| シリアル | BDF1 645 ステップ | 8548 s | − | − | − |
| シリアル (b2y) | BDF2 282 ステップ + パルスで BDF1 リセット | 5923 s | 2624 s（907 回） | 7.6855 µA @429.1 µs | 0.090 / 0.030 µA |
| **MPI 6 ランク (m6a)** | 同上 | **1596 s** | 1200 s（898 回） | 7.6854 µA @429.1 µs | 0.090 / 0.030 µA |
| COMSOL | − | − | − | 7.775 µA @428.0 µs | − |

- m6a と b2y の差は最大 0.0018 µA。系列は `transient_75ms/`。
- MPI+HYPRE が誤った定常（154.67 µA、シリアル 138.676 µA）を出していた原因：Phase24 版 `SolveHypre.c` が、界面行が未合算の局所 RHS／行列再設定を `SetValues`（上書き）で渡していた。上流と同じ「ゼロクリア → `AddToValues`」に戻して解消（6 ランク定常 138.676 µA）。
- WSL の Open MPI は hwloc の OpenCL プラグインが AMD の OpenCL/HSA ランタイムを読み込み CPU を奪う。`HWLOC_COMPONENTS=-opencl,-rsmi,-levelzero`、`OCL_ICD_VENDORS=<空ディレクトリ>` で回避（runner が自動設定）。
- BoomerAMG は 1 反復の残差縮小率が約 0.95 と効きが悪い。緩和法・粗化法・強結合閾値の 5 通りと Residual Mode を試したが、1 回の解き時間は約 3 s（定常）で改善しなかった。
- Elmer 本体の差分は上流 9916c3f05 からの累積パッチ `elmer_native_source.patch`。

## 残課題

- COMSOL の G はプリズム収束値より約 4% 高い。COMSOL 側の離散化誤差（共形 FEM の上界）なのか、物理設定の差なのかは、COMSOL のメッシュ・モデル情報がないと判別できない。
- GPU：HIP 版 HYPRE はホストメモリのままでは実質 CPU で計算しており、デバイスへの移行（PHASE24_HYPRE_MIGRATE）はハングして AMD ドライバの TDR を起こした。
- 行列組み立ての高速化（k(T) 動的プリズム経路・ソース項のみの経路）は反映済み。

## 再現

- G_eff：`scripts/support/run_phase24_conformal_hybrid_fixed_power.py <mesh> --key <k> [--full-conformal]`
- HYPRE の精度：`scripts/support/probe_phase24_hypre_accuracy.py`
- 定常：`scripts/support/run_phase24_gate3_hypre_steady.py --source-project artifacts/phase24_conformal_hybrid/project.json --mesh mesh_hybrid_fullconf_h12 --base-case case_conformal_hybrid_steady_base --linear-system iterative_hypre_pcg_boomeramg --linear-tolerance 1e-12`
- 過渡：`scripts/support/run_phase24_gate4_5_nomortar.py --window 100us --backend hypre --hypre-system iterative_hypre_pcg_boomeramg --linear-tolerance 1e-12 --conformal-tes-stack --base-project artifacts/phase24_conformal_hybrid/project.json --mesh mesh_hybrid_fullconf_h12 ...`
- 比較の数値：`transient_100us/summary.json`
- 75 ms MPI：`scripts/support/run_phase24_gpu_elmergpu.py transient --mesh mesh_hybrid_fullconf_h8 --window 75ms_bdf2 --mpi 6 --elmer cpu/elmer --reference-case <MPI 定常ケース> --linear-system iterative_hypre_pcg_boomeramg --case-option bdf_order=2 ...`（完全な SIF は `transient_75ms/mpi6_bdf2_m6a.sif`、メッシュに `ElmerGrid 2 2 <mesh> -metiskway 6` で partitioning.6 が必要）
- ランタイム：`TES_TOOLCHAIN_BIN=G:\Elmer-Migration\msys64\ucrt64\bin`、runner には `--toolchain-bin` で同じパスを渡す

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
| MPI 6 ランク (m6a) | 同上 | 1596 s | 1200 s（898 回） | 7.6854 µA @429.1 µs | 0.090 / 0.030 µA |
| **MPI 6 ランク + 集中化前処理 (m6b)** | 同上 | **1383 s** | 877 s（883 回、平均 33 反復） | 7.6854 µA @429.1 µs | 0.090 / 0.030 µA |
| COMSOL | − | − | − | 7.775 µA @428.0 µs | − |

- m6a と b2y の差は最大 0.0018 µA。系列は `transient_75ms/`。
- MPI+HYPRE が誤った定常（154.67 µA、シリアル 138.676 µA）を出していた原因：Phase24 版 `SolveHypre.c` が、界面行が未合算の局所 RHS／行列再設定を `SetValues`（上書き）で渡していた。上流と同じ「ゼロクリア → `AddToValues`」に戻して解消（6 ランク定常 138.676 µA）。
- WSL の Open MPI は hwloc の OpenCL プラグインが AMD の OpenCL/HSA ランタイムを読み込み CPU を奪う。`HWLOC_COMPONENTS=-opencl,-rsmi,-levelzero`、`OCL_ICD_VENDORS=<空ディレクトリ>` で回避（runner が自動設定）。
- BoomerAMG は 1 反復の残差縮小率が約 0.95 と効きが悪い。緩和法・粗化法・強結合閾値の 5 通りと Residual Mode を試したが、1 回の解き時間は約 3 s（定常）で改善しなかった。
- 原因：薄いプリズム（TES 層厚 0.16 µm、面内 8 µm）では z 方向項が層内節点を質量行列状の正の値で結び、非対角の約 29% が正になる。遅い誤差モードは TES・メンブレン・Stycast に集中（pyamg で確認）。
- 対策：`Phase24 Prism Lumped Preconditioner = True`。プリズム剛性の z／面内の質量状結合を集中化した行列（M 行列に近い）を `PrecValues` として作り、BoomerAMG だけに使う。PCG は元の整合行列で解くので解は不変（m6a と最大 0.0025 µA）。集中化を離散化そのものに使うと G が 56% 増えるため不可。
- 上流 `SolveHypre.c` は前処理用行列を渡すと PCG もその行列で解いていたため、PCG は真の A、BoomerAMG は前処理用行列でセットアップ・適用するよう修正。
- 定常（6 ランク）の PCG 反復：216 →（集中化前処理 + relax 18）148 →（+ relax 8 / coarsen 8）60。系列と SIF は `transient_75ms/mpi6_bdf2_lumpprec_m6b*`。
- 前処理用行列があると Phase24 の AMG 使い回しが無効になる（残課題）。

### 設定・メッシュによる追加短縮（m6b 基準、6 ランク、2026-09-27）

| 実行 | 変更 | ステップ | 線形解き回数 | 実時間 | m6b との差（最大） | 判定 |
|---|---|---|---|---|---|---|
| m6b | − | 282 | 883 | 1383 s | − | − |
| nl6 | NL 許容 1e-8→1e-6、TES 温度許容 1e-8→1e-7 K | 282 | 588 | 1062 s | 0.0026 µA | 可 |
| c2 | ステップ上限を全区間 2 倍（`75ms_bdf2c`） | 153 | 563 | 1020 s | 0.021 µA | 不可（立ち上がり 0–100 µs） |
| nlh | nl6 + 500 µs 以降 2 倍（`75ms_bdf2h`） | 198 | 457 | 773 s | 0.014 µA | 不可（5–15 ms） |
| **nlh15** | nl6 + 500 µs 以降 1.5 倍（`75ms_bdf2h15`） | 225 | 494 | **961 s** | **0.0041 µA** | **採用** |

- nlh15：ピーク 7.6851 µA @429.1 µs、COMSOL との差 最大 0.091 / RMSE 0.030 µA（m6b と同等）。

### 前処理の使い回し（`Phase24 Precond Matrix Reuse = True`、nlh15 と同じ設定に追加）

HYPRE コンテナ（IJ 行列・PCG・前処理）を保持し、A は毎回値だけ更新、前処理用行列と BoomerAMG は「使用回数が上限（既定 20）」または「PCG 反復が作り直し直後の 1.5 倍かつ +5 回を超えた」ときだけ作り直す（`solvehypre11`）。作り直し判定はランク間で MPI_MAX をとって揃える（反復回数はランク 0 しか知らない）。PCG の初期値は毎回 Elmer の現在値（Elmer は解くたびに対角スケーリングし直すため、HYPRE 内の前回解は初期値として不適）。

| 実行 | 実時間 | AMG 作り直し | PCG 反復合計 | m6b との差（最大） |
|---|---|---|---|---|
| nlh15（使い回しなし） | 961 s | 494 回 | 20,007 | 0.0041 µA |
| pmr（初期値が HYPRE 内の前回解） | 874 s | 46 回 | 25,294 | 0.0042 µA |
| pmr12（同、作り直し閾値 1.2 倍） | 851 s | 88 回 | 25,121 | 0.0042 µA |
| **pmrx（初期値 = Elmer の現在値）** | **736 s** | 45 回 | 20,177 | 0.0043 µA |

- 定常（6 ランク）：47 回の解きで作り直し 3 回、138.678 µA、109 s → 100 s。

### 時間方向の予測器と電力残差判定（2026-09-28）

- 予測器（`phase24_bdf2_predictor`、SolverUtils の可変刻み外挿 + HeatSolve の外挿／TES 電力予測）で値がずれた原因は収束不足ではなく、TES 回路が新しいステップの最初の呼び出しで前ステップを確定するとき、予測器が上書きした外挿温度場を読んでいたこと。系列の電流とインダクタ履歴（PreviousCircuitCurrent）が毎ステップ予測温度で計算され、反復を増やしても消えない系統誤差になっていた（立ち上がりが速すぎる：30 µs で 0.309 対 0.244 µA）。
- 修正：確定に使う TES 温度を前ステップの収束解 `PrevValues(:,1)` から計算（`TESMeanPreviousStep`）。予測器なしでは従来と同一（1 ms 区間で pmrx と差 0）。
- 収束判定 `TES Inner Circuit Power Residual Tolerance`（相対、5e-6）：解いた後の TES 平均温度で回路則を評価し、使った電力との差で判定（増分判定は良い初期値にだまされうる）。予測が十分なら 1 回で抜けられる。
- 前ステップの dT/dP を流用する `TES Inner Circuit Newton Slope Reuse` は、傾きの過大評価で電力修正が停滞する（上限 120 回）ため不採用。
- 注意：`run_phase24_gate4_5_nomortar.py` は BDF2 のとき予測器を既定で有効にする。

| 実行 | 実時間 | 線形解き | m6b との差（最大） | COMSOL との差（最大 / RMSE） |
|---|---|---|---|---|
| pmrx（予測器なし） | 736 s | 494 | 0.0043 µA | 0.091 / 0.030 µA |
| prd（予測器 + 電力残差 5e-6） | 668 s | 466 | 0.0038 µA | 0.090 / 0.030 µA |
| **pre（予測器 + TES 温度誤差判定 1e-7 K + 温度からの電力予測）** | **499 s** | 380 | **0.0041 µA** | **0.090 / 0.030 µA** |

### 1 反復で抜けるための判定（`Phase24 Coupled Temperature Tolerance`、K）

prd では 225 ステップ中 213 が 2 反復だった。1 回目で止めていたのは、電力残差（相対 5e-6、電熱フィードバックを無視した過度に厳しい値）と、温度場全体の増分（吸収体を含む予測誤差を見ていて、非線形性を表していない）。両者を TES 温度の誤差見積もりに置き換えた。

- 回路：|F| / |dP/dT|（負の電熱フィードバック dP/dT < 0 で、1 ステップの熱応答 S を知らずに成り立つ上限）。L/R より十分短いステップ（パルス前後の ns 刻み）は定電流駆動で dP/dT > 0 となるため、反復間の TES 温度変化で判定。
- k(T)：解き始めに使った膜熱伝導率と現在温度での値の最大相対差 × 膜内温度差。
- `TES Inner Circuit Predict From Temperature`：ステップ開始時の電力を、外挿した TES 温度での回路則の解にする（電力の線形外挿の代わり）。1 ms 区間で 1 反復のステップ 43 → 53。
- 75 ms（pre）：1 反復 88 ステップ、2 反復 126、3 反復以上 11。

### 適応時間刻み（`Phase24 Step Controller`、2026-09-28）と予測器上限の緩和

- 実装：HeatSolve が各ステップ後に、直前 3 ステップの TES 温度からの 2 次予測と収束値の差を BDF2 局所誤差（一定刻みで差の 2/9）に換算し、次の刻み dt·0.9·(tol/LTE)^(1/3)（拡大は最大 1.5 倍、1.0〜1.2 倍は据え置き）を提案。ElmerSolver は最後の区間だけその値を使い、終了時刻で `Exit Condition` を立てる。窓 `75ms_adapt500`（500 µs までは 75ms_bdf2h15 の固定刻み）。
- パルス直後から適応させると、TES がまだ反応していない時期に刻みを広げてしまい精度が落ちる（1 ms 区間で pmrx 差 0.0094 µA）。
- 減衰区間のみ適応（ad1）：225 → 184 ステップだが線形解き 380 → 370、459 s、m6b 差 0.0070 µA。刻みが変わるステップで予測器が上限（`Phase24 BDF2 Predictor Max Delta`、既定 1e-5 K）を超えて切れ、3〜4 反復になる。
- 上限を 1e-3 K に緩めると、適応版（ad3）は非線形許容（1e-7 K）が刻み制御の許容（2e-8 K）より大きいため誤差見積もりが収束誤差に埋もれ、刻みが 50〜150 µs に張り付いて 682 ステップ・1026 s に悪化。固定刻みでは有効。
- 結論：適応刻みは手作りの固定刻みに対して正味の利得なし（オプションとして残す）。
- 前ステップ dT/dP の流用は、割線を電力変化が相対 1e-7 超のときだけ採用し、0 < S ≤ T/P に制限する安全策を入れて有効化。

| 実行 | 実時間 | ステップ | 線形解き | m6b との差 |
|---|---|---|---|---|
| pre | 499 s | 225 | 380 | 0.0041 µA |
| ad1（適応、上限 1e-5 K） | 459 s | 184 | 370 | 0.0070 µA |
| ad3（適応、上限 1e-3 K） | 1026 s | 682 | 765 | 0.0030 µA |
| pr2（固定、上限 1e-3 K + 傾き流用） | 464 s | 225 | 363 | 0.0041 µA |

### PCG の許容値（‖b‖ 相対、前処理付き残差ノルム）

1 ms 区間（1e-12 基準）：1e-11 で平均 30.7 → 24.8 反復・差 0.0006 µA、1e-10 で 14.9 反復・0.0032 µA、1e-9 で 7.8 反復・0.033 µA（不可）。いずれも反復 0 回の誤収束はなし。75 ms 全体：

| 線形許容値 | 実時間 | 線形解き | 平均 PCG 反復 | m6b との差 |
|---|---|---|---|---|
| 1e-12（pr2） | 464 s | 363 | 38.5 | 0.0041 µA |
| **3e-11** | **364 s** | 385 | 25.9 | **0.0036 µA** |
| 1e-10 | 1082 s | 1341 | 18.9 | 0.0143 µA |

### BoomerAMG の演算量（`BoomerAMG P Max Elements` / `Truncation Factor` / `Aggressive Levels` / `Aggressive Interpolation Type` / `Print Level` を追加）

現状の階層：7 段（33.5 万 → 13 万 → 4.3 万 → 1.4 万 → 3 千 → 498 → 2）、operator complexity 1.83、補間 平均 2.2 要素/行（HYPRE 既定で P max elements = 4 が既に効いている）。1 ms 区間（線形許容 3e-11）：

| 設定 | 時間 | 平均 PCG 反復 | 線形解き合計 |
|---|---|---|---|
| 現状 | 159 s | 21.1 | 111 s |
| P max 2 | 155 s | 22.9 | 104 s |
| P max 3 + trunc 0.2 | 157 s | 22.4 | 108 s |
| **積極的粗化 1 段（multipass）** | **136 s** | 23.2 | **90 s** |
| 積極的粗化 1 段（ext+i） | 200 s | 39.4 | 152 s |
| 積極的粗化 2 段 | 180 s | 29.5 | 130 s |
| 積極的粗化 1 段 + P max 2 | 142 s | 24.0 | 93 s |

75 ms（積極的粗化 1 段）：364 s → **320 s**、m6b との差 0.0036 → 0.0065 µA（COMSOL 差 最大 0.091 / RMSE 0.030 µA）。

### GPU 切り出し実験（`gpu_probe/`）

h8 の線形系を HIP 版 HYPRE と rocALUTION で解いた。WSL では CPU 6 ランク（0.56 s/解き）より遅い（HYPRE GPU 0.78 s、rocALUTION RS-AMG 0.95 s）。原因は同期遅延（起動＋同期 313 µs、読み戻し 210 µs）で、PCG + AMG は 1 反復に数十回同期する。GPU 路線は打ち切り、CPU で続ける。詳細は `gpu_probe/README.md`。

1e-10 は線形解きの誤差が非線形の TES 温度判定（1e-7 K）を超え、非線形反復が収束しにくくなる（線形解き約 4 倍）。採用は 3e-11（ピーク 7.6843 µA @429.1 µs、COMSOL 差 最大 0.091 / RMSE 0.030 µA）。
- Stycast 16 層メッシュ（`mesh_hybrid_fullconf_h8s16`、節点 −17%）：定常値は同じ（138.677 µA）だが定常の計算時間も変わらず（116 s 対 109 s）、不採用。
- 12 ランク（SMT、`--oversubscribe`）：最初の線形解きから進まず中止、不採用。
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

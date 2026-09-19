# Phase24 COMSOL parity 現状整理

更新日: 2026-09-19 (JST)

この文書は、Phase24 の Gate3 以降について、2026-09-19 時点で確認できている
結果、切り分け、未確定事項をまとめた現行判断資料である。元の設計書に残っている
旧 no-mortar ルートの記述と混在しないよう、今回の測定結果を優先する。

## 結論

現在の判断は次の通り。

1. **Gate3 の定常解は PASS**。32 層 Stycast の非共形メッシュと mortar を使った
   HYPRE は、COMSOL の定常電流差 0.6% 以内に入っている。
2. **HYPRE の過渡 COMSOL parity は FAIL**。40 µs の途中で、COMSOL が示す TES
   応答を再現せず、100 µs まで延長する意味がないことを 20 µs 時点で確認した。
3. **GPU は原因ではない**。CPU-HYPRE と GPU-HYPRE は同じ遅い波形を生成し、GPU
   parity 自体は PASS している。
4. **今回の新しい非共形32層+mortarでは、MUMPSも過渡parityに失敗**した。79 stepを
   正常完走したが、pulse前から電流が定常基準と一致せず、最大差17.826591 µA、RMSE
   12.445186 µAとなった。
5. 旧 CPU MUMPS の良好な100 µs実績は、今回の新メッシュへそのまま移植できない。
   COMSOLに近い波形を得るには、旧 hybrid/mortar ルートのメッシュ・界面処理・初期化を
   同一条件で再現する必要がある。

## Gate 判定の現状

| Gate | 現時点の判定 | 根拠 |
|---|---|---|
| Gate3 定常 HYPRE | **PASS** | 非共形 32 層 + mortar。143.567589 µA、COMSOL 差 +0.358282%。 |
| Gate4 HYPRE 過渡 parity | **FAIL** | 40 µs 部分比較で最大差 0.617301 µA、RMSE 0.217332 µA。目標は各 0.05 / 0.03 µA 以下。 |
| Gate4 MUMPS 過渡 parity | **FAIL** | 非共形32層+mortarで79 step完走。最大差17.826591 µA、RMSE 12.445186 µA。初期baseline自体が不一致。 |
| Gate5 HYPRE 長時間安定性 | **物理parity未達のため採用不可** | 旧 refine20/no-mortar の完走記録はあるが、COMSOL波形parityを満たさないため現行物理ルートの成功根拠にはしない。 |
| Gate6 GPU parity | **solver parity PASS、COMSOL 100 µs parity FAIL** | GPUはCPU-HYPREと一致するが、両者ともCOMSOLの遅い過渡応答を再現しない。 |

## 主要な数値

### Gate3: 非共形 32 層 Stycast + mortar + HYPRE

対象ケースは `case_phase24_g3_s32m_s32m2_10f17b`。

| 項目 | 値 |
|---|---:|
| HYPRE 定常電流 | 143.567589 µA |
| COMSOL 定常電流 | 143.055049 µA |
| HYPRE - COMSOL | +0.358282% |
| 同条件 MUMPS基準 | 143.534427 µA |
| HYPRE - MUMPS | +0.023104% |
| HYPRE最終線形残差 | 1.99717e-6 |
| Elmer非線形残差 | 3.99394e-12 |
| TES回路残差 | 4.85e-14 W |

線形許容誤差は `2e-6`。定常の DC 電流と回路残差については十分な証跡があり、
Gate3 は COMSOL 優先判定で PASS とした。

### HYPRE 過渡: 共形 no-mortar / 共形 mortar

32 層 Stycast の共形メッシュについて、内部境界面を除去した no-mortar と、同じ
共形メッシュに mortar を指定したケースを比較した。両者の波形は実質的に同一で、
この操作だけでは問題は解消しなかった。

40 µs 窓の比較は pulse 後 0--38.75 µs。各モデルの pre-pulse baseline を差し引いて
比較している。

| 指標 | COMSOL | HYPRE |
|---|---:|---:|
| pre-pulse baseline | 143.055049 µA | 143.534426 µA (+0.335099%) |
| 20 µs の電流低下 | 0.0365617 µA | 0.0001277 µA |
| 最大絶対差 | - | 0.617301 µA @ 38.75 µs |
| RMSE | - | 0.217332 µA |
| t10 | 41.7428 µs | 未到達 |
| t50 | 92.8374 µs | 未到達 |

この結果により、40 µs の途中で Gate4 不合格を判断できる。100 µs HYPRE を続けても
parity成功にはつながらないため、長時間計算を打ち切って切り分けを優先した。

### HYPRE 過渡: 非共形 32 層 Stycast + mortar

物理的な界面結合を担保するため、ElmerGrid が `non-conforming` と判定する別メッシュ
を作成し、mortar を有効化した。

- 線形許容誤差 `2e-6` では Gate3 定常は通るが、過渡では 20 µs 付近まで TES 電流が
  ほぼ一定で、COMSOL の微小な応答を失う。
- 非線形許容誤差を `1e-6` に緩めても、TES 応答が平坦なまま 48/79 step 付近まで
  進むだけだった。
- 線形許容誤差を `2e-8` に戻すと、残差が約 `1.8e-6` で停滞し、4000 iteration
  では収束しなかった。
- 共形 no-mortar で `1e-9` まで締めた試験でも TES 応答は回復しなかった。

つまり、単に iteration 数を増やす、GPUを使う、内部境界面を削除する、という対策では
解決していない。

### MUMPS 40 µs: 今回の非共形32層+mortar

対象ケースは `case_phase24_g45_s32m_40us_mumps_mortar_mumps`。ElmerSolverは79 stepを
18:05:40 JSTに正常完走した。

| 指標 | 値 |
|---|---:|
| pre-pulse baseline (COMSOL) | 143.055049 µA |
| MUMPS baseline | 165.706588 µA (+15.834142%) |
| 20 µs時点のCOMSOL電流低下 | +0.0004876 µA |
| 20 µs時点のMUMPS電流低下 | -12.644986 µA |
| 最大絶対波形差 | 17.826591 µA @ 38.75 µs |
| RMSE | 12.445186 µA |
| t10 / t50 | いずれも未到達 |

最初のseries出力は pulse 前の `20.018 ms` で、TES温度 `167.9606 mK`、電流
`165.1575 µA` だった。Gate3定常のTES温度 `168.5691 mK`、電流 `143.5676 µA` と
一致していない。このため、今回のMUMPS波形はCOMSOLとの物理差を直接評価する前に、
初期状態／restart／界面離散化の整合性で不合格と判断する。

### 旧 MUMPS 過渡実績

旧 CPU MUMPS + production-v2 hybrid/mortar 系の 100 µs 実績は次の通り。

| 指標 | 値 |
|---|---:|
| MUMPS baseline | 143.777852 µA (+0.505% vs COMSOL) |
| 最大絶対波形差 | 0.028714 µA |
| RMSE | 0.015865 µA |
| t10 (COMSOL / MUMPS) | 41.7428 / 41.2819 µs |
| t50 (COMSOL / MUMPS) | 92.8374 / 92.8493 µs |

この結果は旧ルートについては「MUMPSの方がCOMSOLに近い」ことを示す。しかし、今回の
非共形32層+mortarでは同じ傾向が再現されなかった。したがって、MUMPSをbackendだけ
差し替えれば解決する、とは判断できない。

## 原因の切り分け

HYPREの平坦な過渡については、**mortarを含む疎行列の条件性と、微小な過渡信号に
対するHYPREの反復残差不足**が最も整合する。一方、今回のMUMPS結果はそれとは別に、
**非共形32層+mortarの過渡初期化または界面離散化が、restart定常状態と一致していない**
ことを示している。

パルス後 20 µs の COMSOL の電流低下は約 0.0366 µA であり、定常電流 143 µA に
対して非常に小さい。一方、HYPRE の定常 Gate3 は最終線形残差約 `2e-6` で成立して
いる。定常 DC 値を合わせるには十分でも、DC値に重なる小さな温度・電流変化を保持する
には不足し得る。

さらに、残差を `2e-8` へ締めると同じ問題が収束せず、約 `1.8e-6` で停滞する。この
挙動は「物理モデルが全く違う」よりも、mortar/熱-回路連成を含む反復線形系で、要求
精度まで解を解像できていないことを示す。MUMPSが近い波形を返すのは、直接法がこの
反復停止・前処理依存性を避けているためと考えるのが妥当である。

次の観測もこの判断を支持する。

- CPU-HYPREとGPU-HYPREの波形が一致するため、GPU演算誤差は主因ではない。
- 共形メッシュで mortar の有無を切り替えても波形は同じであり、内部境界面の重複だけ
  が主因ではない。
- 32 層 Stycast を導入しても HYPRE過渡は平坦なままであり、Stycast の z分割不足
  だけでは説明できない。
- 旧MUMPSの実績では同じ物理イベントの時間スケールがCOMSOLに近いが、今回の新メッシュ
  ではpulse前から22.65 µAの電流オフセットが生じた。したがって、旧実績の再利用には
  メッシュだけでなくrestartと界面処理の完全な再現が必要である。

## Sol に判断してほしい論点

1. 物理結果を最優先する場合、旧MUMPS良好ケースの構成をそのまま再現する専用ルートを
   別途固定し、今回の32層非共形ルートを採用候補から外すか。
2. HYPREを必須とする場合、次の調査対象を「solver許容値」ではなく、mortar連成を
   含む行列の条件数、スケーリング、ブロック前処理、回路自由度の扱いに限定するか。
3. 現行の Gate4/5 の成功条件を、HYPRE専用ではなく「COMSOL parityを満たす採用
   backend」として定義し直すか。

現時点の推奨は、**今回の32層非共形+mortarルートを成功扱いにせず、旧MUMPS良好ケース
との差分を、restart状態・pulse前定常・界面自由度の3点に分けて再現する**、である。
HYPREのsolver調整だけでなく、MUMPSでも初期状態が崩れた事実を優先して扱う。

## 再現可能な成果物

- [Gate3 非共形32層+mortar HYPRE PASS](../artifacts/phase24_gate3_hypre_steady_s32m2/summary.md)
- [共形32層 no-mortar HYPRE 過渡比較](../artifacts/phase24_gate4_5_nomortar/40us/comparison_s32nm_strip/summary.md)
- [共形32層 mortar HYPRE 過渡比較](../artifacts/phase24_gate4_5_nomortar/40us/comparison_stycast32_mortar_partial/summary.md)
- [旧 CPU MUMPS 100 µs 比較](../artifacts/comparison/comsol_cpu_singlepixel_prod_v2_hybrid_100us/summary.md)
- [今回の非共形32層 MUMPS 40 µs 比較](../artifacts/phase24_gate4_5_nomortar/40us/comparison_s32m_mumps/summary.md)
- [Gate6 GPU parity](../artifacts/phase24_gate6_gpu_parity/summary.md)
- [Gate設計書（履歴・基準）](phase24_comsol_parity_gate_design.md)

## MUMPS 40 µs の実行証跡

`case_phase24_g45_s32m_40us_mumps_mortar_mumps` を 2026-09-19 17:34 JST に開始した。
2026-09-19 18:05:40 JSTに ElmerSolver が正常終了し、series.csv と manifest.json が
生成された。比較結果は [summary.md](../artifacts/phase24_gate4_5_nomortar/40us/comparison_s32m_mumps/summary.md)
に保存した。

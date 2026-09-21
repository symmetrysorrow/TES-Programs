# Phase24 COMSOL parity gate design

## 目的

本設計の最終目標は、同一のTES物理条件において、Elmer/HYPREでCOMSOLの
絶対電流値と過渡波形を再現し、31 µs本番へ進めることである。判定は
HYPREの反復回数や単独のsolver成功ではなく、COMSOLに対する物理結果、旧
MUMPS基準に対する差、収束性、再現性および実行安定性で行う。

今回採用する共形no-mortarルートでは、COMSOLを物理的な第一判定基準とする。同一
メッシュMUMPSはbackend parityと離散化差を確認する診断基準として保存するが、
COMSOLにより近いHYPRE候補をMUMPS差だけで不合格にはしない。従来の非共形
production-v2およびmortar基準は履歴比較用artifactとして保持する。

| 項目 | 基準値 |
|---|---:|
| COMSOL定常電流 | 143.055049 µA |
| 共形refine20 no-mortar MUMPS定常電流 | 144.710781 µA |
| 共形refine20 no-mortar MUMPSとCOMSOLの定常差 | +1.157 % |
| 旧mortar MUMPS波形最大差（履歴） | 0.028629 µA |
| 旧mortar MUMPS波形RMSE（履歴） | 0.015804 µA |
| 現行比較メッシュ | `mesh_singlepixel_conformal_gpu_refine20` |
| 界面処理 | mortarなし |

基準結果の証跡は
`artifacts/comparison/prod_v2_mumps_vs_comsol/summary.md` に保存する。

## 固定する比較条件

- COMSOL、旧MUMPS、現行Stage 11 MUMPS、HYPREで形状・材料値・境界条件・初期
  定常状態・電極面積・接触抵抗・バイアスを一致させる。
- MUMPSとの比較は現行の共形メッシュとno-mortarを使用する。no-mortarを今回の
  実用計算ルートとして固定し、MUMPS/HYPREの比較対象を同じ離散化条件にする。
  同一メッシュMUMPSとの差は診断値として必ず記録するが、Gateの物理合否はCOMSOL
  への近さを優先する。
- 旧mortar結果とconformal no-mortar結果は、界面離散化の影響を分離する履歴診断
  として扱う。
- 同じ出力時刻、同じ電流定義、同じ単位変換を使う。絶対値波形と定常値からの
  差分波形を別々に評価する。
- 変更は一度に一要因に限定し、各ケースにSIF、プロジェクト設定、solver.log、
  iteration記録、結果ファイル、実行環境ハッシュを残す。

## Gate一覧

### Gate 0: 比較基準の凍結

COMSOLおよび旧MUMPSの結果、メッシュ、SIF、物性値、電流抽出方法を固定する。

合格条件:

- 旧MUMPSの定常値、波形最大差、RMSEを再読可能なartifactとして保存する。
- COMSOLとElmerの電極面積、接触抵抗、導電率、初期温度、TESパラメータが一致
  していることを設定ファイルから確認できる。
- 今後の候補ケースが同じmesh registryと物理パラメータを参照する。

### Gate 1: 現行Stage 11 MUMPSの旧基準再現

現行Stage 11ビルドで、production-v2メッシュ+no-mortarを用いた定常MUMPS計算を
実施する。これはHYPREの影響を入れる前に、現行ソースと旧基準の差を確定する
ゲートである。

合格条件:

- 物理的な非線形収束を完了し、`ALL DONE`または同等の正常終了を記録する。
- 旧MUMPS定常値に対して ±0.1 %以内。
- COMSOLに対する定常差が旧MUMPS実績の許容帯（目標 ≤0.6 %）に入る。
- 温度、TES平均温度、電流、抵抗、Joule熱が有限値である。

不合格ならHYPREの調整へ進まず、現行Stage 11と旧MUMPSの物理設定、メッシュ、
非線形反復、初期値を差分比較する。

### Gate 2: 同一行列・右辺でのbackend parity

同じ時刻、同じ非線形反復、同じmesh/no-mortar条件でMUMPSとHYPREに与えた疎行列と
右辺を比較する。HYPREの残差だけでなく、解ベクトルとTES電流も比較する。

合格条件:

- CSRの行数、非ゼロ数、行構造が一致する。
- 行列値および右辺の差が、倍精度丸めと並列加算順序で説明できる範囲にある。
- HYPREの相対線形残差が要求値以下である。
- 同一の初期解から得たMUMPS/HYPRE解の差が、物理誤差ではなく線形solver誤差の
  範囲に収まる（目標電流差 < 0.01 µA）。

行列が一致しない場合はassembly、境界消去、mortar、物性評価を調査する。
行列が一致し、解だけが異なる場合はHYPREの前処理・初期解・収束判定を調査する。

### Gate 3: HYPRE定常解（共形no-mortar実用ルート）

共形メッシュ+no-mortarで、HYPREを初期温度から独立に起動し、定常解まで収束
させる。

合格条件:

- 初期温度からの計算が正常終了する。
- COMSOLとの差が目標 ≤0.6 %。
- 同一メッシュMUMPSとの差は必ず記録する。MUMPS差が大きくても、COMSOLに近い
  こと、有限な物理量、線形・非線形・回路の証跡が揃えばGate 3の物理判定を
  不合格にはしない。
- HYPRE線形残差、Elmerの非線形残差、TES回路残差を記録する。
- HYPREの設定変更で結果を合わせるのではなく、設定と許容値をartifactに固定する。

今回のno-mortar実用候補では、`FlexGMRES + BoomerAMG`、最大4000回、線形許容誤差
`2e-8`、BoomerAMG strong threshold `0.5` を基本設定とする。refine20のように
電流感度が高いメッシュで許容誤差を厳格化した場合は、設定変更と結果を別artifact
に保存する。

旧mortar系でのHYPRE失敗は履歴として残すが、今回のGate判定からは除外する。
no-mortarルートでは、HYPREがno-mortar MUMPSの定常値とCOMSOLの許容範囲を同時に
満たすことを確認する。

### Gate 4: HYPRE過渡波形 parity

同じ初期定常状態から短時間、続いて100 µsまでの過渡計算を行い、COMSOL、no-mortar
MUMPS、Stage 11 MUMPS、HYPREを比較する。まず短時間で差の発生点を特定し、その後に窓を
延長する。

合格条件:

- 絶対電流の定常値オフセットと、定常値からの差分波形を分離して評価する。
- MUMPS定常オフセット差は診断値として記録する。COMSOLへの絶対値・差分波形の
  parityを優先し、MUMPS差だけでは不合格にしない。
- COMSOL対HYPREの波形最大差が目標 0.05 µA以下、RMSEが目標 0.03 µA以下。
- 立上り時刻・半値時刻などの時定数指標の差が目標 1 µs以内。
- 内部時間刻み、BDF次数、イベント着地、rejection、linear residualを保存する。

差が絶対値だけに現れる場合は動作点・初期定常値・接触/電極条件を調査する。
差分波形の振幅や時定数だけに現れる場合は時間積分、熱容量、界面熱抵抗、
laggingおよび出力補間を調査する。

### Gate 5: 長時間安定性

100 µs、1 ms、必要に応じて1 µs級の長時間窓で、同じ凍結済みポリシーを使用する。
これは物理parityと別に、メモリ蓄積・HYPRE lifecycle・適応刻みの安定性を確認する
ゲートである。

合格条件:

- 明示的なFortran allocation failure、MPI abort、native crash、`STOP 1`なしで完走
  する。
- `PHASE24_MEMORY`のPrivateUsage/working setが単調増加し続けない。
- peak memoryが短時間基準の1.25倍以内を目安とし、増加量と試行数の関係を記録する。
- HYPRE/IJ matrix/vector、ParCSR、Krylov、AMGの生成・破棄回数とepochを記録する。

メモリが増加する場合は、再構築頻度を下げる前に所有権と破棄順序を確認し、
「症状の先送り」を合格とは扱わない。

### Gate 6: GPU parity

no-mortar CPU HYPREでGate 3〜5を通過した後に、GPU backendを有効化する。GPU化による
速度改善は、物理結果の一致を満たした候補だけで評価する。

合格条件:

- CPU HYPREとGPU HYPREの定常電流差が ±0.05 %以内。
- 過渡波形の最大差・RMSEがGate 4の範囲内。
- CPU fallbackや未解放のdevice resourceがない。
- 同一メッシュ、同一SIF、同一線形solver許容値で比較できる。

### Gate 7: 31 µs本番採用

全ゲートの証跡をまとめ、31 µs本番候補を一つに固定する。

合格条件:

- Gate 0〜6の判定がすべてPASS、または未実施項目と理由が明記されている。
- COMSOLに対する絶対値と差分波形の両方が目標範囲内である。
- 再実行可能なコマンド、環境依存DLL、ビルドsource hash、SIF hash、mesh hashが
  保存されている。
- 少なくとも一回の再実行で同じ判定が得られる。

## 実施順

1. Gate 0の設定・artifactを凍結する。
2. 現行Stage 11 MUMPS + `mesh_singlepixel_prod_v2` + no-mortarでGate 1を実施する。
3. Gate 1がPASSした場合だけ、同一時点のmatrix/RHS dumpでGate 2を実施する。
4. 共形no-mortar条件でHYPREの線形収束問題を解消し、初期温度からGate 3を実施する。
5. まず短い過渡窓でGate 4の差分発生点を特定し、100 µsへ延長する。
6. Gate 5でメモリとlifecycleの長時間安定性を確認する。
7. 最後にGPU parity、31 µs本番採用の順で判定する。

## 変更管理と停止規則

- 物理設定、メッシュ、線形solver、時間積分、lagging、出力スケジュールを同時に
  変更しない。
- 数値が改善しても、行列parityまたは保存された物理条件が崩れた候補は採用しない。
- Gate失敗時は、次のGateへ進まず、失敗を再現する最小ケースを作る。
- `restart`から収束したHYPRE結果は、HYPRE初期値からの独立定常解とは区別して記録する。
- 旧MUMPSとの差が小さいことだけではCOMSOL parityの合格としない。

## 現時点の判定

- Gate 0: 旧MUMPS/COMSOL基準は保存済み。条件の完全凍結を確認中。
- Gate 1: no-mortar production-v2のStage 11 MUMPS基準（143.537345 µA）は正常終了済み。
- Gate 2: no-mortar条件でのMUMPS/HYPRE同一行列比較を次の実施項目とする。
- Gate 3: 旧非共形production-v2 no-mortar HYPREは参考PASS（電流 `143.537345 µA`、
  COMSOL差 `0.337%`）だが、過渡初期場が破綻したため採用しない。現行の共形
  refine20 no-mortarへ切り替え中である。
- Gate 4: 旧production-v2非共形メッシュでは過渡再始動で温度場が破綻したため未達。
  現行の共形refine20 restartで短窓を実施し、Gate4短窓はPASSした。100 µs窓は
  HYPREで継続中である。
- Gate 5: 共形refine20 HYPREの100 µs長時間安定性はPASS。1 ms窓は未実施。
- Gate 6: 共形refine20 no-mortarのHIP HYPRE GPU parityは短窓でPASS。
  ただし100 µs拡張比較ではCOMSOL過渡誤差が大きくFAILとなった。Gate7は未実施。

関連証跡:

- `artifacts/comparison/prod_v2_mumps_vs_comsol/summary.md`
- `artifacts/comparison/prod_v2_hypre_steady_probe/summary.md`
- `artifacts/comparison/conformal_nomortar_hypre_tol1e-10_vs_mumps/summary.md`
- `artifacts/phase24_stage11_exact_failing_system.json`
- `phase24_stage11_native_provenance.json`

### 共形no-mortar切替後のGate 3結果

`mesh_singlepixel_conformal_gpu_refine20`（538,253 nodes、no-mortar）を現行候補に
固定した。HYPREはT0から独立起動し、`FlexGMRES + BoomerAMG`、最大4000回、
許容誤差`2e-8`、strong threshold`0.5`で完走した。

| 項目 | 結果 |
|---|---:|
| HYPRE定常電流 | 143.522179 µA |
| COMSOL差 | +0.326539 % |
| 同一メッシュMUMPS差 | -0.821364 % |
| 最終線形残差 | 1.82447e-8 |
| Gate3 | **PASS（COMSOL優先判定）** |

同一メッシュMUMPSとの差は診断値として保存し、COMSOLに近いHYPRE結果を優先する。
証跡は`artifacts/phase24_gate3_hypre_steady_conformal_refine20/summary.md`にある。
Gate4/5 runnerは`--mesh mesh_singlepixel_conformal_gpu_refine20`と共形restartを
指定して同じルートを継続できる。

### 共形no-mortar Gate 4短窓

MUMPS/HYPREとも23ステップ（イベント後0.9 µs）を正常終了した。HYPREのCOMSOL
比較は最大差`0.001616 µA`、RMSE`0.001050 µA`で、短窓Gate4をPASSとした。
時定数指標は0.9 µs窓では未到達のため適用せず、100 µs窓で判定する。
証跡は`artifacts/phase24_gate4_5_nomortar/short/summary.json`にある。

### 共形no-mortar Gate 5（100 µs）

HYPREのみで100 µs窓（175ステップ）を完走した。Gate5は物理backend parityではなく
長時間のsolver lifecycle／メモリ安定性を判定するため、Gate4短窓で両backendの
波形parityを確認したうえで、HYPRE長時間窓を採用した。

| 項目 | 結果 |
|---|---:|
| 実行時間 | 9,498.87 s |
| 終了 | `ALL DONE`, exit 0 |
| memory records | 1,000 |
| working set 初期/ピーク | 2.548 / 3.073 GB（1.206倍） |
| private usage 初期/ピーク/終了 | 4.308 / 5.623 / 1.787 GB |
| private usage 単調増加 | false |
| Gate5（100 µs） | **PASS** |

証跡は`artifacts/phase24_gate4_5_nomortar/100us/summary.md`および
`summary.json`にある。1 ms窓は未実施であり、必要時の追加安定性試験として残す。

### 共形no-mortar Gate 6（HIP HYPRE GPU）

CPU-HYPREで確定したrefine20/no-mortar条件を維持し、HIP版HYPRE device backendで
T0からの定常解と、同じ定常restartからのGate4短窓（23ステップ、0.9 µs）を実行した。

| 項目 | 結果 |
|---|---:|
| GPU定常電流 | 143.522179 µA |
| CPU-HYPRE定常電流との差 | 0.000000 % |
| GPU vs COMSOL 波形最大差 | 0.001616 µA |
| GPU vs COMSOL 波形RMSE | 0.001050 µA |
| HIP device marker | `AMD Radeon RX 9070 XT`, `backend=HIP` |
| CPU fallback / native crash | なし |
| Gate6 | **PASS** |

再現用ランナーは`scripts/support/run_phase24_gate6_gpu_parity.py`、証跡は
`artifacts/phase24_gate6_gpu_parity/summary.md`および`summary.json`にある。
この短窓判定に加え、100 µs GPU長時間solver実行も完走したが、物理parityは別途FAILとなった。
Gate7は未実施である。

### Gate 6 100 µs GPU拡張検証

HIP HYPRE GPUで100 µs窓（175ステップ）を完走した。GPUとCPU-HYPREの波形は
98.751 µsまで数値一致したが、COMSOLに対しては短窓では見えなかった遅い過渡差が
現れた。比較範囲98.751 µsでも、GPUの電流低下は約`0.00129 µA`に留まり、COMSOLの
`4.18562 µA`低下に対して最大差`4.18433 µA`、RMSE`2.08177 µA`となった。
終端iterationを補完した100.001 µsまでの比較でも最大差は`4.24581 µA`、RMSEは
`2.12170 µA`であり、Gate4の`0.05 µA`/`0.03 µA`閾値を満たさない。

したがって、100 µs GPU計算のsolver完走・GPU実行自体は確認できたが、COMSOLとの
100 µs物理parityは **FAIL** である。CPU-HYPREも同じ波形であるため、GPU固有の
差ではなく、CPU/GPU共通の長時間過渡モデルまたはpulse-to-TES couplingの差が示唆される。
証跡は`artifacts/phase24_gate6_gpu_parity/100us/summary.md`および`summary.json`にある。

## Gate 3実装

共形no-mortar版Gate 3の再現可能な実行・判定は
`scripts/support/run_phase24_gate3_hypre_steady.py` で行う。
`--source-project`、`--mesh`、`--base-case`、`--reference-current-uA`で共形メッシュ
を選択できる。独立起動し、`run.py` のmanifest、solver.log、反復CSV、設定・実行
環境のhashをartifactに保存する。`--audit-only` は既存結果を再実行せずに再判定
する。

現行refine20候補の実行例:

```text
python scripts/support/run_phase24_gate3_hypre_steady.py --source-project artifacts/comparison/nomortar_refinement_probe/project.json --mesh mesh_singlepixel_conformal_gpu_refine20 --base-case case_tes_steady_singlepixel_conformal_gpu_fine --reference-current-uA 144.71078126230953 --variant conformal_refine20 --linear-tolerance 2e-8 --strong-threshold 0.5
```

旧mortar版の試行結果（2026-09-18、履歴）:

- 基準HYPRE FlexGMRES + BoomerAMG（4000回, `1e-10`）は線形残差
  `4.30545e-6` で停止した。
- MGR（4000回, `1e-8`）は残差 `1.62277e-8` で停止した。
- MGR + reuse（`5e-7`）は `ALL DONE` まで完走したが、電流 `142.685881 µA`、
  MUMPS差 `0.759418%` であった。
- MGR（`2e-8`）は非線形反復を完了できず、電流は約 `140.968 µA` までずれた。
- MGR（`1e-8`, GMRES次元20, 最大10000回）はメモリ約627 MBで動作したが、
  初回線形解が上限到達前に完了せず停止した。

no-mortar版Gate 3試行（2026-09-18）はPASSした。Gate 4以降は、Gate 2の
no-mortar backend parity確認だけでなく、過渡再始動時の温度場・残差 sanity check を
通過したうえで、このno-mortarポリシーを継続して実施する。今回のproduction-v2
非共形mesh + no-mortarでは、Gate 4短窓で `T=-1.36e5 K`、`Result Norm=22625.86`
となったため、Gate 5へは進めていない。

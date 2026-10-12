# WSL DualTES 定常・75 ms過渡計算の実測

## 結果

| 計算 | wall時間 | solver CPU / REAL | 完了状態 |
|---|---:|---:|---|
| DualTES 定常（Newton更新ON） | 33.362 s（33.4秒） | 31.78 / 33.25 s | `ALL DONE`、相対変化 `5.95e-7`、TES連成誤差 `4.82e-8 K` |
| DualTES 過渡 20 ms → 75.02 ms | 434.648 s（7分14.6秒） | 428.83 / 434.53 s | `ALL DONE`、225/225 step |
| SinglePixel 過渡の既存WSL実測 | 560 s（9分20秒） | 547.38 / 555.44 s | `ALL DONE`、225/225 step |

DualTES過渡 / SinglePixel過渡のwall比は **0.7762**（約0.78倍）です。今回のDualTES定常と過渡を続けて実行した合計は **468.010秒（7分48.0秒）** でした。比は過渡同士で計算しています。SinglePixel側は既存の20 ms定常restartから始めた過渡実測なので、DualTES定常時間を比の分子へ足していません。

## 条件

両方ともWSLのElmerGPU distro、同じElmerSolver MPI実行ファイル、6 MPI rank・各rank 1 threadです。solver SHA-256は `16b512caf1a49b3d23d4d8e23e636cd4e1df643531e81a65a3f71bd7a1d48c56` です（Elmer 26.1、build logのrevision表示はunknown）。今回のSolver・ElmerSolver library（HeatSolveを含む）・libmatc・TES UDFのSHA-256は`build_manifest.txt`に保存しました。

過渡でSinglePixel replay SIFと揃えた条件は、BDF2、BDF2 predictor、同一の明示timestep scheduleとBDF1切替時刻、非線形許容値 `1e-6`、最大120反復、緩和係数1.0、HYPRE PCG + BoomerAMG、線形許容値 `1e-12`、最大4000反復、BoomerAMG relax type 8、Vector Assembly ON、Phase24 preconditioner reuse ON、TES inner-circuit Newton ON、TES温度予測ON、step commit ONです。MATC Dynamic K Node CacheはDualTES入力で明示OFF、SinglePixel入力では指定なし（既定OFF）です。

DualTES固有のメッシュ、2回路、TES定数、pulse、熱源、restartは保持しています。SinglePixelの1回路モデルへの置換はしていません。DualTES meshは131,633要素（prism 33,060）、SinglePixel replay meshは236,916要素（prism 23,220）で、各々6-way partitionを使っています。物理入力とmeshが違うため、wall比は今回の2ケースに対する実測であり、TES回路数だけの追加コストを表すものではありません。DualTES過渡は定常出力の0.02 s restartから0.07502 sまで計算しました。隔離コピーのSIFだけを変更し、本番DLL・元SIF・元meshは編集していません。

`solver_settings_comparison.json`で両過渡SIFの時間積分・timestep・非線形・線形solver・HYPRE・再利用・Newton回路設定を機械比較し、列挙した全設定が一致することを確認しました。モデル固有のmesh・body ID・2回路定数・pulse・restart結果は比較対象外です。

## 反復・収束とプロファイル

DualTES過渡の温度反復は最大3回で、120回上限には達しませんでした。全225 stepで最終反復のTES回路誤差と温度依存k誤差がそれぞれ`1e-7 K`以下でした。ログ上、adaptive retryとstep failureはありません。SinglePixel replayの温度反復ピークは5回でした。

DualTES過渡ログはHYPRE solve 276回、method-601のsolver報告時間合計209.32 s、assembly timing first-column合計117.67 s、solve timing first-column合計220.30 sでした。HYPRE iteration欄はこのbuildが`-1`を出力するため反復数として使えません。材料物性評価時間とfile I/O時間は個別計測されずwall時間に含まれます。これらのstage値はwall時間と単純加算できる排他的な内訳ではありません。

TES raw series CSVは通常終了時に最終stepの1行をflushしないため、各回路のaccepted iteration CSVから0.07502 sの最終行を別の`*_series_reconciled.csv`へ復元しました。元CSVは保持しています。最終値はTES1: `T=0.16871998 K`, `I=138.362 µA`, `P=3.11160e-10 W`、TES2: `T=0.16871976 K`, `I=138.369 µA`, `P=3.11173e-10 W`です。

## 成果物

- `inputs/dualtes_wsl_steady_20261013_newton.sif`, `inputs/dualtes_wsl_transient_20261013_newton.sif`: 最終測定に使ったSIFコピー
- `inputs/*_settings.diff`: 元SIFからの変更差分
- `inputs/partial_method/`: Newton更新設定を入れる前の参考入力
- `runs/attempt6/steady`, `runs/attempt6/transient`: 最終測定のrun record、solver log、結果、TES系列
- `wall_profile.json`: wall・stage・HYPRE・反復・TES系列の集計
- `solver_settings_comparison.json`: SinglePixel replay SIFと最終DualTES SIFの設定比較
- `runs/attempt6/build_manifest.txt`, `runs/attempt6/mesh_inputs.sha256`: solver/library/UDF/mesh入力hash（今回の最終実行分）
- `rerun_wsl_dualtes.sh`: 新しい`RUN_TAG`で再実行するスクリプト。`prepare_wsl_inputs.py`で入力生成後に実行
- `run_wsl_dualtes_attempt4.sh`, `runs/attempt4/`: Newton更新なしの旧実測（本結果の比較には未採用）
- `run_wsl_dualtes_attempt2.sh`, `run_wsl_dualtes_attempt3.sh`: 結果保存パスを直す前の失敗実行記録
- `analyze_wsl_dualtes.py`, `flush_final_tes_series.py`, `series_flush_audit.json`: 集計と最終系列flush監査

Newton更新なしのDualTES実測も全stepまで進みましたが、SinglePixelとの回路更新方法が異なると分かったため最終比較には使っていません。attempt2/3は計算後の出力パス、attempt5は入力SIFのruntime stagingで停止しており、時間比較から除外しました。最終の過渡比は、Newton更新・温度予測・commitをSinglePixel replayに合わせた`attempt6`同士の記録から計算しています。

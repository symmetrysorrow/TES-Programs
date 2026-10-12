# WSLでの過去単ピクセル過渡計算の再実行

## 判定

**10分以内に完走しました。** WSLのElmerGPU環境に残っていた過去と同じElmer実行ファイルで、名目75 msの単ピクセル過渡ケースを225/225 stepまで実行し、正常終了しました。

- Wall時間: **560秒（9分20秒）**
- Solver報告時間: CPU 547.38秒、REAL **555.44秒（9分15秒）**
- 停止上限: 600秒。上限より40秒早く完走
- MPI: 6 ranks、各1 thread
- Solver SHA-256: `16b512caf1a49b3d23d4d8e23e636cd4e1df643531e81a65a3f71bd7a1d48c56`。過去manifestのSHAと一致
- 終了時刻: **75.020 ms**。完走行を既存のfinal-series flush処理でTES系列に追加

過去のSolver REAL時間499.08秒（8分19秒）より56.36秒、約11.3%遅くなりました。したがって今回確認できたのは「同じsolverで10分以内に完走する」ことで、8分ちょうどの再現ではありません。

## 入力と再現性の範囲

元の生成SIFを隔離コピーし、case identifierと出力名だけを変えました。時間積分、パルス、回路、材料、収束条件、線形solver設定、MPI分割はこのコピー内で変更していません。メッシュと再スタート結果はハードリンクで隔離実行領域に用意し、元ファイルは上書きしていません。

ただし、過去manifestに記録されたSIF SHA-256 `210a753a...feedbf00` と、現在の生成SIF SHA-256 `e202aaa...734a8e88` は一致しません。過去実行フォルダーにはSIF本体でなくhashだけが残っているため、**過去のSIFをバイト単位で再実行した試験ではありません**。今回のrun SIFは現在の生成SIFからcase identifierだけを置換したものです。

## TES系列の比較

過去系列と今回の系列は225行で、時刻列は全行一致しました。過去系列の最終1行は完走後flush済みで、今回も同じ処理を適用してから比較しています。

| 系列 | 最大相対差 | 発生時刻 |
|---|---:|---:|
| TES温度 | 0.00119% | 20.018 ms |
| TES電流 | 0.0425% | 20.019 ms |
| TES抵抗 | 0.179% | 20.018 ms |
| TES電力 | 0.130% | 20.018 ms |

系列の差は小さく、旧応答に近い出力を再現しました。ただし、過去SIF本体のhash不一致が残るため、この比較だけで入力まで完全一致したとは判定しません。

## 主な証拠

- `run_record.json`: 完走・Wall時間・solver hash
- `solver.log`: 225/225 step、`ALL DONE`、Solver CPU/REAL時間
- `historical_manifest.json`, `historical_solver.log`: 過去実行の入力hashと時間
- `historical_series.csv`, `replay_series.csv`, `series_comparison.json`: 数値比較
- `replay.sif`, `preparation.json`: 今回のSIFとsource hash
- `run_wsl_replay.sh`, `flush_final_series.py`, `series_flush_audit.json`: 再実行・最終系列flush手順
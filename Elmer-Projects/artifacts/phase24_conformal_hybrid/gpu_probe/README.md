# GPU 切り出し実験（2026-09-28、RX 9070 XT / WSL ElmerGPU、ROCm 7.2）

Elmer を通さず、h8 の線形系（整合行列 A、集中化した前処理用行列 Al）を GPU で解いて CPU と比較した。
結論：WSL では GPU は CPU 6 ランクより遅い。原因は CPU–GPU 間の同期遅延。GPU 路線は打ち切り、CPU で続ける。

## ファイル

- `export_csr.py`：`/root/linsys/{A,A_lumped}.npz` を raw CSR（`A.csr`、`Al.csr`）と右辺 `b.bin` に書き出す（行列の作り方は `../summary.md` の集中化前処理の節）。
- `hypre_probe.c`：HYPRE PCG + BoomerAMG（AMG は Al、PCG は A）。CPU 版・HIP 版 HYPRE の両方でビルド。IJ はホストで組み立てて `HYPRE_IJMatrixMigrate` で GPU へ移す。
- `roc_probe.cpp`：rocALUTION CG + AMG。CG の Build が前処理を作り直すため、AMG を Al で作るラッパー前処理（`TildePrec`）を使う。
- `lat.hip`：HIP のカーネル起動・同期・読み戻しの遅延測定。

## 結果（ゼロ初期値から相対 3e-11）

| 方式 | 反復 | 1 反復 | 1 回の解き | AMG 構築 |
|---|---|---|---|---|
| CPU 6 ランク HYPRE（relax 8、積極的粗化 1 段） | 42 | 13.4 ms | 0.56 s | 0.19 s |
| GPU HYPRE（relax 18） | 33 | 23.7 ms | 0.78 s | 0.67 s |
| GPU rocALUTION Ruge-Stüben AMG | 28 | 33.6 ms | 0.95 s | 0.36 s |
| GPU rocALUTION 平滑化集約 AMG | 450–700 | 100–200 ms | 54–96 s | – |
| GPU rocALUTION Pairwise AMG | – | – | – | 異常終了 |

遅延（`lat.hip`）：起動のみ 3.7 µs、起動＋同期 313 µs、8 バイト読み戻し 210 µs、GPU 内コピー 201 GB/s。
8000 点のラプラシアンでも HYPRE GPU は 1 反復 15.8 ms かかり、計算量ではなく同期回数（1 反復あたり数十回）が時間を決めている。ネイティブ Linux なら同期は 10–20 µs 程度で、1 反復 2–3 ms が見込まれる。

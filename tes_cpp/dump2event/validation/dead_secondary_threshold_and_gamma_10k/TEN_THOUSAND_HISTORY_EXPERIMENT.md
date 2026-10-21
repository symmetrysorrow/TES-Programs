# 独立seed 2 MeV gamma / 10,000履歴実験

状態: **2026-10-22 01:34 JSTにユーザー指示で対象PHITS PID 41452を停止。PHITS正常終了していない。Dumpall完全性、HDF5 v2、全加算join、ヒストグラム比較は未検証。**

## 条件・入力保全

- 元入力: `H:\TES-Programs-validation\dead_secondary_threshold_and_gamma_10k\gamma2mev_10000\trace_fresh_10k.inp`
- 元入力SHA-256: `EF76A6FED7302CCAC0E35B1D54B130D9A7FB9E64411420A0FACC7C1DF8D8F815`
- run入力は実行driverがASCIIで複製。SHA-256は元入力と同一。
- 条件: `icntl=0`, `maxcas=10000`, `maxbch=1`, `dumpall=-1`, `negs=1`, `e-mode=2`, 2.0 MeV mono-gamma, Pb region 1000/material 2, EGS5 medium 2, unit weight, 1 keV cutoffs, `[t-deposit]` 0–1.5 MeV/150 bins/unit 3/output=deposit。
- seed: `1100101110010010000011011010111011000101100101010011110001101011`。

## 診断ビルド対照

- O2診断実行ファイルSHA-256: `E55A7E0DE0175A7814644F48C8C4DC1E6FF710CBA67EDEA008DFBFE4BE2C517D`。PHITS配布本体は変更していない。
- 1履歴対照runはexit 0。Dumpall SHA-256は既存対照と同じ `9040FF8FE789B1C6DB40C2B36ED6FF7449B4A19491B1619354183FF5CE842163`、`[t-deposit]` SHA-256は `5E18C6155F0D538AFDB2C2014839F903288E1D46486069DBCFD61647D3ED05A4`。

## 10,000履歴run

- O0診断run `run_full_trace` は `2026-10-21 23:11:42` に開始、診断traceの性能上の問題で `23:25:31` に中断(exit `-1`)。部分出力は保存しているが解析に使用しない。
- H:上のO2 runもI/O律速のため23:34:02に中断(exit `-1`, 部分ファイル保持)。有効runは `C:\work\phits_gamma2mev_10k_diag_o2_20261021` に出力し、2026-10-21 23:34:09に開始。実行コマンド:

```powershell
& 'D:\github\TES-Programs-physdeposit-20261020\tes_cpp\dump2event\validation\weighted_cutoff_positron_photon\scripts\run_phits_diagnostic.ps1' `
  -InputFile 'H:\TES-Programs-validation\dead_secondary_threshold_and_gamma_10k\gamma2mev_10000\trace_fresh_10k.inp' `
  -RunDirectory 'C:\work\phits_gamma2mev_10k_diag_o2_20261021' `
  -MaxCas 10000 `
  -Executable 'H:\TES-Programs-validation\esmin_gate_diagnostic\phits_WinGfort_diag_minimal_10k_O2.exe'
```

- 対象PHITS PID 41452と実行ファイルを照合後にそのPIDだけ停止した。driver記録は `finish=2026-10-22T01:34:18.7698007+09:00`, `exit_code=-1`。正常終了を示すPHITS終端記録はなく、途中の `gamma_2mev_10k_phits.out` が残っている。
- 停止後は19ファイル、論理サイズ合計75,980,375,565 bytes。Dumpall 39,103,713,156 bytes、内部trace 23,248,221,345 bytes、tally sequence 7,546,685,328 bytes、Dumpall sequence 5,058,679,248 bytes。これらは破棄・再利用せず、未完了データとして保管。
- C:空き容量の停止直前/直後 `Get-PSDrive` 観測は144,345,436,160 / 144,296,779,776 bytes。後続観測でWMI値と約1.74 GB差があったため、runの物理割当量は断定しない。
- NCOL=2終端、正常終了履歴数、Dumpall/tally/traceの正常runとしてのハッシュ、HDF5 v2、全加算join、t-depositヒストグラム比較は未確認。中断Dumpallを正常PHITS入力として扱わない。

以後、ユーザーから新たな許可があるまで10,000履歴PHITS計算、HDF5化、巨大trace生成を行わない。1000履歴dead-secondaryの再読込照合結果は `STOPPED_10K_AND_1000_REREAD_20261022.md` に記録。

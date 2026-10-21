# 2026-10-22 run停止と2 MeV/1,000履歴の再照合

## 判定

- 10,000履歴run: **ユーザー指示により停止。正常終了していない。** 途中出力は不完全な診断データであり、PHITS完了run、完全Dumpall、tally結果として扱わない。
- dead電子候補: **PHITS 3.341のスコア対象外と確定。** `NOT_SCORED_BY_PHITS` 監査行に残し、depositゼロやHDF5加算に置き換えない。
- 2 MeV/1,000履歴: 通常CLIが既に作成したHDF5 v2を全件再読込し、同じrunのPHITS内部加算と再照合した。新しいPHITS計算や、重複する1.79 GB HDF5生成は行っていない。

## 10,000履歴runを特定して停止

停止前に `Win32_Process` で対象を照合した。

| 項目 | 実測値 |
|---|---|
| Run ID / 出力ディレクトリ | `phits_gamma2mev_10k_diag_o2_20261021` / `C:\work\phits_gamma2mev_10k_diag_o2_20261021` |
| PHITS PID / 実行ファイル | `41452` / `H:\TES-Programs-validation\esmin_gate_diagnostic\phits_WinGfort_diag_minimal_10k_O2.exe` |
| 起動時刻 / driver記録終了時刻 | 2026-10-21 23:34:09.791 JST / 2026-10-22 01:34:18.769 JST |
| 実行条件 | PHITS 3.341、`maxcas=10000`, `maxbch=1`, 2 MeV mono-gamma |
| 終了記録 | driverの `exit_code.txt` は `-1`。PHITS正常終了記録なし |
| 停止方法 | PIDと実行ファイルの一致を再確認後、対象PIDだけに `Stop-Process -Id 41452 -Force`。停止後PID不在を確認 |

実行driverは次のコマンドで起動されていた。

```powershell
& 'D:\github\TES-Programs-physdeposit-20261020\tes_cpp\dump2event\validation\weighted_cutoff_positron_photon\scripts\run_phits_diagnostic.ps1' `
  -InputFile 'H:\TES-Programs-validation\dead_secondary_threshold_and_gamma_10k\gamma2mev_10000\trace_fresh_10k.inp' `
  -RunDirectory 'C:\work\phits_gamma2mev_10k_diag_o2_20261021' `
  -MaxCas 10000 `
  -Executable 'H:\TES-Programs-validation\esmin_gate_diagnostic\phits_WinGfort_diag_minimal_10k_O2.exe'
```

停止対象の初回確認時点では15ファイル、論理サイズ合計25,927,994,368 bytes。C:空き容量は144,345,436,160 bytesだった。その後、PIDと実行ファイルを再確認して停止処理を行うまでrunが出力を続け、停止後の19ファイルは論理サイズ合計75,980,375,565 bytes（70.762 GiB）になっていた。主要ファイルはDumpall 39,103,713,156 bytes、内部trace 23,248,221,345 bytes、tally sequence 7,546,685,328 bytes、Dumpall sequence 5,058,679,248 bytes。途中データを削除・解析して正常run扱いすることはしていない。

| C:の空き容量観測 | 空きbytes |
|---|---:|
| 停止直前の `Get-PSDrive C` | 144,345,436,160 |
| 停止直後の `Get-PSDrive C` | 144,296,779,776 |
| 後続の監査時 `Get-PSDrive C` | 144,205,062,144 |
| 後続の監査時 `Win32_LogicalDisk` | 145,945,341,952 |

PSDriveとWMIの空き容量値に約1.74 GBの差があるため、上記からrunの物理割当量を推定しない。実ファイルの合計は論理長で記録した。PHITS本体、元Dumpall、元入力、PR #6は変更していない。

## NOCAS 398 dead電子の判定

PHITS 3.341の入力仕様では `esmin` は `parz(162)` に対応し、入力で省略した場合の既定値は `0.001 MeV`。同じPHITS 3.341 GNU実行系統・同じ物理設定の隔離runtime probeで、閾値生成後に実測した値は `esmin=1.0000000000000000E-3 MeV`, `eletb(1)=9.9999999999999995E-7 MeV`。probeの追加出力は乱数を呼ばず、PHITSインストール版には変更を加えていない。

3.341ソース経路は `tallsm1.f` の `emins=min(esmin,1.d-6)`、`eletb(1)=emins` と、`talls08.f` のdead-secondary条件 `QCLUSTS(7)>eletb(1)`（厳密なgreater-than）である。`emin(12/13/14)=0.001 MeV` は輸送cutoffであり、このスコアgateとは別の値。

同じ2 MeV/1,000履歴Dumpallからの候補:

| NOCAS | source line | parent ID / ITYP / NCOL / JCOLL | 子index / ITYP / JCLUSTS(4,5) | Q7 (MeV) | Q8 | region/material | 判定 |
|---:|---:|---|---|---:|---:|---|---|
| 398 | 39,594,769 | 180 / 14 / 13 / 13 | 0 / 12 / -1,-1 | 7.5933446597176157e-7 | 1 | 1000 / 2 | `NOT_SCORED_BY_PHITS` |

`Q7-eletb(1)=-2.4066553402823838e-7 MeV` なのでscore関数へ進まない。位置は `(-0.89993733682112143,-0.0080077241081570032,-3.4727350725342681) cm`。親反応に対する別の内部加算は2件（0および0.0010959999999999859 MeV）あるが、dead子index 0のdepositではない。

研究用C++ ledgerはこの判定を `classify_dead_secondary_score(Q7,eletb1)` の厳密比較で適用する。対象外候補は `not_scored_candidates.csv` に記録し、deposit rowを生成しない。閾値が未解決ならfail-closedする。既定converter経路とPR #6は変更しない。

## HDF5 v2の通常CLI結果を全件再読込

通常の `dump2event.exe <Dumpall> <HDF5> --physics pb-egs5 --hdf5-version 2 --physics-config <JSON> --stopping-table <CSV>` が2026-10-21に作成した同一run HDF5 v2を使用。既存CLI記録はexit 0、elapsed 2,349.8716 s、HDF5 1,790,020,441 bytes。追加で同一Dumpall・内部trace・HDF5をストリーム再読込する照合コマンドを実行した。

```powershell
python 'D:\github\TES-Programs-physdeposit-20261020\tes_cpp\dump2event\validation\dead_secondary_threshold_and_gamma_10k\scripts\compare_hdf5_to_phits_trace_streaming.py' `
  --run-dir 'D:\temp\gamma_pb_cli_physics_20261021\seq_n1000_2mev' `
  --dumpall 'D:\temp\gamma_pb_cli_physics_20261021\seq_n1000_2mev\gamma_2mev_thickPb_dumpall.dat' `
  --tracks-csv 'H:\TES-Programs-validation\dead_secondary_threshold_and_gamma_10k\two_mev_1000\join_inputs\raw_tracks.csv' `
  --internal-csv 'D:\temp\gamma_pb_cli_physics_20261021\seq_n1000_2mev\internal_additions.csv' `
  --hdf5 'H:\TES-Programs-validation\dead_secondary_threshold_and_gamma_10k\two_mev_1000\gamma_2mev_1000.h5' `
  --output-json 'H:\TES-Programs-validation\dead_secondary_threshold_and_gamma_10k\two_mev_1000\reread_join_20261022.json'
```

再読込結果: 5,978,912/5,978,912加算一致、欠落・重複・event/particle/region不一致0。5,978,912座標を比較し不一致0、最大位置差0 cm。最大エネルギー差 `4.680023729664029e-14 MeV`、誤差1e-12 MeV超過0、合計deposit差約 `1e-12 MeV`、最大履歴合計差 `6.756817327868703e-13 MeV`。別途再読込したHDF5 completeness auditは1000 event complete、deposit 5,978,912、zero 4,797、negative/invalid/unresolved 0。新しいHDF5ファイルの重複生成はしなかった。

機械可読の再読込結果:

- `TWO_MEV_1000_REREAD_JOIN_20261022.json`
- `TWO_MEV_1000_REREAD_HDF5_AUDIT_20261022.json`

## 検証・範囲

- MSVC Release build: 成功。
- CTest: 1/1成功。
- `python -m pytest tes_cpp\dump2event\tests -q`: 34 passed, 10 skipped。
- 未検証: 中断した10,000履歴runの正常終端、全履歴数、tally数値、HDF5化、分布比較。停止後に新しいPHITS計算を行っていない。
- 統合判断: dead-secondary gateは今回のPHITS 3.341/Pb/2 MeV/1,000履歴条件に限定した研究用実装として検証。PR #6の統合や既定動作への導入は行わない。

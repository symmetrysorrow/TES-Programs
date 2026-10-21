# PHITS 3.341実行時 `esmin` / `eletb(1)` 監査

状態: **実測完了**（2026-10-21 JST）。元のPHITS配布物・入力・Dumpallは変更していない。

## 根拠

- `C:\PHITS\src\read02.f:24940,24947` は入力名 `esmin` を `parz(162)` に対応づける。
- `read02.f:25204-25228` の `data (parz(i),i=101,200)` では、161–170番の行（25217–25218）の第2値が `1.d-3`。したがって入力に `esmin` がないこの入力の既定値は **0.001 MeV**。
- `read00.f:3065-3067` で `esmin=parz(162)` とし、`mntsc` による上限処理を適用する。対象入力は通常EGS5で `mntsc=0`。
- `tallsm1.f:239` は `emins=min(esmin,1.d-6)`、同ファイル:246は `eletb(1)=emins` を生成する。
- `talls08.f:4929-4930` のdead-secondary分岐は `JCLUSTS(4)<0`, `JCLUSTS(5)!=0`, `QCLUSTS(7)>eletb(1)` を要求する。

## 実行時測定

隔離コピーしたシリアルGNU Fortranソースで `tallsm1.f` の閾値生成直後に一度だけ `write(6,*) 'ESMIN_ELETB1_RUNTIME',esmin,eletb(1)` を加え、PHITS 3.341を再ビルドした。診断行は乱数関数を呼ばず、状態を変更しない。

- 入力: `D:\temp\gamma_pb_cli_physics_20261021\seq_n1000_2mev\trace.inp` のコピー、`maxcas=1`、同一seed・同一物理設定。
- 診断実行終了コード: 0。
- stdout実測: `ESMIN_ELETB1_RUNTIME 1.0000000000000000E-003 9.9999999999999995E-007`。
- `esmin`: **0.001 MeV**。
- `eletb(1)`: **9.9999999999999995e-7 MeV**（二進浮動小数点表示）。

同seedの診断なし対照（`phits_WinGfort_lossless_seq.exe`）と診断版で、1履歴のDumpall SHA-256が両方 `9040FF8FE789B1C6DB40C2B36ED6FF7449B4A19491B1619354183FF5CE842163`、`[t-deposit]` SHA-256が両方 `5E18C6155F0D538AFDB2C2014839F903288E1D46486069DBCFD61647D3ED05A4`。比較したファイルは完全一致し、この一回の診断フックによる輸送・タリー差は観測されなかった。

## 影響

設定解決器のPHITS 3.341既定 `esmin` は1.0から0.001 MeVへ修正した。`eletb(1)` は常にPHITS実測/入力に基づき `min(esmin,1e-6)` としてprovenanceへ含める。対象外二次粒子はscore ledgerへ偽のゼロdepositとして出さず、監査候補として別集計する。厳密な内部加算との完全照合は別レポートに記録する。

## 実行環境

- PHITS版: 3.341。
- 診断ビルド: `H:\TES-Programs-validation\esmin_gate_diagnostic\phits_WinGfort_esmin_diag.exe`。
- 診断ビルドの基礎はシリアル `WinGfort` lossless-sequence版。隔離ソースへの1行追加のみ。
- 実測ログと生ファイルはローカル `H:\TES-Programs-validation\esmin_gate_diagnostic\runtime_probe\` にあり、Gitには含めない。

# COMSOLフル区間(パルス+180 ms)比較: 現状

更新日: 2026-07-28

## 目的

既存の `case_tes_mpi_comsol_grid`(パルス+600 µsで打ち切り)を、COMSOLの実計算終端である
パルス+180 ms(絶対時刻200 ms)まで拡張し、MPI(4-rank MUMPS)がCOMSOLの立ち上がり・
ピーク・減衰全体と一致するかを検証した。

## 生成スクリプト

- `scripts/prep/prepare_comsol_timegrid_case.py`
  - 実行: `python scripts/prep/prepare_comsol_timegrid_case.py`
  - 出力: `elmer_project_comsol_timegrid.json`
  - このスクリプトが以下のケース定義をすべて生成する。

## 採択ケース(最終)

**`case_tes_mpi_comsol_grid_full_uniform_continuous`**

- 定常状態(`case_tes_steady_3x_refined`)からの単一連続実行、リスタート区切りなし
- t=0〜パルス+500 µs: COMSOLの実サンプル時刻に1対1で追従(既存`case_tes_mpi_comsol_grid`と同じ細かさ)
- パルス+500 µs〜180 ms: 100 µs一定刻み(1795ステップ)
- 合計205ステージ・1999物理ステップ、絶対時刻0〜200.02 ms
- 実行時間: 約5時間(1回目はMUMPSのメモリ不足`INFO(1)=-9`で失敗、再試行で完走)
- 結果: `results/case_tes_mpi_comsol_grid_full_uniform_continuous/`
  - `tes_mpi_comsol_grid_full_uniform_continuous_series.csv`(TES系列、1999行)
  - `tes_mpi_comsol_grid_full_uniform_continuous_iterations.csv`(非線形反復ログ)
  - `solver.log`, `manifest.json`

実行コマンド:

```powershell
python run.py case_tes_mpi_comsol_grid_full_uniform_continuous --project elmer_project_comsol_timegrid.json --mpi-procs 4 `
  --elmer-solver "D:\Github\TES-Programs\tools\elmer-hypre\install-phase13-step-commit\bin\ElmerSolver.exe" `
  --runtime-bin "C:\msys64\ucrt64\bin"
```

## 実行環境の注意点

1. **`--elmer-solver`は必須。** `run.py`のデフォルト(`C:\Program Files\Elmer 26.1-Release\bin\ElmerSolver.exe`)は
   `libmpi_stubs.dll`にリンクされたシリアル専用ビルドで、`--mpi-procs`を指定してもMUMPSが使えない
   (`ERROR:: CheckLinearSolverOptions: MUMPS solver has not been installed.`)。
   本物のMPI(`msmpi.dll`)にリンクされた自前ビルド`tools/elmer-hypre/install-phase13-step-commit`を
   `--elmer-solver`と`--runtime-bin C:\msys64\ucrt64\bin`で明示指定すること。
2. **空きメモリを4〜5 GB以上確保すること。** 不足すると`OpenBLAS error: Memory allocation still failed`や
   `DMUMPS INFO(1)=-9`(作業領域不足)でMPIジョブが落ちる。落ちた場合は他アプリを閉じて再実行すれば通ることが多い。
3. MS-MPI自体(`mpiexec`/SMPD)は正常なら`ParCommInit: Initialize #PEs: 4`とログに出る。
   `#PEs: 1`のまま4プロセス立ち上がる場合はリンク先MPIライブラリの問題であり、
   `--elmer-solver`の指定漏れ(1.)を疑うこと。

## 判明した2つの数値的アーティファクトとその対処

### (a) BDF1リンギング(パルス+10〜60 ms付近)

パルス+500 µs以降の刻み幅が数百µs〜1.8 ms級になると、COMSOLの実サンプル刻みでも
成長率を1.2倍以内に抑えた滑らかなランプでも、振幅0.3 µA程度(ピーク比約4%)のリンギングが
残ることを確認した(`case_tes_mpi_comsol_grid_full`、`case_tes_mpi_comsol_grid_full_smooth`で検証)。
原因は刻み幅の**急変**ではなく**絶対的な大きさ**で、パルス+500 µsチェックポイントから
50 µs/100 µs一定刻みで再計算する診断ケース(`case_tes_mpi_comsol_grid_ringing_diag[_100us]`、
パルス+500 µs〜25 msのみ)でリンギングがほぼ消える(最大誤差0.3→0.04 µA)ことを確認し、
100 µsを採択した。

### (b) チェックポイント再起動の継ぎ目

パルス+500 µsチェックポイントから`restart_from`で再起動する構成
(`case_tes_mpi_comsol_grid_full_uniform_100us`等)では、継ぎ目(パルス+580〜700 µs付近)に
小さな段差が出た。TES回路UDF(`tes_parallel_circuit`)の内部状態(前ステップ電流の確定タイミング、
Aitken緩和履歴など、FortranモジュールのSAVE変数)がプロセス再起動をまたいで引き継がれないためと
考えられる。単一連続プロセス(`case_tes_mpi_comsol_grid_full_uniform_continuous`)にすることで解消した。
チェックポイント(`Output File`)自体は途中に残しているが、実際にresumeケースを作らない限り
この継ぎ目は発生しない。

## 最終比較結果

`case_tes_mpi_comsol_grid_full_uniform_continuous`とCOMSOL(`docs/Single-Pixel.txt`)の比較:

| 指標 | COMSOL | Elmer MPI | 誤差 |
|---|---:|---:|---:|
| ピーク電流降下 | 7.7748 µA | 7.8330 µA | +0.75% |
| 立ち上がり時間(10-90%) | 161.38 µs | 147.67 µs | -8.50% |
| 立下り時間(90-10%) | 25955.52 µs | 26414.65 µs | +1.77% |

ピーク・立下りは1%台、立ち上がりは-8.5%程度の差が残る(グリッド解像度によらない、
モデル側の物理的な差と考えられる)。

## 比較・プロット生成

`scripts/analysis/plot_comsol_direct_mpi_current.py`(このタスクで拡張):

- `--skip-direct`: 直接法(serial)曲線を省略
- `--linear-x`: 180 ms全体を見る際に既定のsymlogを線形軸に切替
- `--tag <name>`: 同じ`--out`フォルダに複数バリアント(例: `_logx`)を上書きせず並存させる
- Difference(誤差)パネルは廃止し、代わりにpeak/rise/fall時間とCOMSOL比誤差率を
  `pulse_metrics.csv`と標準出力に表形式で出力する

実行例:

```powershell
python scripts/analysis/plot_comsol_direct_mpi_current.py `
  --mpi results/case_tes_mpi_comsol_grid_full_uniform_continuous/tes_mpi_comsol_grid_full_uniform_continuous_series.csv `
  --start-us -20000 --end-us 180000 --linear-x --skip-direct `
  --out artifacts/comparison/comsol_mpi_full_uniform_continuous_final
```

最終成果物: `artifacts/comparison/comsol_mpi_full_uniform_continuous_final/`
- `current_timeseries_comparison.png` / `.svg`(線形軸)
- `current_timeseries_comparison_logx.png` / `.svg`(symlog軸)
- `pulse_metrics.csv` / `pulse_metrics_logx.csv`
- `baseline_corrected_current_smooth.csv` / `_logx.csv`

## 途中生成した非採択ケース(参考・履歴)

生成はされているが最終結果には使わない中間ケース。いずれも
`scripts/prep/prepare_comsol_timegrid_case.py`内に定義が残っている。

| ケース名 | 内容 | 非採択の理由 |
|---|---|---|
| `case_tes_mpi_comsol_grid_coarse_tail` | 0〜500 µs細かい+500-600 µsを5倍間引き、500 µsにチェックポイント | 600 µsで打ち切り、180 msまで届いていない(が、このチェックポイントは他ケースの再利用元として現存) |
| `case_tes_mpi_comsol_grid_coarse_tail_resume` | 上記チェックポイントから600 µsまでの続き | 同上 |
| `case_tes_mpi_comsol_grid_full_tail` | 500 µsチェックポイントから10倍間引きで180 msまで延長 | 継ぎ目にリンギング由来の不自然なこぶ(パルス+100〜150 ms付近) |
| `case_tes_mpi_comsol_grid_full` | 単一連続、COMSOL実サンプルを180 msまでそのまま追従(1262ステップ) | パルス+10〜60 msにBDF1リンギング(振幅0.3 µA) |
| `case_tes_mpi_comsol_grid_full_smooth` | 上記のパルス+9.98-19.98 ms区間だけ滑らかな成長率ランプに置換 | 境界の8倍ジャンプは直ったが、リンギング自体(絶対刻み幅由来)は消えず |
| `case_tes_mpi_comsol_grid_ringing_diag` / `_100us` | 500 µsチェックポイントから50/100 µs一定刻みで25 msまでの短い診断 | 診断専用。100 µsで十分と判断し本番へ反映 |
| `case_tes_mpi_comsol_grid_full_uniform_100us` / `_50us` | 500 µsチェックポイントから100/50 µs一定刻みで180 msまで | リンギングは解消したが、チェックポイント継ぎ目の段差が残る |

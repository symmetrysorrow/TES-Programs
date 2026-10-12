# Windows Release single-pixel 75 ms: nonlinear 1e-6 / preconditioner reuse ON

## 結果

600秒で停止し、75 ms計算は未完走です。225 step中55 stepまで進み、最後に記録された時刻は20.121595 msでした。停止直前のsolver表示は推定残り30.7分でした。したがって完走時間や75 msの性能比は確定できません。

- Wall: 600.06 s（10分上限）
- MPI: 6 rank × 1 thread
- Windows normal Release solver SHA-256: `a4bc7af575c765857b544d97439689320ed5ce84ea40dad39ef776a9d326003a`
- HeatSolve.dll SHA-256: `d817c48ea150fd92137a9fa1221219ecad42ee819561716873593f1fa92570ea`
- 非線形許容値: `1e-6`
- Phase24 preconditioner matrix reuse: SIFでTrueを指定
- MATC node cache: ON、Vector Assembly: OFF
- Linear solver条件は前回Windows Releaseのまま（HYPRE PCG/BoomerAMG、max 2000、tol 1e-11）
- TES系列: 54行。最終時刻20.121595 ms、TES温度0.1694585 K、電流114.17 uA、抵抗20.553 mΩ、電力0.2935 nW

## 再利用フラグについて

preconditioner reuseはHYPREの前処理データの再利用で、MATC node cacheとは別の処理です。両者が直接競合する設定ではありません。OFFは各線形solveで前処理を更新し、変化した行列に対する頑健性を優先する選択です。ONはsetup時間を減らせる一方、非線形反復で行列が変わるとKrylov反復が増える可能性があります。線形solverが指定残差まで収束すれば線形解自体は許容誤差内ですが、収束性と実時間は変わり得ます。

今回のSIFではONを要求しましたが、ログには `Phase24 precond reuse:` の実行時診断が0件で、HYPRE setup行が637件ありました。したがって、このビルド・実行で再利用が実際に有効になったとは確認できません。OFFとの速度差をこのrunから結論づけません。MATC cacheとの競合を示すログもありません。

今回のrunは10分時点で停止し、`ALL DONE`はありません。ログとTES/反復系列は保存済みです。10分後にsolver/MPIプロセスが残っていないことを確認しました。元の75 ms SIF、restart、DLLは変更していません。

## ファイル

- `run_windows_prior_flags_full.py`: 10分timeout付き再実行スクリプト
- `runs/windows-prior-flags-full/singlepixel/case.sif`: 実行用隔離SIF
- `runs/windows-prior-flags-full/singlepixel/run_record.json`: hashesと実行条件
- `runs/windows-prior-flags-full/singlepixel/solver.log`: MPI solver log
- `runs/windows-prior-flags-full/singlepixel/singlepixel75_full_nl1e6_precondreuseon_cacheon_1_series.csv`: TES系列
- `runs/windows-prior-flags-full/singlepixel/singlepixel75_full_nl1e6_precondreuseon_cacheon_1_iterations.csv`: nonlinear iteration履歴

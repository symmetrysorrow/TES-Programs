# DualTES過渡解析の短時間プロファイル調査

調査日: 2026-10-11

対象: Elmer 26.1、MPI 6 ranks、短時間4-step診断と既存75 ms完走ログ

## 結論

現状の測定からは「DualTESの要素数が少ないのにassemblyが遅い」原因を特定できない。外側のAssembly値は `CPUTime()`、Phase24の内訳は `RealTime()` で計測されており、差分を未計測の実時間として扱えない。さらに同じ実時間時計で記録されたbulkイベント合計はDualTES 37.204 s、single-pixel 83.340 sで、DualTES側が短い。assembly全体の差を説明する内訳はまだない。

安全に既定値へ入れられる高速化は今回確認できなかった。Phase24 vector assemblyを有効化した短時間試験は、DualTESとsingle-pixelの両方で最初のAssembly開始直後にSIGSEGVとなった。したがってvector assemblyおよび、それを前提とするstatic matrix reuseは採用しない。75 msの再実行はしていない。

次に必要なのは、現在のassembly区間内で `DefaultInitialize()`、element/bulkループ、`DefaultDirichletBCs()`、任意のpreconditioner生成、後処理を同じ `RealTime()` 時計で測ること。rankごとの値とrank最大値を集めるが、測定のためのMPI barrierは挿入しない。HYPREについては既存の診断出力から初回の約12秒の遅延は再現せず、報告されたsetup時間では説明できない。

## 基準ログ

| 実行 | steps / nonlinear assembly | 経過時間 | Assembly CPUTime合計 | Solve CPUTime合計 | 外側Phase24 solve RealTime合計 | HYPRE solution time合計 |
|---|---:|---:|---:|---:|---:|---:|
| DualTES short | 4 / 6 | 121.81 s | 103.80 s | 14.29 s | 14.279 s | 0.916 s |
| single-pixel short | 4 / 6 | 90.26 s | 83.58 s | 2.38 s | 2.381 s | 0.299 s |
| DualTES既存full pulse | 225 saved steps / 507 calls | 8653.66 s | 8279.85 s | 369.67 s | — | 238.948 s |

full pulseは既存ログを解析した。再実行していない。短い2ケースは同じ時間刻み列と6 nonlinear callsだが、meshとrestart fieldが異なる診断比較であり、完全なscaling試験ではない。

### Assembly時計の読み方

DualTESのPhase24詳細RealTimeイベント合計は element traversal 14.836 s、material-property evaluation 2.108 s、local element calculation 18.709 s、global sparse insertion 1.308 s、boundary assembly 0.243 s、合計37.204 s。single-pixelの同イベント合計は83.340 s。ユーザー提供の初期メモにあるbulk時間推定とは一致しないため、本レポートはsolver log上の全6 callの生値を採用した。

HeatSolveの外側Assembly計測は `CPUTime()` で `DefaultInitialize()` より前に開始する。Phase24 bulkのRealTime計測は `DefaultInitialize()` より後に始まる。外側区間にはbulkのほか、`DefaultDirichletBCs()`、任意のlumped-preconditioner生成、operator bookkeeping等が含まれる。CPU時計の103.80 sからRealTime合計37.204 sを引いた66.596 sは、有効な実時間内訳ではない。原因の確定には各区間を同一時計で再計測する必要がある。

## HYPREの初回遅延

DualTES shortの外側solve RealTimeは6回合計14.279 s。そのうち2回目が12.119 sだが、対応する `SolveHypre: Solution time` は0.051 sだった。このoutlierは同一条件の再実行で再現しなかった。再実行の外側solveは `[1.045, 0.319, 0.288, 0.277, 0.242, 0.227]` s、合計2.398 s。inner solution timeは `[0.725, 0.052, 0.029, 0.028, 0.045, 0.028]` s、合計0.907 s。

`PHASE24_HYPRE_REPORT=1` は既存の診断出力のみを有効にする。6 rank分のsetup報告は各callで平均約0.137 s（個々の報告0.094–0.171 s）、setup timerを含むentry以来時間の平均は約0.196 s（0.144–0.269 s）。この計測値は12.119 sの遅延を説明しない。MPI待機、OSスケジューリング、wrapper周辺などは調査候補だが、ログだけで原因とは断定できない。

診断環境変数ありのrepeatは6 nonlinear iterationを完走し、6個すべての `Result Norm` が保存baselineと一致した。series CSVの確認行も一致した（例: t=0.020018、TES温度0.16868560274410294、電流1.3954045342313423e-4）。この範囲では診断出力による数値変更は見られない。実行時間自体は環境変動があり、repeatは107.44 s、baselineは121.81 sだった。

full pulseの外側Solve CPUTime合計369.67 sとinner HYPRE solution合計238.948 sの差も、それぞれ異なるtimerおよび範囲を含むため、特定の未計測処理に割り当てられない。

## vector assembly / reuse試験

既存restart状態を相対パスで指定し、元の状態ファイルSHA256を維持してDualTESとsingle-pixelをそれぞれ試験した。両方ともmesh上のbulk要素は100% fast pathに入り、generic要素数は0だった。DualTESは164,693 bulk elements（tetra 131,633、prism 33,060）、single-pixelは260,136（tetra 236,916、prism 23,220）。それでも両ケースとも “Starting Assembly” の直後に終了コード `4294967293`（Windows表示上のSIGSEGV相当）で停止し、solveや数値結果に到達しなかった。

過去のPhase24 optimized SIFには `Phase24 Vector Assembly = Logical True` が設定されているが、現在の短時間baseline SIFでは指定されず、現行sourceのdefaultはfalse。高速化例の設定をそのまま使う根拠にはならない。`Phase24 Static Matrix Reuse=True` は今回実行していない。vector pathが最初のAssemblyでクラッシュするため、reuseだけを組み合わせて試す前提がない。原因究明にはdebug buildとnative stack traceが必要。

## source / binaryの対応

今回の実行manifestのruntime hashは ElmerSolver_mpi.exe `64922c6d98390e7edebb16db87df8118b46fda4b68ed9a483eb7b8e556a26448`、loaded HeatSolve.dll `be54d5d6087c38ce14bb02241e292ac9e3f2baf2f846a70543fbc39ca2b4002f`。無視対象の外部build tree内にある出力artifactも同じDLL hashで、実行元runtime artifactと一致する。確認したHeatSolve.F90 source hashは `965fa34d9d051467d821e87c62806c8cea866e1478e3fc83f6bc5d20d2fe3f50`、SolveHypre.cは `5ef6ee098e34128d0e433cedde1b744e02121e0cc63b5cecf91f53eefca7677f`。したがってsource、build output、loaded runtimeのartifact同一性は以前のhandoffより強く確認できた。ただしこれは完全なcompiler/link provenance（build commandや全依存関係の再現証明）ではない。

歴史的な `install-stage11` の最適化HeatSolve DLLは `851efd...` であり、現在ロードされたDLLとは別物。過去の最適化ビルドの効果を今回のruntimeに帰属させない。

## 推奨する次の作業

1. isolated diagnostic buildに、`DefaultInitialize()`、bulk loop、`DefaultDirichletBCs()`、preconditioner生成、bookkeeping各区間の同じRealTime計測を追加する。
2. 各MPI rankの時間を記録し、rank平均・最小・最大を集約する。新しいMPI barrierは追加しない。
3. まず4-step短時間ケースでDualTESとsingle-pixelを再測定し、assembly全体も同じRealTimeで囲んで区間合計と照合する。
4. vector assemblyのSIGSEGVをdebug build / stack traceで直してから、その短いケースで数値norm・series・iteration数をbaselineと比較する。数値一致と速度改善が確認できた場合だけ設定変更を検討する。
5. 75 msの本計算は再実行しない。短時間計測と数値ゲートを通過した最適化ができてから必要性を再判断する。

## 再現用artifact

同じディレクトリの `run_short_ab.py` が短時間の診断runner、`manifest.json` がruntime・入力状態・生成物のhashと実行statusを記録する。`dualtes_hypre_report.*` は診断出力repeat、`dualtes_vector.*` と `singlepixel_vector.*` はfast pathのcrash試験である。solver log、生成SIF、series、iteration CSV、state snapshotも同梱する。

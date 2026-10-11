# Phase24 MATC Node Cache: Release検証と統合判断

## 結果の要約

同一の隔離Release build、MPI 6 ranks、Vector Assembly OFFで4-step DualTESとsingle-pixelを比較した。MATC node cacheは明示的にONにした試験で機能し、温度場・TES出力・収束履歴は機械精度内で一致した。Release全timestepの外側`DefaultSolve`直前・直後で保存した線形系も、全6 nonlinear solveについて行列係数とsparsity metadataがbyte一致した。RHSはstep 4の最終桁だけに差が残った。

DualTESの3回測定中央値はOFF 71.004秒、ON 25.615秒で2.77倍。single-pixelは各1回の計測でOFF 53.666秒、ON 28.700秒で1.87倍だった。Single-pixelは反復試験が1回ずつなので、単独性能値のばらつきは評価していない。

キャッシュは既定値Falseのままにした。`TES-Programs`にはElmerの`HeatSolve.F90`製品ソースが含まれず、実装差分は隔離用sourceに適用した既存patchとして保存されている。この結果は統合候補の検証であり、本番DLL・設定への採用ではない。

## Release buildと実験条件

- CMake/Ninja、GNU Fortran/C/C++ 16.1、Fortran `-O3`、C/C++ `-O3 -DNDEBUG`、OpenMP、MS-MPI、HYPRE。
- `HeatSolve.dll`と`ElmerSolver_mpi.exe`はTemp配下の隔離Release buildから実行。ビルドは1771/1771 targetで成功。
- MPI 6 ranks、`OMP_NUM_THREADS=1`、`OMP_DYNAMIC=FALSE`、Vector Assembly OFF。
- 各試験前にTES stateを同一の初期stateへ戻した。DualTESではOFF/ONの各runで2つのstateをコピーし直した。Single-pixelも同じstateファイルを各run前に戻した。
- DualTESはOFF/ON warm-up後、OFF1→ON1→ON2→OFF2→OFF3→ON3の順で実行。Wall計測はsolver開始から終了まで。
- Single-pixelは4-stepの4 timestepを実行。75 ms本計算は行っていない。

## 実測性能

| ケース／指標 | Cache OFF | Cache ON | 比較 |
|---|---:|---:|---:|
| DualTES 全体Wall中央値（3回） | 71.004 s | 25.615 s | 2.77倍、64.0%短縮 |
| DualTES 全体Wall平均 ± 標準偏差 | 70.921 ± 0.766 s | 25.249 ± 0.646 s | 同一Release build |
| DualTES rank 5・Material 7 `material_eval` 合計中央値 | 52.258 s | 0.149 s | 要素内評価の計測値 |
| Single-pixel 全体Wall（各1回） | 53.666 s | 28.700 s | 1.87倍、46.5%短縮 |

DualTES rank 5のMaterial 7評価合計はOFFで51.382–52.307秒、ONで0.135–0.171秒だった。ONではrank 5に22,024 unique nodeがあり、各nonlinear assemblyでmaterial単位の`ListGetReal`呼出しは1回となる。検証用Release runではcache evaluation自体の計測値がrank 5で約0.88–0.96秒／assemblyだった。この時間は要素ループ内の`material_eval`欄に含まれないため、要素内値だけではなく全体Wallも併記した。

Single-pixelのRelease verification runではMaterial 7のrank-local unique-node数はrank 0–5順に10,407、6,969、8,441、10,079、10,133、11,849。各assemblyでmaterial単位1回の`ListGetReal`を使い、前assemblyで再評価を避けたnode数はrank 5で18,852。cache構築時間はrankごと0.029–0.042秒だった。

## 数値検証

### DualTES

- Release温度場: 6 rank × 4 timestepの全nodeを比較。最大絶対差 `5.55e-17 K`、最大相対差 `3.29e-16`。
- TES 1/2のseries CSVとiteration CSVは完全一致。
- `HeatSolve: Result Norm` 6件は完全一致。
- HYPRE PCGのiteration/residualログ36件はrank出力を含めて完全一致。
- 全timestepの線形系: step 1/4で76ファイル、step 2/3で38ファイル、step 4で76ファイルをそれぞれ比較。matrix係数ファイルとsparsity/size metadataは全ファイルbyte一致。RHSはstep 1–3で完全一致。step 4は1,101,672値中238値に最終桁差があり、最大絶対差 `1.02e-20`、最大相対差 `6.37e-16`。行対応は全て一致。
- Release cache verification runのMaterial 7 sampled node照合は8件、mismatch 0。

### Single-pixel

- Release温度場: 6 rank × 4 timestepを全nodeで比較。最大絶対差 `8.33e-17 K`、最大相対差 `4.94e-16`。
- TES series 3行、iteration CSV 6行は完全一致。
- `HeatSolve: Result Norm` 6件とHYPRE PCGのiteration/residualログ36件は完全一致。
- cache verificationは6 rankでcacheを構築し、旧`ListGetRealArray`経路とのsample照合12件、mismatch 0。

Result shardはElmerが開始時刻を記録するためSHA256 byte一致の指標にはせず、出力値をrank・timestep・nodeで読み出して比較した。Release結果のshard hash、差分統計、CSV hashはJSON manifestに保存している。

## MATC cacheの処理範囲と安全性

既存patchの適用条件は、熱伝導率がscalar MATC string、依存変数がTemperatureのみ、Lua/procedureなし、Heat Conductivity定義が単一であること。material IDとrank-local node indexで値を分離し、現在のTemperatureから各nonlinear assemblyで再評価する。条件を満たさないmaterialは既存Generic経路に戻る。行列/RHS組立て、CSR挿入、時間積分、TES circuit演算は変更していない。

実測ではDualTES Material 7とsingle-pixel Material 7を検証した。Temperature以外の依存式、異方性、Lua/procedure、複数conductivity entryを使う負例テストは今回実施していない。判定条件はMATC式の依存メタデータに依存するため、未検証表現形式を有効対象に拡げないこと。

キャッシュは`SAVE`状態を持つ。今回のrun間ではsolverを別プロセス起動し、同一mesh/stateから再現したが、同一プロセス中に別Modelへ切替えるケースの寿命・破棄は検証していない。MPI再partitionやmesh再初期化後に別Modelの状態を再利用しないことを型・Model identityで保証する実装も確認できていない。このため既定OFFを維持し、汎用統合前にcache cleanupとModel identity guardを追加検討する。

コード上では既存`Phase24EvaluateDynK()`とGeneric cacheがそれぞれ動的熱伝導率を評価する。前者は既存のPhase24 prism/assembly pathに束縛され、Generic prism routeからその配列を安全に共有できるか証明できなかったため、patchはrank-local Generic mapを別に構築する。保守統合時には評価ルーチン共通化を検討するが、配列shape・material/node対応の数値検証を先に維持する。

## 以前の75 ms single-pixel最適化との関係

75 msは物理時間幅を表す過去の計算ケースであり、今回実行していない。過去のcommit `a0dc3dfbc9fe2d3be996d3abc12e0e2ae0ee8139`の主な高速化は、scalar P1 tetraに適用されるVector Assembly経路、CSR destination再利用、静的matrix reuse等である。hybrid prism/tetraのGeneric経路と温度依存MATCは適用外または従来経路に残る。今回のcacheはVector Assembly OFFのGeneric prism経路で、shared nodeごとに重複するMaterial 7のMATC評価を削減する。

| 機構 | 主に削減する処理 | 対象と既知の制限 |
|---|---|---|
| Vector Assembly / CSR reuse | P1 scalar tetraの組立て・CSR位置探索、静的matrix再利用 | 現在のVector Assembly SIGSEGVは別問題として未解決。今回ONにしていない |
| Generic MATC node cache | Material/node単位のTemperature-only scalar MATC重複評価 | Generic prism経路。温度依存外のMATC等はfallback。既定OFF |

2機構は理論上別経路を短縮し、将来の併用余地があるが、同時有効化は今回検証していない。Vector Assembly SIGSEGVは未解決のまま。

## 文字化け修正と成果物

以前の3つのMarkdownに含まれた`?`はUTF-8デコード時点で置換済みだったため、原文そのものの復元ではなく、既存manifest・JSON・ログ・patchに基づいて日本語を再構成した。対象ファイルの文字化けを検査し、UTF-8で保存した。

- `report_ja.md`: Release測定・数値結果・安全性・未解決事項
- `release_runs.json`, `singlepixel_release_runs.json`: solver wall、rank-local profile、状態とSIF hash
- `release_all_assembly_capture.json`, `release_all_assembly_comparison.json`: 全nonlinear solveの行列/RHS捕捉と比較値
- `release_dualtes_field_comparison.json`, `singlepixel_temperature_comparison.json`: Release温度場全node比較
- `dualtes_numerical_logs.json`, `singlepixel_numerical_logs.json`: TES、Result Norm、HYPRE記録の比較
- `run_*.py`, `*_6rank.sif`, `.solver.log`, `*.csv`: 再実行入力・生ログ・出力
- `HeatSolve_generic_matc_node_cache.patch`: 隔離Elmer sourceに適用した候補差分
- `artifact_inventory.json`: 保存済み成果物のバイト数とSHA-256

数値比較のため一時行列/結果ファイルはTemp内に置き、全capture shardとresult shardのhashをmanifestに保存した。既存の本番DLL・本番SIF・物理材料値・時間刻みは変更していない。

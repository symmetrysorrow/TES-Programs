# 4-step追加切り分け：MPI待機とVector Assembly

調査日: 2026-10-11
対象コミット: `6294c3893462e880dcdd43c3f2089165b799abf0`
実験範囲: 4ステップ、6 MPI ranks。Vector OFFのDualTES / single-pixel。診断ビルドは GNU Fortran 16.1.0、GNU OpenMP 5.2、MS-MPI、Debug `-O0 -g3 -fbacktrace -fno-omit-frame-pointer`、`OMP_NUM_THREADS=1`。75 ms解析は実行していない。

## 1. MPI待機とローカル計算

`DefaultDirichletBCs()`の内部を同じ `RealTime()` 時計で分割し、各rankのローカル境界条件数の集計、整数 `ParallelReduction`、Dirichlet値の隣接rank通信、行列処理を記録した。計測用MPI barrierは追加していない。DualTESとsingle-pixelは同じ診断ビルドと4-step設定で実行し、各々36 rank/assemblyサンプルを得た。両方とも正常終了した。

| ケース / rank | bulk (s) | Dirichlet (s) | うち最初の整数 reduction (s) | Dirichlet値通信 `CommunicateDirichletBCs` (s) | assembly合計 (s) |
|---|---:|---:|---:|---:|---:|
| DualTES 0 | 4.851 | 9.678 | 9.662 | 0.0002 | 14.566 |
| DualTES 1 | 14.168 | 0.389 | 0.380 | <0.001 | 14.565 |
| DualTES 2 | 13.519 | 1.038 | 1.029 | <0.001 | 14.565 |
| DualTES 3 | 4.862 | 9.667 | 9.651 | 0.0002 | 14.567 |
| DualTES 4 | 12.866 | 1.692 | 1.682 | <0.001 | 14.565 |
| DualTES 5 | 14.547 | 0.010 | 0.000 | <0.001 | 14.565 |
| single-pixel 0 | 12.145 | 0.013 | 0.000 | <0.001 | 12.167 |
| single-pixel 1 | 7.119 | 4.975 | 4.947 | <0.001 | 12.171 |
| single-pixel 2 | 7.915 | 4.185 | 4.159 | <0.001 | 12.171 |
| single-pixel 3 | 9.246 | 2.910 | 2.899 | <0.001 | 12.168 |
| single-pixel 4 | 9.318 | 2.841 | 2.830 | <0.001 | 12.168 |
| single-pixel 5 | 10.111 | 2.048 | 2.036 | <0.001 | 12.167 |

表は36呼び出しのrank別平均。DualTESはrank 0/3、single-pixelはrank 1–5のDirichlet時間のほぼ全てが最初の `ParallelReduction(n)` にあります。局所カウント処理は計測精度以下で、Dirichlet値通信のsetup/pack/send/recv-wait/unpackも各rankで概ね0–0.001秒でした。従って、観測されたDirichlet時間は境界値データ転送やローカルDirichlet行列処理では説明できず、collective内での経過時間です。

この計測だけではcollectiveのネットワーク処理時間と、rank到着差による待ち時間をさらに分解できません。よって「9.66秒全てが待ち」とは断定しません。ただし、Dirichletが短いrankほどbulkが長いという逆相関と、assembly合計がrank間で揃うことから、早くbulkを終えたrankが遅いrankを待つ到着差が主因という解釈を強く支持します。collective内部の通信時間と待ち時間を厳密に割るには、barrierを入れずPMPI/MPI実装の計測を追加する必要があります。

`HeatSolve.F90`のbulk要素ループと `DefaultFinishBulkAssembly()` の該当経路にはMPI collective / point-to-point callがありません。したがってrank間のbulk差はMPI通信待ちより、ローカルな要素処理量・要素種別・材料演算の差です。partitionファイルの内訳は次の通りです（rank 0 = `part.1`として対応付け）。

### DualTES partition

| rank | local nodes / equation rows | bulk elems | tetra | prism | boundary elems | target BC tag 56 elems / nodes | shared nodes |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 45,632 | 164,693 | 131,633 | 33,060 | 38,422 | 4,146 / 2,177 | 1,993 |
| 1 | 46,395 | 131,064 | 70,716 | 60,348 | 42,677 | 0 / 0 | 3,004 |
| 2 | 46,234 | 141,770 | 87,444 | 54,326 | 40,710 | 0 / 0 | 2,718 |
| 3 | 45,848 | 166,517 | 133,433 | 33,084 | 38,428 | 4,150 / 2,179 | 2,131 |
| 4 | 45,581 | 148,150 | 99,804 | 48,346 | 38,806 | 0 / 0 | 2,900 |
| 5 | 45,728 | 116,541 | 50,513 | 66,028 | 43,815 | 0 / 0 | 3,214 |

rank 0/3はtetra主体で総要素数が多くてもbulkが約4.85秒、rank 5は要素総数が最少でもprism数が最多でbulkが約14.55秒でした。この2点はDualTESの遅いbulkが単純な要素数負荷ではなく、prism側の要素当たり演算コストの高さを含むことを示唆します。single-pixelでもrank 1/2は約99k prismのみで7.1–7.9秒、rank 0は約237k tetraと23k prismで12.1秒でした。要素種別に加え、材料・係数評価を含めた要素当たり仕事量を確認する必要があります。

### single-pixel partition

| rank | local nodes / equation rows | bulk elems | tetra | prism | boundary elems | target BC tag 30 elems / nodes | shared nodes |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 57,890 | 260,136 | 236,916 | 23,220 | 29,541 | 19 / 19 | 4,024 |
| 1 | 57,021 | 99,108 | 0 | 99,108 | 58,196 | 8,401 / 4,348 | 2,555 |
| 2 | 56,951 | 99,196 | 0 | 99,196 | 57,860 | 8,072 / 4,177 | 2,343 |
| 3 | 59,482 | 124,834 | 25,290 | 99,544 | 18,460 | 0 / 0 | 7,289 |
| 4 | 59,394 | 122,270 | 20,817 | 101,453 | 18,130 | 0 / 0 | 7,023 |
| 5 | 59,372 | 120,561 | 18,884 | 101,677 | 19,977 | 0 / 0 | 6,873 |

BC数はmeshの指定境界tagに一致する要素/節点数。ログの `nrows` は上表のlocal nodesと一致し、このP1温度問題のpartition自由度数として記録した。値0のrankを含むrank別カウント集計は `local_count` 区間で計測され、ほぼ0秒だった。

## 2. OpenMP CRITICAL SIGSEGV

- 元の6-rank GDB記録は `CRS_PreparePositionCache()` の `PHASE24_POSITION_CACHE_SETUP` critical入口で停止し、Prepare側criticalを隔離ソース上で外した追試では `CRS_GlueLocalMatrixCached()` の同名critical入口でも停止しました。critical構文だけを消して本番化する変更は行っていません。
- Debug build cacheはGNU Fortran 16.1.0 / GNU OpenMP 5.2、`-fopenmp`、`gomp;mingwthrd`。solver DLL群の依存は同じ `libgomp-1.dll` を指し、実行パス上で確認したlibgompは1個、SHA256 `F53CB4DA77CD3EC8E9EC15C11A8FD5CFCFF448E33B37DA1E1EB225C455028037`。GDBのロードモジュールにも同パスのruntimeが記録されています。明白な異なるlibgompの重複ロードは確認できません。
- 同じGNU Fortran/libgompでnamed criticalの最小EXEテストは1/6 threadsで成功し、DLL内のparallel loopも成功しました。EXEと別DLLで同じcritical名を重複生成する人工テストでは更新競合が起きましたが、CRSの該当named criticalは `CRSMatrix.F90` にのみ存在し、HeatSolve側など別DLLに重複定義はありません。この人工例は本件の直接原因とは言えません。
- `-fcheck=all`で `Phase24ElementNodes` / `Phase24ElementEq` のshape mismatch（6行配列へ4要素全体を代入）を検出しました。診断用複製では代入先を `1:4` に直しましたがVector ONのSIGSEGVは解消していません。この既知不具合は修正対象ですが、今回のcritical障害の根本原因ではありません。
- 従ってABI/runtime重複、named criticalそのもの、配列shape mismatchのいずれも根本原因として確定していません。実行時点のロック記憶域破壊、先行する範囲外書き込み、キャッシュindex/position配列の不整合は未解決です。CRITICALの削除やlock回避を含む本番変更はしていません。根因未確定のため隔離ビルドへの本番候補修正もまだありません。

## 3. Vector ON/OFF数値・速度比較

Vector OFFの4-step DualTESとsingle-pixelは完走しました。Vector ONは現在も最初のassemblyでSIGSEGVとなるため、ON時の温度場、TES回路量、行列、RHS、残差、反復数を比較できず、速度改善も評価していません。数値的同等性が確認できない状態で性能主張は行いません。既存OFF/ON診断はDebug `-O0`でもあり、本番速度比較の証拠には使えません。

## 4. 実験物と再現情報

この追試は同一コミットの証拠用worktree、`%TEMP%\\phase24-debug-src/build/runtime`の隔離Elmer複製で行いました。ユーザーの本番設定・通常Elmer source tree・既存のstage1–3レポートは変更していません。既存の4-step SIF、stdout、series/iteration CSV、時刻・hashを `dirichlet_split_manifest.json` に保存しています。新規の計測ログとpartition値はこのレポートと隣接する `*_dirichlet_6rank.solver.log` にあります。

次に必要な作業は、SIGSEGV発生時にGDBからcritical lock変数・キャッシュ配列のアドレスと前後のheap integrityを採取し、Bounds/Address sanitizer相当または配列境界検査で最初の破壊命令を捕まえることです。MPI側はPMPIで `ParallelReduction` の入退時刻とrank到着時刻を採り、collective内通信量と待ち時間を分離します。Vector ONが完走した後に限り、同一最適化設定の4-stepケースで数値等価性を比較し、その後に性能を測ります。
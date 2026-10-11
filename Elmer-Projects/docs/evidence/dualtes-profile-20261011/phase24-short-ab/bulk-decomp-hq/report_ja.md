# DualTES Bulk律速・Vector Assemblyクラッシュ調査（4-step短時間実験）

実施日: 2026-10-11  
基準コミット: `8bc3e403988450ebb86d9d60654f8c42adf00f21`  
ビルド: GNU Fortran 16.1.0 / GNU libgomp 5.2 / MS-MPI、Debug `-O0 -g3 -fbacktrace -fno-omit-frame-pointer`、`OMP_NUM_THREADS=1`  
ケース: DualTESとsingle-pixel、各4ステップ、各6 MPI rank、Vector Assembly OFF。75 ms本計算は実行していない。計測用barrierは追加していない。

## 結論

DualTESでrank 5が遅い直接要因は、prism一般の演算量ではなく、body 14/16に割り当てられたMaterial 7の温度依存`Heat Conductivity`を、各prismのGeneric assembly経路で繰り返しMATC評価していることだった。rank 5のMATC評価時間は1回のBulk assembly当たり平均約10.59秒で、Bulk平均14.84秒の約71%を占める。MATCの要素当たり評価時間は約0.30 msで、他rankと同程度だった。したがって、rank間の差は特定rankに高コスト材料の要素が偏っていることが主で、rank 5内部の要素演算が異常に遅い証拠はない。

温度依存MATCのnode値を材料ごとのunique node集合で一度だけ評価し、既存の`Phase24EvaluateDynK()`用キャッシュをVector OFFでも再利用する案が有力である。ただし隔離試作はビルドの依存情報生成で止まり、数値比較に到達していない。高速化としては未採用・未検証である。

Vector AssemblyのSIGSEGV原因は未確定。named criticalの実ロック記憶域、実引数をGDBで採取する試みは、Windows GDBの共有ライブラリ・シンボル解決に失敗し、attach時に`0xc0000374`例外を出した。既存のスタックはcritical入口での障害を示すが、ロック破損か先行メモリ破壊かは判別できない。CRITICALは削除しておらず、本番Vector ONも行っていない。

## 計測結果

rank-local wallの6回のassembly平均（秒）:

| ケース | r0 | r1 | r2 | r3 | r4 | r5 | 1回のassembly中の最大値 |
|---|---:|---:|---:|---:|---:|---:|---:|
| DualTES Bulk | 5.019 | 14.405 | 13.908 | 5.293 | 13.178 | 14.844 | 16.884 |
| single-pixel Bulk | 12.108 | 6.987 | 7.757 | 9.188 | 9.197 | 10.018 | 13.183 |

このrank wallは各rankが個別に出力するローカル値であり、MPI最大値ではない。上表の最大値は生ログを後処理して得た。

rank別の1 assembly当たり処理要素数（generic profileで数えたclass 1=element code 504 tetra、class 2=code 706 prism）:

| ケース / rank | tetra | prism | 主な高コスト材料処理 |
|---|---:|---:|---|
| DualTES r0 | 131,633 | 33,060 | Material 7の少ないpartition。Bulk 5.02 s |
| DualTES r1 | 70,716 | 60,348 | body 5/7・Material 7のprism約33,413個。Bulk 14.41 s |
| DualTES r2 | 87,444 | 54,326 | body 5/7・Material 7のprism約30,657個。Bulk 13.91 s |
| DualTES r3 | 133,433 | 33,084 | Material 7の少ないpartition。Bulk 5.29 s |
| DualTES r4 | 99,804 | 48,346 | body 14/16・Material 7のprism約27,807個。Bulk 13.18 s |
| DualTES r5 | 50,513 | 66,028 | body 14/16・Material 7のprism約36,023個。Bulk 14.84 s |
| single-pixel r5 | 18,884 | 101,677 | body 5/7・Material 7のprism約18,852個。Bulk 10.02 s |

要素数は各rankの6回分を6で割った平均。single-pixel r5はDualTES r5より全prism数が大きいがBulkは短い。一方、Material 7のprismはDualTES r5の方が多い。よってprism数単独ではBulk時間を説明できない。

DualTES rank 5のbody/material別内訳（6回平均、秒/assembly）:

| body / material / class | 要素数 | body区間wall | 材料評価 | 局所計算 | global update | 残差 |
|---|---:|---:|---:|---:|---:|---:|
| 14 / 7 / prism | 28,822 | 9.522 | 8.512 | 0.804 | 0.058 | 0.149 |
| 16 / 7 / prism | 7,201 | 2.332 | 2.082 | 0.194 | 0.015 | 0.042 |
| Material 7合計 | 36,023 | 11.854 | 10.594 | 0.998 | 0.073 | 0.191 |

Material 7以外の例として、rank 3のbody 9 / Material 3は1回当たり約76,108 prism、MATC/material評価約0.025秒、局所計算約0.319秒だった。rank 5 body 19 / Material 1のtetra約133,433個も材料評価約0.229秒、局所計算約1.64秒/assembly。処理種別・要素数だけでなく材料評価経路が支配要因である。

single-pixelでもrank 5 body 5/7・Material 7で約18,852 prism/assemblyの材料評価が約5.85秒/assemblyを占めた。同じMATC材料の高い要素当たりコストが再現した。

Material 7はSIFで温度依存MATC式として定義される。現行generic pathは`HeatSolve.F90`で各要素ごとに`ListGetRealArray(Material,'Heat Conductivity',Hwrk,n,Element % NodeIndexes)`を呼び、共有nodeの物性値を要素ごとに再評価する。対して`Phase24EvaluateDynK()`は材料ごとのunique node列をまとめて評価し、動的prism fast pathで使う既存機構である。Vector OFFではこの機構はgeneric assemblyで使われない。

## 計時の意味と信頼性

* `PHASE24_BULK_LOCAL`はrank-local。各body/class区間のwallと、既存の材料評価・局所計算・global update stage timerを同じwall時計で出した。rankごとに6回の計測がある。
* `PHASE24 rank wall`もrank-local。rank 0だけの値でも全rank最大値でもない。
* 旧`PHASE24_EVENT`はElmer `Info()`経由のため、設定既定のOutput PE=0に従うrank 0の値であり、MPI max reduction値ではない。
* CSR/RHSの個別サンプルは1/32呼出しで`omp_get_wtime()`を採り、call/sample比で展開する方法を実装した。しかし生ログでCSRに0や約0.032秒刻みの推定値が混在し、RHSはほぼ0となるため、定量内訳として採用しない。global update合計の方を参照し、この個別分離は未確定とする。
* `omp_get_wtime()`自体は100万回で0.044秒（約44 ns/call）。body/class区間の切替時だけ新たなwall時刻を読むため、タイマー呼出しそのものは処理時間に対して十分小さい。
* すべてDebug `-O0`なので、これらのwall時間をRelease本番性能へ外挿しない。

## 高速化案の試行と状態

隔離ソースだけに`Phase24 Dynamic K Node Cache` opt-inを足し、Vector AssemblyをOFFのままgeneric dynamic prismのk(T)を既存unique-node cacheから読む試作を行った。物性はassemblyごとに現在温度から再評価し、局所剛性/RHS/CSR経路は変えない想定だった。しかし診断ビルドはCMake/NinjaのFortran依存情報生成で停止し、リンク・4-step実行に至らなかった。安全性・数値同等性・速度改善はいずれも未確認で、試作はコミット対象外。

次の安全な検証は、MATC材料を`LIST_TYPE_VARIABLE_SCALAR`に限定する型チェックを追加した隔離ビルドで、OFF baseline対cache ONの同一4-step比較を実施すること。比較項目はtemperature, TES circuit current/power, outer matrix, RHS, residual, iteration history。これに加えてRelease最適化設定で測定を繰り返す必要がある。MPI partition負荷分散は材料別計算コストを重みに加えた実験候補だが、meshと初期状態を固定し、同じ数値比較を完了するまで採用しない。

## Vector Assembly SIGSEGV

既存の4-step ON試験は1 rank / 6 rankとも失敗。初回は`CRS_PreparePositionCache()`のnamed critical、同criticalを隔離ソースでコメント化した別試行は`CRS_GlueLocalMatrixCached()`のnamed criticalでSIGSEGVした。CRITICALを無条件削除して本番化していない。前回の最小named critical再現は成功し、重複libgompロードも未検出であり、今回それらの再試験はしていない。

対象Debug DLLでは`.gomp_critical_user_phase24_position_cache_build`がRVA `0x4d018`、`.gomp_critical_user_phase24_position_cache_setup`がRVA `0x4d020`にあることをobjdumpで再確認した。実行時のmodule baseとロックbytes、critical直前の`Gmtr/Indices/Lmtr`の実アドレスはGDBから採取できていない。GDB attach試行は共有ライブラリsource symbol解決失敗、実行中プロセスへのattach時`0xc0000374`例外で終了した。これは破損原因の確定証拠ではない。配列bounds-checkビルドとCRSキャッシュ配列/引数の監視も未実施。

Vector ONが解決して正常終了するまでは数値A/B比較を行えない。OFFの温度/TES数値結果は各ログ・series/iterations CSVに保存されている。

## 成果物

* `dualtes_6rank.sif`, `singlepixel_6rank.sif`: 再実行用4-step入力。
* `run_bulkdecomp_hq.py`: 実験実行スクリプト（既存runtime workspaceを使う。出力先はこのフォルダ）。
* `manifest.json`: compiler/runtime/build、SIF/output SHA256、return code、時間範囲。
* `dualtes_6rank.solver.log`, `singlepixel_6rank.solver.log`: rank/body/material/class計測の生ログ。
* `*_series.csv`, `*_iterations.csv`: 4-step solver outputs。
* `gdb_attach_critical.gdb`, `gdb_attach_critical.log`: 実ロック調査試行と失敗内容。既存`../debug_manifest.json`とstackログも参照。
* `HeatSolve_diagnostic.patch` ? `SolverUtils_diagnostic_snapshot.F90`: この計測を出すために一時ビルドへ加えたdiagnostic-only差分。製品ソースへの適用はしていない。

## 未解決

1. Material 7 MATCをunique nodeで評価する試作の動作・数値精度・速度は未検証。
2. 6-rank時の局所計算・MATC・collective間同期待ちの因果分離は未実施。bulk計測からMPI待ちを直接差し引いていない。
3. rank-local body/materialの詳細timingはpartitionの物理的element attributionであり、複数rank間のghost/haloコスト分析ではない。
4. Vector Assembly SIGSEGVのcritical lock記憶域、引数、先行する範囲外アクセス/メモリ破壊原因は未確定。
5. Debug `-O0`計測のみで、Releaseでの速度優位性は未評価。

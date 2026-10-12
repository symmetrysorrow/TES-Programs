# Phase24 AST導入後 DualTES MATC Cache rank別監査

## 結論

「no safe temperature-only MATC prisms」は、DualTES全体でAST判定が失敗したことを示すログではありません。AST導入後のReleaseでCache Verifyを有効にして再実行した結果、rank 1・2・4・5ではMaterial 7のMATC Node Cacheが構築・使用され、rank 0・3ではローカルのMaterial 7 prismが0個でした。したがってrank 0・3は対象要素がないため空キャッシュとなった正常経路です。対象prismがあるrankでAST安全性条件に拒否された証拠はありません。

この結果により、前回の [`matc-ast-release/report_ja.md`](../matc-ast-release/report_ja.md) にある「DualTESはGeneric cacheを利用しなかった」という解釈を訂正します。同レポート作成時のVerify OFFログはrank別のBuild/Evalを出しておらず、rank 0/3の空キャッシュ時に出るrankなしのInfo行だけから全rankを判断していました。

## Rank別のVerify結果

同一のAST導入後Release DLL、DualTES 4-step SIF、同じrestart、MS-MPI 6 ranks、Vector Assembly OFFでVerifyを有効にして実行しました。solver終了コードは0で4 stepを完了しています。Verifyは各対象rankからMaterial 7の要素を2個ずつ標本比較し、8個すべて旧評価との差が0でした。

| Rank | Material 7 prism数 | MATC metadata | AST | unique node数 | MATC再評価イベント | MATC呼出数 | cache-hit prism（再利用1 assembly） | 4-step中の合計hit数 | fallback対象element数 |
|---:|---:|---|---|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 対象要素なし、未評価 | 対象要素なし、未評価 | 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | 33,413 | PASS | PASS | 20,441 | 6 | 6 | 33,413 | 167,065 | 0 |
| 2 | 30,657 | PASS | PASS | 18,780 | 6 | 6 | 30,657 | 153,285 | 0 |
| 3 | 0 | 対象要素なし、未評価 | 対象要素なし、未評価 | 0 | 0 | 0 | 0 | 0 | 0 |
| 4 | 27,807 | PASS | PASS | 17,088 | 6 | 6 | 27,807 | 139,035 | 0 |
| 5 | 36,023 | PASS | PASS | 22,024 | 6 | 6 | 36,023 | 180,115 | 0 |

Material 7はSIF上で `Heat Conductivity = Variable Temperature` と単一の `Real MATC` 式として定義されています。既存metadata条件（scalar variable-string、依存名Temperature、procedure/Luaなし、熱伝導率定義が単一）に加え、同じ式は保存済みAST単体ログで受理されています。rank 1・2・4・5のBuild記録はこのmetadata gateとAST gateを含む実際の統合判定が通ったことを示します。rank 0・3では対象element走査が空のため、材料式のrank別判定自体が実行されていません。

Verifyコードの `no safe temperature-only MATC prisms` 行にはrank番号が付いていません。隣接するBuild記録はrank 3とrank 0で `materials` / `unique_nodes` が出ない空キャッシュを示し、builderの当該Info呼出しはローカル適格material数が0の分岐内にあります。したがって該当分岐に入ったrankは0と3です。マージ済みログ上の単一Info行だけから、どちらのMPI processが文字列を出力したかまでは区別できません。

4 timestepの計算ではHeatSolve assembly評価イベントは6回あり、cache対象4 rankで各イベント1回ずつMATC一括評価（計6回/rank）を行いました。最初の評価後の5 assemblyではMaterial 7の全prismがnode cacheをhitしました。総hit数は6 rank合計639,500 element評価です。Verify標本の8要素はすべて `max_abs=0` でした。

## AST導入後Releaseでの反復性能

上記のように対象全rankでcache利用を確認したため、VerifyをOFFに戻し、同一solver/HeatSolve/MATC Release DLL、同じ4-step入力/restart、6 ranks、Vector OFFで各3回実行しました。実行順はOFF1→ON1→ON2→OFF2→OFF3→ON3です。時間は外側wall時計です。

| Cache | 3回のWall秒 | 平均 ± 標本標準偏差 | 中央値 |
|---|---|---:|---:|
| OFF | 74.841 / 74.278 / 75.835 | 74.985 ± 0.788 s | 74.841 s |
| ON | 25.615 / 25.401 / 27.983 | 26.333 ± 1.433 s | 25.615 s |

平均Wall比は2.848倍、平均時間短縮は64.88%でした。6実行すべて終了コード0です。各モードのResult Norm列、HYPRE PCG反復・残差列、TES 1/2 series CSVおよびiteration CSVのSHA-256は3回を通じて同じです。Verify OFFの性能runのSIF、stdout/stderr、生CSV、各rank result shard hashは `performance_repeats.json` に記録しています。

Release DLL SHA-256はSolver `70e1d4cd12b6177b44961f7afbb85379b4faa568b7e44fc432fc1ef4cba683a5`、HeatSolve `9f9a1ddb7e5737defa9508233024880a5d7083e765e69d3685673ba4fa8dc6fd`、libmatc `3fb50c8464cc171ac7b3227394e45d688f1eb2418848481a5badec9325b7fc7a` です。

全node温度場の既存OFF/ON数値比較、TES履歴、assembly matrix/RHS比較はAST Releaseの [`runtime_runs.json`](../matc-ast-release/runtime_runs.json) と [`matrix_rhs_comparison.json`](../matc-ast-release/matrix_rhs_comparison.json) にあります。今回のVerify runでは追加で対象8要素の直接値比較を行いました。性能反復runではresult shardのhashを保存していますが、別個に温度場差分を再計算してはいません。

## 制約・未解決事項

- キャッシュ既定OFF、Vector Assembly OFFのままです。本番DLL/SIFは変更していません。
- 75 ms本計算は実行していません。
- `PHASE24_HYPRE_REPORT=1` による既存capture診断が各MPI rankからstderrへ出力されるため、MS-MPIの統合時に文字単位の交錯が見られます。solver stdoutの完了、終了コード、CSV、result shard出力は正常でした。stdoutとstderrは分離保存しています。
- AST導入済みでもcacheのModel/Mesh寿命管理、solver再初期化を含む別途の安全性ブロッカーはこの実験では解決していません。

## 再実行と成果物

`run_dualtes_rankaudit_perf.py` は性能の3組を同じ順で実行します。`analyze_dualtes_rankaudit.py` はVerifyログからrank別prism数、cache構築・評価・hit、性能の平均・ばらつきを抽出し、`rank_audit.json` を生成します。全ファイルのSHA-256は `evidence_manifest.json` に記録します。AST単体の正例・負例ログも `matc_ast_purity_test.log` として同梱しています。

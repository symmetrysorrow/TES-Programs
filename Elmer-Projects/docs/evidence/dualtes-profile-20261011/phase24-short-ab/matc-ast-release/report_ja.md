# Phase24 MATC Node Cache：再構成Release検証とAST安全性判定

## 結果

公式Elmer `5a8de867068be0568f09af40fb90fee300dfbede` から既存のPhase24 patch 01を適用した16ファイルは、patch 01 manifestのGit blobと再照合しました。続けて独立patch 02を4ファイルに適用し、新しい隔離Releaseをクリーンビルドしました。第三の一時source worktreeでpatch 01、02を順番に `git apply --check`・適用し、両段階のblob検証が成功しました。patch 01は変更していません。

構文木を評価せず検査するMATC APIとHeatSolveの安全ゲートを追加しました。C単体正例・負例試験、MATC既存評価との数値一致、DualTES / single-pixelの4-step OFF/ON実行は完了しました。AST変更後もキャッシュ既定OFF、Vector Assembly OFFを維持しています。Model/Mesh寿命管理や既定ON化は行っていません。

## 実測した4-step runtime

全実行はMS-MPI 6 ranks、同じ再構成Release構成、同じrestart状態、4 steps、Vector Assembly OFFで終了コード0でした。各モードのOFF/ON実行は同じsolver/HeatSolve/MATC DLL群で行いました。最終候補SHA-256はsolver `70e1d4cd…c683a5`、HeatSolve `9f9a1ddb…8dc6fd`、libmatc `3fb50c84…7fc7a` です。完全hashとtoolchainは `evidence_manifest.json` を参照してください。

| ケース | Cache | Wall秒 | 全node温度比較 | Result Norm / HYPRE / TES履歴 |
|---|---:|---:|---:|---|
| DualTES | OFF | 73.078 | 対AST前最大差 `1.39e-16 K` | AST前と一致 |
| DualTES | ON | 31.778 | 対AST前最大差 `5.55e-17 K` | AST前と一致 |
| single-pixel | OFF | 68.389 | 対AST前最大差 `1.39e-16 K` | AST前と一致 |
| single-pixel | ON | 29.879 | 対AST前最大差 `1.11e-16 K` | AST前と一致 |

OFF/ON内でも温度場は丸め誤差範囲です。DualTES最大差は `1.39e-16 K`、single-pixel最大差は `5.55e-17 K`。TES series/iteration CSVのSHAはOFF/ONおよびAST前結果で同一です。Result Norm列とHYPREの反復数・残差も各ケースのOFF/ONと保存済みAST前結果で一致しています。温度場全node・全4時刻の詳細は `runtime_runs.json`、AST前保存結果との比較は `pre_ast_runtime_comparison.json` にあります。

DualTES時刻1でassembly captureを追加して再実行しました。OFF/ONとも76 captureファイルが揃い、sparsity/サイズ情報・行列係数はbyte一致。24 RHS shardの1,101,672値をnode-row keyで比較し、異なる値0、最大絶対/相対差0でした。これは時刻1のcapture結果で、4時刻全ての行列captureではありません。結果は `matrix_rhs_comparison.json` と生capture manifestに保存しています。

## AST API試験

`matc_ast_purity_test_final.log` は失敗0件です。数値・`tx`、基本算術、べき乗、`abs()`およびDualTES Material 7実式を受理し、`tx+ty`、時間状態名、未知関数、ユーザー関数、代入、文字列、制御構文、複数文、括弧不正、末尾token、空/nullを拒否しました。MATCで組み込み `abs` のユーザー関数shadowingがブロックされることも確認しました。受理した式は通常評価と同一値でした。正常なparseと不正parseの後でも出力buffer・`guard`変数・後続評価が維持されました。

single-pixelのONログではMATC cacheが実際に有効になり、Material 7、rank別unique node、build/eval計測行が出ています。DualTESのONログでは既存ゲートが「safe temperature-only MATC prismsなし」としてGeneric cacheを利用しませんでした。ここは安全なfallback動作として数値回帰を確認しましたが、DualTESで今回のAST導入後キャッシュが使われたとは報告しません。DualTES側でどの従来条件が不成立かは次の調査事項です。

## 性能の扱い

同じ最終Releaseバイナリ内の単発OFF/ON Wall値を記録しましたが、反復試験をしておらず、AST APIそのものの性能改善とは評価していません。single-pixelはキャッシュ有効化ログを伴いますが、ばらつき未評価です。DualTESはキャッシュfallbackでした。AST前の保存済みDLL/solverはhashが異なるため、前回報告の速度比を今回の結果として流用していません。

## 実装・ソース再現

- `02_matc_ast_purity_validator.patch` はparser AST検査、C API wrapper/header、Fortran安全ゲートの独立patchです。既存 `01_phase24_source_from_release_26_1.patch` は不変更です。
- MATCの許可文法、builtin identity確認、fail-closed解析と状態復元は `matc_ast_api_design_ja.md` に記録しました。
- `replay_ast_source_and_build.ps1` は公式revision、patch hash、16個のpatch 01 postimage、patch 02 pre/postimageを検証します。`-Build` で同じRelease configure optionとclean buildを実行します。
- Release clean buildは2,237 Ninjaタスク成功。parser.c最終修正後のMATC DLL再リンクも成功しました。

## 未実施・未解決

- Cache lifecycle、Model/Mesh差替え、同一プロセス再初期化、ASTメモリ故障注入は未実施です。
- DualTES Generic cacheが安全条件を満たさずfallbackした理由の各条件別診断は未確定です。
- キャッシュ寿命、異方性や他のMATC表現、未知builtinの安全性は保証しません。
- 4-step以外、75 ms本計算、Vector Assembly ON、既定ON移行、性能反復測定は未実施です。

## 再実行

`run_final_release_ab.py` はDualTES/single-pixel各4-step OFF/ONを実行します。`run_final_matrix_rhs_capture.py` はDualTES時刻1の行列/RHS captureを行い、`compare_matrix_rhs.py` が比較します。`matc_ast_purity_test.c` はMATC DLLにリンクして単体試験します。実行時のツールパスはスクリプト先頭の定数とmanifestに記載しています。

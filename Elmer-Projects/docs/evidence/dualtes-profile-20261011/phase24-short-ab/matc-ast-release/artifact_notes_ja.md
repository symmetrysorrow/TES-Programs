# 証拠ファイルの読み方

主要な結果は `report_ja.md` と `evidence_manifest.json` にまとめています。最終版の runtime 生ログ・SIF・TES CSVは `final_` prefix、確定した行列/RHS captureは `final_matcap_` prefixです。AST前の同じPhase24再構成Release結果は `recon_` prefixと `runtime_baseline_pre_ast.json` です。

次のファイルは作業中に修正・置換した試行の記録です。最終PASS/FAIL判定には使っていません。

- `build_matc_ast_failure.log`: 最初のunit-test compileでinclude pathを1階層誤ったため失敗。その後正しいinclude pathで再ビルドし、`matc_ast_purity_test_final.log` が全件PASS。
- `ast_*.solver.log` / `ast_*.sif`: 最終purity条件のコードレビュー前に実行した中間run。最終検証は `final_` prefix。
- `ast_matcap_*.solver.log`: 最初のcapture試行はcapture出力dir未作成で失敗。ディレクトリ作成を追加し、最終DLL/最終patchで `final_matcap_*.solver.log` を終了コード0で再取得。
- `*_pre_review.json`, `runtime_*pre_purity_review.json`: purity条件のレビュー前に取得した中間値。最終比較は `matrix_rhs_comparison.json` と `pre_ast_runtime_comparison.json`。

この中間ログは実験の履歴として保持しており、成功した最終結果と混同しないよう明示しています。

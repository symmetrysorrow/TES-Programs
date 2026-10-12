# MATC AST Temperature-only 判定 API

## API

`mtc_expression_is_temperature_only(const char *)` はMATC式を実行せず、既存の `parser.c` の `parse()` で構文木を作り、その木だけを検査します。評価器 `evalclause()` とMATC関数呼出しは通りません。

受理する構文は数値リテラル、`tx`、括弧を含む式、単項マイナス、二項の `+ - * / ^`、および監査済みの単項 `abs()` です。AST node形状も検査し、余分なリンク、文、代入、未知node、関数引数リストは拒否します。式長4095文字以下、深さ256以下、node数4096以下の制限もあります。

`abs` は名前だけで判断しません。MATC command registryの登録先が実装関数 `func_abs` と同一であり、引数数が1で、同名のユーザー関数がないことを確認します。`func_abs` の実装は `return abs(arg)` だけです。`CMDFLAG_PW` は純粋性を示すフラグではなく「element-by-element operation」なので、純粋性の根拠には使いません。今回のMATCではユーザー関数が組み込みcommandを上書きできないことを実行試験で確認しました。

## fail-closedと状態復元

- 空/null式、長すぎる式、未許可識別子、関数、演算子、AST形状、末尾token、parse error、内部エラーはすべて拒否します。
- parserの文字列、現在token、lexer bufferを退避して復元します。解析中は独立したMATC allocation listとlocal `setjmp` を使い、ASTを含む一時割当てを解放してから元の状態に戻します。
- C API wrapperはMATC出力bufferとそのcount/capacityを退避し、parse-only中の値を分離して復元します。テストは成功解析・構文エラーの両方で既存出力と変数 `guard=41`、後続のMATC評価が維持されることを確認しています。
- APIはMATC contextの初期化やユーザー関数評価を行いません。ASTに現れない外部状態にはアクセスしません。

## Fortran側の適用

既存 `Phase24IsTemperatureOnlyMatc()` の追加必要条件としてC APIを呼びます。既存のscalar-string型、`Temperature`依存名、Lua/procedure不使用、異方性/追加conductivityなし、単一conductivity条件を維持します。AST判定で拒否した値は既存Generic経路に残します。キャッシュ既定値OFF、Vector Assembly OFFです。

## 確認範囲と限界

単体試験では実際のDualTES Material 7式、算術正例、`tx`以外の変数、未知関数、ユーザー関数、代入、文字列、制御構文、複数文、構文不正、末尾token、空/nullを確認しました。受理された `tx*tx`、`(tx+1)^2`、`abs(tx)`、Material 7式は通常のMATC評価と比較し、一致しました。

このAPIは「Temperature以外の入力を持たない式」という狭い文法だけを安全対象にします。モデルやMeshの寿命管理を保証するものではありません。AST不正・内部エラー試験とNode Cacheの寿命試験は今回の対象外で、未実施のままです。

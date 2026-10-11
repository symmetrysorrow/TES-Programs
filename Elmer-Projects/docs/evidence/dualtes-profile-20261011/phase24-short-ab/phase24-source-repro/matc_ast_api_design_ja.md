# MATC pure-expression AST validator設計案

## 目的と境界

MATC `CValue` の依存性を式本体からfail-closedで判定し、温度だけに依存するscalar式だけをNode Cache候補にする。今回の作業ではAPIも評価処理も変更していない。

## 既存実装の観察

公式Elmer 26.1のMATC `parser.c` は `TREE`/`CLAUSE` を作り、`eval.c` が実行する。Elmer `Lists.F90` は `ValueListEntry_t%CValue` をMATCに渡すが、Fortran側には式本体のpure/dependency判定APIがない。MATC parser stateはグローバル/THREADPRIVATEなscanner状態を含むため、既存 `doit()` を呼ぶ解析判定は式を実行しうる上、評価stateを汚染しうる。これは安全判定として使用しない。

## API案

MATCライブラリ内に、例えば `MatcValidatePureScalarExpression(expr, result, diagnostic)` を追加する。APIは式を一切評価せず、成功時だけ `PURE_SCALAR_TX` を返す。一般式を評価する `Matc()` / `doit()` の動作は変更しない。

1. 既存scanner/parserと同じ構文規則を使うparse-only entryを作る。入力bufferとscanner位置は専用parser contextに分ける。暫定的に既存threadprivate globalsを使う設計なら、全scanner状態を保存・復元し、parse errorのlongjmp/cleanupを安全に処理することを先に証明する。
2. 式全体をASTにparseし、EOF到達・scalar expression単一rootを確認する。文字列の部分一致やregexで許可を決めない。
3. 許可nodeを数値リテラル、識別子 `tx`、監査済みの単項符号/`+ - * / ^`などの算術opだけに限定する。比較、論理、代入、matrix/vector、index、transpose、apply、制御文、system call、文字列、comment/statementは拒否する。
4. `abs`を含む関数nodeは、AST node名だけで許可せず、builtin registryで実装identityを照合し、MATC user-function registryに同名定義がないことを同時に確認する。明示登録のない関数はすべて拒否する。registryがshadowingやruntime差替えを許すなら、その状態を評価時にも固定できない限りキャッシュ候補から外す。
5. `tx`以外の全identifier、global variable、constant alias、user functionは拒否する。未知nodeやvalidator内部エラー、parse errorは「不明」ではなく必ずfalseで返し、ElmerはGeneric経路へfallbackする。
6. Elmer側では既存 `DependName`/`LIST_TYPE_VARIABLE_SCALAR_STR`/scalar-isotropic条件に追加の必要条件としてこのAPIを呼ぶ。AST validator trueだけでMATC dependency metadataのチェックを置き換えない。

## Material 7の正例候補

現DualTES SIFの式:

```text
(7.854e-06)*(0.5*((tx+1.0e-12)+abs(tx-1.0e-12)))^(4.252-1)*0.4*((0.0007)-(0.0005))/(8*((1e-06)+(1.5e-05)+(1e-06))*(0.0005))
```

これは目視ではnumeric/tx/arithmetic/power/absだけだが、AST判定とMATC runtime比較を通過するまでは正例としてPASSにしない。`abs` builtinのshadowing不可と純粋性も確かめる。

## 独立テスト設計

APIを実装したらMATCライブラリ単体のテストで以下を実行し、その後Elmer integrationを試す。

| 式 | 期待 |
|---|---|
| `tx*tx`, Material 7実式 | allow候補。MATC値を通常評価経路と比較 |
| `tx+ty`, `tx+x`, `tx+temperature`, `pi*tx` | reject（`tx`以外の識別子） |
| `tx+unlisted_time_state` | reject（隠れ状態依存） |
| `foo(tx)`, user-defined `abs`, builtin registry外関数 | reject |
| `tx=tx+1`, `begin ... end`, `if`, `$system` / `apply(...)` | reject（代入/制御/副作用） |
| unmatched parenthesis, trailing token, empty string | parse failureとしてreject |
| `tx/0`, `tx^large` | pure判定と数値domain validationを分離。例外/NaNの扱いをGeneric評価と同一にするまでreject候補 |

AST許可はpure/dependencyだけの証明であり、数値domainや物理範囲の証明ではない。候補式の同じtx値で通常のMATC評価と比較し、roundingを含む一致基準を別途定義する。

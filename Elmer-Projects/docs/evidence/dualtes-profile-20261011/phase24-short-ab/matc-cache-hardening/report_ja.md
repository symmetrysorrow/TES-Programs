# MATC Node Cache 汎用統合ブロッカー調査

## 判定

**安全性ゲートは未通過です。キャッシュ既定OFFを維持し、今回の成果を既定OFFの隔離実験候補として扱います。** 本番DLL/SIF、Vector Assembly設定、既存成功済み実験証拠は変更していません。75 ms本計算は実行していません。

今回、MATC式の構文と式保存形式、Elmer 26.1の元ソース、既存キャッシュ候補の適用範囲を調べました。式本体を判定するAPI、寿命ガード、再現可能なPhase24ソース差分はいずれも完成していません。したがって修正DLLも作成せず、未実施テストをPASSにしていません。

## 実測および既存テストの整合

前回の既存Release OFF/ON測定は既知の実験候補パッチの測定値です。DualTESで71.004 s対25.615 s（約2.77倍）、single-pixelで53.666 s対28.700 s（各1試行）でした。これらは安全性修正後のバイナリの結果ではなく、今回の統合可否を支える証拠にはしません。

前回のTime依存runtime負例は、6 rank・4 step・既存DualTES入力の検査用コピーで実行済みです。return code 0、generic fallbackログあり、実測95.783 sでした。この1ケースを既存 `safety_manifest.json` の `runtime_negative_tests_executed` と件数に反映しました。これはTime依存のfallbackだけのPASSであり、式内隠れ依存や寿命試験のPASSを意味しません。

今回、対象キャッシュパッチをクリーンな公式Elmer 26.1 sourceに対して `git apply --check` しました。`HeatSolve.F90` の最初のhunk（line 62）が適用できず失敗しました。パッチは既存Phase24変更を前提とするため、クリーンなElmer 26.1からの再現ビルドにはそのまま利用できません。

## MATCの式表現と判定可能性

公式Elmer `release-26.1` の参照commitは `5a8de867068be0568f09af40fb90fee300dfbede` です。MATCの `ValueListEntry_t%CValue` は `Lists.F90` からMATC文字列として `Matc()` に渡されます。MATC本体は `parser.c` で構文木（`TREE` / `CLAUSE`）を生成し、`eval.c` が評価します。式保存自体は文字列で、現在Elmerの公開Fortran側に純粋性や識別子を検査する関数は確認できませんでした。

DualTES Material 7の実式は次の形です。

```text
(7.854e-06)*(0.5*((tx+1.0e-12)+abs(tx-1.0e-12)))^(4.252-1)*0.4*((0.0007)-(0.0005))/(8*((1e-06)+(1.5e-05)+(1e-06))*(0.0005))
```

字句の目視上は数値定数、`tx`、括弧、基本算術、べき乗、`abs`のみです。ただしこれは監査者の目視であり、MATC構文木を走査した証明ではありません。式本体を調べず`DependName=Temperature`だけに依存する現行ゲートでは隠れたグローバル変数、ユーザー関数、未知構文を除外できません。文字列検索の追加ではこの問題は解決しません。MATC parser内でparse-onlyのAST検査を実装し、数値/tx/許可演算子/実装とshadowing条件を確認できる純粋なbuiltin関数だけを許可する変更が必要です。parse error時は常にGeneric fallbackにする必要があります。

この判定にはMATC CライブラリとElmer Fortran側の両方の変更、および対象の実際のMATC DLL/SOを使った負例runtime試験が要ります。今回はその変更・ビルド・runtime確認はしていません。

## 寿命管理とソース再現性

確認できるElmer標準の元revisionは上記26.1 tagです。しかし隔離中のPhase24 `HeatSolve.F90` はGit管理外で、公式26.1ソースから大きく差分があります。隔離ソースSHA-256は `bde260c7b801349b6c21beacc2772cadc0672fb373811e17d8bbb244fea6193d`、追跡済み証拠コピーは `965FA34D9D051467D821E87C62806C8CEA866E1478E3FC83F6BC5D20D2FE3F50` です。前者の完全な差分履歴・対応する累積パッチ一式を再構成できなかったため、この状態を「元revisionから再現可能」とは扱いません。

MATC cache配列は `HeatSolver` 内部状態に `SAVE` を付けて保持されます。clear routineはあるものの、solver破棄時に必ず呼ばれるfinalize hook、Model/material世代照合、要素接続の世代照合を確認できませんでした。solver呼び出しスコープのallocatableにする案は検討できますが、キャッシュだけでなくPhase24由来のSAVE状態やmesh再利用を含めた再初期化テストをしていないため、今回はソースパッチを作って安全と主張しません。

## テスト状況

| 対象 | 状態 | 根拠 |
|---|---|---|
| Time依存宣言のgeneric fallback | PASS_RUNTIME（既存試験） | 6-rank 4-step、return code 0、fallbackログあり |
| DualTES Material 7式のAST純粋性 | NOT RUN | AST validator未実装 |
| 隠れたMATC依存、未承認関数、構文負例 | NOT RUN | 安全化DLLなし |
| 座標依存、異方性、複数material共有node | NOT RUN | 新DLL/runtime fixtureなし |
| 同一プロセスModel/Mesh/material変更 | NOT RUN | identity/generation管理なし |
| solver再初期化、OFF→ON→OFF、finalize解放 | NOT RUN | lifecycle修正・試験なし |
| DualTES/single-pixel 4-step回帰 | NOT RUN | 安全化ビルドなし。前回の数値結果は別候補の既存証拠 |
| 安全化Release速度 | NOT RUN | 数値安全性ゲート未通過 |

## 次に必要な作業

1. 公式26.1 commitを固定し、Phase24差分を構成する全てのソース変更を出所別に整理して、hash事前条件付きの再現パッチを作る。現在の候補パッチをクリーンソースへ直接適用しない。
2. MATC parserに評価を行わないAST pure-expression判定APIを追加し、Material 7式を正例、隠れた名前・関数・代入・制御構文・未知構文を負例として検査する。
3. cache storageをHeatSolver呼出しまたは明示solver instanceの寿命に結び付ける。Mesh/material/body接続変更時は安全に破棄・再構築し、solver終了時の実際の解放を確認する。
4. 隔離DLLで負例、共有node、同一プロセスModel/Mesh変更、再初期化、OFF→ON→OFFを実行し、全fallbackログと生ログを保存する。
5. その後に同じRelease設定でDualTES/single-pixel 4-step数値回帰を行い、一致した場合のみ速度比較する。

**今回の結論：** MATC本体には構文木を使った安全判定を実装できる余地がありますが、現候補にその判定はなく、キャッシュ寿命とPhase24ソース由来も証明できていません。汎用統合の安全性ゲート通過とは報告できません。キャッシュ既定OFFを維持してください。

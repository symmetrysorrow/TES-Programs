# DualTES / Single-pixel MATC Node Cache：安全性監査レポート

## 判定

**現状のパッチは「既定OFFの実験候補」としてのみ保持し、汎用の安全な統合候補とは判定しません。** キャッシュOFFは既定値のままです。今回の監査では本番ソース、DLL、SIF、材料式を変更していません。既存のRelease性能・数値結果は再利用できる範囲を明示し、キャッシュ判定と寿命管理の不足をコード監査で検出しました。

この判定は「計算結果が既に比較されたDualTES / single-pixelケースで誤差が見つかった」という意味ではありません。既存比較は良好ですが、任意のTemperature-only MATC式や同一プロセス内のModel/Mesh変更に対する安全性を立証できていません。

## 既存Release実測の再利用

同一の隔離Release `-O3`、MPI 6 ranks、Vector Assembly OFF、4-step入力で行われた前回計測を参照しています。

| ケース | Cache OFF | Cache ON | 試行数 | 解釈 |
|---|---:|---:|---:|---|
| DualTES Wall中央値 | 71.004 s | 25.615 s | 各3回 | 約2.77倍。既存候補パッチの実測値 |
| single-pixel Wall | 53.666 s | 28.700 s | 各1回 | 約1.87倍。ただしばらつき評価は未実施 |

DualTESの全node温度場の最大絶対差は `5.55e-17 K`、single-pixelは `8.33e-17 K`。DualTESの全6 nonlinear solveで行列係数とsparsity metadataはbyte一致し、RHSは最終ステップだけ最大 `1.02e-20` の丸め差でした。TES series、iteration履歴、Result Norm、HYPRE反復数・残差は保存済み比較で一致しています。値と条件は `../matc-cache-release/` のReleaseレポート、manifest、比較JSON、生ログを参照してください。

既存4-step入力はrestart time `20 ms` から開始し、pulse start `20.02 ms`、pulse duration `1 ns`、刻み `18 us, 1 us, 1 us, 1 ns` の短い区間です。パルス時刻はSIF記載値です。温度・TES状態を含む前回のOFF/ON数値比較は有効な既存証拠ですが、BDF次数切替、restart後の同一プロセス初回assembly、MeshChanged、Model差替えを網羅した証拠ではありません。

## 今回のコード監査

対象は `bulk-matc-cache/HeatSolve_generic_matc_node_cache.patch`（SHA-256 `53f97c44e82215fb3d230c0b90b024b126cd1dfde716a3ff3ec56c171c02c045`）です。実装ファイルはElmer本体の非Git管理ソースにあり、Git管理対象リポジトリにはありません。既存パッチを新しい検証フォルダにそのまま保存しました。パッチ自体は変更していません。

監査スクリプトは、適用条件関数、cache clear/build/evaluate処理、初期化分岐、既存Release比較manifestを読みます。結果は `safety_audit_results.json` に保存しました。

確認できたガードは次のとおりです。

- MATC scalar string型、宣言依存名が `Temperature` のみ、Lua/procedureなし、Heat Conductivity定義が1件という判定。
- Material IDごとにrank-local node mapを分け、同じnodeを共有するmaterialの値を同じslotへ混ぜない構造。
- nonlinear assemblyごとに `ListGetReal` で現在のTemperatureからキャッシュ値を再評価する経路。
- option既定値OFF、Vector Assembly OFFで独立する構成。

不足している点は安全性の統合を妨げます。

1. 判定はValueListの依存メタデータを見ますが、MATC式本体 `CValue` の内容や純粋性を検査しません。従って「MATC式に隠れた依存や副作用がない」ことは機械的に立証されていません。`Temperature` と宣言されていても式を読まない判定は、要求された負例を保証できません。
2. cacheはFortran `SAVE` 配列で、clearの主な条件は `AllocationsDone` / `MeshChanged` です。Model/material/body-to-material構成の世代・署名を保存し照合する処理が見当たりません。同じプロセスでModelが入れ替わるケースやMesh接続変更の寿命テストは未実施です。
3. clearルーチンはありますが、solver終了・破棄時に必ず呼ぶfinalize hookを今回確認できていません。同一プロセスでの再初期化やメモリ解放は未検証です。
4. 診断用ElmerソースのGit revision / pristine source hashを既存Release manifestから特定できません。既存Release DLLはハッシュ記録済みですが、同一の元ソースに新しい修正を再適用して再構築した結果ではありません。

## ケース別テスト状況

`matc_dependency_cases.json` に再実行対象となる材料fragmentと期待動作を記録しました。`audit_cache_safety.py` はコード上の静的なガードを検査し、解決できない項目をPASS扱いにしません。

| ケース | 今回の結果 | 根拠 |
|---|---|---|
| 宣言Temperature-only scalar MATC | 既存DualTES / single-pixelで数値一致を確認 | Release A/B証拠。任意のMATC式一般の証明ではない |
| 宣言Time依存 | **runtime PASS** | 6-rank 4-step検査用コピーでreturn code 0、generic fallbackログを確認。 |
| 座標依存 | 静的ガードPASS、runtime SKIP | 依存名完全一致の条件で候補から外れる。runtime入力は未実施 |
| 異方性 / 複数Heat Conductivity | 静的ガードPASS、runtime SKIP | scalar型と単一定義のコード条件を確認。SIF実行は未実施 |
| Lua / 外部procedure | 静的ガードPASS、runtime SKIP | `PROCEDURE` / `LuaFun` の条件を確認。SIF実行は未実施 |
| Temperature宣言と式内の隠れた依存 | **静的監査FAIL** | 式本体を検査していない。MATC runtimeでのfixture挙動は未確認 |
| 同nodeを共有する複数dynamic material | SKIP | per-material mapはコード上分離。Elmer runtimeの直接比較なし |
| 同一プロセスでModel/Meshを差替え | **寿命ガード不足、runtime SKIP** | model/material signatureなし。旧値不再利用を実証できない |
| OFF → ON → OFF、restart初回assembly | SKIP | 同一プロセス制御テストなし |
| 大きな過渡 / BDF切替 | SKIP | 既存短区間はpulse付近を含むがBDF次数切替と状態再初期化を網羅しない |
| Release性能 | 既存値を再利用 | 今回は安全性修正と新DLLがないため再測定なし |

`SKIP` は合格ではありません。明示Time依存の負例は既存の隔離Release DLLで追加実行し、6-rank 4-stepをreturn code 0で完了、ログに `no safe temperature-only MATC prisms` が出てgeneric fallbackを確認しました。入力は既存DualTES SIFの検査用コピーで、Material 7の依存宣言だけをTemperatureからTimeへ変更しています。実測Wall 95.783 sはfallback検証の値で、性能比較や物理数値の証拠には使いません。入力SIF・生ログ・JSONをこのフォルダに保存しました。

MATC runtimeの座標依存、異方性、Lua/procedure、複数material共有node、隠れ依存、Model/Mesh寿命、OFF→ON→OFF、restart初回assemblyは未実行です。安全性修正済みの隔離DLLは作成していません。診断用の既存ソース/DLLを本番へコピーしておらず、製品DLL・設定も変更していません。

## 次の作業と統合条件

1. 元Elmer revisionと全対象source hashを固定し、そのrevisionに対するパッチ適用前チェックを追加する。hash不一致なら適用を止める。
2. MATC依存判定を式本体を含めてfail-closedにする。安全性を構文・関数集合から証明できない式はGenericへ戻す。既存Material 7式を含む正例を維持する。
3. cacheをsolver初期化・MeshChanged・Model/material変更・終了のライフサイクルに結び付け、識別不能なら再構築または無効化する。単なるポインタアドレス一致に依存しない。
4. 独立した小規模Elmer入力で負例runtime fallback、共有nodeの複数material、temperature変化、restart、同一プロセスmodel/mesh置換を実行する。
5. 修正DLLでDualTES/single-pixelの4-step pulse区間を数値比較し、数値合格後に同一Release構成で速度を再測定する。

**現時点の答え：既存ケースでは大きな高速化と浮動小数点丸め範囲の一致が実測されています。しかし、物理モデルやpixel数が変わっても誤った物性値を使わない汎用機構としては、現パッチの依存式検査と寿命管理の証拠が足りず、安全に提供できる状態ではありません。既定OFFを維持してください。**

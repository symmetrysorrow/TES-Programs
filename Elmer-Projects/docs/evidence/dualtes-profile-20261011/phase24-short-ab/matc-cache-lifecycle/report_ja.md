# Phase24 MATC Node Cache：activation局所寿命管理の検証

## 判定

Generic MATC Node Cacheの配列・mapを、`HeatSolver` 1回のactivationに属する非`SAVE` allocatableローカル状態へ移した。cache defaultはOFFのままで、Vector AssemblyもOFF。本番DLL/SIFと75 ms解析は使用していない。

隔離ReleaseでAST負例・正例、DualTES/single-pixel 4-step OFF/ON、rank別cache診断、DualTES timestep 1の行列/RHS capture、DualTES各3回の性能比較を実施した。実施済みの数値ゲートは通過した。**明示的ONでのopt-in候補（A）**と判定する。未知のModel/Mesh変更・同一プロセス再初期化を含むruntime寿命試験とmemory sanitizerは未実施のため、**既定ON候補（B）ではない**。default OFFを維持する。

## ソース変更・寿命の根拠

対象ソースはElmer公式 `5a8de867068be0568f09af40fb90fee300dfbede`。再構成済み01、AST 02は変更せず、独立した `03_matc_cache_lifecycle.patch` をその後に適用した。`evidence_manifest.json` にパッチSHA、03前後のGit blobとSHAを記録した。01→02→03適用・preimage検証・Release buildの再実行手順は `replay_01_02_03_and_build.ps1` に保存した。

03は`HeatSolve.F90`のGeneric MATC cache状態だけを変更する。

- `Phase24GenericKNodeSlot/Nodes/Mats/Start/Count/Nodal`、ready/count/hit/eval/verify状態から`SAVE`と宣言時初期値を除去した。
- `HeatSolver`の最初の処理としてclearを行い、早期return前でも当該activationの状態を残さない。
- 既存clearは`ALLOCATED`確認後にのみ解放し、状態カウンタも初期化する。allocatableローカルは通常return時にFortranが自動解放する。
- cache buildは現在の`CurrentModel`のBodies→Material対応と現在のrank-local `Solver % Mesh`からmapを作る。slotはmaterialごとに別に割り当てるため、異なるmaterial間でnode番号を共有しても評価値を誤共有しない。
- 同一activation内ではnode mapを使い回す一方、各非線形Assemblyの直前にmaterial/node値を再評価する。確認ログではactive rankごとに`ListGetReal`呼び出し1回/material/Assemblyだった。
- 関連するElmer呼び出し経路を調べ、HeatSolver activation中にModel/Mesh/material/connectivityを書き換える呼び出しは確認されなかった。同期した通常のElmer lifecycleでは、次activationの冒頭で現在値から再構築する。外部スレッドによるactivation中の並行変更は保証対象外。
- 他のPhase24 `SAVE`状態、行列/RHS組立て、TES回路、時間積分は変更していない。

### ローカル状態の静的寿命監査

|項目|結果|
|---|---|
|早期return|clearをHeatSolver先頭に置いたため、既存activationからの状態持ち越しはない。build前のreturnも未割当状態。|
|再帰|`HeatSolver`のrecursive activationごとにローカルallocatable状態が分離される。|
|正常return/early return|allocatableローカルはactivation終了時に自動解放。明示解放はALLOCATED guard付き。|
|OpenMP|この作業の確認範囲ではcache map/値は構築・評価後にassembly中read-only。HeatSolver自体の一般的なreentrancyを保証する変更ではない。|
|MPI|各rankのprocess-localローカル配列。rank間共有なし。|
|in-flight Model/Mesh mutation|同一activation内の変更経路は監査範囲で見つからず。並行外部変更の安全性は未保証。|

## 4-step寿命・数値試験

同一の隔離Release DLL、同一restart state、MPI 6 ranks、同一4-step SIF条件でDualTESとsingle-pixelをOFF/ON実行した。キャッシュONのVerify診断は別の同一条件DualTES runで有効化した。各通常OFF/ON runのreturn codeはすべて0。

|比較|DualTES|single-pixel|
|---|---:|---:|
|温度node数/step|rank shards全件を比較、全4 step|rank shards全件を比較、全4 step|
|最大絶対温度差|`1.388e-16 K`|`8.327e-17 K`|
|最大相対差|`8.227e-16`|`4.935e-16`|
|TES series/iteration CSV|4ファイル、byte一致|2ファイル、byte一致|
|Result Norm列|完全一致|完全一致|
|HYPRE反復数・残差列|完全一致|完全一致|
|solver終了|OFF/ONとも0|OFF/ONとも0|

DualTESの行列/RHS captureはtimestep 1の全実行非線形assembly（capture対象2回）で行った。OFF/ON各76 shard/fileが全て対応し、matrix/sparsity関連の非RHSファイルはbyte一致。RHSは24 rank-shard、合計1,101,672値でrow map一致、差分0（max absolute/relative差0）。**timestep 2–4の全matrix/RHSは取得していない。**

AST単体テストを03 Releaseの新しい`libmatc.dll`で再実行し、正例・負例すべてPASS（failures=0）。ログは`matc_ast_unit.log`。

全node温度にはOFF/ON間の最大数ulp程度の差があるが、Result Norm、HYPRE履歴、TES CSV、capture範囲内のmatrix/RHSは完全一致。cache計算経路の数値同等性は確認した。

## 全rankの寿命診断

DualTES Verify ONの再起動4-step実行はreturn code 0。次表は`lifecycle_dualtes_verify.solver.log`から抽出した。cache mapは各rankで4 activationに各1回buildされた。

|rank|prism対象material|unique nodes|map build回数|MATC評価Assembly|観測された前Assembly hit数/対象prism数|
|---:|---:|---:|---:|---:|---:|
|0|なし|0|4|0|0|
|1|1 (Material 7)|20,441|4|6|33,413|
|2|1 (Material 7)|18,780|4|6|30,657|
|3|なし|0|4|0|0|
|4|1 (Material 7)|17,088|4|6|27,807|
|5|1 (Material 7)|22,024|4|6|36,023|

rank 0/3の「no safe temperature-only MATC prisms」は対象prismを持たないrankの空cacheであり、安全性拒否を示すものではない。この実行では警告が計4回出た。rank 1/2/4/5はMaterial 7が1つずつacceptedとなり、1回/material/assemblyの`ListGetReal`評価と過去assembly hitが確認できた。cache利用ゼロのactive rank、安全性条件による拒否は見つからなかった。

rank別map構築時間はactive rank平均約0.303–0.388秒/activation、空rank平均約0.0048–0.0050秒。評価時間は各assemblyでrank 1/2/4/5ごとに概ね0.9–1.3秒（Verify計測、Assembly内の値評価）。Verify計測を性能比較に含めていない。

## Release性能比較

同じ03 Release binary、MPI 6、4 step、Vector OFF、Verify OFFでDualTESをOFF/ON各3回、交互に測定した。

|mode|各run (秒)|平均 ± 標本SD|中央値|
|---|---|---:|---:|
|Cache OFF|91.009, 90.357, 87.076|89.481 ± 2.108|90.357|
|Cache ON|33.417, 32.493, 31.239|32.383 ± 1.093|32.493|

平均Wall比は**2.763倍**（同一build内）。旧02段階の平均2.848倍は別binaryの過去証拠であり、今回の成果値とは混同していない。03のVerify build診断でactivationごとにmapを再構築する追加時間を計測した。`Phase24GenericKNodeSlot(6,nBulk)`等からの静的割当概算は、最大rankのbulk要素数を用いておよそ16 MiB以下/active process（allocator overheadを除く、build時NodeMarkは一時追加）。RSS増分は直接測定していない。

## テスト一覧

|テスト|状態|根拠/制限|
|---|---|---|
|01→02→03 hash検証・隔離Release build|PASS|`release_build.log`末尾に`PATCH_SEQUENCE_01_02_03_PASS`と`PATCH_03_RELEASE_BUILD_PASS`。|
|MATC AST正例/負例|PASS|新libmatcで全テスト、failures=0。|
|同じModelの複数step・nonlinear再評価|PASS|DualTES 4-step、rank別map build 4回、各active rankで6回のMATC評価。|
|restart後の初回Assembly|PASS|既存restart stateから両モデルを4-step実行。|
|異なるmaterialの共有node|静的監査PASS / 専用runtime SKIP|slot割当をmaterialごとに独立化するコード。現在の入力で専用合成共有node試験は作っていない。|
|同一プロセスCache OFF→ON→OFF|SKIP|SIFのsolver parameterを同一Elmer process内で切替える既存経路がなく、この4-step runは別process。activation内stateは非SAVEで毎回初期化。|
|同一プロセスsolver再初期化|SKIP|専用driver未用意。通常の別activation/restart runは確認。|
|同一プロセスmaterial変更|SKIP|専用mutation driver未用意。各activationで現Modelから再構築する静的根拠あり。|
|node数一定でconnectivity変更|SKIP|専用mesh mutation driver未用意。各activationで接続から再構築。|
|bounds/memory sanitizer|SKIP|Release buildにbounds/memory sanitizerを付けていない。二重解放回避はALLOCATED guardと自動解放の構造で監査。|
|75 ms本計算|NOT RUN|制約に従い未実施。|

## 統合判定

**A — 明示的Cache ONのopt-in候補：PASS。** AST fail-closedは変わらず、Generic cacheの状態はHeatSolver activationを越えず、現行restart DualTES/single-pixelの数値比較とrank別実行で動作を確認した。

**B — 既定ONの汎用統合候補：未達。** Mesh/material/connectivity変更、同一プロセス再初期化、境界ケースの専用runtime試験とmemory sanitizerがSKIP。既定値はOFFのままとする。

## 再実行・成果物

`replay_01_02_03_and_build.ps1`で上流revision固定・01→02→03 hash検証・Release buildを再現できる。runtime入力と状態、SIF、ログ、比較JSON、ソース/ DLL hashは本ディレクトリと隣接する既存`phase24-source-repro` / `matc-ast-release` / `bulk-matc-cache` evidenceに保存している。行列/RHS captureの大きなraw shardは隔離Temp runtime内に置き、Gitには全shard SHAと比較結果を保存した。再生成には`run_lifecycle_matrix_rhs_capture.py`を使う。

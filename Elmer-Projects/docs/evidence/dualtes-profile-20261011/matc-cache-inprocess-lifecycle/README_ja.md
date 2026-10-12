# MATC Node Cache 同一プロセス寿命調査

## 結果

**B: 実行ハーネス未完成。PASSは報告しない。**

公式Elmer 26.1では `Mesh Mode` が同一 `ElmerSolver` プロセス内で `LoadModel` とsolver初期化を繰り返す経路であることを確認した。`ReloadInput` は安全なModel/Mesh差し替えAPIとして使わない。

今回は実際の同一PIDログとcacheの解放ログを得る04パッチおよびランナーを完成できていない。既存01・02・03の再監査や4-step結果を今回の寿命検証として扱わない。

## 公式ソースとパッチ適用済み構造

公式tag `release-26.1`、commit `5a8de867068be0568f09af40fb90fee300dfbede` を確認した。既存01→02→03は隔離ソースに順番通り適用可能で、MATC cache stateは `HeatSolver` の非 `SAVE` local allocatableになっている。

- `fem/src/ElmerSolver.F90:388-413`: `Mesh Mode` の繰り返しで `FreeModel(CurrentModel)`、`MeshIndex` 増分、同じ `ModelName` で `LoadModel(...)` を実行する。1つのElmerSolverプロセスの実行ループ内にあるため、ここが安全な複数activation入口。
- `fem/src/ElmerSolver.F90:429`: 読んだ `Mesh Mode` を以後の繰り返し制御に使う。
- `fem/src/ModelDescription.F90:279-302`: `ReloadInputFile` は既存Modelに対して `LoadInputFile` を呼ぶだけで、`FreeModel` と `LoadModel` のmesh reload lifecycleを実行しない。
- `fem/src/modules/ReloadInput.F90:36-68`: `ReloadInput` のコメント自体が「carefully」「not all information may be read」と制限を明記。Model/Mesh変更を保証できないため本ハーネスには使わない。
- `fem/src/modules/HeatSolve.F90:112`: `HeatSolver` entry。
- `fem/src/modules/HeatSolve.F90:2536`: cache helperを含む内部 `CONTAINS` 節。
- `fem/src/modules/HeatSolve.F90:3606-3727`: `Phase24ClearGenericDynKCache`、`Phase24BuildGenericDynKCache`。Build/clear/evaluate routineはHeatSolver内部手続きであり、cache配列と `CurrentModel`、`Solver` をホスト関連で参照する。独立driverからは直接呼べない。
- `fem/src/modules/HeatSolve.F90:3718` 以降: Generic評価は `ListGetReal` に現在の `Temperature` nodal valuesを渡す。各assembly評価の同値確認にもこのElmer MATC経路が必要。

## 実装を止めた具体的な障害

Harness入口に使えるElmer APIは見つかったが、既存MATC cache実装は独立した型/モジュール/APIではない。手続きは `HeatSolver` 内部で、動的cache stateを非 `SAVE` ローカル配列に持つ。よって通常のテストdriverから「本番と同じ実装」を呼ぶ手段がない。ここで似た配列処理を複製したテストを書くと模擬実装となり、要求された本番cacheの証拠にならない。

Mesh Modeで2つのfixture meshを読み込む実行試験、PID/activation idの付与、activation終了時の解放配列数記録をつなぐ独立04変更まで今回完成していない。従って同一PID・複数activationの実行ログ、build/release数、cache/Generic差分の証拠は未取得である。

## 最小設計変更案

次の独立04で本番実装とtestが同じコードを通るようにする。

1. 新規 `fem/src/modules/Phase24GenericDynKCache.F90` に `TYPE Phase24GenericDynKCache_t` を置き、`NodeSlot`, `Nodes`, `Materials`, `Start`, `Count`, `Nodal`, `Ready`, cache hit/eval/build counterを保持する。
2. `BuildGenericDynKCache(Cache, Model, Solver, ActivationId, Pid)`、`EvaluateGenericDynKCache(Cache, Model, Temperature, ...)`、`ClearGenericDynKCache(Cache, FreedArrays)`を同モジュールに実装する。既存HeatSolver内部手続きのマップ構築、 MATC `ListGetReal` 評価をこのモジュールへ移す。コピーした別実装は作らない。
3. `HeatSolver.F90` はローカル `Phase24GenericDynKCache_t` を所有し、現状のBuild/Evaluate/Clear呼び出しを同モジュールに置き換える。cache optionの既定値は `.FALSE.` のまま。通常終了時には明示clearし、Fortran自動解放も維持する。
4. 新規test-only modeは `ElmerSolver.F90` の既存 `Mesh Mode` を使う。短いSIFと `mesh` / `mesh1` fixtureをこの成果物ディレクトリに置く。Mesh AはMaterial 7の六節点prism、Mesh Bは別node connectivityまたは別material assignmentにし、両方にtemperature-only MATC式を置く。
5. 各HeatSolver activation entry/exit、cache build/clear/evalで `getpid()`、activation ID、材料数、node数、build count、解放配列数、material/node/slot対応、hit数、cached値と `ListGetReal` Generic値のmax abs/rel差分を記録する。assemblyごとに温度を変え、MATC値を再評価した回数と差分を記録する。
6. 実行スクリプトは公式26.1 clean checkoutに01→02→03→04を適用してbuildし、Mesh Modeの1 rank・短時間SIFを起動する。ログ検査は同一PIDかつactivation Aの終了行より後にB開始行があること、A/B map build各1回以上、release数、Bのmaterial/node mapがfixtureに一致、assembly evalとGeneric差分が許容誤差内であることを必須にする。

### 04 patchで変更するファイル

- `fem/src/modules/Phase24GenericDynKCache.F90` (追加)
- `fem/src/modules/HeatSolve.F90` (cache ownershipとmodule呼び出しへ置換。01/02/03 patchファイル自体は変更しない)
- `fem/src/modules/CMakeLists.txt` または実際のmodule source list (新規moduleの依存順を追加)
- test-onlyのSIF、2つのmesh fixture、runner、ログ検査器はこのディレクトリに追加

## 制約と未実施

- 01・02・03パッチファイルは変更していない。
- production DLL/SIF、物理モデル、Vector Assembly設定には触れていない。
- 75 ms本計算は実行していない。
- 同一PID実行やcache寿命のPASSは報告しない。

## 再開手順

1. Elmer公式 `release-26.1` (`5a8de867068be0568f09af40fb90fee300dfbede`) のclean checkoutを用意する。
2. 既存 `phase24-source-repro/replay_source_and_build.ps1` と `matc-cache-lifecycle/replay_01_02_03_and_build.ps1` のmanifest検証を使い01→02→03を適用する。
3. 上記のmodule抽出とtest modeを独立04として実装し、clean Release buildする。
4. fixtureでMesh Modeを2回通し、PIDが同じの生ログと数値比較JSONを保存する。
5. ログ判定スクリプトを再実行し、全条件を満たした場合だけAへ更新する。


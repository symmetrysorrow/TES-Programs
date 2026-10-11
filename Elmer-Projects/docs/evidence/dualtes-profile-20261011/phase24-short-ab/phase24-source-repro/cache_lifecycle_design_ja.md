# MATC Node Cache寿命監査と設計案

## 現行コードの確認事項

- `HeatSolver_init(Model,Solver,...)` はsolver初期化時に呼ばれるが、cache clearを保証するfinalize対称hookは確認できない。
- `HeatSolver` は `RECURSIVE SUBROUTINE` で、assemblyとNewton反復を実行する。MATC Node Cacheの配列とready/counter stateはローカル宣言に `SAVE` を付けており、呼び出し間もプロセス内に残る。
- 初期化分岐は `AllocationsDone` / `Solver % MeshChanged` とcache optionを見てclear/buildする。Model pointer・material list content・body-to-material mapping・connectivity generationの完全なidentity署名照合は確認できない。
- current sourceのdefault option `Phase24 Dynamic K Node Cache` は `.FALSE.`。今回も変更なし。

## 安全な所有方法の候補

1. **solver instance ownershipが使えるなら最優先:** cacheをSolver_tに紐付く専用 stateに置き、Model/Mesh/material generationと連携する。solver destructor/finalizerで確実に解放する。generationを持たないModel/Meshではidentityを推測せず、assembly開始時にinvalidateする。
2. **その仕組みがない場合:** cache allocatableを`HeatSolver`呼び出しactivationのlocal stateにしてautomatic deallocationを利用する。ただしassembly loopがactivation内にあること、再帰/reentrancyとOpenMPで共有しないこと、各callで再構築した場合も性能上許容することを検証する。Fortran declaration initializationは暗黙SAVEを導く場合があるため、宣言初期化を避け、procedure entryで明示初期化する。
3. `Solver % MeshChanged` がtrue、Model/material/body mapが更新、node connectivityやmaterial property list世代が不明、cache optionがOFF、または材料式validatorがfalseになったら、lookup前にそのcacheを無効化/解放する。
4. `Temperature`が非線形反復ごとに変わるため、node→property slot mapだけを保持し、k値を各assembly前に現在温度で計算する。値そのものをassembly間で再利用しない。
5. Cache OFFの経路がcacheの有無にかかわらず従来Generic経路へ進むことを同一プロセスOFF→ON→OFF試験で確認する。

## 解放・再初期化検証ゲート

ソース変更後、別々のElmer solver instanceを同一プロセスで生成/破棄し、MemorySanitizer相当またはFortran runtime checkingを使い、clear/finalizeの呼び出しとallocatable解放を記録する。Model材料の差替え、mesh node connectivity変更（同一node数でも接続だけ変更）、MeshChanged、異なるmaterialが共有nodeを持つ場合、restart初回assembly、OFF→ON→OFFを実行する。各assemblyでgeneric基準値とproperty値/RHS/行列を比較し、古いslotが一度でも読まれたら失敗とする。

この設計案はまだ実装・実行されていない。runtime finalization経路の証明が完了するまでcache既定値はOFF。

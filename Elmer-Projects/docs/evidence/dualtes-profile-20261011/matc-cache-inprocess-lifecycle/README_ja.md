# MATC Node Cache 同一プロセス寿命検証

## 結果

**PASS。公式Elmer 26.1の実`HeatSolver`とGeneric MATC Node Cacheを使い、Mesh Modeで異なる2 meshを同じOSプロセスから順に処理した。**

- Cache ON: PID `39424`、HeatSolver activation `1, 2` が`meshA`、続く`3, 4` が`meshA1`。
- Cache OFF: PID `34760`、同じmesh順と4 activationをGeneric経路で実行。
- ONの各activationでmaterial 1種・6 nodeのcache mapを1回構築し、assembly中に2回評価。2回目にはcache hitを記録。
- 各activation終了時に割り当て済み6配列をすべて実際に`DEALLOCATE`し、全配列が未割当てであることを記録。
- `meshA`と`meshA1`はnode数をともに6に保ちながらprism connectivity/node対応を変更。map fingerprintはA=`12ac3491e6d8d13969f914ae24d733b7f08b8723d48baf50a0bf0f39dde406a8`、B=`ab36c89d350ac43dea515933f0ead051e3fdd8340c4f9afe59cd208d66b33a3b`で異なり、Elmer mesh fixtureの期待node順と一致した。
- 各assemblyのcached MATC値とGeneric MATC評価の`max_abs`は0。ON/OFFの全24 nodal temperature比較も`max_abs=0`、`max_rel=0`だった。許容差は`atol=1e-13`、`rtol=1e-12`。
- Mesh A/BそれぞれのON/OFF Elmer `.result`ファイルもSHA-256が一致し、byte単位で同一だった。
- OFF→ON→OFFの同一プロセス切替は**SKIP**。Solver optionは各ModelのSolver ValueListから読むため、この同一設定SIFによるMesh Mode試験で実行時にoptionを書き換える安全なAPIを使わなかった。代替として、同じmesh/物理入力のOFFとONを別プロセスで比較した。

判定詳細は[`pass_fail_skip.json`](pass_fail_skip.json)、通常Releaseの無診断ビルド確認は[`normal_release_build_validation.json`](normal_release_build_validation.json)に保存した。実行ログは`logs/cache_on_raw.log`、`logs/cache_off_raw.log`と読みやすいUTF-8版にある。

## Elmer 26.1 Mesh Mode経路

調査した公式ソースはtag `release-26.1`、commit `5a8de867068be0568f09af40fb90fee300dfbede`。既存パッチ01→02→03適用後の隔離ソースに独立patch 04を作った。`ReloadInput`は使用していない。

- [`ElmerSolver.F90`のMesh Mode loop](https://github.com/ElmerCSC/elmerfem/blob/release-26.1/fem/src/ElmerSolver.F90)は`MeshIndex=0`から始まり、loopごとに`FreeModel(CurrentModel)`、MeshIndex更新、`LoadModel(...,MeshIndex)`を通る。これにより同一ElmerSolverプロセスでMesh A→Mesh BのModel loadとsolver activationが実際に発生する。
- `HeatSolver`はsolver activationごとに呼び出される。diagnostic buildではmodule内のSAVE activation counterとOS `getpid()`で各呼び出しを識別した。
- ElmerのMesh DB形式は`Mesh DB "fixtures/meshes" "meshA"`。Mesh Modeは連番mesh名を順に読み、`meshA`の次に`meshA1`を読む。各mesh directoryにはElmerGrid形式の`mesh.header`, `mesh.nodes`, `mesh.elements`, `mesh.boundary`を置いた。
- ElmerSolverは次のmesh番号を要求し続ける仕様なので、test-only compile macro時に2 meshを処理した後loopを抜けるガードを追加した。macroなしReleaseのpreprocessed sourceに診断markerが残らず、`HeatSolve.dll`と`ElmerSolver_mpi.exe`をビルドできたことを記録した。

今回の実行では、2 meshそれぞれで2 timestepを計算した。そのためHeatSolverは各meshにつき2 activation発生し、A activation終了後にB activationが始まった。

## Fixtureと設定

- [`meshmode_cache_on.sif`](meshmode_cache_on.sif)と[`meshmode_cache_off.sif`](meshmode_cache_off.sif)はMesh Modeを有効にした短い過渡問題。initial temperature=1、timestep=0.001、2 step、密度/比熱/ソースは1。
- Body 1はMaterial 7を参照し、temperature-only MATC conductivity `1.0 + tx*tx`を使う。Cache ON/OFF以外の物理設定を揃えた。
- Vector Assemblyは両SIFでOFF。ONのSIFのみNode CacheとGeneric verificationをONにする。既定Cache-OFF設定やProduction DLL/SIFには変更を加えていない。
- Mesh A/Bはいずれも6-node prismを1要素持つ。A connectivityは`1 2 3 4 5 6`、Bは`3 1 2 6 4 5`。material assignmentはどちらもmaterial 7。

## ログと自動判定

`PHASE24_LIFECYCLE`診断はtest compile時のみ有効で、各activationのPID/rank/id/model/material/node、mapのmaterial/element/local-node/node/slot、build/eval/hit、cache clearの割当数/実解放数、cached-vs-Generic max差を記録する。

診断ログ上、cache clearはactivation開始時にも呼ばれるが、寿命判定に使うのは終了時の`allocated_before=6`, `freed_arrays=6`, `all_unallocated=T`。終了ログと次activation開始ログの順もcheckerが検査する。BのmapはB自身のmesh connectivityに一致する。BではAの前回slot tableはすでに未割当てであり、MATC値はBのModel/温度を用いてassemblyごとに再評価される。

`check_lifecycle.py`はPID一致、activation順、Mesh Mode完了marker、fixture connectivityから計算したmap fingerprint、build/clear、cache hit、Generic差分、ON/OFF nodal temperature許容差を自動判定する。fallbackだけで全期間を計算した試験はPASSにならない。

## 再実行

必要なもの: 公式Elmer 26.1 source checkout、そこへ適用済みの既存01→02→03、CMake/Ninja、GNU Fortran/C/C++、MS-MPI、Elmerの依存ライブラリ、Python 3。実行環境に応じて compiler、source/build/install pathを指定する。

1. 01→03適用済みで未変更のElmer sourceに診断buildを作成する。スクリプトがpatch 04を検査・適用し、test-only macro付きReleaseをbuild/installする。

   ```powershell
   .\build_diagnostic_solver.ps1 `
     -ElmerSource C:\src\elmerfem `
     -BuildDir C:\build\elmer-cache-diag `
     -InstallDir C:\install\elmer-cache-diag `
     -FortranCompiler G:\Elmer-Migration\msys64\ucrt64\bin\gfortran.exe `
     -CCompiler G:\Elmer-Migration\msys64\ucrt64\bin\cc.exe `
     -CxxCompiler G:\Elmer-Migration\msys64\ucrt64\bin\c++.exe
   ```

2. `ElmerGrid 14 2 fixtures/meshes/meshA.msh -out fixtures/meshes/meshA`と同様にA/Bを変換する。保存済みのElmerGrid出力も同梱している。
3. 同一のdiagnostic `ElmerSolver_mpi.exe`をそれぞれ1 MPI rankで実行し、自動判定を行う。

   ```powershell
   .\run_meshmode.ps1 -ElmerHome C:\install\elmer-cache-diag
   ```

`run_meshmode.ps1`は同じSIF/mesh入力でONとOFFを順に実行して生ログを保存し、最後に`check_lifecycle.py`を実行する。通常Releaseの既定OFF確認には診断macroなしで`HeatSolve`と`ElmerSolver_mpi.exe`を再buildする。今回の検証ではこの2 targetが成功し、Fortran flagsは空、preprocessed sourceからdiagnostic markerが除外されることを確認した。

## Patchと制約

- 独立変更は[`04_matc_cache_meshmode_diagnostics.patch`](04_matc_cache_meshmode_diagnostics.patch)のみ。変更対象は`fem/src/ElmerSolver.F90`と`fem/src/modules/HeatSolve.F90`。
- Generic cache Build/Evaluate/Clear本体を移動・複製していない。診断は既存routineの周囲に追加した。
- パッチ01/02/03ファイル、Production DLL/SIF、TESモデル、production Vector Assembly設定は変更していない。
- 75 ms本計算、6-rank実行、長時間性能測定は行っていない。
- ソースpatch pre/post SHA-256は[`patch04_source_manifest.json`](patch04_source_manifest.json)、最終artifact一覧は[`SHA256SUMS`](SHA256SUMS)を参照。

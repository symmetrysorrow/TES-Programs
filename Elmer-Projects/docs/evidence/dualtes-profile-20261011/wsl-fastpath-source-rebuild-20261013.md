# WSL高速solver source snapshotのRelease rebuild確認

## 結果

Elmer公式base `5a8de867068be0568f09af40fb90fee300dfbede`を新規checkoutし、Git管理下の16 source filesをhash検証してoverlay後、Elmer Release buildとinstallが正常終了した。solver、`libelmersolver.so`、`libmatc.so`を生成した。

実行したscriptは[`build_elmer_tes_fastpath_wsl.sh`](../../../src/support/build_elmer_tes_fastpath_wsl.sh)。build設定はRelease、Ninja、6 jobs、GCC/GFortran 13.3.0、CMake 3.28.3、Ninja 1.11.1、OpenMPI 4.1.6、HYPRE SONAME 3.1.0。build済みbinaryのSHA-256:

- ElmerSolver_mpi: `ceca2f8484cdeeb2d446ee527fdc702f7782f38b4698dde2a166bd152bae5dca`
- libelmersolver.so: `4835b772e0fbaeb666d54323e75e5ceb271a244237b147620133250bbd5c12af`
- libmatc.so: `65420cb3a8b64b862cac3ab4bb9fdd0d0d102f818cef2fdadcf766abb2e36c8d`
- HYPRE libHYPRE.so: `48700d996a8de2df9bdc409f2c9f2f8b8b4882f301c5d6deadd13f1f35a71bf0`

生成solver hashは、75 ms実測solverの`16b512...`とは異なる。これは別build binaryであり、bit-identical buildや75 ms wall-time一致を示す結果ではない。この検証では75 ms SIFを実行していない。

## 過去solverへの来歴上の制限

過去runでsource hashを採取していなかったため、現在保存した`/root/src/elmer`の16ファイルが、当時solverをコンパイルした瞬間の全source内容と完全一致していたことは遡及証明できない。今回確立した成果は、現在のbuild source snapshotのGit保存と、新しい環境から同sourceをRelease build/installできる手順である。

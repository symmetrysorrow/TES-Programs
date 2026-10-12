# WSL/Linux DualTES高速solverのソース保存とRelease再構築

## 75 ms実測solverのsource追跡

2026-10-13のWSL計測solver SHA-256は`16b512caf1a49b3d23d4d8e23e636cd4e1df643531e81a65a3f71bd7a1d48c56`です。既存`/root/build-cpu/elmer/CMakeCache.txt`の`CMAKE_HOME_DIRECTORY`は`/root/src/elmer`を指していました。このsource treeから、solver/HYPRE/MPI/assembly/heatの変更が入った16ファイルをこのリポジトリの[`measured-wsl-source`](../patches/elmer-26.1/measured-wsl-source/)に保存しました。`SHA256SUMS`はsource snapshotの全ファイルを照合します。`.gitattributes`は改行を自動変換せず、実測source bytesを保持します。

注意点として、別の`/mnt/d/github/TES-Programs/tools/elmer-hypre/src` checkoutのHEADは`9916c3f`ですが、build directoryが参照していた`/root/src/elmer/fem/src/modules/HeatSolve.F90`はそのcommit treeと異なります。また以前作成した`01_phase24_source_from_release_26_1.patch`によるsource replayも、ここに保存した`/root/src/elmer`の内容とは一致しません。したがって、どちらかのcommit/patchを過去バイナリの厳密なprovenanceとして扱わず、build directoryが参照したsource files自体を保存しました。

過去のbuild時点で全source hashをmanifestへ出していなかったため、現存solver binaryが当時このsnapshotと同一のsourceから作られたかを遡及的に完全証明することはできません。今回のsnapshotにより、以後の再buildではsource bytesとbase revisionを照合して記録できます。

## 新しい隔離Release build

必要なものはGit、CMake、Ninja、GCC/GFortran、MPI開発環境、HYPREのMPI共有ライブラリ、MUMPS、METIS、ParMETIS、およびElmerの有効機能が要求する開発ライブラリです。build scriptはElmer公式`5a8de867068be0568f09af40fb90fee300dfbede`をclean checkoutし、上記16ファイルだけをhash検証後にoverlayし、新規build/install directoryへRelease buildします。既存source/build/installは上書きしません。

実測したElmerGPU WSL環境はUbuntu 24.04系、GCC/GFortran 13.3.0、CMake 3.28.3、Ninja 1.11.1、OpenMPI 4.1.6、HYPRE shared library SONAME 3.1.0でした。構成はCPU solverです。HYPRE prefixには`include/HYPRE.h`と`lib/libHYPRE.so`が必要です。

```bash
BUILD_JOBS=6 \
  bash Elmer-Projects/src/support/build_elmer_tes_fastpath_wsl.sh \
  "$HOME/build/elmer-tes-fastpath" \
  "$HOME/opt/elmer-tes-fastpath" \
  "/path/to/hypre-prefix"

source "$HOME/build/elmer-tes-fastpath/runtime_env.sh"
mpiexec -n 6 "$HOME/opt/elmer-tes-fastpath/bin/ElmerSolver_mpi"
```

MUMPS/METIS/ParMETISのprefixが`/usr`以外なら、`MUMPS_ROOT`、`METIS_INCLUDE_DIR`、`METIS_LIBRARY`、`PARMETIS_INCLUDE_DIR`、`PARMETIS_LIBRARY`を環境変数で指定します。成功時にElmer base、source snapshot manifest、toolchain version、solver・`libelmersolver.so`・`libmatc.so`のSHA-256を`BUILD_ROOT/build_manifest.txt`へ出力します。実行環境は`BUILD_ROOT/runtime_env.sh`から読み込みます。この環境設定は別prefixのHYPRE共有ライブラリもloaderへ登録します。

同じsourceとCMake設定から再buildできますが、compiler、MPI、HYPRE、CPUアーキテクチャが異なる環境ではバイナリhashや計算bit列が同一になるとは限りません。wall timeを比較する際は同じMPI rank/thread数、SIF、UDF、mesh、restartを使い、それらのhashと収束ログを記録します。基準計測のruntime artifactは[`wsl-dualtes-time-20261013`](evidence/dualtes-profile-20261011/wsl-dualtes-time-20261013/)にあります。

## 本計算時の設定

- Elmer Release build、6 MPI ranks、rankあたり1 thread
- Vector Assembly ON
- MATC Dynamic K Node Cache OFF（SinglePixelは既定OFF）
- DualTES/SinglePixel固有の物理モデル、timestep、BDF、収束許容値、HYPRE設定を維持
- Cache ONの数値試験は別の同一モデル短区間比較として管理

この再build scriptを使っても、他環境で75 msのwall時間や温度系列が自動的に同一になるわけではありません。source再現性と実機での数値・性能再確認は別の検証段階です。

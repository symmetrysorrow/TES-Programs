#!/usr/bin/env bash
# CPU-only HYPRE + Elmer in the ElmerGPU distro (no HIP runtime in the solve
# path).  Same sources as build_elmer_hypre_hip_elmergpu.sh; installs to
# /opt/elmer-gpu/cpu.  Run as root inside ElmerGPU.
set -euo pipefail
SRC=/root/src
BUILD=/root/build-cpu
PREFIX=/opt/elmer-gpu/cpu
export PATH=/usr/sbin:/usr/bin:/sbin:/bin
mkdir -p "$BUILD" "$PREFIX"
rsync -a --exclude .git /mnt/d/github/TES-Programs/tools/elmer-hypre/src/ "$SRC/elmer/"
cmake -S "$SRC/hypre/src" -B "$BUILD/hypre" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX/hypre" \
  -DBUILD_SHARED_LIBS=ON -DHYPRE_ENABLE_MPI=ON -DHYPRE_BUILD_EXAMPLES=OFF \
  -DHYPRE_ENABLE_HIP=OFF -DHYPRE_ENABLE_CUDA=OFF -DHYPRE_ENABLE_OPENMP=OFF
cmake --build "$BUILD/hypre" --parallel
cmake --install "$BUILD/hypre"
cmake -S "$SRC/elmer" -B "$BUILD/elmer" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX/elmer" \
  -DWITH_MPI=ON -DWITH_Mumps=ON -DWITH_Hypre=ON -DWITH_AMGX=OFF -DBUILD_TESTING=OFF \
  -DCMAKE_Fortran_FLAGS=-ffree-line-length-none \
  -DMUMPSROOT=/usr -DPARMETISROOT=/usr \
  -DMetis_INCLUDE_DIR=/usr/include -DMetis_LIBRARIES=/usr/lib/x86_64-linux-gnu/libmetis.so \
  -DParMetis_INCLUDE_DIR=/usr/include -DParMetis_LIBRARIES=/usr/lib/libparmetis.so \
  -DHypre_INCLUDE_DIR="$PREFIX/hypre/include" -DHypre_LIBRARIES="$PREFIX/hypre/lib/libHYPRE.so"
cd "$BUILD/elmer"
ninja MagnetoDynamics || true
ninja
cmake --install .
echo BUILD_OK

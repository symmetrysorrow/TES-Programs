#!/usr/bin/env bash
# Build HIP HYPRE + Elmer inside the dedicated "ElmerGPU" WSL distro
# (Ubuntu 24.04 on G:, ROCm 7.2.4 + librocdxg).  Sources come from the
# repository's tools/ tree; build trees and installs stay on the distro's
# ext4 disk.  Run as root:
#   wsl -d ElmerGPU -u root -- bash /mnt/d/github/TES-Programs/Elmer-Projects/scripts/support/build_elmer_hypre_hip_elmergpu.sh
set -euo pipefail
TOOLS=/mnt/d/github/TES-Programs/tools
PREFIX=/opt/elmer-gpu
SRC=/root/src
BUILD=/root/build
export ROCM_PATH=/opt/rocm HIP_PATH=/opt/rocm
export PATH=/opt/rocm/bin:/opt/rocm/lib/llvm/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export LD_LIBRARY_PATH=/opt/rocm/lib
mkdir -p "$SRC" "$BUILD" "$PREFIX"

# Copies keep the D: source trees read-only and make the build IO local.
rsync -a --delete --exclude .git "$TOOLS/hypre-hip-v3-1-0-src/" "$SRC/hypre/"
rsync -a --delete --exclude .git "$TOOLS/elmer-hypre/src/" "$SRC/elmer/"

cmake -S "$SRC/hypre/src" -B "$BUILD/hypre" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX/hypre" \
  -DBUILD_SHARED_LIBS=ON -DHYPRE_ENABLE_MPI=ON -DHYPRE_BUILD_EXAMPLES=OFF \
  -DHYPRE_ENABLE_UNIFIED_MEMORY=OFF -DHYPRE_ENABLE_GPU_AWARE_MPI=OFF -DHYPRE_ENABLE_UMPIRE=OFF \
  -DHYPRE_ENABLE_HIP=ON -DCMAKE_HIP_ARCHITECTURES=gfx1201 \
  -DHYPRE_ENABLE_ROCSPARSE=ON -DHYPRE_ENABLE_ROCRAND=ON -DHYPRE_ENABLE_ROCBLAS=ON \
  -DHYPRE_ENABLE_ROCSOLVER=ON -DHYPRE_ENABLE_ROCTHRUST=ON \
  -DHYPRE_ENABLE_DEVICE_MALLOC_ASYNC=ON -DHYPRE_ENABLE_CXX=ON
cmake --build "$BUILD/hypre" --parallel
cmake --install "$BUILD/hypre"

cmake -S "$SRC/elmer" -B "$BUILD/elmer" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX/elmer" \
  -DWITH_MPI=ON -DWITH_Mumps=ON -DWITH_Hypre=ON -DWITH_AMGX=OFF -DBUILD_TESTING=OFF \
  -DMUMPSROOT=/usr -DPARMETISROOT=/usr \
  -DMetis_INCLUDE_DIR=/usr/include -DMetis_LIBRARIES=/usr/lib/x86_64-linux-gnu/libmetis.so \
  -DParMetis_INCLUDE_DIR=/usr/include -DParMetis_LIBRARIES=/usr/lib/libparmetis.so \
  -DHypre_INCLUDE_DIR="$PREFIX/hypre/include" -DHypre_LIBRARIES="$PREFIX/hypre/lib/libHYPRE.so"
cmake --build "$BUILD/elmer" --parallel
cmake --install "$BUILD/elmer"
echo BUILD_OK

#!/usr/bin/env bash
# Build HIP HYPRE + Elmer inside the dedicated "ElmerGPU" WSL distro
# (Ubuntu 24.04, ROCm 7.2.4 + librocdxg). Source commits and patches are
# pinned in docs/elmergpu-sources.lock; build trees/install stay on ext4.
# Install dependencies first with install_elmergpu_build_deps.sh. Run as root:
#   wsl -d ElmerGPU -u root -- bash /mnt/d/github/TES-Programs/Elmer-Projects/src/support/build_elmer_hypre_hip_elmergpu.sh
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(realpath "$SCRIPT_DIR/../../..")"
# shellcheck disable=SC1091
source "$REPO_ROOT/Elmer-Projects/docs/elmergpu-sources.lock"
PREFIX=/opt/elmer-gpu
SRC="${ELMERGPU_SOURCE_ROOT:-/root/elmergpu-src}"
BUILD="${ELMERGPU_BUILD_ROOT:-/root/elmergpu-build}"
export ROCM_PATH=/opt/rocm HIP_PATH=/opt/rocm
export PATH=/opt/rocm/bin:/opt/rocm/lib/llvm/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export LD_LIBRARY_PATH=/opt/rocm/lib
mkdir -p "$SRC" "$BUILD" "$PREFIX"

bash "$SCRIPT_DIR/prepare_elmergpu_sources.sh" "$SRC"

if ! command -v hipcc >/dev/null || ! command -v ninja >/dev/null; then
  echo "Missing HIP or Ninja build dependencies; run install_elmergpu_build_deps.sh" >&2
  exit 2
fi
rocm_version="$(dpkg-query -W -f='${Version}' rocm-core 2>/dev/null || true)"
if [[ "$rocm_version" != "$rocm"* ]]; then
  echo "Expected ROCm/HIP $rocm; install the pinned runtime described in docs/elmergpu_build.md" >&2
  exit 2
fi

cmake -S "$SRC/hypre/src" -B "$BUILD/hypre" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX/hypre" \
  -DBUILD_SHARED_LIBS=ON -DHYPRE_ENABLE_MPI=ON -DHYPRE_BUILD_EXAMPLES=OFF \
  -DHYPRE_ENABLE_UNIFIED_MEMORY=OFF -DHYPRE_ENABLE_GPU_AWARE_MPI=OFF -DHYPRE_ENABLE_UMPIRE=OFF \
  -DHYPRE_ENABLE_HIP=ON -DCMAKE_HIP_ARCHITECTURES="$hip_architecture" \
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

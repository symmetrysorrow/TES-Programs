#!/usr/bin/env bash
# Rebuild the Elmer source snapshot used by the WSL DualTES measurements.
# Usage: build_elmer_tes_fastpath_wsl.sh BUILD_ROOT INSTALL_PREFIX HYPRE_ROOT
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 BUILD_ROOT INSTALL_PREFIX HYPRE_ROOT" >&2
  exit 2
fi

BUILD_ROOT=$(realpath -m "$1")
PREFIX=$(realpath -m "$2")
HYPRE_ROOT=$(realpath -e "$3")
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../../.." && pwd)
SOURCE_SNAPSHOT="$REPO_ROOT/Elmer-Projects/patches/elmer-26.1/measured-wsl-source"
SOURCE_HASHES="$SOURCE_SNAPSHOT/SHA256SUMS"
BASE_COMMIT=5a8de867068be0568f09af40fb90fee300dfbede

for tool in git cmake ninja sha256sum gcc g++ gfortran mpicc mpicxx mpif90; do
  command -v "$tool" >/dev/null || { echo "Missing build tool: $tool" >&2; exit 1; }
done
[[ -f "$SOURCE_HASHES" ]] || { echo "Missing source snapshot manifest: $SOURCE_HASHES" >&2; exit 1; }
[[ -f "$HYPRE_ROOT/include/HYPRE.h" ]] || { echo "HYPRE headers not found under $HYPRE_ROOT/include" >&2; exit 1; }
[[ -f "$HYPRE_ROOT/lib/libHYPRE.so" ]] || { echo "HYPRE shared library not found under $HYPRE_ROOT/lib" >&2; exit 1; }
[[ ! -e "$BUILD_ROOT" ]] || { echo "Refusing existing build root: $BUILD_ROOT" >&2; exit 1; }
[[ ! -e "$PREFIX" ]] || { echo "Refusing existing install prefix: $PREFIX" >&2; exit 1; }

mkdir -p "$BUILD_ROOT" "$PREFIX"
SRC="$BUILD_ROOT/elmer-source"
BUILD="$BUILD_ROOT/elmer-build"
git -c core.autocrlf=false clone --filter=blob:none https://github.com/ElmerCSC/elmerfem.git "$SRC"
git -C "$SRC" -c core.autocrlf=false checkout --detach "$BASE_COMMIT"
[[ "$(git -C "$SRC" rev-parse HEAD)" == "$BASE_COMMIT" ]] || { echo "Unexpected upstream source revision" >&2; exit 1; }
(cd "$SOURCE_SNAPSHOT" && sha256sum -c SHA256SUMS)
while read -r expected rel; do
  rel=${rel%$'\r'}
  [[ -n "$expected" && -n "$rel" ]] || continue
  install -D -m 0644 "$SOURCE_SNAPSHOT/$rel" "$SRC/$rel"
  actual=$(sha256sum "$SRC/$rel" | awk '{print $1}')
  [[ "$actual" == "$expected" ]] || { echo "Source snapshot copy mismatch: $rel" >&2; exit 1; }
done < "$SOURCE_HASHES"

cmake -S "$SRC" -B "$BUILD" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DELMER_SOLVER_HOME="$PREFIX/share/elmersolver" \
  -DELMER_INSTALL_LIB_DIR:STRING=lib/elmersolver \
  -DWITH_MPI=ON \
  -DWITH_Mumps=ON \
  -DWITH_Hypre=ON \
  -DWITH_AMGX=OFF \
  -DBUILD_TESTING=OFF \
  -DCMAKE_Fortran_FLAGS=-ffree-line-length-none \
  -DMUMPSROOT="${MUMPS_ROOT:-/usr}" \
  -DMetis_INCLUDE_DIR="${METIS_INCLUDE_DIR:-/usr/include}" \
  -DMetis_LIBRARIES="${METIS_LIBRARY:-/usr/lib/x86_64-linux-gnu/libmetis.so}" \
  -DParMetis_INCLUDE_DIR="${PARMETIS_INCLUDE_DIR:-/usr/include}" \
  -DParMetis_LIBRARIES="${PARMETIS_LIBRARY:-/usr/lib/libparmetis.so}" \
  -DHypre_INCLUDE_DIR="$HYPRE_ROOT/include" \
  -DHypre_LIBRARIES="$HYPRE_ROOT/lib/libHYPRE.so"
# This tree has an upstream Fortran module dependency that Ninja does not
# always order correctly on a clean parallel build. Build its producer first.
cmake --build "$BUILD" --target MagnetoDynamics --parallel "${BUILD_JOBS:-$(nproc)}"
cmake --build "$BUILD" --parallel "${BUILD_JOBS:-$(nproc)}"
cmake --install "$BUILD"

SOLVER="$PREFIX/bin/ElmerSolver_mpi"
ELMERSOLVER_LIB="$PREFIX/lib/elmersolver/libelmersolver.so"
MATC_LIB="$PREFIX/lib/elmersolver/libmatc.so"
for artifact in "$SOLVER" "$ELMERSOLVER_LIB" "$MATC_LIB"; do
  [[ -f "$artifact" ]] || { echo "Expected build artifact missing: $artifact" >&2; exit 1; }
done

RUNTIME_ENV="$BUILD_ROOT/runtime_env.sh"
cat > "$RUNTIME_ENV" <<EOF
export ELMER_HOME="$PREFIX/share/elmersolver"
export ELMER_SOLVER_HOME="$PREFIX/share/elmersolver"
export LD_LIBRARY_PATH="$HYPRE_ROOT/lib:$PREFIX/lib/elmersolver\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
EOF

cat > "$BUILD_ROOT/build_manifest.txt" <<EOF
upstream_url=https://github.com/ElmerCSC/elmerfem.git
upstream_base_commit=$BASE_COMMIT
custom_source_snapshot=$SOURCE_SNAPSHOT
source_snapshot_manifest_sha256=$(sha256sum "$SOURCE_HASHES" | awk '{print $1}')
source_snapshot_file_count=$(wc -l < "$SOURCE_HASHES")
build_type=Release
generator=Ninja
build_jobs=${BUILD_JOBS:-$(nproc)}
gcc=$(gcc -dumpfullversion -dumpversion)
gfortran=$(gfortran -dumpfullversion -dumpversion)
cmake=$(cmake --version | sed -n '1p')
ninja=$(ninja --version)
mpi=$(mpiexec --version | head -n 1)
hypre_root=$HYPRE_ROOT
hypre_lib_sha256=$(sha256sum "$HYPRE_ROOT/lib/libHYPRE.so" | awk '{print $1}')
runtime_env=$RUNTIME_ENV
solver_sha256=$(sha256sum "$SOLVER" | awk '{print $1}')
libelmersolver_sha256=$(sha256sum "$ELMERSOLVER_LIB" | awk '{print $1}')
libmatc_sha256=$(sha256sum "$MATC_LIB" | awk '{print $1}')
EOF
echo "BUILD_PASS"
echo "Solver: $SOLVER"
echo "Manifest: $BUILD_ROOT/build_manifest.txt"
echo "Runtime environment: source $RUNTIME_ENV"

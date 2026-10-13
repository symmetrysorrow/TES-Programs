#!/usr/bin/env bash
# Install CPU-only Elmer/HYPRE with the repository's pinned TES/Phase24 patches.
# Runs on Ubuntu 24.04 x86_64 (including WSL). No ROCm, HIP, CUDA or GPU is needed.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(realpath "$SCRIPT_DIR/../..")"
if [[ "$(id -u)" -ne 0 ]]; then
  echo "Run as root: sudo bash $0" >&2
  exit 2
fi
# shellcheck disable=SC1091
source /etc/os-release
if [[ "${ID:-}" != ubuntu || "${VERSION_ID:-}" != 24.04 || "$(uname -m)" != x86_64 ]]; then
  echo "CPU bootstrap supports Ubuntu 24.04 x86_64, found ${PRETTY_NAME:-unknown} ($(uname -m))" >&2
  exit 2
fi

PREFIX="${TES_CPU_PREFIX:-/opt/tes-elmer-cpu}"
SRC="${TES_CPU_SOURCE_ROOT:-/opt/tes-elmer-cpu-sources}"
BUILD="${TES_CPU_BUILD_ROOT:-/opt/tes-elmer-cpu-build}"
export DEBIAN_FRONTEND=noninteractive

apt-get update
apt-get install -y software-properties-common
add-apt-repository -y multiverse
apt-get update

# The package manifest contains validated versions; accept current Ubuntu 24.04
# revisions for installation and record their exact versions after setup.
mapfile -t packages < <(
  sed -E '/^[[:space:]]*(#|$)/d; s/=.*$//' "$PROJECT_ROOT/docs/elmergpu-ubuntu-24.04-packages.txt"
)
apt-get install -y "${packages[@]}" python3 python3-venv python3-pip libglu1-mesa libgl1

mkdir -p "$PREFIX" "$SRC" "$BUILD"
bash "$SCRIPT_DIR/prepare_elmergpu_sources.sh" "$SRC"

cmake -S "$SRC/hypre/src" -B "$BUILD/hypre" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX/hypre" \
  -DBUILD_SHARED_LIBS=ON -DHYPRE_ENABLE_MPI=ON -DHYPRE_BUILD_EXAMPLES=OFF \
  -DHYPRE_ENABLE_HIP=OFF -DHYPRE_ENABLE_CUDA=OFF -DHYPRE_ENABLE_OPENMP=OFF
cmake --build "$BUILD/hypre" --parallel
cmake --install "$BUILD/hypre"

cmake -S "$SRC/elmer" -B "$BUILD/elmer" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX/elmer" \
  -DWITH_MPI=ON -DWITH_Mumps=ON -DWITH_Hypre=ON -DWITH_AMGX=OFF \
  -DBUILD_TESTING=OFF -DCMAKE_Fortran_FLAGS=-ffree-line-length-none \
  -DMUMPSROOT=/usr -DPARMETISROOT=/usr \
  -DMetis_INCLUDE_DIR=/usr/include \
  -DMetis_LIBRARIES=/usr/lib/x86_64-linux-gnu/libmetis.so \
  -DParMetis_INCLUDE_DIR=/usr/include \
  -DParMetis_LIBRARIES=/usr/lib/x86_64-linux-gnu/libparmetis.so \
  -DHypre_INCLUDE_DIR="$PREFIX/hypre/include" \
  -DHypre_LIBRARIES="$PREFIX/hypre/lib/libHYPRE.so"
cmake --build "$BUILD/elmer" --parallel
cmake --install "$BUILD/elmer"

for executable in ElmerSolver ElmerGrid elmerf90; do
  test -x "$PREFIX/elmer/bin/$executable" || {
    echo "CPU build did not install $executable" >&2
    exit 2
  }
done

python3 -m venv "$PREFIX/venv"
"$PREFIX/venv/bin/python" -m pip install -r "$PROJECT_ROOT/requirements-cpu.txt"
mkdir -p "$PREFIX/udf"
export ELMER_HOME="$PREFIX/elmer"
export PATH="$ELMER_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/hypre/lib:$ELMER_HOME/lib/elmersolver:${LD_LIBRARY_PATH:-}"
elmerf90 "$PROJECT_ROOT/src/fortran/tes_parallel_circuit.f90" -o "$PREFIX/udf/tes_parallel_circuit.so"
test -s "$PREFIX/udf/tes_parallel_circuit.so"

{
  echo "CPU installed at $(date -u +%FT%TZ)"
  echo "Elmer base $(git -C "$SRC/elmer" rev-parse HEAD)"
  echo "HYPRE base $(git -C "$SRC/hypre" rev-parse HEAD)"
  echo "Elmer patch $(sha256sum "$PROJECT_ROOT/docs/elmer_phase24_custom.patch" | cut -d ' ' -f1)"
  echo "HYPRE patch $(sha256sum "$PROJECT_ROOT/docs/hypre_hip_rocm714_compat.patch" | cut -d ' ' -f1)"
  dpkg-query -W -f='${binary:Package}=${Version}\n' "${packages[@]}"
} > "$PREFIX/build-provenance.txt"
"$PREFIX/venv/bin/python" -m pip freeze > "$PREFIX/python-packages.txt"
echo "CPU_BOOTSTRAP_OK: $PREFIX"

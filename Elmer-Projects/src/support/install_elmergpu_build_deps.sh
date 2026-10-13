#!/usr/bin/env bash
# Install Ubuntu packages required to build the pinned Elmer + HYPRE sources.
set -euo pipefail

if [[ "$(id -u)" -ne 0 ]]; then
  echo "Run as root inside Ubuntu 24.04 WSL" >&2
  exit 2
fi
. /etc/os-release
if [[ "${ID:-}" != ubuntu || "${VERSION_ID:-}" != 24.04 ]]; then
  echo "This dependency list targets Ubuntu 24.04; found ${PRETTY_NAME:-unknown}" >&2
  exit 2
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_FILE="$SCRIPT_DIR/../../docs/elmergpu-ubuntu-24.04-packages.txt"
ROCM_PACKAGE_FILE="$SCRIPT_DIR/../../docs/elmergpu-rocm-7.2.4-packages.txt"
apt-get update
mapfile -t packages < <(cat "$PACKAGE_FILE" "$ROCM_PACKAGE_FILE" | grep -Ev '^[[:space:]]*(#|$)')
apt-get install -y "${packages[@]}"

LOCK="$SCRIPT_DIR/../../docs/elmergpu-sources.lock"
# shellcheck disable=SC1090
source "$LOCK"
rocm_version="$(dpkg-query -W -f='${Version}' rocm-core 2>/dev/null || true)"
if ! command -v hipcc >/dev/null || [[ "$rocm_version" != "$rocm"* ]]; then
  cat >&2 <<EOF
ROCm/HIP $rocm was not found. Install the AMD ROCm $rocm WSL runtime and
librocdxg component as documented in docs/elmergpu_build.md, then rerun this
script. ROCm packages are pinned in the separate AMD package list.
EOF
  exit 2
fi

#!/usr/bin/env bash
set -euo pipefail
export MSYSTEM=UCRT64
export ELMER_HOME=/d/Github/TES-Programs/tools/elmer-hypre/install-outer-capture-v2
export ELMER_LIB="$ELMER_HOME/share/elmersolver/lib"
export PATH="$ELMER_HOME/bin:/ucrt64/bin:/usr/bin:/c/Program Files/CMake/bin:/c/Windows/System32"
cd /d/Github/TES-Programs/Elmer-Projects
exec "$ELMER_HOME/bin/ElmerSolver_mpi.exe" artifacts/phase24_outer_capture_build_audit/run/case_phase24_outer_capture.sif

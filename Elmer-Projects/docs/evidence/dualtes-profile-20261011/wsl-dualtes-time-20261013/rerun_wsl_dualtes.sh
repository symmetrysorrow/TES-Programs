#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
SRC=/mnt/d/Github/TES-Programs/Elmer-Projects/workspace/work/meshes/mesh_toml_95333801e1a3
SOLVER=/opt/elmer-gpu/cpu/elmer/bin/ElmerSolver_mpi
EXPECTED_SOLVER=16b512caf1a49b3d23d4d8e23e636cd4e1df643531e81a65a3f71bd7a1d48c56
RUN_TAG="${RUN_TAG:-rerun-$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_DIR="$ROOT/runs/$RUN_TAG"
RUNTIME="$ROOT/runtime_$RUN_TAG"
MESH="$RUNTIME/workspace/work/meshes/mesh_toml_95333801e1a3"

[[ ! -e "$RUN_DIR" && ! -e "$RUNTIME" ]] || { echo "Refusing to overwrite $RUN_TAG"; exit 2; }
ACTUAL_SOLVER=$(sha256sum "$SOLVER" | awk '{print $1}')
[[ "$ACTUAL_SOLVER" == "$EXPECTED_SOLVER" ]] || { echo "Unexpected solver hash $ACTUAL_SOLVER"; exit 2; }

export ELMER_HOME=/opt/elmer-gpu/cpu/elmer
export LD_LIBRARY_PATH=/opt/elmer-gpu/udf:/opt/elmer-gpu/cpu/hypre/lib:/opt/elmer-gpu/cpu/elmer/lib/elmersolver
export PATH=/opt/elmer-gpu/cpu/elmer/bin:$PATH
export OMP_NUM_THREADS=1
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mkdir -p "$MESH" "$RUNTIME/workspace/generated/cases" "$RUN_DIR/steady" "$RUN_DIR/transient"
ln -s "$MESH" "$RUNTIME/mesh_toml_95333801e1a3"
cp -a "$SRC"/mesh.* "$SRC/entities.sif" "$MESH/"
cp -a "$SRC/partitioning.6" "$MESH/"
cp -a "$SRC/case_toml_2cc534365acf_steady_1.state" "$SRC/case_toml_2cc534365acf_steady_2.state" "$MESH/"
cp "$ROOT/inputs/dualtes_wsl_steady_20261013_newton.sif" "$RUNTIME/workspace/generated/cases/"
cp "$ROOT/inputs/dualtes_wsl_transient_20261013_newton.sif" "$RUNTIME/workspace/generated/cases/"
(cd "$SRC" && sha256sum mesh.* entities.sif case_toml_2cc534365acf_steady_1.state case_toml_2cc534365acf_steady_2.state > "$RUN_DIR/mesh_inputs.sha256")
find "$SRC/partitioning.6" -type f -print0 | sort -z | xargs -0 sha256sum >> "$RUN_DIR/mesh_inputs.sha256"
printf 'solver_sha256=%s\n' "$ACTUAL_SOLVER" > "$RUN_DIR/build_manifest.txt"
sha256sum /opt/elmer-gpu/cpu/elmer/lib/elmersolver/libelmersolver.so /opt/elmer-gpu/cpu/elmer/lib/elmersolver/libmatc.so /opt/elmer-gpu/udf/tes_parallel_circuit.so /opt/elmer-gpu/udf/tes_transient_heat_source_t0.so >> "$RUN_DIR/build_manifest.txt"

run_case() {
  local label=$1 sif=$2 log=$3 record=$4
  local start_utc start_ns end_utc end_ns elapsed rc
  start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ); start_ns=$(date +%s%N)
  echo "START $label $start_utc"
  set +e
  (cd "$RUNTIME" && mpiexec -n 6 "$SOLVER" "workspace/generated/cases/$sif") > "$log" 2>&1
  rc=$?
  set -e
  end_ns=$(date +%s%N); end_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  elapsed=$(awk -v a="$start_ns" -v b="$end_ns" 'BEGIN {printf "%.3f", (b-a)/1000000000}')
  RUN_LABEL="$label" RUN_SIF="$ROOT/inputs/$sif" RUN_LOG="$log" RUN_RECORD="$record" RUN_START="$start_utc" RUN_END="$end_utc" RUN_ELAPSED="$elapsed" RUN_RC="$rc" RUN_HASH="$ACTUAL_SOLVER" python3 - <<'PY'
import hashlib,json,os,re
from pathlib import Path
log=Path(os.environ['RUN_LOG']); text=log.read_text(encoding='utf-8',errors='replace')
total=re.findall(r'SOLVER TOTAL TIME\(CPU,REAL\):\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)',text)
steps=re.findall(r'MAIN:\s*Time:\s*(\d+)\s*/\s*(\d+):\s*([0-9.eE+-]+)',text)
sif=Path(os.environ['RUN_SIF'])
r={'case':os.environ['RUN_LABEL'],'start_utc':os.environ['RUN_START'],'end_utc':os.environ['RUN_END'],'wall_seconds':float(os.environ['RUN_ELAPSED']),'exit_code':int(os.environ['RUN_RC']),'normal_completion':int(os.environ['RUN_RC'])==0 and 'ALL DONE' in text,'last_logged_step':int(steps[-1][0]) if steps else None,'step_count':int(steps[-1][1]) if steps else None,'last_logged_time_s':float(steps[-1][2]) if steps else None,'solver_cpu_real_seconds':list(map(float,total[-1])) if total else None,'mpi_ranks':6,'omp_threads_per_rank':1,'solver':'/opt/elmer-gpu/cpu/elmer/bin/ElmerSolver_mpi','solver_sha256':os.environ['RUN_HASH'],'sif':str(sif),'sif_sha256':hashlib.sha256(sif.read_bytes()).hexdigest(),'log':os.environ['RUN_LOG']}
Path(os.environ['RUN_RECORD']).write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
print(json.dumps(r,indent=2))
PY
  if [[ "$rc" -ne 0 ]] || ! grep -q 'ALL DONE' "$log"; then echo "FAILED $label; see $log"; return 1; fi
}

run_case dualtes_steady dualtes_wsl_steady_20261013_newton.sif "$RUN_DIR/steady/solver.log" "$RUN_DIR/steady/run_record.json"
cp "$MESH/dualtes_wsl_steady_20261013_newton.result"* "$RUN_DIR/steady/"
cp "$RUNTIME"/dualtes_wsl_steady_20261013_newton*.csv "$RUN_DIR/steady/" 2>/dev/null || true
run_case dualtes_transient_75ms dualtes_wsl_transient_20261013_newton.sif "$RUN_DIR/transient/solver.log" "$RUN_DIR/transient/run_record.json"
cp "$RUNTIME"/dualtes_wsl_transient_20261013_newton*.csv "$RUN_DIR/transient/" 2>/dev/null || true
echo "ALL RUNS COMPLETE: $RUN_DIR"

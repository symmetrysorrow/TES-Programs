#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/d/Github/TES-Programs-dualtes-profile-evidence-20261011/Elmer-Projects/docs/evidence/dualtes-profile-20261011/wsl-dualtes-time-20261013
SRC=/mnt/d/Github/TES-Programs/Elmer-Projects/workspace/work/meshes/mesh_toml_95333801e1a3
RUNTIME=$ROOT/runtime_attempt3
MESH=$RUNTIME/workspace/work/meshes/mesh_toml_95333801e1a3
SOLVER=/opt/elmer-gpu/cpu/elmer/bin/ElmerSolver_mpi
EXPECTED_SOLVER=16b512caf1a49b3d23d4d8e23e636cd4e1df643531e81a65a3f71bd7a1d48c56
ACTUAL_SOLVER=$(sha256sum "$SOLVER" | awk '{print $1}')
[[ "$ACTUAL_SOLVER" == "$EXPECTED_SOLVER" ]] || { echo "Unexpected solver hash $ACTUAL_SOLVER"; exit 2; }
[[ ! -e "$RUNTIME" ]] || { echo "Refusing to overwrite existing runtime: $RUNTIME"; exit 2; }
export ELMER_HOME=/opt/elmer-gpu/cpu/elmer
export LD_LIBRARY_PATH=/opt/elmer-gpu/udf:/opt/elmer-gpu/cpu/hypre/lib:/opt/elmer-gpu/cpu/elmer/lib/elmersolver
export PATH=/opt/elmer-gpu/cpu/elmer/bin:$PATH
export OMP_NUM_THREADS=1
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
mkdir -p "$MESH" "$RUNTIME/workspace/generated/cases" "$ROOT/runs/attempt3/steady" "$ROOT/runs/attempt3/transient"
cp -a "$SRC"/mesh.* "$SRC/entities.sif" "$MESH/"
cp -a "$SRC/partitioning.6" "$MESH/"
cp -a "$SRC/case_toml_2cc534365acf_steady_1.state" "$SRC/case_toml_2cc534365acf_steady_2.state" "$MESH/"
cp "$ROOT/inputs/dualtes_wsl_steady_20261013.sif" "$RUNTIME/workspace/generated/cases/"
cp "$ROOT/inputs/dualtes_wsl_transient_20261013.sif" "$RUNTIME/workspace/generated/cases/"
(cd "$SRC" && sha256sum mesh.* entities.sif case_toml_2cc534365acf_steady_1.state case_toml_2cc534365acf_steady_2.state > "$ROOT/mesh_inputs.sha256")
find "$SRC/partitioning.6" -type f -print0 | sort -z | xargs -0 sha256sum >> "$ROOT/mesh_inputs.sha256"
printf 'solver_sha256=%s\n' "$ACTUAL_SOLVER" > "$ROOT/build_manifest.txt"
sha256sum /opt/elmer-gpu/cpu/elmer/lib/elmersolver/libelmersolver.so /opt/elmer-gpu/cpu/elmer/lib/elmersolver/libmatc.so /opt/elmer-gpu/udf/tes_parallel_circuit.so /opt/elmer-gpu/udf/tes_transient_heat_source_t0.so >> "$ROOT/build_manifest.txt"
cd "$RUNTIME"
run_case() {
  local label=$1 sif=$2 log=$3 record=$4
  local start_utc start_ns end_utc end_ns elapsed rc
  start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ); start_ns=$(date +%s%N)
  echo "START $label $start_utc"
  set +e
  mpiexec -n 6 "$SOLVER" "workspace/generated/cases/$sif" > "$log" 2>&1
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
run_case dualtes_steady dualtes_wsl_steady_20261013.sif "$ROOT/runs/attempt3/steady/solver.log" "$ROOT/runs/attempt3/steady/run_record.json"
cp "$MESH/dualtes_wsl_steady_20261013.result"* "$ROOT/runs/attempt3/steady/" 2>/dev/null || true
cp "$RUNTIME"/dualtes_wsl_steady_20261013*.csv "$ROOT/runs/attempt3/steady/" 2>/dev/null || true
run_case dualtes_transient_75ms dualtes_wsl_transient_20261013.sif "$ROOT/runs/attempt3/transient/solver.log" "$ROOT/runs/attempt3/transient/run_record.json"
cp "$RUNTIME"/dualtes_wsl_transient_20261013*.csv "$ROOT/runs/attempt3/transient/" 2>/dev/null || true
printf 'ALL RUNS COMPLETE\n'



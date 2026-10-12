#!/usr/bin/env bash
set -uo pipefail
RUN_ROOT="/mnt/d/Github/TES-Programs-dualtes-profile-evidence-20261011/Elmer-Projects/outputs/runs/wsl-replay2/singlepixel/runtime"
EVIDENCE="/mnt/d/Github/TES-Programs-dualtes-profile-evidence-20261011/Elmer-Projects/outputs/runs/wsl-replay2/singlepixel"
export ELMER_HOME=/opt/elmer-gpu/cpu/elmer
export LD_LIBRARY_PATH=/opt/elmer-gpu/udf:/opt/elmer-gpu/cpu/hypre/lib:/opt/elmer-gpu/cpu/elmer/lib/elmersolver
export PATH=/opt/elmer-gpu/cpu/elmer/bin:$PATH
export OMP_NUM_THREADS=1
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
cd "$RUN_ROOT"
if [[ -L mesh_hybrid_fullconf_h8 ]]; then unlink mesh_hybrid_fullconf_h8; fi
mkdir -p mesh_hybrid_fullconf_h8
start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
start_epoch=$(date +%s)
set +e
timeout --signal=INT --kill-after=20s 600s mpiexec -n 6 ElmerSolver_mpi workspace/generated/cases/case_phase24_g45_pre_75ms_bdf2h15_hypre_gpu_wsl_replay.sif > "$EVIDENCE/solver.log" 2>&1
rc=$?
set -e
end_epoch=$(date +%s)
end_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
elapsed=$((end_epoch-start_epoch))
solver_hash=$(sha256sum /opt/elmer-gpu/cpu/elmer/bin/ElmerSolver_mpi | awk '{print $1}')
sif_hash=$(sha256sum workspace/generated/cases/case_phase24_g45_pre_75ms_bdf2h15_hypre_gpu_wsl_replay.sif | awk '{print $1}')
status=INCOMPLETE
if [[ "$rc" -eq 0 ]] && grep -q 'MAIN: \*\*\* Elmer Solver: ALL DONE \*\*\*' "$EVIDENCE/solver.log"; then status=COMPLETED; elif [[ "$rc" -eq 124 ]] || [[ "$rc" -eq 137 ]]; then status=STOPPED_AT_10_MINUTES; fi
printf '{\n  "case": "singlepixel_75ms_wsl_replay",\n  "status": "%s",\n  "exit_code": %s,\n  "wall_seconds": %s,\n  "start_utc": "%s",\n  "end_utc": "%s",\n  "ranks": 6,\n  "omp_threads_per_rank": 1,\n  "solver": "/opt/elmer-gpu/cpu/elmer/bin/ElmerSolver_mpi",\n  "solver_sha256": "%s",\n  "run_sif": "runtime/workspace/generated/cases/case_phase24_g45_pre_75ms_bdf2h15_hypre_gpu_wsl_replay.sif",\n  "run_sif_sha256": "%s",\n  "time_limit_seconds": 600,\n  "log": "solver.log"\n}\n' "$status" "$rc" "$elapsed" "$start_utc" "$end_utc" "$solver_hash" "$sif_hash" > "$EVIDENCE/run_record.json"
exit 0
"""Run Phase24 steady/transient cases on the RX 9070 XT via the ElmerGPU WSL distro.

The case is generated and synced on Windows with the existing Gate 3 / Gate 4
runners, so the physics and SIF generation are identical to the CPU route.
Only ElmerSolver runs inside WSL (ElmerGPU: Ubuntu 24.04 on G:, ROCm 7.2.4 +
librocdxg, HIP HYPRE and Elmer under /opt/elmer-gpu; see
scripts/support/build_elmer_hypre_hip_elmergpu.sh).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1].parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DISTRO = "ElmerGPU"
PREFIX = "/opt/elmer-gpu"


def wsl_path(path: Path) -> str:
    resolved = path.resolve()
    return f"/mnt/{resolved.drive[0].lower()}{resolved.as_posix()[2:]}"


def sync(project: Path) -> None:
    result = subprocess.run([sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project)],
                            cwd=ROOT, capture_output=True, text=True)
    (project.parent / "sync.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"sync failed; see {project.parent / 'sync.log'}")


def run_wsl(project: Path, case: str, mesh: str, threads: int, elmer: str = "elmer", extra_env: list[str] | None = None, mpi: int = 1) -> dict:
    # Linux Elmer resolves "<mesh>/../work/..." textually; keep the anchor.
    (ROOT / mesh).mkdir(parents=True, exist_ok=True)
    bash = f"""set -euo pipefail
export ROCM_PATH=/opt/rocm HIP_PATH=/opt/rocm HIP_VISIBLE_DEVICES=0
export HSA_ENABLE_DXG_DETECTION=1
export PATH=/opt/rocm/bin:/usr/sbin:/usr/bin:/sbin:/bin
export ELMER_HOME={PREFIX}/{elmer}
export LD_LIBRARY_PATH={PREFIX}/udf:{PREFIX}/hypre/lib:{PREFIX}/{elmer}/lib/elmersolver:/opt/rocm/lib:/usr/lib/wsl/lib
export OMP_NUM_THREADS={threads}
export PHASE24_MEMORY_TRACE=1
# hwloc (inside Open MPI, also for 1 rank) would otherwise load the AMD
# OpenCL ICD and the HSA runtime, whose polling threads steal the cores.
export HWLOC_COMPONENTS=-opencl,-rsmi,-levelzero OCL_ICD_VENDORS=/opt/elmer-gpu/no-ocl-icd
{"".join(f"export {kv}" + chr(10) for kv in (extra_env or []))}cd '{wsl_path(ROOT)}'
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
python3 run.py '{case}' --project '{wsl_path(project)}' --skip-sync --mpi-procs {mpi} --elmer-solver {PREFIX}/{elmer}/bin/{'ElmerSolver_mpi' if mpi > 1 else 'ElmerSolver'} --runtime-bin ''
"""
    log = project.parent / f"{case}.gpu_launcher.log"
    start = time.monotonic()
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(["wsl.exe", "-d", DISTRO, "-u", "root", "--exec", "bash", "-lc", bash],
                                cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    return {"case": case, "exit_code": result.returncode, "elapsed_s": time.monotonic() - start, "launcher_log": str(log)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("steady", "transient"))
    parser.add_argument("--project", type=Path, default=ROOT / "artifacts/phase24_conformal_hybrid/project.json")
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--base-case", default="case_conformal_hybrid_steady_base")
    parser.add_argument("--window", default="100us")
    parser.add_argument("--reference-case", default=None)
    parser.add_argument("--linear-system", default="iterative_hypre_pcg_boomeramg_gpu")
    parser.add_argument("--linear-tolerance", type=float, default=1e-12)
    parser.add_argument("--case-option", action="append", default=[], metavar="KEY=JSON")
    parser.add_argument("--solver-option", action="append", default=[], metavar="KEY=JSON",
                        help="merged into the generated case's solver block")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--mpi", type=int, default=1, help="MPI ranks (needs mesh/partitioning.N)")
    parser.add_argument("--env", action="append", default=[], metavar="KEY=VALUE",
                        help="extra environment variable exported inside WSL")
    parser.add_argument("--elmer", default="elmer", help="install dir under /opt/elmer-gpu (e.g. elmer-dev)")
    parser.add_argument("--keep-result", action="store_true",
                        help="transient: also write the per-step .result field history (several GB)")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    options = {k: json.loads(v) for k, _, v in (item.partition("=") for item in args.case_option)}
    solver_options = {k: json.loads(v) for k, _, v in (item.partition("=") for item in args.solver_option)}

    if args.mode == "steady":
        import scripts.support.run_phase24_gate3_hypre_steady as g3
        g3.APPLY_MORTAR_BCS = False
        g3.configure_mesh(source_project=args.project, mesh=args.mesh, base_case=args.base_case,
                          reference_current_uA=None, comsol_current_uA=None)
        g3.configure_variant(args.tag)
        g3.LINEAR_SYSTEM = args.linear_system
        g3.LINEAR_TOLERANCE = args.linear_tolerance
        g3.BOOMER_AMG_STRONG_THRESHOLD = 0.5
        project = g3.build_project()
        data = json.loads(project.read_text(encoding="utf-8"))
        data["cases"][g3.CASE_NAME].update(options)
        project.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        case = g3.CASE_NAME
    else:
        import scripts.support.run_phase24_gate4_5_nomortar as g45
        g45.configure_mesh(base_project=args.project, mesh=args.mesh, mesh_tag=args.tag,
                           reference_case=args.reference_case, reference_current_uA=None)
        g45.APPLY_MORTAR_BCS = False
        g45.CONFORMAL_TES_STACK = True
        g45.HYPRE_SYSTEM = args.linear_system
        g45.HYPRE_TOLERANCE = args.linear_tolerance
        g45.RUN_TAG = "gpu"
        g45.EXTRA_CASE_OPTIONS.update(options)
        project, case = g45.prepare_project(args.window, "hypre")
        if not args.keep_result:
            # The per-step field history is several GB per 75 ms run and is
            # not needed for the series comparison.
            data = json.loads(project.read_text(encoding="utf-8"))
            data["cases"][case]["output_result"] = False
            project.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if solver_options:
        data = json.loads(project.read_text(encoding="utf-8"))
        data["cases"][case].setdefault("solver", {}).update(solver_options)
        project.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sync(project)
    print(json.dumps({"project": str(project), "case": case}))
    if args.dry_run:
        return 0
    record = run_wsl(project, case, args.mesh, args.threads, args.elmer, args.env, args.mpi)
    print(json.dumps(record))
    return record["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())

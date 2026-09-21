"""Run the exact production-v2 short backend parity test.

This is a diagnostic runner, not a production qualification.  It creates two
distinct five-stage hybrid-grid cases (about 0.9 us after the pulse):

* direct MUMPS
* native HYPRE/BoomerAMG

The mesh, restart, pulse, BDF order, nonlinear tolerances, and UDF source are
otherwise generated from the same production-v2 case definition.  It then
compares both traces with COMSOL and compares HYPRE directly with MUMPS.

Run from the repository root::

    python scripts/analysis/run_production_backend_parity_short.py

Use ``--dry-run`` to print the commands without launching Elmer.
"""
from __future__ import annotations

import argparse
import copy
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "elmer_project_singlepixel_prod_v2_original_timegrid.json"
COMSOL = ROOT / "docs" / "Single-Pixel.txt"
OUT = ROOT / "artifacts" / "production_backend_parity_short"
CASE_MUMPS = "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid_cpu_smoke_5step"
CASE_HYPRE = "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid_amgx_smoke_5step"
PARITY_STAGES = [
    ["18[us]", 1],
    ["1[us]", 2],
    ["1[ns]", 1],
    ["10[ns]", 10],
    ["100[ns]", 9],
]


def run(command: list[str], log: Path, dry_run: bool) -> dict:
    record = {
        "command": command,
        "log": str(log),
        "started": datetime.now(timezone.utc).isoformat(),
    }
    if dry_run:
        record["exit_code"] = None
        record["finished"] = record["started"]
        return record
    start = time.monotonic()
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    record.update(
        {
            "exit_code": result.returncode,
            "finished": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic() - start,
        }
    )
    return record


def prepare_project(
    backend: str,
    case_name: str,
    hypre_tolerance: float = 5e-7,
    hypre_max_iterations: int = 2000,
) -> None:
    """Create one short case with a stage-truncated, not step-truncated grid."""
    sys.path.insert(0, str(ROOT))
    from scripts.prep import run_singlepixel_prod_v2_original_timegrid as prep

    full_case = prep.make_project(
        prep.BASE_MESH,
        backend,
        False,
        True,
        "hybrid",
        None,
        "amgx",
    )
    project = json.loads(PROJECT.read_text(encoding="utf-8"))
    candidate = copy.deepcopy(project["cases"].pop(full_case))
    candidate["timesteps"] = copy.deepcopy(PARITY_STAGES)
    candidate["output_intervals"] = [999999] * len(PARITY_STAGES)
    candidate["output_intervals"][-1] = 1
    candidate["series_file"] = f"{case_name}_series.csv"
    candidate["iteration_series_file"] = f"{case_name}_iterations.csv"
    candidate["output_file_path"] = (
        f"../work/meshes/{prep.BASE_MESH}/{case_name}.result"
    )
    candidate["state_file"] = f"work/meshes/{prep.BASE_MESH}/{case_name}.state"
    candidate["solver"] = copy.deepcopy(candidate["solver"])
    candidate["solver"]["linear_system"] = backend
    if backend == "iterative_hypre_flexgmres_boomeramg":
        # The native HYPRE build has a measured residual floor near 1e-7 on
        # this matrix.  5e-7 is the existing production acceptance tolerance;
        # 1e-10 would test an already-known impossible convergence target.
        candidate["solver"]["linear_system_max_iterations"] = hypre_max_iterations
        candidate["solver"]["linear_system_convergence_tolerance"] = hypre_tolerance
    candidate["solver_comment"] = (
        f"Production-v2 backend parity diagnostic: {backend}, 0.9 us window"
    )
    project["cases"][case_name] = candidate
    PROJECT.write_text(
        json.dumps(project, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    prep.write_known_steady_state(ROOT / candidate["state_file"])


def solver_command(case: str, args: argparse.Namespace) -> list[str]:
    return [
        sys.executable,
        str(ROOT / "run.py"),
        case,
        "--project",
        str(PROJECT),
        "--mpi-procs",
        "1",
        "--elmer-solver",
        str(args.solver),
        "--runtime-bin",
        str(args.runtime_bin),
        "--toolchain-bin",
        str(args.toolchain_bin),
    ]


def series(case: str) -> Path:
    return ROOT / "results" / case / f"{case}_series.csv"


def compare_comsol(case: str, label: str, out: Path, dry_run: bool) -> dict:
    command = [
        sys.executable,
        str(ROOT / "scripts" / "analysis" / "compare_singlepixel_amgx_comsol.py"),
        "--elmer",
        str(series(case)),
        "--comsol",
        str(COMSOL),
        "--out",
        str(out),
        "--end-us",
        "0.9",
        "--solver-label",
        label,
    ]
    return run(command, out / "comparison.log", dry_run)


def compare_backends(out: Path, dry_run: bool) -> dict:
    command = [
        sys.executable,
        str(ROOT / "scripts" / "analysis" / "phase24_difference_evidence.py"),
        "--candidate",
        str(series(CASE_HYPRE)),
        "--reference",
        str(series(CASE_MUMPS)),
        "--candidate-label",
        "production-v2 native HYPRE",
        "--reference-label",
        "production-v2 direct MUMPS",
        "--checkpoint-us",
        "0.1",
        "0.2",
        "0.3",
        "0.4",
        "0.5",
        "0.6",
        "0.7",
        "0.8",
        "0.9",
        "--out",
        str(out),
    ]
    return run(command, out / "comparison.log", dry_run)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--solver",
        type=Path,
        default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin\ElmerSolver_mpi.exe"),
    )
    parser.add_argument(
        "--runtime-bin",
        type=Path,
        default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\lib\elmersolver"),
    )
    parser.add_argument(
        "--toolchain-bin",
        type=Path,
        default=Path(r"C:\msys64\ucrt64\bin"),
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--hypre-tolerance",
        type=float,
        default=5e-7,
        help="HYPRE convergence tolerance for the backend-only follow-up",
    )
    parser.add_argument(
        "--hypre-max-iterations",
        type=int,
        default=2000,
        help="HYPRE maximum iterations for the backend-only follow-up",
    )
    parser.add_argument(
        "--hypre-only",
        action="store_true",
        help="reuse the existing MUMPS result and execute only the HYPRE case",
    )
    args = parser.parse_args()
    for path, label in ((args.solver, "solver"), (args.runtime_bin, "runtime-bin"), (args.toolchain_bin, "toolchain-bin")):
        if not args.dry_run and not path.exists():
            raise SystemExit(f"{label} not found: {path}")

    OUT.mkdir(parents=True, exist_ok=True)
    summary: dict = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(PROJECT),
        "window": "0..0.9 us after the 20.02 ms pulse",
        "conditions": {
            "mesh": "mesh_singlepixel_prod_v2",
            "restart": "validated serial production-v2 steady result",
        "bdf_order": 1,
        "time_grid_stages": PARITY_STAGES,
            "nonlinear_max_iterations": 120,
            "nonlinear_convergence_tolerance": 1e-8,
            "hypre_tolerance": args.hypre_tolerance,
            "hypre_max_iterations": args.hypre_max_iterations,
        },
        "runs": {},
        "comparisons": {},
    }

    cases = (
        (("iterative_hypre_flexgmres_boomeramg", "hypre", CASE_HYPRE),)
        if args.hypre_only
        else (
            ("mumps", "mumps", CASE_MUMPS),
            ("iterative_hypre_flexgmres_boomeramg", "hypre", CASE_HYPRE),
        )
    )
    for backend, label, case in cases:
        prepare_project(
            backend,
            case,
            args.hypre_tolerance,
            args.hypre_max_iterations,
        )
        execution = run(solver_command(case, args), OUT / f"run_{label}.log", args.dry_run)
        summary["runs"][case] = execution

    if args.dry_run or series(CASE_MUMPS).is_file():
        summary["comparisons"]["mumps_vs_comsol"] = compare_comsol(
            CASE_MUMPS, "production-v2 direct MUMPS", OUT / "mumps_vs_comsol", args.dry_run
        )
    if args.dry_run or series(CASE_HYPRE).is_file():
        summary["comparisons"]["hypre_vs_comsol"] = compare_comsol(
            CASE_HYPRE, "production-v2 native HYPRE", OUT / "hypre_vs_comsol", args.dry_run
        )
    if (args.dry_run or series(CASE_HYPRE).is_file()) and (args.dry_run or series(CASE_MUMPS).is_file()):
        summary["comparisons"]["hypre_vs_mumps"] = compare_backends(
            OUT / "hypre_vs_mumps", args.dry_run
        )
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"summary={OUT / 'summary.json'}")
    print(f"mumps_vs_comsol={OUT / 'mumps_vs_comsol' / 'summary.md'}")
    print(f"hypre_vs_comsol={OUT / 'hypre_vs_comsol' / 'summary.md'}")
    print(f"hypre_vs_mumps={OUT / 'hypre_vs_mumps' / 'summary.md'}")
    return 0 if all(
        record.get("exit_code") in (0, None)
        for record in summary["runs"].values()
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())

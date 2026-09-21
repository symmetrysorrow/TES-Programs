"""Test a fresh outer FlexGMRES setup after Phase24 matrix refreshes.

This is a short, non-production diagnostic.  It clones the existing 1-us
Phase24 HYPRE 5e-7 case into a new case name, runs it with the separately
installed rebuilt solver, and enables ``PHASE24_HYPRE_FORCE_KRYLOV_SETUP``.
The existing result directory is never reused.
"""
from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PROJECT = ROOT / "artifacts/phase24_hypre_tolerance_diagnostic_1us/phase24_hypre_tolerance_diagnostic.json"
ARTIFACT_DIR = ROOT / "artifacts/phase24_hypre_refresh_setup_diagnostic_1us"
PROJECT = ARTIFACT_DIR / "phase24_hypre_refresh_setup_diagnostic.json"
CASE = "case_phase24_short_hypre_refresh_setup_1us"
RESULT_DIR = ROOT / "results" / CASE
LAUNCHER_LOG = RESULT_DIR / "refresh_setup_launcher.log"
SOLVER_LOG = RESULT_DIR / "solver.log"


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def make_project() -> None:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    base = copy.deepcopy(project["cases"]["case_phase24_short_hypre_tol5e7_1us"])
    base["series_file"] = f"{CASE}_series.csv"
    base["iteration_series_file"] = f"{CASE}_iterations.csv"
    base["solver_comment"] = "Diagnostic: force fresh outer FlexGMRES setup after matrix refresh"
    base["phase24_smoke"] = {
        "purpose": "short matrix-refresh correctness diagnostic",
        "variant": "force fresh outer Krylov setup",
        "no_matrix_dump": True,
        "no_vtu": True,
    }
    base["output_result_path"] = None
    base["output_file_path"] = None
    project["cases"] = {CASE: base}
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11-refresh-test\bin\ElmerSolver_mpi.exe"))
    parser.add_argument("--runtime-bin", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11-refresh-test\lib\elmersolver"))
    parser.add_argument("--toolchain-bin", type=Path, default=Path(r"C:\msys64\ucrt64\bin"))
    args = parser.parse_args()
    solver, runtime_bin, toolchain_bin = map(absolute, (args.solver, args.runtime_bin, args.toolchain_bin))
    for path, label in ((solver, "solver"), (runtime_bin, "runtime DLL directory"), (toolchain_bin, "toolchain DLL directory")):
        if not path.exists():
            raise SystemExit(f"{label} not found: {path}")
    make_project()
    sync = subprocess.run(
        [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(PROJECT)],
        cwd=ROOT, capture_output=True, text=True,
    )
    (ARTIFACT_DIR / "sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
    if sync.returncode != 0:
        return sync.returncode
    command = [
        sys.executable, str(ROOT / "run.py"), CASE,
        "--project", str(PROJECT), "--skip-sync",
        "--elmer-solver", str(solver),
        "--runtime-bin", str(runtime_bin),
        "--toolchain-bin", str(toolchain_bin),
        "--mpi-procs", "1",
    ]
    env = os.environ.copy()
    env["PHASE24_HYPRE_FORCE_KRYLOV_SETUP"] = "1"
    env["PHASE24_HYPRE_FORCE_AMG_SETUP"] = "1"
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    t0 = time.monotonic()
    with LAUNCHER_LOG.open("w", encoding="utf-8") as handle:
        process = subprocess.run(command, cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "started": started,
        "finished": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - t0,
        "case": CASE,
        "project": str(PROJECT),
        "solver": str(solver),
        "runtime_bin": str(runtime_bin),
        "environment": {
            "PHASE24_HYPRE_FORCE_KRYLOV_SETUP": "1",
            "PHASE24_HYPRE_FORCE_AMG_SETUP": "1",
        },
        "command": command,
        "exit_code": process.returncode,
        "solver_log": str(SOLVER_LOG),
        "all_done": SOLVER_LOG.is_file() and "ALL DONE" in SOLVER_LOG.read_text(encoding="utf-8", errors="replace"),
    }
    (ARTIFACT_DIR / "run.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return process.returncode


if __name__ == "__main__":
    raise SystemExit(main())

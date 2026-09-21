"""Test legacy HYPRE BiCGStab/BoomerAMG on the pre-pulse operating point."""
from __future__ import annotations

import copy
import csv
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/phase24_hypre_tolerance_diagnostic_1us/phase24_hypre_tolerance_diagnostic.json"
ARTIFACT = ROOT / "artifacts/phase24_hypre_bicgstab_2step_diagnostic"
CASE = "case_phase24_hypre_bicgstab_1e8_2step"
SOLVER = Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11-refresh-test\bin\ElmerSolver_mpi.exe")
RUNTIME = Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11-refresh-test\lib\elmersolver")
TOOLCHAIN = Path(r"C:\msys64\ucrt64\bin")


def main() -> int:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    candidate = copy.deepcopy(project["cases"]["case_phase24_short_hypre_tol5e7_1us"])
    candidate.update({
        "timesteps": [["18[us]", 1], ["1[us]", 1]],
        "output_intervals": [999999, 999999],
        "series_file": f"{CASE}_series.csv",
        "iteration_series_file": f"{CASE}_iterations.csv",
        "solver_comment": "Diagnostic: HYPRE BiCGStab/BoomerAMG, tolerance 1e-8, two pre-pulse steps",
    })
    candidate["solver"] = copy.deepcopy(candidate["solver"])
    candidate["solver"]["linear_system"] = "iterative_hypre_boomeramg"
    candidate["solver"]["linear_system_convergence_tolerance"] = 1.0e-8
    candidate["solver"]["linear_system_max_iterations"] = 2000
    candidate["phase24_smoke"] = {"purpose": "HYPRE algorithm isolation", "variant": "BiCGStab + BoomerAMG", "no_matrix_dump": True, "no_vtu": True}
    project["cases"] = {CASE: candidate}
    ARTIFACT.mkdir(parents=True, exist_ok=True)
    project_path = ARTIFACT / "phase24_hypre_bicgstab_2step.json"
    project_path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sync = subprocess.run([sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project_path)], cwd=ROOT, capture_output=True, text=True)
    (ARTIFACT / "sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
    if sync.returncode:
        return sync.returncode
    solver, runtime, toolchain = (p if p.is_absolute() else ROOT / p for p in (SOLVER, RUNTIME, TOOLCHAIN))
    result_dir = ROOT / "results" / CASE
    result_dir.mkdir(parents=True, exist_ok=True)
    log = result_dir / "solver.log"
    command = [sys.executable, str(ROOT / "run.py"), CASE, "--project", str(project_path), "--skip-sync", "--elmer-solver", str(solver), "--runtime-bin", str(runtime), "--toolchain-bin", str(toolchain), "--mpi-procs", "1"]
    env = os.environ.copy()
    env["PHASE24_DISABLE_NATIVE_CAPTURE"] = "1"
    started = time.monotonic()
    with (result_dir / "bicgstab_2step_launcher.log").open("w", encoding="utf-8") as handle:
        process = subprocess.run(command, cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
    series = result_dir / f"{CASE}_series.csv"
    first = None
    if series.is_file():
        with series.open(newline="", encoding="utf-8") as handle:
            row = next(csv.DictReader(handle), None)
        if row:
            first = {"time_s": float(row["time_s"]), "temperature_K": float(row["tes_temperature_K"]), "current_uA": float(row["tes_current_A"]) * 1e6}
    payload = {"generated_utc": datetime.now(timezone.utc).isoformat(), "case": CASE, "exit_code": process.returncode, "elapsed_seconds": time.monotonic() - started, "all_done": log.is_file() and "ALL DONE" in log.read_text(encoding="utf-8", errors="replace"), "first_series": first, "solver_log": str(log), "references": {"COMSOL_uA": 143.05504932879472, "Phase24_MUMPS_uA": 144.268506298}}
    (ARTIFACT / "summary.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (ARTIFACT / "summary.md").write_text("# Phase24 HYPRE BiCGStab diagnostic\n\n" + json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return process.returncode


if __name__ == "__main__":
    raise SystemExit(main())

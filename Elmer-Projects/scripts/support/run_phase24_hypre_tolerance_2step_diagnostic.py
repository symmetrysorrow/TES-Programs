"""Run a cheap pre-pulse HYPRE tolerance diagnostic.

The full 1-us comparison takes several minutes because the mesh is large.  A
two-step case reaches the same pre-pulse operating point and is enough to
measure whether the HYPRE linear residual is responsible for the absolute
current offset.  Native capture is disabled so the test measures the solver
path rather than diagnostic vector reads.
"""
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
SOURCE_PROJECT = ROOT / "artifacts/phase24_hypre_tolerance_diagnostic_1us/phase24_hypre_tolerance_diagnostic.json"
ARTIFACT_DIR = ROOT / "artifacts/phase24_hypre_tolerance_2step_diagnostic"
SOLVER = Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11-refresh-test\bin\ElmerSolver_mpi.exe")
RUNTIME_BIN = Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11-refresh-test\lib\elmersolver")
TOOLCHAIN_BIN = Path(r"C:\msys64\ucrt64\bin")
BASE_CASE = "case_phase24_short_hypre_tol5e7_1us"
VARIANTS = {
    "1e-8": 1.0e-8,
    "1e-9": 1.0e-9,
}


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def make_project() -> tuple[Path, dict[str, str]]:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    base = project["cases"][BASE_CASE]
    cases: dict[str, str] = {}
    for label, tolerance in VARIANTS.items():
        case = f"case_phase24_hypre_tol_{label.replace('-', 'm').replace('.', 'p')}_2step"
        candidate = copy.deepcopy(base)
        candidate["timesteps"] = [["18[us]", 1], ["1[us]", 1]]
        candidate["output_intervals"] = [999999, 999999]
        candidate["series_file"] = f"{case}_series.csv"
        candidate["iteration_series_file"] = f"{case}_iterations.csv"
        candidate["solver_comment"] = f"Diagnostic: HYPRE linear tolerance {tolerance:.1e}, two pre-pulse steps"
        candidate["solver"] = copy.deepcopy(candidate["solver"])
        candidate["solver"]["linear_system_convergence_tolerance"] = tolerance
        candidate["phase24_smoke"] = {
            "purpose": "pre-pulse linear tolerance isolation",
            "variant": f"HYPRE tolerance {tolerance:.1e}",
            "no_matrix_dump": True,
            "no_vtu": True,
        }
        project["cases"][case] = candidate
        cases[label] = case
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    path = ARTIFACT_DIR / "phase24_hypre_tolerance_2step.json"
    path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path, cases


def read_first_series(path: Path) -> dict[str, float] | None:
    if not path.is_file():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle), None)
    if row is None:
        return None
    return {
        "time_s": float(row["time_s"]),
        "temperature_K": float(row["tes_temperature_K"]),
        "current_uA": float(row["tes_current_A"]) * 1.0e6,
        "resistance_ohm": float(row["tes_resistance_ohm"]),
        "power_W": float(row["tes_power_W"]),
    }


def main() -> int:
    solver = absolute(SOLVER)
    runtime_bin = absolute(RUNTIME_BIN)
    toolchain_bin = absolute(TOOLCHAIN_BIN)
    for path, label in ((SOURCE_PROJECT, "source project"), (solver, "solver"),
                        (runtime_bin, "runtime DLL directory"),
                        (toolchain_bin, "toolchain DLL directory")):
        if not path.exists():
            raise SystemExit(f"{label} not found: {path}")
    project_path, cases = make_project()
    sync = subprocess.run(
        [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project_path)],
        cwd=ROOT, capture_output=True, text=True,
    )
    (ARTIFACT_DIR / "sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
    if sync.returncode != 0:
        return sync.returncode

    results: dict[str, dict] = {}
    for label, case in cases.items():
        result_dir = ROOT / "results" / case
        log = result_dir / "solver.log"
        series = result_dir / f"{case}_series.csv"
        command = [
            sys.executable, str(ROOT / "run.py"), case,
            "--project", str(project_path), "--skip-sync",
            "--elmer-solver", str(solver), "--runtime-bin", str(runtime_bin),
            "--toolchain-bin", str(toolchain_bin), "--mpi-procs", "1",
        ]
        env = os.environ.copy()
        env["PHASE24_DISABLE_NATIVE_CAPTURE"] = "1"
        launcher = result_dir / "tolerance_2step_launcher.log"
        result_dir.mkdir(parents=True, exist_ok=True)
        started = time.monotonic()
        with launcher.open("w", encoding="utf-8") as handle:
            process = subprocess.run(command, cwd=ROOT, env=env,
                                     stdout=handle, stderr=subprocess.STDOUT)
        results[label] = {
            "case": case,
            "tolerance": VARIANTS[label],
            "exit_code": process.returncode,
            "elapsed_seconds": time.monotonic() - started,
            "all_done": log.is_file() and "ALL DONE" in log.read_text(encoding="utf-8", errors="replace"),
            "first_series": read_first_series(series),
            "solver_log": str(log),
        }

    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(project_path),
        "environment": {"PHASE24_DISABLE_NATIVE_CAPTURE": "1"},
        "results": results,
        "reference": {
            "comsol_baseline_uA": 143.05504932879472,
            "same_phase24_mumps_baseline_uA": 144.268506298,
        },
    }
    (ARTIFACT_DIR / "summary.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = [
        "# Phase24 HYPRE linear tolerance two-step diagnostic", "",
        "Only the two pre-pulse steps are run; native capture is disabled.", "",
        "| tolerance | exit | ALL DONE | first current [uA] | first T [K] | elapsed [s] |",
        "|---:|---:|:---:|---:|---:|---:|",
    ]
    for label, result in results.items():
        first = result["first_series"] or {}
        lines.append(f"| {result['tolerance']:.1e} | {result['exit_code']} | {result['all_done']} | {first.get('current_uA', 'n/a')} | {first.get('temperature_K', 'n/a')} | {result['elapsed_seconds']:.1f} |")
    lines += ["", "References: COMSOL baseline 143.055049 uA; same-Phase24 direct MUMPS baseline 144.268506 uA.", ""]
    (ARTIFACT_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0 if all(item["exit_code"] == 0 and item["all_done"] for item in results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())

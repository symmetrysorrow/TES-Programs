"""Test whether the Phase24 HYPRE 1e-8 failure is only an iteration limit.

The one-microsecond Phase24 case is reproduced with the same matrix assembly,
BDF2, nonlinear settings, HYPRE reuse, and adaptive lagging.  The linear
tolerance remains ``1e-8``; only ``Linear System Max Iterations`` changes from
2000 to 10000.  Native failure telemetry is collected from the rebuilt
solver, and completed series are compared with CPU/MUMPS and the existing
5e-7 HYPRE result.

Run from the repository root::

    python scripts/support/run_phase24_hypre_max_iterations_diagnostic.py
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PROJECT = ROOT / "elmer_project_phase24_production.json"
BASE_CASE = "case_phase24_hypre_cpu_pulse_1ms_5us"
CASE = "case_phase24_short_hypre_tol1e8_max10000_1us"
TOLERANCE = 1.0e-8
MAX_ITERATIONS = 10000
REFERENCE_SERIES = ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv"
LOOSE_SERIES = ROOT / "results/case_phase24_short_hypre_tol5e7_1us/case_phase24_short_hypre_tol5e7_1us_series.csv"
DIAGNOSTIC_DIR = ROOT / "artifacts/phase24_hypre_max_iterations_diagnostic_1us"
PROJECT_PATH = DIAGNOSTIC_DIR / "phase24_hypre_max_iterations_diagnostic.json"
COMPARISON_ROOT = ROOT / "artifacts/comparison"
RESULT_DIR = ROOT / "results" / CASE
LAUNCHER_LOG = RESULT_DIR / "max_iterations_diagnostic_launcher.log"
SOLVER_LOG = RESULT_DIR / "solver.log"
SERIES = RESULT_DIR / f"{CASE}_series.csv"


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def build_project() -> Path:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    base = project["cases"][BASE_CASE]
    candidate = copy.deepcopy(base)
    candidate["timesteps"] = base["timesteps"][:5]
    candidate["output_intervals"] = base["output_intervals"][:5]
    candidate["series_file"] = f"{CASE}_series.csv"
    candidate["iteration_series_file"] = f"{CASE}_iterations.csv"
    candidate["output_result_path"] = None
    candidate["output_file_path"] = None
    candidate["solver_comment"] = (
        "Diagnostic: same Phase24 HYPRE path; tolerance=1.0e-8; max iterations=10000"
    )
    candidate["comparison_time_grid"] = {
        "mode": "Phase24 HYPRE maximum-iteration diagnostic",
        "purpose": "distinguish iteration-limit failure from residual stagnation",
        "post_pulse_end": "1[us]",
    }
    candidate["solver"]["linear_system_convergence_tolerance"] = TOLERANCE
    candidate["solver"]["linear_system_max_iterations"] = MAX_ITERATIONS
    candidate["phase24_hypre_reuse"] = True
    candidate["phase24_preconditioner_lagging"] = "adaptive"
    candidate["phase24_smoke"] = {
        "purpose": "short HYPRE max-iteration diagnostic; approximately 1 us after pulse",
        "variant": "HYPRE tolerance 1e-8, max iterations 10000",
        "no_matrix_dump": True,
        "no_vtu": True,
    }
    project["cases"][CASE] = candidate
    DIAGNOSTIC_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_PATH.write_text(
        json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return PROJECT_PATH


def command_for(solver: Path, runtime_bin: Path, toolchain_bin: Path) -> list[str]:
    return [
        sys.executable, str(ROOT / "run.py"), CASE,
        "--project", str(PROJECT_PATH), "--skip-sync",
        "--elmer-solver", str(solver),
        "--runtime-bin", str(runtime_bin),
        "--toolchain-bin", str(toolchain_bin),
        "--mpi-procs", "1",
    ]


def audit_log(path: Path) -> dict:
    if not path.is_file():
        return {"solver_log_exists": False, "telemetry_records": []}
    text = path.read_text(encoding="utf-8", errors="replace")
    records = []
    for match in re.finditer(r"PHASE24_HYPRE_SOLVE_FAILURE\s+([^\n\r]+)", text):
        record = {}
        for item in match.group(1).split():
            if "=" not in item:
                continue
            key, value = item.split("=", 1)
            try:
                record[key] = int(value)
            except ValueError:
                try:
                    record[key] = float(value.replace("D", "E").replace("d", "e"))
                except ValueError:
                    record[key] = value
        records.append(record)
    return {
        "solver_log_exists": True,
        "all_done": "ALL DONE" in text,
        "stop_1": text.count("STOP 1"),
        "telemetry_records": records,
        "telemetry_missing": "PHASE24_HYPRE_SOLVE_FAILURE" not in text,
        "required_iterations": [int(v) for v in re.findall(r"Required iterations\s+([0-9]+)", text)],
        "required_residuals": [float(v.replace("D", "E").replace("d", "e")) for v in re.findall(r"to norm\s+([-+0-9.EeDd]+)", text)],
    }


def compare(candidate: Path, reference: Path, candidate_label: str,
            reference_label: str, output: Path) -> dict:
    if not candidate.is_file() or not reference.is_file():
        return {"status": "skipped", "reason": "input series missing", "output": str(output)}
    command = [
        sys.executable, str(ROOT / "scripts/analysis/phase24_difference_evidence.py"),
        "--candidate", str(candidate), "--reference", str(reference),
        "--candidate-label", candidate_label, "--reference-label", reference_label,
        "--out", str(output),
        "--checkpoint-us", "0.1", "0.2", "0.3", "0.4", "0.5",
        "0.6", "0.7", "0.8", "0.9", "1.0",
    ]
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    return {
        "status": "completed" if process.returncode == 0 else "failed",
        "exit_code": process.returncode, "stdout": process.stdout,
        "stderr": process.stderr, "output": str(output), "command": command,
    }


def write_summary(summary: dict) -> None:
    DIAGNOSTIC_DIR.mkdir(parents=True, exist_ok=True)
    (DIAGNOSTIC_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    audit = summary["audit"]
    records = audit.get("telemetry_records", [])
    lines = [
        "# Phase24 HYPRE maximum-iteration diagnostic",
        "",
        f"Generated: `{summary['generated_utc']}`",
        "",
        "- Same Phase24/BDF2/HYPRE policy as the 1e-8 diagnostic",
        "- Only `Linear System Max Iterations` changed to `10000`; tolerance remains `1e-8`",
        f"- Exit code: `{summary['run'].get('exit_code')}`",
        f"- ALL DONE: `{audit.get('all_done')}`",
        f"- Successful solve records: `{len(audit.get('required_iterations', []))}`",
        f"- Failure telemetry records: `{len(records)}`",
        "",
        "## Decision",
        "",
    ]
    if audit.get("all_done"):
        lines.append("The requested tolerance was reached in at least one complete run. Compare the generated absolute and baseline-corrected waveforms before changing production settings.")
    elif records and all(item.get("iterations") == MAX_ITERATIONS for item in records):
        lines.append("The solver still reaches the 10000-iteration limit without satisfying 1e-8; this is residual stagnation or an effective convergence barrier, not merely the previous 2000-iteration cap.")
    else:
        lines.append("The run did not complete the intended convergence test; inspect the telemetry records and solver log.")
    lines.extend(["", "## Failure telemetry", ""])
    if records:
        lines.extend(["```json", json.dumps(records, indent=2, ensure_ascii=False), "```"])
    else:
        lines.append("No native failure telemetry record was found.")
    lines.extend([
        "", f"Solver log: `{SOLVER_LOG}`", f"Launcher log: `{LAUNCHER_LOG}`",
        "", f"1e-8 vs CPU/MUMPS: `{summary['comparisons']['tol1e8_vs_cpu']['output']}`",
        f"1e-8 vs HYPRE 5e-7: `{summary['comparisons']['tol1e8_vs_tol5e7']['output']}`",
        "",
    ])
    (DIAGNOSTIC_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin\ElmerSolver_mpi.exe"))
    parser.add_argument("--runtime-bin", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\lib\elmersolver"))
    parser.add_argument("--toolchain-bin", type=Path, default=Path(r"C:\msys64\ucrt64\bin"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if not REFERENCE_SERIES.is_file():
        raise SystemExit(f"reference series not found: {REFERENCE_SERIES}")
    solver, runtime_bin, toolchain_bin = map(absolute, (args.solver, args.runtime_bin, args.toolchain_bin))
    if not args.dry_run:
        for path, label in ((solver, "solver"), (runtime_bin, "runtime DLL directory"), (toolchain_bin, "toolchain DLL directory")):
            if not path.exists():
                raise SystemExit(f"{label} not found: {path}")

    project = build_project()
    sync_command = [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project)]
    if args.dry_run:
        print("[dry-run] " + " ".join(sync_command))
    else:
        sync = subprocess.run(sync_command, cwd=ROOT, capture_output=True, text=True)
        (DIAGNOSTIC_DIR / "sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
        if sync.returncode != 0:
            raise SystemExit(f"sync_elmer_parameters.py failed with exit code {sync.returncode}")
        generated_sif = ROOT / "generated/cases" / f"{CASE}.sif"
        sif_text = generated_sif.read_text(encoding="utf-8", errors="replace")
        expected_line = f"Linear System Max Iterations = {MAX_ITERATIONS}"
        if expected_line not in sif_text:
            raise SystemExit(
                f"generated SIF did not contain the requested iteration limit: "
                f"{expected_line}; inspect {generated_sif}"
            )

    command = command_for(solver, runtime_bin, toolchain_bin)
    print("[max-iterations] command: " + " ".join(command))
    if args.dry_run:
        run = {"command": command, "exit_code": None, "log": str(LAUNCHER_LOG)}
        audit = {}
        comparisons = {
            "tol1e8_vs_cpu": {"status": "dry-run", "output": str(COMPARISON_ROOT / "phase24_hypre_tol1e8_max10000_vs_cpu_1us")},
            "tol1e8_vs_tol5e7": {"status": "dry-run", "output": str(COMPARISON_ROOT / "phase24_hypre_tol1e8_max10000_vs_tol5e7_1us")},
        }
    else:
        start = time.monotonic()
        started = datetime.now(timezone.utc).isoformat()
        env = os.environ.copy()
        env["PHASE24_DISABLE_NATIVE_CAPTURE"] = "1"
        RESULT_DIR.mkdir(parents=True, exist_ok=True)
        with LAUNCHER_LOG.open("w", encoding="utf-8") as handle:
            process = subprocess.run(command, cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
        run = {
            "command": command, "started": started,
            "finished": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic() - start,
            "exit_code": process.returncode, "log": str(LAUNCHER_LOG),
        }
        audit = audit_log(SOLVER_LOG)
        comparisons = {
            "tol1e8_vs_cpu": compare(SERIES, REFERENCE_SERIES, "HYPRE 1e-8, max iterations 10000", "CPU/MUMPS Phase23", COMPARISON_ROOT / "phase24_hypre_tol1e8_max10000_vs_cpu_1us"),
            "tol1e8_vs_tol5e7": compare(SERIES, LOOSE_SERIES, "HYPRE 1e-8, max iterations 10000", "HYPRE tolerance 5e-7", COMPARISON_ROOT / "phase24_hypre_tol1e8_max10000_vs_tol5e7_1us"),
        }

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(project), "tolerance": TOLERANCE,
        "max_iterations": MAX_ITERATIONS, "run": run,
        "solver_log": str(SOLVER_LOG), "series": str(SERIES),
        "audit": audit, "comparisons": comparisons,
    }
    write_summary(summary)
    print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
    return 1 if run.get("exit_code") not in (0, None) else 0


if __name__ == "__main__":
    raise SystemExit(main())

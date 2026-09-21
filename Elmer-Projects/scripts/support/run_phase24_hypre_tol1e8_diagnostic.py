"""Run the next Phase24 HYPRE tolerance diagnostic at ``1e-8``.

The case is identical to the previously tested one-microsecond HYPRE case;
only ``Linear System Convergence Tolerance`` is changed to ``1e-8``.  The
completed result is compared with the CPU/MUMPS reference and the existing
``5e-7`` HYPRE result, including both absolute current and baseline-corrected
waveform evidence.

Run from the repository root::

    python scripts/support/run_phase24_hypre_tol1e8_diagnostic.py

Use ``--dry-run`` to generate the project and print the command without
starting Elmer.
"""
from __future__ import annotations

import argparse
import copy
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PROJECT = ROOT / "elmer_project_phase24_production.json"
BASE_CASE = "case_phase24_hypre_cpu_pulse_1ms_5us"
CASE = "case_phase24_short_hypre_tol1e8_1us"
TOLERANCE = 1.0e-8
REFERENCE_SERIES = ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv"
LOOSE_SERIES = ROOT / "results/case_phase24_short_hypre_tol5e7_1us/case_phase24_short_hypre_tol5e7_1us_series.csv"
DIAGNOSTIC_DIR = ROOT / "artifacts/phase24_hypre_tol1e8_diagnostic_1us"
PROJECT_PATH = DIAGNOSTIC_DIR / "phase24_hypre_tol1e8_diagnostic.json"
COMPARISON_ROOT = ROOT / "artifacts/comparison"


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
        "Diagnostic: same Phase24 HYPRE path; linear tolerance=1.0e-8"
    )
    candidate["comparison_time_grid"] = {
        "mode": "Phase24 HYPRE intermediate linear-tolerance diagnostic",
        "purpose": "approximately 1 us post-pulse absolute and baseline-corrected comparison",
        "post_pulse_end": "1[us]",
    }
    candidate["solver"]["linear_system_convergence_tolerance"] = TOLERANCE
    candidate["phase24_hypre_reuse"] = True
    candidate["phase24_preconditioner_lagging"] = "adaptive"
    candidate["phase24_smoke"] = {
        "purpose": "short HYPRE tolerance diagnostic; approximately 1 us after pulse",
        "variant": "HYPRE tolerance 1e-8",
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
        sys.executable,
        str(ROOT / "run.py"),
        CASE,
        "--project", str(PROJECT_PATH),
        "--skip-sync",
        "--elmer-solver", str(solver),
        "--runtime-bin", str(runtime_bin),
        "--toolchain-bin", str(toolchain_bin),
        "--mpi-procs", "1",
    ]


def run_case(command: list[str], log_path: Path, dry_run: bool) -> dict:
    started = datetime.now(timezone.utc).isoformat()
    if dry_run:
        return {"command": command, "started": started, "finished": started,
                "exit_code": None, "log": str(log_path)}
    start = time.monotonic()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        process = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    return {
        "command": command,
        "started": started,
        "finished": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - start,
        "exit_code": process.returncode,
        "log": str(log_path),
    }


def audit_log(path: Path) -> dict:
    if not path.is_file():
        return {"log_exists": False}
    text = path.read_text(encoding="utf-8", errors="replace")
    iterations = [int(value) for value in re.findall(r"Required iterations\s+([0-9]+)", text)]
    relative_changes = []
    for value in re.findall(r"Relative Change\s*:\s*([-+0-9.EeDd]+)", text):
        try:
            relative_changes.append(float(value.replace("D", "E").replace("d", "e")))
        except ValueError:
            pass
    return {
        "log_exists": True,
        "all_done": "ALL DONE" in text,
        "abort_markers": sum(text.count(marker) for marker in (
            "STOP 1", "job aborted", "application aborted", "Could not allocate memory"
        )),
        "hypre_solve_calls": len(iterations),
        "required_iterations_zero": iterations.count(0),
        "required_iterations_max": max(iterations, default=None),
        "required_iterations_sum": sum(iterations),
        "relative_change_samples": len(relative_changes),
        "relative_change_zero": sum(value == 0.0 for value in relative_changes),
        "relative_change_max": max(relative_changes, default=None),
        "tolerance_failure_marker": "failed to satisfy the production linear tolerance" in text,
        "native_xvec_selector_warning": "NATIVE_XVEC_CAPTURE_CANDIDATE_SEEN_BUT_SELECTOR_MISMATCH" in text,
    }


def compare(candidate: Path, reference: Path, candidate_label: str,
            reference_label: str, output: Path) -> dict:
    command = [
        sys.executable,
        str(ROOT / "scripts/analysis/phase24_difference_evidence.py"),
        "--candidate", str(candidate),
        "--reference", str(reference),
        "--candidate-label", candidate_label,
        "--reference-label", reference_label,
        "--out", str(output),
        "--checkpoint-us", "0.1", "0.2", "0.3", "0.4", "0.5",
        "0.6", "0.7", "0.8", "0.9", "1.0",
    ]
    if not candidate.is_file() or not reference.is_file():
        return {"status": "skipped", "reason": "input series missing", "output": str(output)}
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    return {
        "status": "completed" if process.returncode == 0 else "failed",
        "command": command,
        "exit_code": process.returncode,
        "stdout": process.stdout,
        "stderr": process.stderr,
        "output": str(output),
    }


def write_summary(summary: dict) -> None:
    (DIAGNOSTIC_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    run = summary["run"]
    audit = summary["audit"]
    lines = [
        "# Phase24 HYPRE intermediate tolerance diagnostic",
        "",
        f"Generated: `{summary['generated_utc']}`",
        "",
        "## Fixed conditions",
        "",
        "- Same mesh, restart, pulse, inner circuit, Phase24 assembly, BDF2, nonlinear controls, HYPRE reuse, and adaptive preconditioner lagging",
        "- Same approximately `1 us` post-pulse window",
        "- Only `Linear System Convergence Tolerance` changes from the previous `5e-7` case",
        "",
        "## Result",
        "",
        "| condition | exit | ALL DONE | HYPRE solves | zero-iteration solves | max iterations | tolerance failure |",
        "|---|---:|---:|---:|---:|---:|---:|",
        f"| HYPRE tolerance 1e-8 | {run.get('exit_code', 'n/a')} | {audit.get('all_done', 'n/a')} | "
        f"{audit.get('hypre_solve_calls', 'n/a')} | {audit.get('required_iterations_zero', 'n/a')} | "
        f"{audit.get('required_iterations_max', 'n/a')} | {audit.get('tolerance_failure_marker', 'n/a')} |",
        "",
        "## Comparisons",
        "",
        f"- 1e-8 vs CPU/MUMPS: `{summary['comparisons']['tol1e8_vs_cpu']['output']}`",
        f"- 1e-8 vs HYPRE 5e-7: `{summary['comparisons']['tol1e8_vs_tol5e7']['output']}`",
        f"- Existing HYPRE 5e-7 vs CPU/MUMPS: `{summary['comparisons']['tol5e7_vs_cpu']['output']}`",
        "",
        "Interpretation: compare both absolute current and baseline-corrected waveform. A successful 1e-8 run that moves both toward the reference supports insufficient linear convergence as a cause; unchanged results point to another backend or scaling issue. The selector-mismatch warning, if present, is recorded separately and is not treated as the root cause.",
        "",
        f"Generated project: `{PROJECT_PATH}`",
        "",
    ]
    (DIAGNOSTIC_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin\ElmerSolver_mpi.exe"))
    parser.add_argument("--runtime-bin", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\lib\elmersolver"))
    parser.add_argument("--toolchain-bin", type=Path, default=Path(r"C:\msys64\ucrt64\bin"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    solver = absolute(args.solver)
    runtime_bin = absolute(args.runtime_bin)
    toolchain_bin = absolute(args.toolchain_bin)
    required = [(solver, "solver"), (runtime_bin, "runtime DLL directory"), (toolchain_bin, "toolchain DLL directory")]
    if not REFERENCE_SERIES.is_file():
        raise SystemExit(f"reference series not found: {REFERENCE_SERIES}")
    if not args.dry_run:
        for path, label in required:
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

    launcher_log = ROOT / "results" / CASE / "tol1e8_diagnostic_launcher.log"
    solver_log = ROOT / "results" / CASE / "solver.log"
    series = ROOT / "results" / CASE / f"{CASE}_series.csv"
    command = command_for(solver, runtime_bin, toolchain_bin)
    print("[tol1e8] command: " + " ".join(command))
    run = run_case(command, launcher_log, args.dry_run)
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tolerance": TOLERANCE,
        "project": str(project),
        "run": run,
        "solver_log": str(solver_log),
        "series": str(series),
        "audit": audit_log(solver_log) if not args.dry_run else {},
        "comparisons": {},
    }
    if args.dry_run:
        summary["comparisons"] = {
            "tol1e8_vs_cpu": {"status": "dry-run", "output": str(COMPARISON_ROOT / "phase24_hypre_tol1e8_vs_cpu_1us")},
            "tol1e8_vs_tol5e7": {"status": "dry-run", "output": str(COMPARISON_ROOT / "phase24_hypre_tol1e8_vs_tol5e7_1us")},
            "tol5e7_vs_cpu": {"status": "dry-run", "output": str(COMPARISON_ROOT / "phase24_hypre_tol5e7_vs_cpu_1us")},
        }
    else:
        summary["comparisons"]["tol1e8_vs_cpu"] = compare(
            series, REFERENCE_SERIES, "HYPRE tolerance 1e-8", "CPU/MUMPS Phase23",
            COMPARISON_ROOT / "phase24_hypre_tol1e8_vs_cpu_1us"
        )
        summary["comparisons"]["tol1e8_vs_tol5e7"] = compare(
            series, LOOSE_SERIES, "HYPRE tolerance 1e-8", "HYPRE tolerance 5e-7",
            COMPARISON_ROOT / "phase24_hypre_tol1e8_vs_tol5e7_1us"
        )
        summary["comparisons"]["tol5e7_vs_cpu"] = compare(
            LOOSE_SERIES, REFERENCE_SERIES, "HYPRE tolerance 5e-7", "CPU/MUMPS Phase23",
            COMPARISON_ROOT / "phase24_hypre_tol5e7_vs_cpu_1us"
        )
    write_summary(summary)
    print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
    return 1 if run.get("exit_code") not in (0, None) else 0


if __name__ == "__main__":
    raise SystemExit(main())

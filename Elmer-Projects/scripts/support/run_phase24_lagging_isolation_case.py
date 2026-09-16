"""Run the isolated Phase24 preconditioner-lagging diagnostic.

This is the next cause-isolation case after the short difference campaign:

* BDF2 is kept unchanged;
* ``Phase24 HYPRE Reuse`` remains enabled;
* only ``Phase24 Preconditioner Lagging`` is changed from ``adaptive`` to
  ``disabled``.

Unlike the previous ``reuse_off`` cases, this keeps the HYPRE lifecycle and
container reuse enabled.  It therefore isolates the AMG/preconditioner
lagging policy instead of changing several lifecycle controls at once.

From the repository root::

    python scripts/support/run_phase24_lagging_isolation_case.py

Use ``--dry-run`` to generate the project and print the commands without
starting Elmer.  The run is intentionally limited to the same approximately
1 us post-pulse window as the preceding campaign.
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
REFERENCE_SERIES = ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv"
DIAGNOSTIC_DIR = ROOT / "artifacts/phase24_lagging_isolation_1us"
PROJECT_PATH = DIAGNOSTIC_DIR / "phase24_lagging_isolation.json"
CASE_NAME = "case_phase24_short_bdf2_reuse_on_lagging_disabled_1us"
COMPARISON_DIR = ROOT / "artifacts/comparison/phase24_bdf2_reuse_on_lagging_disabled_1us"
BASE_CASE = "case_phase24_hypre_cpu_pulse_1ms_5us"


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def build_project() -> Path:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    base = project["cases"][BASE_CASE]
    candidate = copy.deepcopy(base)
    candidate["timesteps"] = base["timesteps"][:5]
    candidate["output_intervals"] = base["output_intervals"][:5]
    candidate["series_file"] = f"{CASE_NAME}_series.csv"
    candidate["iteration_series_file"] = f"{CASE_NAME}_iterations.csv"
    candidate["output_result_path"] = None
    candidate["output_file_path"] = None
    candidate["comparison_time_grid"] = {
        "mode": "Phase24 preconditioner-lagging isolation diagnostic",
        "purpose": "approximately 1 us post-pulse comparison",
        "post_pulse_end": "1[us]",
    }
    candidate["solver_comment"] = (
        "Diagnostic: BDF2, HYPRE lifecycle reuse ON, preconditioner lagging disabled"
    )
    candidate["bdf_order"] = 2
    candidate["phase24_hypre_reuse"] = True
    candidate["phase24_preconditioner_lagging"] = "disabled"
    candidate["phase24_smoke"] = {
        "purpose": "short isolation run; approximately 1 us after pulse",
        "variant": "BDF2 / HYPRE reuse ON / preconditioner lagging disabled",
        "no_matrix_dump": True,
        "no_vtu": True,
    }
    project["cases"][CASE_NAME] = candidate
    DIAGNOSTIC_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_PATH.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    return PROJECT_PATH


def command_for(solver: Path, runtime_bin: Path, toolchain_bin: Path) -> list[str]:
    return [
        sys.executable,
        str(ROOT / "run.py"),
        CASE_NAME,
        "--project", str(PROJECT_PATH),
        "--skip-sync",
        "--elmer-solver", str(solver),
        "--runtime-bin", str(runtime_bin),
        "--toolchain-bin", str(toolchain_bin),
        "--mpi-procs", "1",
    ]


def run_process(command: list[str], log_path: Path, dry_run: bool) -> dict:
    started = datetime.now(timezone.utc).isoformat()
    if dry_run:
        return {
            "command": command,
            "started": started,
            "finished": started,
            "elapsed_seconds": None,
            "exit_code": None,
            "log": str(log_path),
        }
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


def audit_solver_log(log_path: Path) -> dict:
    """Collect convergence/lifecycle evidence without declaring ALL DONE valid."""
    if not log_path.is_file():
        return {"log_exists": False}
    text = log_path.read_text(encoding="utf-8", errors="replace")
    required_iterations = [
        int(value) for value in re.findall(r"Required iterations\s+([0-9]+)", text)
    ]
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
            "job aborted", "application aborted", "STOP 1", "Could not allocate memory",
        )),
        "hypre_lifecycle_lines": text.count("Phase24 HYPRE lifecycle:"),
        "preconditioner_decision_lines": text.count("Phase24 preconditioner decision:"),
        "hypre_solve_calls": len(required_iterations),
        "hypre_required_iterations_zero": required_iterations.count(0),
        "hypre_required_iterations_max": max(required_iterations, default=None),
        "relative_change_samples": len(relative_changes),
        "relative_change_zero": sum(value == 0.0 for value in relative_changes),
        "relative_change_max": max(relative_changes, default=None),
        "nonlinear_iteration_limit_markers": text.count("maximum number of nonlinear iterations"),
    }


def run_comparison(series: Path, dry_run: bool) -> dict:
    command = [
        sys.executable,
        str(ROOT / "scripts/analysis/phase24_difference_evidence.py"),
        "--candidate", str(series),
        "--reference", str(REFERENCE_SERIES),
        "--candidate-label", "native HYPRE / reuse ON / lagging disabled",
        "--reference-label", "CPU/MUMPS Phase23",
        "--out", str(COMPARISON_DIR),
        "--checkpoint-us", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0",
    ]
    if dry_run:
        return {"command": command, "exit_code": None}
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    return {
        "command": command,
        "exit_code": process.returncode,
        "stdout": process.stdout,
        "stderr": process.stderr,
        "output": str(COMPARISON_DIR),
    }


def write_summary(summary: dict) -> None:
    (DIAGNOSTIC_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    run = summary.get("run", {})
    audit = summary.get("log_audit", {})
    comparison = summary.get("comparison", {})
    lines = [
        "# Phase24 preconditioner-lagging isolation",
        "",
        f"Generated: `{summary['generated_utc']}`",
        "",
        "## Controlled change",
        "",
        "- BDF order: `2` (unchanged)",
        "- `Phase24 HYPRE Reuse`: `True` (unchanged/enabled)",
        "- `Phase24 Preconditioner Lagging`: `disabled` (only intended change)",
        "- Window: approximately `1 us` after the pulse",
        "",
        "## Run result",
        "",
        f"- Process exit code: `{run.get('exit_code', 'n/a')}`",
        f"- `ALL DONE` marker: `{audit.get('all_done', 'n/a')}`",
        f"- Abort/allocation markers: `{audit.get('abort_markers', 'n/a')}`",
        f"- HYPRE solves: `{audit.get('hypre_solve_calls', 'n/a')}`; zero-iteration solves: `{audit.get('hypre_required_iterations_zero', 'n/a')}`",
        f"- Relative-change samples: `{audit.get('relative_change_samples', 'n/a')}`; zero values: `{audit.get('relative_change_zero', 'n/a')}`",
        "",
        "`ALL DONE` is a process-completion marker, not proof that every nonlinear solve converged. Review the solver log and the zero-iteration count together.",
        "",
        f"Comparison: `{comparison.get('output', 'not run')}` (exit `{comparison.get('exit_code', 'n/a')}`)",
        "",
        f"Solver log: `{summary.get('solver_log', run.get('log', 'n/a'))}`",
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
    parser.add_argument(
        "--audit-existing",
        action="store_true",
        help="re-audit the existing solver.log and comparison without rerunning Elmer",
    )
    args = parser.parse_args()

    solver = absolute(args.solver)
    runtime_bin = absolute(args.runtime_bin)
    toolchain_bin = absolute(args.toolchain_bin)
    if not REFERENCE_SERIES.is_file():
        raise SystemExit(f"reference series not found: {REFERENCE_SERIES}")
    if not args.dry_run:
        for path, label in ((solver, "solver"), (runtime_bin, "runtime DLL directory"), (toolchain_bin, "toolchain DLL directory")):
            if not path.exists():
                raise SystemExit(f"{label} not found: {path}")

    if args.audit_existing:
        summary_path = DIAGNOSTIC_DIR / "summary.json"
        if not summary_path.is_file():
            raise SystemExit(f"existing diagnostic summary not found: {summary_path}")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        solver_log_path = ROOT / "results" / CASE_NAME / "solver.log"
        series = ROOT / "results" / CASE_NAME / f"{CASE_NAME}_series.csv"
        summary["solver_log"] = str(solver_log_path)
        summary["log_audit"] = audit_solver_log(solver_log_path)
        if series.is_file():
            summary["comparison"] = run_comparison(series, False)
        else:
            summary["comparison"] = {"error": f"series not found: {series}"}
        write_summary(summary)
        print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
        audit_ok = summary["log_audit"].get("all_done") and not summary["log_audit"].get("abort_markers")
        comparison_ok = summary["comparison"].get("exit_code") == 0
        return 0 if audit_ok and comparison_ok else 1

    project = build_project()
    result_dir = ROOT / "results" / CASE_NAME
    launcher_log_path = result_dir / "lagging_isolation_launcher.log"
    solver_log_path = result_dir / "solver.log"
    series = ROOT / "results" / CASE_NAME / f"{CASE_NAME}_series.csv"
    command = command_for(solver, runtime_bin, toolchain_bin)
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(project),
        "case": CASE_NAME,
        "solver": str(solver),
        "runtime_bin": str(runtime_bin),
        "toolchain_bin": str(toolchain_bin),
        "run": {},
        "log_audit": {},
        "comparison": {},
    }

    sync_command = [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(PROJECT_PATH)]
    if args.dry_run:
        print("[dry-run] " + " ".join(sync_command))
    else:
        sync = subprocess.run(sync_command, cwd=ROOT, capture_output=True, text=True)
        (DIAGNOSTIC_DIR / "sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
        if sync.returncode != 0:
            summary["sync_exit_code"] = sync.returncode
            write_summary(summary)
            return sync.returncode

    print("[lagging-isolation] command: " + " ".join(command))
    summary["run"] = run_process(command, launcher_log_path, args.dry_run)
    if not args.dry_run:
        summary["solver_log"] = str(solver_log_path)
        summary["log_audit"] = audit_solver_log(solver_log_path)
        if series.is_file():
            summary["comparison"] = run_comparison(series, False)
        else:
            summary["comparison"] = {"error": f"series not found: {series}"}
    else:
        summary["comparison"] = run_comparison(series, True)
    write_summary(summary)
    print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
    return 0 if summary["run"].get("exit_code") in (0, None) and "error" not in summary["comparison"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

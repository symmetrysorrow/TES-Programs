"""Run the same-Phase24 direct-MUMPS control case.

The case is copied from the Phase24 native-HYPRE production case and keeps
the mesh, restart, pulse, inner-circuit update, BDF2 order, timestep grid,
nonlinear controls, and Phase24 assembly settings unchanged.  Only the
linear-system backend changes from native HYPRE FlexGMRES/BoomerAMG to
Elmer's direct MUMPS solver.

Run from the repository root::

    python scripts/support/run_phase24_same_path_mumps_case.py

Use ``--dry-run`` to generate the project and print commands without starting
Elmer.  The diagnostic stops at the same approximately 1 us post-pulse
window as the preceding short cases.
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
CASE_NAME = "case_phase24_short_same_path_mumps_1us"
DIAGNOSTIC_DIR = ROOT / "artifacts/phase24_same_path_mumps_1us"
PROJECT_PATH = DIAGNOSTIC_DIR / "phase24_same_path_mumps.json"
RESULT_DIR = ROOT / "results" / CASE_NAME
SERIES = RESULT_DIR / f"{CASE_NAME}_series.csv"
CPU_REFERENCE = ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv"
HYPRE_REFERENCE = ROOT / "results/case_phase24_short_bdf2_reuse_on_1us/case_phase24_short_bdf2_reuse_on_1us_series.csv"
CPU_COMPARISON = ROOT / "artifacts/comparison/phase24_same_path_mumps_vs_cpu_1us"
HYPRE_COMPARISON = ROOT / "artifacts/comparison/phase24_same_path_mumps_vs_hypre_1us"


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
    candidate["solver_comment"] = (
        "Diagnostic: same Phase24 path, BDF2, HYPRE settings retained; direct MUMPS backend"
    )
    candidate["comparison_time_grid"] = {
        "mode": "Phase24 same-path direct-MUMPS isolation diagnostic",
        "purpose": "approximately 1 us post-pulse comparison against native HYPRE and CPU/MUMPS",
        "post_pulse_end": "1[us]",
    }
    candidate["solver"]["linear_system"] = "mumps"
    candidate["phase24_smoke"] = {
        "purpose": "short backend isolation run; approximately 1 us after pulse",
        "variant": "Phase24 BDF2 / direct MUMPS",
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


def run_process(command: list[str], launcher_log: Path, dry_run: bool) -> dict:
    started = datetime.now(timezone.utc).isoformat()
    if dry_run:
        return {"command": command, "started": started, "finished": started, "exit_code": None, "log": str(launcher_log)}
    start = time.monotonic()
    launcher_log.parent.mkdir(parents=True, exist_ok=True)
    with launcher_log.open("w", encoding="utf-8") as handle:
        process = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    return {
        "command": command,
        "started": started,
        "finished": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - start,
        "exit_code": process.returncode,
        "log": str(launcher_log),
    }


def audit_solver_log(path: Path) -> dict:
    if not path.is_file():
        return {"log_exists": False}
    text = path.read_text(encoding="utf-8", errors="replace")
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
        "mumps_markers": text.count("MUMPS"),
        "relative_change_samples": len(relative_changes),
        "relative_change_zero": sum(value == 0.0 for value in relative_changes),
        "relative_change_max": max(relative_changes, default=None),
        "nonlinear_iteration_limit_markers": text.count("maximum number of nonlinear iterations"),
    }


def compare(candidate: Path, reference: Path, label: str, output: Path, dry_run: bool) -> dict:
    command = [
        sys.executable,
        str(ROOT / "scripts/analysis/phase24_difference_evidence.py"),
        "--candidate", str(candidate),
        "--reference", str(reference),
        "--candidate-label", label,
        "--reference-label", "reference",
        "--out", str(output),
        "--checkpoint-us", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0",
    ]
    if dry_run:
        return {"command": command, "exit_code": None, "output": str(output)}
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    return {
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
    run = summary.get("run", {})
    audit = summary.get("log_audit", {})
    lines = [
        "# Phase24 same-path direct-MUMPS isolation",
        "",
        f"Generated: `{summary['generated_utc']}`",
        "",
        "## Controlled change",
        "",
        "- Mesh, restart, pulse, inner circuit, BDF2, timestep grid, nonlinear controls: unchanged from Phase24 HYPRE case",
        "- Phase24 assembly/operator settings: unchanged",
        "- Linear backend: `iterative_hypre_flexgmres_boomeramg` -> `mumps`",
        "",
        "## Run result",
        "",
        f"- Process exit code: `{run.get('exit_code', 'n/a')}`",
        f"- `ALL DONE`: `{audit.get('all_done', 'n/a')}`",
        f"- Abort/allocation markers: `{audit.get('abort_markers', 'n/a')}`",
        f"- MUMPS markers: `{audit.get('mumps_markers', 'n/a')}`",
        f"- Relative-change samples: `{audit.get('relative_change_samples', 'n/a')}`; zero values: `{audit.get('relative_change_zero', 'n/a')}`",
        "",
        "The two comparison bundles separate backend mismatch from the older CPU/MUMPS reference's broader Phase23 differences.",
        "",
        f"Direct MUMPS vs CPU/MUMPS: `{summary.get('cpu_comparison', {}).get('output', 'not run')}`",
        f"Direct MUMPS vs native HYPRE: `{summary.get('hypre_comparison', {}).get('output', 'not run')}`",
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
    args = parser.parse_args()

    solver = absolute(args.solver)
    runtime_bin = absolute(args.runtime_bin)
    toolchain_bin = absolute(args.toolchain_bin)
    for path, label in ((CPU_REFERENCE, "CPU reference"), (HYPRE_REFERENCE, "HYPRE reference")):
        if not path.is_file():
            raise SystemExit(f"{label} series not found: {path}")
    if not args.dry_run:
        for path, label in ((solver, "solver"), (runtime_bin, "runtime DLL directory"), (toolchain_bin, "toolchain DLL directory")):
            if not path.exists():
                raise SystemExit(f"{label} not found: {path}")

    project = build_project()
    launcher_log = RESULT_DIR / "same_path_mumps_launcher.log"
    solver_log = RESULT_DIR / "solver.log"
    command = command_for(solver, runtime_bin, toolchain_bin)
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(project),
        "case": CASE_NAME,
        "solver": str(solver),
        "runtime_bin": str(runtime_bin),
        "toolchain_bin": str(toolchain_bin),
        "run": {},
        "solver_log": str(solver_log),
        "log_audit": {},
        "cpu_comparison": {},
        "hypre_comparison": {},
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

    print("[same-path-mumps] command: " + " ".join(command))
    summary["run"] = run_process(command, launcher_log, args.dry_run)
    if not args.dry_run:
        summary["log_audit"] = audit_solver_log(solver_log)
        if SERIES.is_file():
            summary["cpu_comparison"] = compare(SERIES, CPU_REFERENCE, "same Phase24 / direct MUMPS", CPU_COMPARISON, False)
            summary["hypre_comparison"] = compare(SERIES, HYPRE_REFERENCE, "same Phase24 / direct MUMPS", HYPRE_COMPARISON, False)
        else:
            summary["cpu_comparison"] = {"error": f"series not found: {SERIES}"}
            summary["hypre_comparison"] = {"error": f"series not found: {SERIES}"}
    else:
        summary["cpu_comparison"] = compare(SERIES, CPU_REFERENCE, "same Phase24 / direct MUMPS", CPU_COMPARISON, True)
        summary["hypre_comparison"] = compare(SERIES, HYPRE_REFERENCE, "same Phase24 / direct MUMPS", HYPRE_COMPARISON, True)
    write_summary(summary)
    print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
    errors = any("error" in summary[key] for key in ("cpu_comparison", "hypre_comparison"))
    return 0 if summary["run"].get("exit_code") in (0, None) and not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())

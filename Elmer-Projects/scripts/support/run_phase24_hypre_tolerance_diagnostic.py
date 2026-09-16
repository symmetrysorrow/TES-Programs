"""Compare two HYPRE linear convergence tolerances on the same Phase24 case.

Only ``Linear System Convergence Tolerance`` changes between the variants:

* ``5e-7``: current production diagnostic setting;
* ``1e-10``: strict linear-solve diagnostic setting.

Mesh, restart, pulse, inner circuit, Phase24 assembly, BDF2, HYPRE reuse,
adaptive preconditioner lagging, nonlinear controls, and the approximately
1-us post-pulse window are kept identical.  Each result is compared with the
historical CPU/MUMPS reference and the strict result is also compared directly
with the loose result.

Run from the repository root::

    python scripts/support/run_phase24_hypre_tolerance_diagnostic.py

Use ``--dry-run`` to generate the cases and print commands without running
Elmer.
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
REFERENCE_SERIES = ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv"
DIAGNOSTIC_DIR = ROOT / "artifacts/phase24_hypre_tolerance_diagnostic_1us"
PROJECT_PATH = DIAGNOSTIC_DIR / "phase24_hypre_tolerance_diagnostic.json"
COMPARISON_ROOT = ROOT / "artifacts/comparison"
VARIANTS = {
    "tol5e7": {
        "case": "case_phase24_short_hypre_tol5e7_1us",
        "tolerance": 5.0e-7,
        "label": "HYPRE tolerance 5e-7",
    },
    "tol1e10": {
        "case": "case_phase24_short_hypre_tol1e10_1us",
        "tolerance": 1.0e-10,
        "label": "HYPRE tolerance 1e-10",
    },
}


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def build_project() -> Path:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    base = project["cases"][BASE_CASE]
    for key, options in VARIANTS.items():
        candidate = copy.deepcopy(base)
        candidate["timesteps"] = base["timesteps"][:5]
        candidate["output_intervals"] = base["output_intervals"][:5]
        candidate["series_file"] = f"{options['case']}_series.csv"
        candidate["iteration_series_file"] = f"{options['case']}_iterations.csv"
        candidate["output_result_path"] = None
        candidate["output_file_path"] = None
        candidate["solver_comment"] = (
            f"Diagnostic: same Phase24 HYPRE path; linear tolerance={options['tolerance']:.1e}"
        )
        candidate["comparison_time_grid"] = {
            "mode": "Phase24 HYPRE linear-tolerance diagnostic",
            "purpose": "approximately 1 us post-pulse absolute and baseline-corrected comparison",
            "post_pulse_end": "1[us]",
        }
        candidate["solver"]["linear_system_convergence_tolerance"] = options["tolerance"]
        candidate["phase24_hypre_reuse"] = True
        candidate["phase24_preconditioner_lagging"] = "adaptive"
        candidate["phase24_smoke"] = {
            "purpose": "short HYPRE tolerance diagnostic; approximately 1 us after pulse",
            "variant": options["label"],
            "no_matrix_dump": True,
            "no_vtu": True,
        }
        project["cases"][options["case"]] = candidate
    DIAGNOSTIC_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_PATH.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    return PROJECT_PATH


def command_for(case: str, solver: Path, runtime_bin: Path, toolchain_bin: Path) -> list[str]:
    return [
        sys.executable,
        str(ROOT / "run.py"),
        case,
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
        return {"command": command, "started": started, "finished": started, "exit_code": None, "log": str(log_path)}
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
        "abort_markers": sum(text.count(marker) for marker in ("STOP 1", "job aborted", "application aborted", "Could not allocate memory")),
        "hypre_solve_calls": len(iterations),
        "required_iterations_zero": iterations.count(0),
        "required_iterations_max": max(iterations, default=None),
        "required_iterations_sum": sum(iterations),
        "relative_change_samples": len(relative_changes),
        "relative_change_zero": sum(value == 0.0 for value in relative_changes),
        "relative_change_max": max(relative_changes, default=None),
    }


def compare(candidate: Path, reference: Path, candidate_label: str, reference_label: str, output: Path, dry_run: bool) -> dict:
    command = [
        sys.executable,
        str(ROOT / "scripts/analysis/phase24_difference_evidence.py"),
        "--candidate", str(candidate),
        "--reference", str(reference),
        "--candidate-label", candidate_label,
        "--reference-label", reference_label,
        "--out", str(output),
        "--checkpoint-us", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0",
    ]
    if dry_run:
        return {"command": command, "exit_code": None, "output": str(output)}
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    return {"command": command, "exit_code": process.returncode, "stdout": process.stdout, "stderr": process.stderr, "output": str(output)}


def write_summary(summary: dict) -> None:
    (DIAGNOSTIC_DIR / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = [
        "# Phase24 HYPRE linear-tolerance diagnostic",
        "",
        f"Generated: `{summary['generated_utc']}`",
        "",
        "## Fixed conditions",
        "",
        "- Same mesh, restart, pulse, inner circuit, Phase24 assembly, BDF2, nonlinear controls, HYPRE reuse, and adaptive preconditioner lagging",
        "- Same approximately `1 us` post-pulse window",
        "- Only `Linear System Convergence Tolerance` changes",
        "",
        "## Results",
        "",
        "| condition | exit | ALL DONE | HYPRE solves | zero-iteration solves | max iterations |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for key, options in VARIANTS.items():
        item = summary["variants"].get(key, {})
        audit = item.get("audit", {})
        lines.append(
            f"| {options['label']} | {item.get('run', {}).get('exit_code', 'n/a')} | "
            f"{audit.get('all_done', 'n/a')} | {audit.get('hypre_solve_calls', 'n/a')} | "
            f"{audit.get('required_iterations_zero', 'n/a')} | {audit.get('required_iterations_max', 'n/a')} |"
        )
    lines.extend([
        "",
        f"5e-7 vs CPU/MUMPS: `{summary.get('comparisons', {}).get('tol5e7_vs_cpu', {}).get('output', 'not run')}`",
        f"1e-10 vs CPU/MUMPS: `{summary.get('comparisons', {}).get('tol1e10_vs_cpu', {}).get('output', 'not run')}`",
        f"1e-10 vs 5e-7: `{summary.get('comparisons', {}).get('tol1e10_vs_tol5e7', {}).get('output', 'not run')}`",
        "",
        "Interpretation: if the strict tolerance moves the absolute current toward the direct/COMSOL value and changes the baseline-corrected waveform, linear convergence was insufficient. If it does not, investigate HYPRE vector transfer/scaling or the solver backend implementation.",
        "",
        f"Generated project: `{PROJECT_PATH}`",
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

    solver = absolute(args.solver)
    runtime_bin = absolute(args.runtime_bin)
    toolchain_bin = absolute(args.toolchain_bin)
    if not REFERENCE_SERIES.is_file():
        raise SystemExit(f"reference series not found: {REFERENCE_SERIES}")
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

    summary = {"generated_utc": datetime.now(timezone.utc).isoformat(), "project": str(project), "variants": {}, "comparisons": {}}
    for key, options in VARIANTS.items():
        case = options["case"]
        launcher_log = ROOT / "results" / case / "tolerance_diagnostic_launcher.log"
        solver_log = ROOT / "results" / case / "solver.log"
        series = ROOT / "results" / case / f"{case}_series.csv"
        command = command_for(case, solver, runtime_bin, toolchain_bin)
        print(f"[{key}] command: " + " ".join(command))
        run = run_case(command, launcher_log, args.dry_run)
        summary["variants"][key] = {
            "tolerance": options["tolerance"],
            "run": run,
            "solver_log": str(solver_log),
            "series": str(series),
            "audit": audit_log(solver_log) if not args.dry_run else {},
        }
        write_summary(summary)

    if not args.dry_run:
        tol5_series = Path(summary["variants"]["tol5e7"]["series"])
        tol10_series = Path(summary["variants"]["tol1e10"]["series"])
        if tol5_series.is_file() and tol10_series.is_file():
            summary["comparisons"]["tol5e7_vs_cpu"] = compare(tol5_series, REFERENCE_SERIES, "HYPRE tolerance 5e-7", "CPU/MUMPS Phase23", COMPARISON_ROOT / "phase24_hypre_tol5e7_vs_cpu_1us", False)
            summary["comparisons"]["tol1e10_vs_cpu"] = compare(tol10_series, REFERENCE_SERIES, "HYPRE tolerance 1e-10", "CPU/MUMPS Phase23", COMPARISON_ROOT / "phase24_hypre_tol1e10_vs_cpu_1us", False)
            summary["comparisons"]["tol1e10_vs_tol5e7"] = compare(tol10_series, tol5_series, "HYPRE tolerance 1e-10", "HYPRE tolerance 5e-7", COMPARISON_ROOT / "phase24_hypre_tol1e10_vs_tol5e7_1us", False)
        else:
            summary["comparisons"] = {"error": "one or both candidate series are missing"}
    else:
        summary["comparisons"] = {"dry_run": True}
    write_summary(summary)
    print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
    failures = [item for item in summary["variants"].values() if item["run"].get("exit_code") not in (0, None)]
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

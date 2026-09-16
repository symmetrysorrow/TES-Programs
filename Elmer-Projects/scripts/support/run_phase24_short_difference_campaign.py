"""Run the minimal Phase24 cause-isolation campaign and compare all outputs.

This campaign deliberately stops at approximately one microsecond after the
20.02 ms pulse.  It runs three otherwise identical native-HYPRE cases:

* BDF2 with the production HYPRE/preconditioner reuse policy;
* BDF2 with HYPRE and preconditioner reuse disabled;
* BDF1 with HYPRE and preconditioner reuse disabled.

Each case is generated from the checked-in Phase24 production case in a
separate campaign project.  The existing CPU/MUMPS Phase23 result is used as
the comparison reference.  No files are deleted; the generated project,
solver logs, manifests, CSV comparisons, and summary are retained under
``artifacts/phase24_short_difference_campaign`` and ``artifacts/comparison``.

The campaign is intentionally not run on import.  From the repository root:

    python scripts/support/run_phase24_short_difference_campaign.py

Use ``--dry-run`` to generate the cases and print the commands without
starting Elmer.  A single campaign normally takes minutes rather than the
roughly one-hour 1 ms run, but the exact time depends on HYPRE convergence.
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
SOURCE_PROJECT = ROOT / "elmer_project_phase24_production.json"
REFERENCE_SERIES = ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv"
CAMPAIGN_DIR = ROOT / "artifacts/phase24_short_difference_campaign"
PROJECT_PATH = CAMPAIGN_DIR / "phase24_short_campaign.json"
BASE_CASE = "case_phase24_hypre_cpu_pulse_1ms_5us"
CASES = {
    "bdf2_reuse_on": {
        "bdf_order": 2,
        "phase24_hypre_reuse": True,
        "phase24_preconditioner_lagging": "adaptive",
        "label": "BDF2 / HYPRE reuse ON",
    },
    "bdf2_reuse_off": {
        "bdf_order": 2,
        "phase24_hypre_reuse": False,
        "phase24_preconditioner_lagging": "disabled",
        "label": "BDF2 / reuse OFF",
    },
    "bdf1_reuse_off": {
        "bdf_order": 1,
        "phase24_hypre_reuse": False,
        "phase24_preconditioner_lagging": "disabled",
        "label": "BDF1 / reuse OFF",
    },
}


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def build_project() -> tuple[Path, dict[str, str]]:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    base = project["cases"][BASE_CASE]
    case_names: dict[str, str] = {}
    for suffix, options in CASES.items():
        name = f"case_phase24_short_{suffix}_1us"
        case_names[suffix] = name
        candidate = copy.deepcopy(base)
        candidate["timesteps"] = base["timesteps"][:5]
        candidate["output_intervals"] = base["output_intervals"][:5]
        candidate["series_file"] = f"{name}_series.csv"
        candidate["iteration_series_file"] = f"{name}_iterations.csv"
        candidate["output_result_path"] = None
        candidate["output_file_path"] = None
        candidate["comparison_time_grid"] = {
            "mode": "Phase24 short cause-isolation diagnostic",
            "purpose": "approximately 1 us post-pulse comparison",
            "post_pulse_end": "1[us]",
        }
        candidate["solver_comment"] = f"Short cause-isolation diagnostic: {options['label']}"
        candidate["phase24_hypre_reuse"] = options["phase24_hypre_reuse"]
        candidate["phase24_preconditioner_lagging"] = options["phase24_preconditioner_lagging"]
        candidate["bdf_order"] = options["bdf_order"]
        candidate["phase24_smoke"] = {
            "purpose": "short cause-isolation run; approximately 1 us after pulse",
            "variant": options["label"],
            "no_matrix_dump": True,
            "no_vtu": True,
        }
        project["cases"][name] = candidate
    CAMPAIGN_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_PATH.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    return PROJECT_PATH, case_names


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
    finished = datetime.now(timezone.utc).isoformat()
    return {
        "command": command,
        "started": started,
        "finished": finished,
        "elapsed_seconds": time.monotonic() - start,
        "exit_code": process.returncode,
        "log": str(log_path),
    }


def run_comparison(candidate: Path, label: str, output: Path, dry_run: bool) -> dict:
    command = [
        sys.executable,
        str(ROOT / "scripts/analysis/phase24_difference_evidence.py"),
        "--candidate", str(candidate),
        "--reference", str(REFERENCE_SERIES),
        "--candidate-label", label,
        "--reference-label", "CPU/MUMPS Phase23",
        "--out", str(output),
        "--checkpoint-us", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0",
    ]
    if dry_run:
        return {"command": command, "exit_code": None}
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    return {"command": command, "exit_code": process.returncode, "stdout": process.stdout, "stderr": process.stderr}


def write_summary(summary: dict) -> None:
    (CAMPAIGN_DIR / "campaign_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = [
        "# Phase24 short difference campaign",
        "",
        f"Generated: `{summary['generated_utc']}`",
        f"Project: `{PROJECT_PATH}`",
        "",
        "## Decision rule",
        "",
        "- If BDF2 reuse-OFF moves materially toward BDF1 reuse-OFF, reuse/lifecycle is implicated.",
        "- If BDF1 reuse-OFF matches the CPU reference while BDF2 reuse-OFF does not, BDF2/history handling is implicated.",
        "- If both HYPRE variants differ similarly from CPU/MUMPS, investigate Phase24 assembly, HYPRE tolerance/conditioning, or the runtime UDF before changing the time integrator.",
        "",
        "The CPU/MUMPS reference is not a bit-identical solver stack: it uses the historical Phase23 binary and BDF1. The campaign therefore localizes the cause; it is not itself a final numerical qualification.",
        "",
        "## Runs",
        "",
        "| variant | exit code | comparison |",
        "|---|---:|---|",
    ]
    for suffix, item in summary["variants"].items():
        comparison = item.get("comparison", {})
        comparison_dir = comparison.get("output", "n/a")
        lines.append(f"| {CASES[suffix]['label']} | {item.get('exit_code', 'n/a')} | `{comparison_dir}` |")
    lines.extend([
        "",
        "All solver output is retained in the corresponding `results/case_phase24_short_*_1us/` directory. Each comparison directory contains `summary.md`, `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, and `comparison.png` when plotting is available.",
        "",
    ])
    (CAMPAIGN_DIR / "campaign_summary.md").write_text("\n".join(lines), encoding="utf-8")


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

    project, case_names = build_project()
    generated_at = datetime.now(timezone.utc).isoformat()
    summary = {
        "generated_utc": generated_at,
        "project": str(project),
        "reference": str(REFERENCE_SERIES),
        "solver": str(solver),
        "runtime_bin": str(runtime_bin),
        "toolchain_bin": str(toolchain_bin),
        "dry_run": args.dry_run,
        "variants": {},
    }
    sync_command = [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(PROJECT_PATH)]
    if args.dry_run:
        print("[dry-run] " + " ".join(sync_command))
    else:
        sync = subprocess.run(sync_command, cwd=ROOT, capture_output=True, text=True)
        (CAMPAIGN_DIR / "sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
        if sync.returncode != 0:
            summary["sync_exit_code"] = sync.returncode
            write_summary(summary)
            return sync.returncode
    for suffix, options in CASES.items():
        case = case_names[suffix]
        log_path = ROOT / "results" / case / "campaign_launcher.log"
        command = command_for(case, solver, runtime_bin, toolchain_bin)
        print(f"[{suffix}] starting: {' '.join(command)}")
        run = run_case(command, log_path, args.dry_run)
        item = dict(run)
        if not args.dry_run:
            series = ROOT / "results" / case / f"{case}_series.csv"
            comparison_dir = ROOT / "artifacts" / "comparison" / f"phase24_{suffix}_1us"
            if series.is_file():
                comparison = run_comparison(series, options["label"], comparison_dir, args.dry_run)
                comparison["output"] = str(comparison_dir)
                item["comparison"] = comparison
            else:
                item["comparison"] = {"error": f"series not found: {series}"}
        else:
            item["comparison"] = run_comparison(ROOT / "results" / case / f"{case}_series.csv", options["label"], ROOT / "artifacts" / "comparison" / f"phase24_{suffix}_1us", True)
        summary["variants"][suffix] = item
        write_summary(summary)
        print(f"[{suffix}] exit_code={run['exit_code']}")

    failures = [item for item in summary["variants"].values() if item.get("exit_code") not in (0, None)]
    comparison_failures = [
        item for item in summary["variants"].values()
        if item.get("comparison", {}).get("exit_code") not in (0, None)
    ]
    missing = [item for item in summary["variants"].values() if "error" in item.get("comparison", {})]
    summary["campaign_exit_code"] = 1 if failures or comparison_failures or missing else 0
    write_summary(summary)
    print(f"campaign_summary={CAMPAIGN_DIR / 'campaign_summary.md'}")
    return summary["campaign_exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())

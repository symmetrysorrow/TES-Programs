"""Run the first pre-pulse Phase24 solve with HYPRE and direct MUMPS.

This is a short backend-parity diagnostic, not a production qualification.
It runs only the first ``18 us`` step before the pulse and limits each case to
one nonlinear iteration.  Both cases retain the Phase24 assembly, restart,
mesh, inner-circuit, BDF2, and timestep settings.  Only the linear backend
changes.  The assembled matrix and RHS are saved at the linear-solve entry
point, then compared automatically.

From the repository root::

    python scripts/support/run_phase24_baseline_backend_parity.py

Use ``--dry-run`` to generate the project and print commands without running
Elmer.  Existing output is never deleted; the two cases use dedicated names
and dump prefixes.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PROJECT = ROOT / "elmer_project_phase24_production.json"
BASE_CASE = "case_phase24_hypre_cpu_pulse_1ms_5us"
DIAGNOSTIC_DIR = ROOT / "artifacts/phase24_baseline_backend_parity_18us"
PROJECT_PATH = DIAGNOSTIC_DIR / "phase24_baseline_backend_parity.json"
DUMP_DIR = DIAGNOSTIC_DIR / "dumps"
VARIANTS = {
    "hypre": {
        "case": "case_phase24_baseline_firstsolve_hypre_18us",
        "linear_system": "iterative_hypre_flexgmres_boomeramg",
        "label": "Phase24 / HYPRE FlexGMRES + BoomerAMG",
    },
    "mumps": {
        "case": "case_phase24_baseline_firstsolve_mumps_18us",
        "linear_system": "mumps",
        "label": "Phase24 / direct MUMPS",
    },
}


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def build_project() -> Path:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    base = project["cases"][BASE_CASE]
    DUMP_DIR.mkdir(parents=True, exist_ok=True)
    for key, options in VARIANTS.items():
        candidate = copy.deepcopy(base)
        candidate["timesteps"] = base["timesteps"][:1]
        candidate["output_intervals"] = [1]
        candidate["series_file"] = f"{options['case']}_series.csv"
        candidate["iteration_series_file"] = f"{options['case']}_iterations.csv"
        candidate["output_result_path"] = None
        candidate["output_file_path"] = None
        candidate["solver_comment"] = (
            f"Diagnostic first pre-pulse solve: {options['label']}; one nonlinear iteration"
        )
        candidate["comparison_time_grid"] = {
            "mode": "Phase24 first pre-pulse backend parity",
            "purpose": "compare the first assembled matrix/RHS before the pulse",
            "step": "18[us]",
        }
        candidate["solver"]["linear_system"] = options["linear_system"]
        candidate["solver"]["nonlinear_max_iterations"] = 1
        candidate["solver"]["matrix_dump_prefix"] = str(
            (DUMP_DIR / f"{key}_firstsolve").resolve()
        )
        candidate["solver"]["matrix_dump_solution"] = False
        candidate["phase24_smoke"] = {
            "purpose": "first pre-pulse backend parity; one 18-us solve",
            "variant": options["label"],
            "no_matrix_dump": False,
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


def dump_path(prefix: Path, suffix: str) -> Path | None:
    direct = prefix.with_name(prefix.name + suffix)
    if direct.is_file():
        return direct
    rank_zero = prefix.with_name(prefix.name + suffix + ".0")
    return rank_zero if rank_zero.is_file() else None


def read_dump(path: Path, columns: int) -> dict[object, float]:
    values: dict[object, float] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.split()
            if len(fields) < columns:
                raise ValueError(f"{path}:{line_number}: expected {columns} fields")
            if columns == 3:
                key: object = (int(fields[0]), int(fields[1]))
                value = float(fields[2])
            else:
                key = int(fields[0])
                value = float(fields[1])
            values[key] = value
    return values


def compare_dump(left: Path | None, right: Path | None, columns: int) -> dict[str, object]:
    if left is None or right is None:
        return {"available": False, "left": str(left) if left else None, "right": str(right) if right else None}
    a = read_dump(left, columns)
    b = read_dump(right, columns)
    keys = set(a) | set(b)
    differences = [abs(a.get(key, 0.0) - b.get(key, 0.0)) for key in keys]
    relative = [
        abs(a.get(key, 0.0) - b.get(key, 0.0)) / max(abs(a.get(key, 0.0)), abs(b.get(key, 0.0)), 1.0e-300)
        for key in keys
    ]
    return {
        "available": True,
        "left": str(left),
        "right": str(right),
        "left_records": len(a),
        "right_records": len(b),
        "union_records": len(keys),
        "nonzero_difference_records": sum(value != 0.0 for value in differences),
        "max_absolute_difference": max(differences, default=0.0),
        "max_relative_difference": max(relative, default=0.0),
    }


def sha256(path: Path | None) -> str | None:
    if path is None:
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def first_series_row(case: str) -> dict[str, float] | None:
    path = ROOT / "results" / case / f"{case}_series.csv"
    if not path.is_file():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle), None)
    if row is None:
        return None
    return {
        "time_s": float(row["time_s"]),
        "tes_current_uA": float(row["tes_current_A"]) * 1.0e6,
    }


def audit_log(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {"log_exists": False}
    text = path.read_text(encoding="utf-8", errors="replace")
    return {
        "log_exists": True,
        "all_done": "ALL DONE" in text,
        "abort_markers": sum(text.count(marker) for marker in ("STOP 1", "job aborted", "application aborted", "Could not allocate memory")),
        "mumps_markers": text.count("MUMPS"),
        "hypre_markers": text.count("SolveHypre:"),
        "linear_system_save_markers": text.count("Saving matrix to:"),
    }


def write_summary(summary: dict) -> None:
    (DIAGNOSTIC_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = [
        "# Phase24 first pre-pulse backend parity",
        "",
        f"Generated: `{summary['generated_utc']}`",
        "",
        "## Fixed conditions",
        "",
        "- Same mesh, restart, pulse definition, inner circuit, Phase24 assembly, BDF2, and first `18 us` timestep",
        "- One nonlinear iteration per case; this is a matrix/RHS diagnostic, not a convergence qualification",
        "- Only linear backend differs: native HYPRE versus direct MUMPS",
        "",
        "## Cases",
        "",
        "| backend | exit | ALL DONE | matrix dump | first TES current [µA] |",
        "|---|---:|---:|---|---:|",
    ]
    for key, options in VARIANTS.items():
        item = summary["variants"].get(key, {})
        audit = item.get("audit", {})
        row = item.get("first_series_row") or {}
        prefix = item.get("dump_prefix", "")
        lines.append(
            f"| {options['label']} | {item.get('run', {}).get('exit_code', 'n/a')} | "
            f"{audit.get('all_done', 'n/a')} | `{prefix}` | {row.get('tes_current_uA', 'n/a')} |"
        )
    matrix = summary.get("matrix_comparison", {})
    lines.extend([
        "",
        "## Matrix/RHS comparison",
        "",
        f"- Matrix (`_a.dat`): `{json.dumps(matrix.get('matrix', {}), ensure_ascii=False)}`",
        f"- RHS (`_b.dat`): `{json.dumps(matrix.get('rhs', {}), ensure_ascii=False)}`",
        "",
        "Interpretation: a material matrix/RHS difference means the Phase24 assembly or constraint path differs. If they match but the first solved current differs, inspect backend solve/scaling and RHS-to-solution handling.",
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

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(project),
        "variants": {},
        "matrix_comparison": {},
    }
    for key, options in VARIANTS.items():
        case = options["case"]
        launcher_log = ROOT / "results" / case / "backend_parity_launcher.log"
        solver_log = ROOT / "results" / case / "solver.log"
        prefix = DUMP_DIR / f"{key}_firstsolve"
        command = command_for(case, solver, runtime_bin, toolchain_bin)
        print(f"[{key}] command: " + " ".join(command))
        run = run_case(command, launcher_log, args.dry_run)
        item = {
            "run": run,
            "solver_log": str(solver_log),
            "dump_prefix": str(prefix),
            "audit": audit_log(solver_log) if not args.dry_run else {},
            "first_series_row": first_series_row(case) if not args.dry_run else None,
            "dump_hashes": {},
        }
        for suffix in ("_a.dat", "_b.dat", "_sizes.dat"):
            item["dump_hashes"][suffix] = sha256(dump_path(prefix, suffix)) if not args.dry_run else None
        summary["variants"][key] = item
        write_summary(summary)

    if not args.dry_run:
        left = DUMP_DIR / "hypre_firstsolve"
        right = DUMP_DIR / "mumps_firstsolve"
        summary["matrix_comparison"] = {
            "matrix": compare_dump(dump_path(left, "_a.dat"), dump_path(right, "_a.dat"), 3),
            "rhs": compare_dump(dump_path(left, "_b.dat"), dump_path(right, "_b.dat"), 2),
        }
    else:
        summary["matrix_comparison"] = {"dry_run": True}
    write_summary(summary)
    print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

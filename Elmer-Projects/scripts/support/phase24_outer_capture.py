"""Prepare and materialize native outer-HeatSolve linear-system captures.

The native opt-in writes ``outer_before_*`` through the normal DefaultSolve
save hook and ``outer_after_*`` immediately after that *same* outer solve.
This tool only packages those raw files; it never treats a Schur, restricted,
or preconditioner dump as a HeatSolve system.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from pathlib import Path
from typing import Any


SIF_KEYS = (
    "Apply Mortar BCs",
    "Phase24 Vector Assembly",
    "Phase24 Static Matrix Reuse",
    "Phase24 Operator Lagging",
)


def read_sif(path: Path | None) -> dict[str, str]:
    if path is None or not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    values: dict[str, str] = {}
    for key in SIF_KEYS:
        match = re.search(rf'^\s*"?{re.escape(key)}"?\s*=\s*(.+?)\s*$', text, re.M | re.I)
        if match:
            values[key] = match.group(1).strip()
    return values


def matrix_stats(path: Path) -> dict[str, int]:
    """Return exact A-derived shape and primal/constraint row counts."""
    rows: set[int] = set()
    diagonal_rows: set[int] = set()
    max_index = 0
    records = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            try:
                row, column = int(fields[0]), int(fields[1])
            except ValueError:
                continue
            rows.add(row)
            if row == column:
                diagonal_rows.add(row)
            max_index = max(max_index, row, column)
            records += 1
    return {
        "total_rows": max_index,
        "matrix_records": records,
        "primal_rows": len(diagonal_rows),
        "constraint_rows": len(rows - diagonal_rows),
    }


def circuit_row(path: Path | None, timestep: int, nonlinear_iteration: int) -> dict[str, float]:
    if path is None or not path.is_file():
        return {}
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle):
            try:
                if int(float(row.get("time_step", "nan"))) != timestep:
                    continue
                if int(float(row.get("nonlinear_iter", "nan"))) != nonlinear_iteration:
                    continue
            except ValueError:
                continue
            result: dict[str, float] = {}
            for key in ("tes_temperature_K", "previous_current_A", "raw_current_A", "relaxed_power_W"):
                try:
                    result[key] = float(row[key])
                except (KeyError, TypeError, ValueError):
                    pass
            return result
    return {}


def capture_log_metadata(path: Path | None, timestep: int, nonlinear_iteration: int) -> dict[str, Any]:
    if path is None or not path.is_file():
        return {}
    wanted = f"PHASE24_OUTER_CAPTURE stage=before, timestep={timestep}, nonlinear_iteration={nonlinear_iteration},"
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if wanted not in line:
            continue
        out: dict[str, Any] = {}
        for key, value in re.findall(r"(time|dt|bdf_order|matrix_epoch|rhs_epoch)=([^,\s]+)", line):
            try:
                out[key] = int(value) if key.endswith("epoch") or key == "bdf_order" else float(value)
            except ValueError:
                out[key] = value
        return out
    return {}


def required_raw(directory: Path) -> dict[str, Path]:
    return {
        "A": directory / "outer_before_a.dat",
        "b": directory / "outer_before_b.dat",
        "x_before": directory / "outer_before_sol.dat",
        "x_after": directory / "outer_after_sol.dat",
    }


def prepare(root: Path, timestep: int, iterations: int) -> None:
    for nonlinear_iteration in range(1, iterations + 1):
        (root / f"ts{timestep:04d}_nl{nonlinear_iteration:04d}").mkdir(parents=True, exist_ok=True)
    print('  "Phase24 Outer HeatSolve Capture" = Logical True')
    print(f'  "Phase24 Outer HeatSolve Capture Prefix" = String "{root}"')
    print(f'  "Phase24 Outer HeatSolve Capture Timestep" = Integer {timestep}')
    print(f'  "Phase24 Outer HeatSolve Capture Max Iterations" = Integer {iterations}')


def materialize(root: Path, timestep: int, iterations: int, sif: Path | None, log: Path | None,
                iteration_csv: Path | None) -> int:
    failed = False
    for nonlinear_iteration in range(1, iterations + 1):
        directory = root / f"ts{timestep:04d}_nl{nonlinear_iteration:04d}"
        raw = required_raw(directory)
        missing = [name for name, path in raw.items() if not path.is_file()]
        if missing:
            print(f"missing {directory}: {', '.join(missing)}")
            failed = True
            continue
        for name, source in raw.items():
            shutil.copy2(source, directory / f"{name}.dat")
        metadata: dict[str, Any] = {
            "capture_kind": "outer_heat_solve",
            "timestep": timestep,
            "nonlinear_iteration": nonlinear_iteration,
            "raw_files": {name: source.name for name, source in raw.items()},
            "matrix": matrix_stats(raw["A"]),
            "mesh_interface": {"mortar": read_sif(sif).get("Apply Mortar BCs")},
            "phase24": {key: value for key, value in read_sif(sif).items() if key != "Apply Mortar BCs"},
            "solver": capture_log_metadata(log, timestep, nonlinear_iteration),
            "circuit": circuit_row(iteration_csv, timestep, nonlinear_iteration),
        }
        (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        print(directory / "metadata.json")
    return 1 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--root", type=Path, required=True)
    common.add_argument("--timestep", type=int, default=1)
    common.add_argument("--iterations", type=int, default=3)
    sub.add_parser("prepare", parents=[common])
    materialize_parser = sub.add_parser("materialize", parents=[common])
    materialize_parser.add_argument("--sif", type=Path)
    materialize_parser.add_argument("--solver-log", type=Path)
    materialize_parser.add_argument("--iteration-csv", type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.root, args.timestep, args.iterations)
        return 0
    return materialize(args.root, args.timestep, args.iterations, args.sif, args.solver_log, args.iteration_csv)


if __name__ == "__main__":
    raise SystemExit(main())

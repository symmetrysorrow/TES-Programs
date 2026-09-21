"""Materialize native ``SolveWithLinearRestriction`` saddle-point captures.

The native hook is intentionally tiny and writes only raw matrix/vector files.
This helper adds provenance, file hashes, and the runtime row split without
guessing the system dimension from filenames or iteration order.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any


SIF_KEYS = (
    "Linear System Solver",
    "Linear System Direct Method",
    "Eliminate Linear Constraints",
    "Penalty Linear Constraints",
    "Apply Mortar BCs",
    "Phase24 Full Restriction Capture",
    "Phase24 Full Restriction Capture Prefix",
    "Phase24 Vector Assembly",
    "Phase24 Static Matrix Reuse",
    "Phase24 Operator Lagging",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def parse_sizes(path: Path) -> dict[str, int]:
    values = [int(line.split()[0]) for line in path.read_text().splitlines() if line.split()]
    if len(values) < 4:
        raise ValueError(f"expected four size records in {path}")
    return {
        "total_rows": values[0],
        "matrix_nnz": values[1],
        "primal_rows": values[2],
        "constraint_matrix_rows": values[3],
        "constraint_rows": max(values[0] - values[2], 0),
    }


def matrix_stats(path: Path) -> dict[str, int]:
    records = 0
    max_index = 0
    with path.open(encoding="utf-8", errors="replace") as source:
        for line in source:
            fields = line.split()
            if len(fields) < 3:
                continue
            try:
                row, column = int(fields[0]), int(fields[1])
            except ValueError:
                continue
            records += 1
            max_index = max(max_index, row, column)
    return {"matrix_records": records, "max_index": max_index}


def capture_log_metadata(path: Path | None, timestep: int, iteration: int) -> dict[str, Any]:
    if path is None or not path.is_file():
        return {}
    prefix = (
        "PHASE24_FULL_RESTRICTION_CAPTURE stage=before, "
        f"timestep={timestep}, nonlinear_iteration={iteration},"
    )
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if prefix not in line:
            continue
        result: dict[str, Any] = {}
        for key, value in re.findall(r"(time|dt|bdf_order|solver|direct_method|total_rows|primal_rows|constraint_rows|matrix_nnz|constraint_matrix_rows|eliminate|penalty|restriction_active)=\s*([^,\s]+)", line):
            try:
                result[key] = int(value) if key in {
                    "bdf_order", "total_rows", "primal_rows", "constraint_rows",
                    "matrix_nnz", "constraint_matrix_rows", "eliminate", "penalty",
                    "restriction_active",
                } else float(value) if key in {"time", "dt"} else value
            except ValueError:
                result[key] = value
        return result
    return {}


def files_for(directory: Path) -> dict[str, Path]:
    return {
        "full_A_before": directory / "full_A_before.dat",
        "full_b_before": directory / "full_b_before.dat",
        "full_x_before": directory / "full_x_before.dat",
        "full_A_after": directory / "full_A_after.dat",
        "full_b_after": directory / "full_b_after.dat",
        "full_x_after": directory / "full_x_after.dat",
    }


def prepare(root: Path, timestep: int, iterations: int) -> None:
    for iteration in range(1, iterations + 1):
        (root / f"ts{timestep:04d}_nl{iteration:04d}").mkdir(parents=True, exist_ok=True)
    print('  "Phase24 Full Restriction Capture" = Logical True')
    print(f'  "Phase24 Full Restriction Capture Prefix" = String "{root}"')
    print(f'  "Phase24 Full Restriction Capture Timestep" = Integer {timestep}')
    print(f'  "Phase24 Full Restriction Capture Max Iterations" = Integer {iterations}')


def materialize(root: Path, timestep: int, iterations: int, sif: Path | None, log: Path | None) -> int:
    failed = False
    sif_values = read_sif(sif)
    for iteration in range(1, iterations + 1):
        directory = root / f"ts{timestep:04d}_nl{iteration:04d}"
        files = files_for(directory)
        missing = [name for name, path in files.items() if not path.is_file()]
        sizes_path = directory / "full_sizes_before.dat"
        if not sizes_path.is_file():
            missing.append("full_sizes_before")
        if missing:
            print(f"missing {directory}: {', '.join(missing)}")
            failed = True
            continue
        sizes = parse_sizes(sizes_path)
        before_stats = matrix_stats(files["full_A_before"])
        after_stats = matrix_stats(files["full_A_after"])
        event = capture_log_metadata(log, timestep, iteration)
        metadata: dict[str, Any] = {
            "capture_kind": "solve_with_linear_restriction_full_system",
            "capture_position": "SolveWithLinearRestriction after CollectionMatrix construction and immediately around SolveLinearSystem",
            "timestep": timestep,
            "nonlinear_iteration": iteration,
            "raw_files": {name: path.name for name, path in files.items()},
            "runtime": {
                **sizes,
                "before_matrix_records": before_stats["matrix_records"],
                "after_matrix_records": after_stats["matrix_records"],
            },
            "provenance": {
                "physical_time": event.get("time"),
                "dt": event.get("dt"),
                "bdf_order": event.get("bdf_order"),
                "solver_name": event.get("solver", sif_values.get("Linear System Solver")),
                "direct_method": event.get("direct_method", sif_values.get("Linear System Direct Method")),
                "mumps_backend": event.get("direct_method", sif_values.get("Linear System Direct Method")),
                "collection_matrix_construction_mode": "MATRIX_LIST assembled, then converted to CRS",
                "eliminate_linear_constraints": event.get("eliminate", 0),
                "penalty_linear_constraints": event.get("penalty", 0),
                "apply_mortar_bcs": sif_values.get("Apply Mortar BCs"),
                "restriction_mortar_active": event.get("restriction_active"),
                "log_event": event,
            },
            "sha256": {
                "matrix_before": sha256(files["full_A_before"]),
                "matrix_after": sha256(files["full_A_after"]),
                "rhs_before": sha256(files["full_b_before"]),
                "rhs_after": sha256(files["full_b_after"]),
            },
        }
        metadata["matrix_sha_before"] = metadata["sha256"]["matrix_before"]
        metadata["matrix_sha_after"] = metadata["sha256"]["matrix_after"]
        metadata["rhs_sha_before"] = metadata["sha256"]["rhs_before"]
        metadata["rhs_sha_after"] = metadata["sha256"]["rhs_after"]
        (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        print(directory / "metadata.json")
    return int(failed)


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
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.root, args.timestep, args.iterations)
        return 0
    return materialize(args.root, args.timestep, args.iterations, args.sif, args.solver_log)


if __name__ == "__main__":
    raise SystemExit(main())

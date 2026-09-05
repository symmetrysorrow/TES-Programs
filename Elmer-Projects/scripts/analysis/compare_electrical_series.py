"""Compare canonical TES electrical series artifacts."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


FIELDS = (
    "tes_temperature_K",
    "tes_current_A",
    "tes_resistance_ohm",
    "tes_power_W",
)


def read_series(path: Path | None) -> list[dict[str, float]] | None:
    if path is None or not path.exists():
        return None
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not set(FIELDS).issubset(rows[0]):
        return None
    return [
        {
            key: float(row[key])
            for key in ("time_s", "time_step", "nonlinear_iter", *FIELDS)
            if key in row and row[key] not in (None, "")
        }
        for row in rows
    ]


def compare(left_path: Path | None, right_path: Path | None) -> dict[str, object]:
    left = read_series(left_path)
    right = read_series(right_path)
    if left is None or right is None:
        return {"status": "NOT_AVAILABLE", "left": str(left_path) if left_path else None, "right": str(right_path) if right_path else None}
    by_step_left = {int(row.get("time_step", index)): row for index, row in enumerate(left)}
    by_step_right = {int(row.get("time_step", index)): row for index, row in enumerate(right)}
    steps = sorted(set(by_step_left) & set(by_step_right))
    if not steps:
        return {"status": "INCOMPLETE", "reason": "no common timestep rows"}
    report: dict[str, object] = {"status": "PASS", "common_timesteps": len(steps), "observables": {}}
    for field in FIELDS:
        errors = []
        relative = []
        for step in steps:
            a = by_step_left[step][field]
            b = by_step_right[step][field]
            errors.append((abs(a - b), step, by_step_left[step].get("time_s", float(step))))
            relative.append(abs(a - b) / max(abs(b), 1.0e-300))
        max_abs, max_step, max_time = max(errors, key=lambda item: item[0])
        report["observables"][field] = {
            "max_absolute_difference": max_abs,
            "max_relative_difference": max(relative),
            "rmse": math.sqrt(sum(error[0] ** 2 for error in errors) / len(errors)),
            "time_of_max_difference_s": max_time,
            "time_step_of_max_difference": max_step,
        }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", type=Path, required=True)
    parser.add_argument("--right", type=Path, required=True)
    parser.add_argument("--label", default="comparison")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {args.label: compare(args.left, args.right)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + chr(10), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

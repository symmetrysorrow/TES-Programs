"""Summarize nonlinear iteration counts per timestep for TES CSV logs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "time_step" not in rows[0] or "nonlinear_iter" not in rows[0]:
        raise ValueError(f"{path}: expected time_step and nonlinear_iter columns")
    per_step = Counter(int(row["time_step"]) for row in rows)
    histogram = Counter(per_step.values())
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "rows": len(rows),
        "timesteps": len(per_step),
        "iterations_per_step_min": min(per_step.values()),
        "iterations_per_step_max": max(per_step.values()),
        "iterations_per_step_histogram": {
            str(count): steps for count, steps in sorted(histogram.items())
        },
        "iterations_by_timestep": {
            str(step): count for step, count in sorted(per_step.items())
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("series", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summaries = {path.stem: summarize(path) for path in args.series}
    result: dict = {"series": summaries}
    if len(summaries) == 2:
        first, second = summaries.values()
        first_steps = first["iterations_by_timestep"]
        second_steps = second["iterations_by_timestep"]
        common = sorted(set(first_steps) & set(second_steps), key=int)
        result["steps_with_different_counts"] = [
            {
                "time_step": int(step),
                "reference": first_steps[step],
                "candidate": second_steps[step],
            }
            for step in common
            if first_steps[step] != second_steps[step]
        ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

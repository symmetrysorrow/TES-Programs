"""Summarize machine-readable Phase20 block-Schur probe output.

The native Phase20 source instrumentation writes one CSV row per outer block
preconditioner application and one row per inner Schur solve.  This utility
keeps the aggregation outside the hot loop and also accepts a log-only run,
which is useful with older Phase19 binaries.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any


SCHUR_RE = re.compile(r"Matrix-free Schur GMRES iterations:\s*(\d+)\s+residual:\s*([0-9.Ee+-]+)")
WARN_RE = re.compile(r"Inner Schur GMRES reached its limit")


def _float(row: dict[str, str], *names: str) -> float | None:
    for name in names:
        if row.get(name, "") != "":
            return float(row[name])
    return None


def summarize(outer_path: Path | None = None, schur_path: Path | None = None,
              log_path: Path | None = None) -> dict[str, Any]:
    outer: list[dict[str, str]] = []
    schur: list[dict[str, str]] = []
    if outer_path:
        with outer_path.open(newline="", encoding="utf-8") as handle:
            outer = list(csv.DictReader(handle))
    if schur_path:
        with schur_path.open(newline="", encoding="utf-8") as handle:
            schur = list(csv.DictReader(handle))

    log_solves: list[dict[str, Any]] = []
    log_warnings = 0
    if log_path:
        for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
            match = SCHUR_RE.search(line)
            if match:
                log_solves.append({"iterations": int(match.group(1)), "final_residual": float(match.group(2))})
            if WARN_RE.search(line):
                log_warnings += 1

    schur_rows = schur or log_solves
    k_actions = [int(float(row["k_actions"])) for row in schur if row.get("k_actions")]
    iterations = [int(float(row["iterations"])) for row in schur_rows if row.get("iterations") is not None]
    finals = [_float(row, "final_residual", "schur_final_residual") for row in schur_rows]
    finals = [value for value in finals if value is not None]
    outer_residuals = [_float(row, "outer_residual", "solver_reported_residual") for row in outer]
    outer_residuals = [value for value in outer_residuals if value is not None]
    wall = [_float(row, "elapsed_wall_seconds", "wall_seconds") for row in outer]
    wall = [value for value in wall if value is not None]

    first_residual = outer_residuals[0] if outer_residuals else None
    last_residual = outer_residuals[-1] if outer_residuals else None
    elapsed = max(wall) if wall else None
    reduction = (last_residual / first_residual) if first_residual and last_residual is not None else None
    return {
        "outer_rows": len(outer),
        "schur_solves": len(schur_rows),
        "schur_iteration_counts": iterations,
        "schur_reached_tolerance": sum(str(row.get("reached_tolerance", "")).lower() in {"true", "t", "1"} for row in schur_rows),
        "schur_hit_maxiter": sum(str(row.get("hit_maxiter", "")).lower() in {"true", "t", "1"} for row in schur_rows) or log_warnings,
        "schur_final_residual_min": min(finals) if finals else None,
        "schur_final_residual_max": max(finals) if finals else None,
        "k_actions_total": sum(k_actions) if k_actions else None,
        "k_actions_per_schur_solve": (sum(k_actions) / len(k_actions)) if k_actions else None,
        "outer_first_residual": first_residual,
        "outer_last_residual": last_residual,
        "outer_residual_reduction": reduction,
        "elapsed_wall_seconds": elapsed,
        "residual_reduction_per_second": ((1.0 - reduction) / elapsed if reduction is not None and elapsed else None),
        "residual_reduction_per_k_action": ((1.0 - reduction) / sum(k_actions) if reduction is not None and k_actions else None),
        "source": {"outer_csv": str(outer_path) if outer_path else None,
                   "schur_csv": str(schur_path) if schur_path else None,
                   "log": str(log_path) if log_path else None},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-csv", type=Path)
    parser.add_argument("--schur-csv", type=Path)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not any((args.outer_csv, args.schur_csv, args.log)):
        parser.error("at least one of --outer-csv, --schur-csv, or --log is required")
    result = summarize(args.outer_csv, args.schur_csv, args.log)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

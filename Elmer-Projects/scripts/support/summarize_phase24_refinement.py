"""Summarize the corrected native MUMPS refinement capture."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.support.analyze_phase24_saddle_sensitivity import indexed_vector, load_matrix


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    meta = json.loads((args.capture / "metadata.json").read_text(encoding="utf-8"))
    total = int(meta["runtime"]["total_rows"])
    primal = int(meta["runtime"]["primal_rows"])
    matrix = load_matrix(args.capture / "full_A_after.dat", (total, total))
    rhs = indexed_vector(args.capture / "full_b_after.dat", total)
    state = indexed_vector(args.capture / "full_x_after.dat", total)
    residual = np.asarray(matrix @ state - rhs).ravel()
    rows = list(csv.DictReader(args.trajectory.open(encoding="utf-8")))
    final = rows[-1]
    result = {
        "trajectory_rows": len(rows),
        "trajectory_first": rows[0],
        "trajectory_last": final,
        "final_refined_tes_temperature_K": float(final["tes_temperature_K"]),
        "final_refined_current_A": float(final["raw_current_A"]),
        "final_refined_current_uA": float(final["raw_current_A"]) * 1.0e6,
        "full_constrained_residual": {"l2": float(np.linalg.norm(residual)), "max_abs": float(np.max(np.abs(residual)))},
        "primal_residual": {"l2": float(np.linalg.norm(residual[:primal])), "max_abs": float(np.max(np.abs(residual[:primal])))},
        "constraint_residual": {"l2": float(np.linalg.norm(residual[primal:])), "max_abs": float(np.max(np.abs(residual[primal:])))},
        "dimensions": {"total_rows": total, "primal_rows": primal, "constraint_rows": total - primal},
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(rows), "current_uA": result["final_refined_current_uA"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

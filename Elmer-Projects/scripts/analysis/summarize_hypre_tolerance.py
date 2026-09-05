"""Summarize the current HYPRE tolerance sweep from solver logs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from summarize_conformal_runs import log_summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", required=True, help="label=solver.log")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    runs = {}
    for item in args.case:
        label, raw_path = item.split("=", 1)
        path = Path(raw_path)
        runs[label] = {"log": str(path.resolve()), **log_summary(path)}
    report = {
        "study_route": "conformal shared-node + external circuit_parallel UDF",
        "runs": runs,
        "decision": "REPORT_ONLY",
        "note": "The current 1e-8 result is not comparable to the earlier inactive-inner-circuit run; changing to circuit_parallel also changed the thermal source path.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Extract HYPRE FlexGMRES residual trajectory from print-level-3 output."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


# HYPRE print level 3 columns are:
# iteration, residual norm, convergence rate, relative residual norm.
LINE = re.compile(
    r"^\s*(\d+)\s+"
    r"([0-9.+-]+(?:[Ee][+-]?\d+)?)\s+"
    r"([0-9.+-]+(?:[Ee][+-]?\d+)?)\s+"
    r"([0-9.+-]+(?:[Ee][+-]?\d+)?)\s*$"
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("csv", type=Path)
    args = parser.parse_args()
    rows: list[tuple[int, float]] = []
    for line in args.log.read_text(encoding="utf-8", errors="replace").splitlines():
        match = LINE.match(line)
        if match:
            rows.append((int(match.group(1)), float(match.group(4))))
    if not rows:
        raise SystemExit("no HYPRE residual lines found")
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("iteration", "relative_residual"))
        writer.writerow((0, 1.0))
        writer.writerows(rows)
    print(f"extracted {len(rows)} iterations; final={rows[-1][1]:.17g}; csv={args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

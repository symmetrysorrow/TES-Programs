#!/usr/bin/env python3
"""Aggregate same-run PHITS PRE/POST additions by particle and JCOLL."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    groups: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    with args.csv.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            key = (row["ityp"], row["jcoll"], row["ncol"])
            groups[key].append(float(row["tally_increment_MeV"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("ityp", "jcoll", "ncol", "count", "positive", "zero", "negative", "sum_MeV", "max_MeV"))
        for (ityp, jcoll, ncol), values in sorted(groups.items()):
            writer.writerow((ityp, jcoll, ncol, len(values), sum(x > 0 for x in values),
                             sum(x == 0 for x in values), sum(x < 0 for x in values),
                             format(sum(values), ".17g"), format(max(values, default=0), ".17g")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

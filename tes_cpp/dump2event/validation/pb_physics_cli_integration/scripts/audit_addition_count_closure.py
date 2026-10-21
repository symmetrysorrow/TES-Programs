#!/usr/bin/env python3
"""Check the known 10k ledger's exclusive process counts and total."""
from __future__ import annotations
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

EXPECTED = {
    "egs5_electron_continuous": 1_901_831,
    "egs5_electron_boundary": 83,
    "egs5_electron_cutoff": 38_266,
    "dead_electron_secondary": 1_310,
    "egs5_positron_continuous": 1_795,
    "dead_positron_secondary": 9,
    "exact_zero": 11_448,
}
TOTAL = 1_954_742


def audit(path: Path) -> dict:
    counts: Counter[str] = Counter()
    total = negative = unresolved = 0
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            total += 1
            counts[row["process"]] += 1
            value = row.get("E_deposit_MeV", "")
            negative += int(bool(value) and float(value) < 0)
            unresolved += int(row.get("status") != "reconstructed")
    result = {"total": total, "counts": dict(sorted(counts.items())),
              "negative": negative, "unresolved": unresolved}
    if total != TOTAL or dict(counts) != EXPECTED or negative or unresolved:
        raise ValueError(f"addition closure failed: {json.dumps(result, sort_keys=True)}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.ledger), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

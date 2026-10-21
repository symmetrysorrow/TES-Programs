#!/usr/bin/env python3
"""Import the EXACT previously verified Pb HATCH CSV into the repo data folder.

This utility NEVER downloads, approximates, fabricates, changes, or commits
PHITS/EGS5 coefficients. Its distribution-permission switch is an explicit
user attestation, not a legal determination or license grant.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import shutil

EXPECTED_SHA = "20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9"
EXPECTED_SIZE = 130754
HEADER = [
    "medium", "medname", "lelke", "eke0", "eke1", "ededx0", "ededx1",
    "pdedx0", "pdedx1", "phits_density", "egs_density",
]
DEST = Path(__file__).resolve().parents[1] / "data" / "pb_egs5_phits3341.csv"


def validate_original(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"verified historical HATCH CSV was not found: {path}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA:
        raise ValueError(f"wrong HATCH table SHA-256: {digest}; expected {EXPECTED_SHA}")
    if path.stat().st_size != EXPECTED_SIZE:
        raise ValueError(f"verified HATCH CSV byte size mismatch: {path.stat().st_size}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != HEADER:
            raise ValueError("HATCH CSV columns do not match PHITS 3.341 Pb probe")
        rows = list(reader)
    lead = [row for row in rows if int(row["medium"]) == 2]
    if len(rows) != 600 or len(lead) != 150:
        raise ValueError("expected 600 coefficient rows total and 150 medium-2 Pb rows")
    indices = {int(row["lelke"]) for row in lead}
    if indices != set(range(1, 151)) or any(
        row["medname"].upper() != "PB" or float(row["phits_density"]) != 9.78
        or float(row["egs_density"]) != 9.78 for row in lead
    ):
        raise ValueError("medium 2 is not the validated Pb/HATCH density and index setup")
    return digest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="full historical 600-row HATCH CSV")
    parser.add_argument("--confirm-public-redistribution-permitted", action="store_true",
                        help="confirm that the source data can be committed to this PUBLIC repository")
    opts = parser.parse_args()
    digest = validate_original(opts.source)
    print(f"Valid Pb HATCH table, SHA-256={digest}")
    if not opts.confirm_public_redistribution_permitted:
        print("NOT copying: permission to redistribute PHITS/EGS5 coefficients publicly was not confirmed.")
        print("If permission is granted, repeat with --confirm-public-redistribution-permitted.")
        return
    if DEST.exists():
        if validate_original(DEST) != digest:
            raise ValueError("destination contains different coefficient bytes")
        print(f"Already present: {DEST}")
        return
    DEST.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(opts.source, DEST)
    validate_original(DEST)
    print(f"Staged locally (not git added): {DEST}")
    print("Review redistribution permission and inspect git diff before git add / commit.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Import the EXACT previously verified Pb HATCH CSV into the repo data folder.

This utility NEVER downloads, approximates, fabricates, changes, or commits
PHITS/EGS5 coefficients. Private-repository mode checks the live GitHub
repository visibility via gh before copying. Private visibility is NOT itself
a redistribution license; review applicable PHITS/EGS5 terms separately.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess

EXPECTED_SHA = "20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9"
EXPECTED_SIZE = 130754
HEADER = [
    "medium", "medname", "lelke", "eke0", "eke1", "ededx0", "ededx1",
    "pdedx0", "pdedx1", "phits_density", "egs_density",
]
REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = (REPO_ROOT / "dump2event_physics_validation" /
                  "egs5_stopping_power_reconstruction" / "runs" /
                  "n100_stopping_probe_v2" / "egs5_stopping_power_coefficients.csv")
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



def require_private_github_repository() -> None:
    """Refuse private-only import until GitHub confirms THIS checkout is private.

    A private repository is an access-control measure, not a legal license
    to redistribute restricted source/derived data.
    """
    try:
        check = subprocess.run(
            ["gh", "repo", "view", "--json", "isPrivate", "--jq", ".isPrivate"],
            cwd=REPO_ROOT, capture_output=True, text=True, check=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise RuntimeError(
            "Cannot confirm GitHub privacy. Authenticate GitHub CLI (gh auth login), "
            "then verify 'gh repo view --json isPrivate --jq .isPrivate'. "
            "Refusing to stage HATCH data."
        ) from error
    if check.stdout.strip().lower() != "true":
        raise RuntimeError(
            "GitHub reports that this checkout's repository is NOT private. "
            "Refusing to stage HATCH data in an accessible repository."
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source", nargs="?", type=Path, default=DEFAULT_SOURCE,
        help="original 600-row HATCH CSV (default: archived main checkout validation file)",
    )
    choices = parser.add_mutually_exclusive_group()
    choices.add_argument(
        "--confirm-private-repository", action="store_true",
        help="verify private visibility of the current GitHub checkout before importing",
    )
    choices.add_argument(
        "--confirm-public-redistribution-permitted", action="store_true",
        help="attest independently established permission for public redistribution",
    )
    opts = parser.parse_args()
    digest = validate_original(opts.source)
    print(f"Valid historical Pb HATCH table, SHA-256={digest}")
    if opts.confirm_private_repository:
        require_private_github_repository()
    elif not opts.confirm_public_redistribution_permitted:
        print("NOT copying: GitHub privacy or public redistribution permission is unverified.")
        print("For your PRIVATE GitHub repository, rerun with --confirm-private-repository.")
        print("The importer does not determine PHITS/EGS5 redistribution rights.")
        return
    if DEST.exists():
        if validate_original(DEST) != digest:
            raise ValueError("destination contains different coefficient bytes")
        print(f"Already present: {DEST}")
        return
    DEST.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(opts.source, DEST)
    validate_original(DEST)
    print(f"Verified CSV copied locally (not committed or pushed): {DEST}")
    print("Review license constraints, git diff, and remote visibility before git add/commit.")


if __name__ == "__main__":
    main()

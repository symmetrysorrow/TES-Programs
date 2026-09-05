"""Wait for the tail1000 1.25 us run, then join it to the validated tail."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CASE = "case_tail1000_fine_time_1p25_alg5qg1_pulse"
RESULT_PREFIX = ROOT / "results" / CASE / f"{CASE}_series.csv"
TAIL = (
    ROOT
    / "artifacts"
    / "comparison"
    / "comsol_cpu_singlepixel_prod_v2_hybrid_full"
    / "elmer_series_cleaned.csv"
)
MERGE = ROOT / "scripts" / "analysis" / "merge_singlepixel_tail.py"
COMPARE = ROOT / "scripts" / "analysis" / "compare_singlepixel_amgx_comsol.py"
OUT_DIR = ROOT / "artifacts" / "comparison" / "comsol_tail1000_fine_time_1p25_alg5qg1"
MERGED = OUT_DIR / f"{CASE}_merged_series.csv"
COMSOL = ROOT / "docs" / "Single-Pixel.txt"
PULSE_S = 20.020e-3


def completed_prefix(path: Path) -> bool:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            return False
        last_time = max(float(row["time_s"]) for row in rows)
        # Elmer writes the last accepted output one 1.25 us step before the
        # nominal 1000 us endpoint for this interval layout.
        return (last_time - PULSE_S) * 1.0e6 >= 997.0
    except (OSError, KeyError, ValueError):
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poll-seconds", type=float, default=60.0)
    parser.add_argument("--no-wait", action="store_true")
    args = parser.parse_args()
    if args.poll_seconds <= 0.0:
        parser.error("--poll-seconds must be positive")

    while not completed_prefix(RESULT_PREFIX):
        if args.no_wait:
            raise RuntimeError(f"completed result not found: {RESULT_PREFIX}")
        print("waiting for completed result:", RESULT_PREFIX, flush=True)
        time.sleep(args.poll_seconds)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    merge_command = [
        sys.executable,
        str(MERGE),
        "--prefix",
        str(RESULT_PREFIX),
        "--tail",
        str(TAIL),
        "--out",
        str(MERGED),
        "--seam-us",
        "1000",
    ]
    subprocess.run(merge_command, cwd=ROOT, check=True)

    compare_command = [
        sys.executable,
        str(COMPARE),
        "--elmer",
        str(MERGED),
        "--comsol",
        str(COMSOL),
        "--out",
        str(OUT_DIR),
        "--end-us",
        "179980",
        "--solver-label",
        "CPU MUMPS (1.25 us prefix + validated tail)",
        "--timestep-description",
        "0.625 us prefix + 1.25 us through 998.751 us + validated tail",
    ]
    subprocess.run(compare_command, cwd=ROOT, check=True)
    print("comparison complete:", OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

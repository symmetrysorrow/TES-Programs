"""Join a completed fine-grid SinglePixel pulse trace to a validated tail.

The tail is joined in CSV space because both inputs use the same mesh and
physical parameters.  The prefix must contain samples through ``--seam-us``;
tail samples at or before the seam are discarded to avoid duplicate time
stamps from the coarse tail grid.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PREFIX = (
    ROOT
    / "results"
    / "case_tail1000_fine_time_1p25_alg5qg1_pulse"
    / "case_tail1000_fine_time_1p25_alg5qg1_pulse_series.csv"
)
DEFAULT_TAIL = (
    ROOT
    / "artifacts"
    / "comparison"
    / "comsol_cpu_singlepixel_prod_v2_hybrid_full"
    / "elmer_series_cleaned.csv"
)
DEFAULT_OUT = (
    ROOT
    / "artifacts"
    / "comparison"
    / "comsol_tail1000_fine_time_1p25_alg5qg1"
    / "case_tail1000_fine_time_1p25_alg5qg1_merged_series.csv"
)
PULSE_S = 20.020e-3
TIME_TOLERANCE_S = 1.0e-10


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "time_s" not in reader.fieldnames:
            raise ValueError(f"{path} has no time_s CSV column")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} is empty")
    return list(reader.fieldnames), rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, default=DEFAULT_PREFIX)
    parser.add_argument("--tail", type=Path, default=DEFAULT_TAIL)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--seam-us",
        type=float,
        default=1000.0,
        help="absolute time from the pulse at which the fine prefix ends",
    )
    args = parser.parse_args()
    if args.seam_us <= 0.0:
        parser.error("--seam-us must be positive")
    for path in (args.prefix, args.tail):
        if not path.is_file():
            raise FileNotFoundError(path)

    prefix_fields, prefix_rows = read_csv(args.prefix)
    tail_fields, tail_rows = read_csv(args.tail)
    if prefix_fields != tail_fields:
        raise ValueError(
            "prefix and tail CSV columns differ:\n"
            f"prefix={prefix_fields}\n"
            f"tail={tail_fields}"
        )

    seam_s = PULSE_S + args.seam_us * 1.0e-6
    prefix_rows = sorted(prefix_rows, key=lambda row: float(row["time_s"]))
    tail_rows = sorted(tail_rows, key=lambda row: float(row["time_s"]))
    prefix_before_seam = [
        row for row in prefix_rows if float(row["time_s"]) <= seam_s + TIME_TOLERANCE_S
    ]
    tail_after_seam = [
        row for row in tail_rows if float(row["time_s"]) > seam_s + TIME_TOLERANCE_S
    ]
    if not prefix_before_seam:
        raise ValueError("prefix has no samples at or before the requested seam")
    prefix_end = float(prefix_before_seam[-1]["time_s"])
    if prefix_end < seam_s - 2.0e-6:
        raise ValueError(
            f"prefix ends at {(prefix_end - PULSE_S) * 1e6:.6f} us after pulse; "
            f"requested seam is {args.seam_us:.6f} us"
        )
    if not tail_after_seam:
        raise ValueError("tail has no samples after the requested seam")

    merged_rows = prefix_before_seam + tail_after_seam
    merged_times = [float(row["time_s"]) for row in merged_rows]
    if any(right <= left for left, right in zip(merged_times, merged_times[1:])):
        raise ValueError("merged time column is not strictly increasing")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=prefix_fields)
        writer.writeheader()
        writer.writerows(merged_rows)

    metadata = {
        "prefix": str(args.prefix),
        "prefix_sha256": sha256(args.prefix),
        "tail": str(args.tail),
        "tail_sha256": sha256(args.tail),
        "pulse_time_s": PULSE_S,
        "seam_after_pulse_us": args.seam_us,
        "prefix_end_after_pulse_us": (prefix_end - PULSE_S) * 1.0e6,
        "tail_first_after_pulse_us": (float(tail_after_seam[0]["time_s"]) - PULSE_S)
        * 1.0e6,
        "prefix_rows_used": len(prefix_before_seam),
        "tail_rows_used": len(tail_after_seam),
        "merged_rows": len(merged_rows),
        "merged_start_time_s": merged_times[0],
        "merged_end_time_s": merged_times[-1],
    }
    args.out.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

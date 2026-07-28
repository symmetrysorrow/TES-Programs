"""Compare same-grid serial and MPI TES current series.

The comparison is intentionally stricter than peak-only regression: it checks
the operating point, pulse height and timing, plus the largest pointwise
current difference with and without removal of the baseline offset.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


TIME = "time_s"
CURRENT = "tes_current_A"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_series(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise ValueError(f"{path}: series is empty")
    missing = {TIME, CURRENT} - rows[0].keys()
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if any(rows[i][TIME] >= rows[i + 1][TIME] for i in range(len(rows) - 1)):
        raise ValueError(f"{path}: time_s must be strictly increasing")
    return rows


def mean(values: list[float]) -> float:
    if not values:
        raise ValueError("comparison window contains no samples")
    return math.fsum(values) / len(values)


def compare(
    reference_path: Path,
    candidate_path: Path,
    *,
    pulse_time: float = 0.020020,
    baseline_start: float = 0.019500,
    baseline_end: float = 0.020020,
    baseline_limit: float = 0.01,
    height_limit: float = 0.02,
    peak_time_limit: float = 10.0e-6,
    max_difference_limit: float = 0.02,
) -> dict:
    reference = read_series(reference_path)
    candidate = read_series(candidate_path)
    time_key = lambda value: round(value, 12)
    candidate_by_time = {time_key(row[TIME]): row for row in candidate}
    paired = [
        (row, candidate_by_time[time_key(row[TIME])])
        for row in reference
        if time_key(row[TIME]) in candidate_by_time
    ]
    if not paired:
        raise ValueError("series have no matching timestamps")
    first_common = paired[0][0][TIME]
    last_common = paired[-1][0][TIME]
    expected_times = [
        time_key(row[TIME])
        for row in reference
        if first_common <= row[TIME] <= last_common
    ]
    if [time_key(pair[0][TIME]) for pair in paired] != expected_times:
        raise ValueError("candidate is missing reference timestamps in the common window")

    baseline_pairs = [
        pair
        for pair in paired
        if baseline_start <= pair[0][TIME] <= baseline_end
    ]
    post_pairs = [pair for pair in paired if pair[0][TIME] >= pulse_time]
    if not baseline_pairs or len(post_pairs) < 3:
        raise ValueError("series does not cover the baseline and post-pulse windows")

    ref_baseline = mean([ref[CURRENT] for ref, _ in baseline_pairs])
    cand_baseline = mean([cand[CURRENT] for _, cand in baseline_pairs])
    ref_peak = min(post_pairs, key=lambda pair: pair[0][CURRENT])
    cand_peak = min(post_pairs, key=lambda pair: pair[1][CURRENT])
    ref_height = ref_baseline - ref_peak[0][CURRENT]
    cand_height = cand_baseline - cand_peak[1][CURRENT]
    if ref_height <= 0.0:
        raise ValueError("reference has no positive pulse height")
    if ref_peak[0][TIME] == last_common or cand_peak[1][TIME] == last_common:
        raise ValueError("comparison ends at a current minimum; extend the run")

    raw_pair = max(
        paired, key=lambda pair: abs(pair[1][CURRENT] - pair[0][CURRENT])
    )
    corrected_pair = max(
        paired,
        key=lambda pair: abs(
            (pair[1][CURRENT] - cand_baseline)
            - (pair[0][CURRENT] - ref_baseline)
        ),
    )
    values = {
        "reference_baseline_A": ref_baseline,
        "candidate_baseline_A": cand_baseline,
        "reference_height_A": ref_height,
        "candidate_height_A": cand_height,
        "reference_peak_time_s": ref_peak[0][TIME],
        "candidate_peak_time_s": cand_peak[1][TIME],
        "raw_max_time_s": raw_pair[0][TIME],
        "raw_max_reference_A": raw_pair[0][CURRENT],
        "raw_max_candidate_A": raw_pair[1][CURRENT],
        "baseline_corrected_max_time_s": corrected_pair[0][TIME],
    }
    metric_values = {
        "baseline_current": abs(cand_baseline - ref_baseline) / abs(ref_baseline),
        "pulse_height": abs(cand_height - ref_height) / ref_height,
        "peak_time": abs(cand_peak[1][TIME] - ref_peak[0][TIME]),
        "raw_max_current": abs(
            raw_pair[1][CURRENT] - raw_pair[0][CURRENT]
        ) / ref_height,
        "baseline_corrected_max_current": abs(
            (corrected_pair[1][CURRENT] - cand_baseline)
            - (corrected_pair[0][CURRENT] - ref_baseline)
        ) / ref_height,
    }
    limits = {
        "baseline_current": baseline_limit,
        "pulse_height": height_limit,
        "peak_time": peak_time_limit,
        "raw_max_current": max_difference_limit,
        "baseline_corrected_max_current": max_difference_limit,
    }
    metrics = {
        name: {
            "value": value,
            "limit": limits[name],
            "passed": value <= limits[name],
        }
        for name, value in metric_values.items()
    }
    return {
        "inputs": {
            "reference": {
                "path": str(reference_path.resolve()),
                "sha256": sha256(reference_path),
            },
            "candidate": {
                "path": str(candidate_path.resolve()),
                "sha256": sha256(candidate_path),
            },
        },
        "comparison_window_s": [first_common, last_common],
        "matched_samples": len(paired),
        "baseline_samples": len(baseline_pairs),
        "values": values,
        "metrics": metrics,
        "overall_passed": all(metric["passed"] for metric in metrics.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="serial reference series CSV")
    parser.add_argument("candidate", type=Path, help="MPI candidate series CSV")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pulse-time", type=float, default=0.020020)
    parser.add_argument("--baseline-start", type=float, default=0.019500)
    parser.add_argument("--baseline-end", type=float, default=0.020020)
    args = parser.parse_args()
    result = compare(
        args.reference,
        args.candidate,
        pulse_time=args.pulse_time,
        baseline_start=args.baseline_start,
        baseline_end=args.baseline_end,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0 if result["overall_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

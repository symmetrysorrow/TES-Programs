"""Check an MPI TES transient series against the frozen direct-solver series."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


CURRENT = "tes_current_A"
TIME = "time_s"


def read_series(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise ValueError(f"{path}: series is empty")
    if any(rows[i][TIME] >= rows[i + 1][TIME] for i in range(len(rows) - 1)):
        raise ValueError(f"{path}: time_s must be strictly increasing")
    return rows


def mean(values: list[float]) -> float:
    if not values:
        raise ValueError("comparison window contains no samples")
    return math.fsum(values) / len(values)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "candidate",
        type=Path,
        help="MPI series CSV",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("artifacts/series/tes_pulse_20ms_3x_series.csv"),
    )
    parser.add_argument("--pulse-time", type=float, default=0.020020)
    parser.add_argument("--baseline-start", type=float, default=0.019500)
    parser.add_argument("--baseline-end", type=float, default=0.020020)
    args = parser.parse_args()

    reference = read_series(args.reference)
    candidate = read_series(args.candidate)
    # Elmer's repeated floating-point timestep accumulation can differ in the
    # last printed digits between serial and MPI runs.  A picosecond key keeps
    # all configured timestep stages distinct while ignoring that noise.
    time_key = lambda value: round(value, 12)
    candidate_by_time = {time_key(row[TIME]): row for row in candidate}
    paired = [
        (row, candidate_by_time[time_key(row[TIME])])
        for row in reference
        if time_key(row[TIME]) in candidate_by_time
    ]
    if not paired:
        raise ValueError("reference and MPI series have no exactly matching times")

    first_common = paired[0][0][TIME]
    last_common = paired[-1][0][TIME]
    expected_reference_times = [
        row[TIME]
        for row in reference
        if first_common <= row[TIME] <= last_common
    ]
    common_times = [time_key(ref[TIME]) for ref, _ in paired]
    if common_times != [time_key(value) for value in expected_reference_times]:
        raise ValueError("MPI series is missing reference timestamps inside the comparison window")

    baseline_pairs = [
        pair
        for pair in paired
        if args.baseline_start <= pair[0][TIME] <= args.baseline_end
    ]
    post_pairs = [pair for pair in paired if pair[0][TIME] >= args.pulse_time]
    if not baseline_pairs or len(post_pairs) < 3:
        raise ValueError("series does not cover the baseline and post-pulse windows")

    reference_baseline = mean([ref[CURRENT] for ref, _ in baseline_pairs])
    candidate_baseline = mean([mpi[CURRENT] for _, mpi in baseline_pairs])
    reference_peak_pair = min(post_pairs, key=lambda pair: pair[0][CURRENT])
    candidate_peak_pair = min(post_pairs, key=lambda pair: pair[1][CURRENT])
    reference_peak_time = reference_peak_pair[0][TIME]
    candidate_peak_time = candidate_peak_pair[1][TIME]
    if reference_peak_time == last_common or candidate_peak_time == last_common:
        raise ValueError(
            "comparison ends at a current minimum; extend the MPI run beyond the peak"
        )

    reference_height = reference_baseline - reference_peak_pair[0][CURRENT]
    candidate_height = candidate_baseline - candidate_peak_pair[1][CURRENT]
    if reference_height <= 0.0:
        raise ValueError("frozen reference has no positive pulse height")

    baseline_error = abs(candidate_baseline - reference_baseline) / abs(reference_baseline)
    height_error = abs(candidate_height - reference_height) / reference_height
    peak_time_error = abs(candidate_peak_time - reference_peak_time)
    max_current_error = max(
        abs(mpi[CURRENT] - ref[CURRENT]) for ref, mpi in paired
    )
    normalized_max_error = max_current_error / reference_height

    checks = [
        ("baseline current", baseline_error, 0.01, "%"),
        ("pulse height", height_error, 0.02, "%"),
        ("peak time", peak_time_error, 10.0e-6, "s"),
        ("maximum current difference", normalized_max_error, 0.02, "% of reference height"),
    ]
    passed = True
    print(f"comparison window: {first_common:.12g} - {last_common:.12g} s")
    for label, value, limit, unit in checks:
        ok = value <= limit
        passed = passed and ok
        shown = 100.0 * value if unit.startswith("%") else value
        shown_limit = 100.0 * limit if unit.startswith("%") else limit
        print(
            f"{'PASS' if ok else 'FAIL'}  {label}: "
            f"{shown:.9g} {unit} (limit {shown_limit:.9g} {unit})"
        )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

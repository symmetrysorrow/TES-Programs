#!/usr/bin/env python3
"""Compare reconstructed per-history bins with PHITS output=deposit histogram."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def phits_histogram(path: Path) -> list[tuple[float, float, float]]:
    bins: list[tuple[float, float, float]] = []
    with path.open(encoding="ascii", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) != 4:
                continue
            try:
                low, high, score = map(float, fields[:3])
            except ValueError:
                continue
            if high > low:
                bins.append((low, high, score))
    return bins


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("phits_histogram", type=Path)
    ap.add_argument("history_comparison", type=Path,
                    help="CSV from compare_internal_history_totals.py")
    ap.add_argument("output_csv", type=Path)
    ap.add_argument("--source-normalization", type=float, required=True,
                    help="PHITS resc2 total source weight; 1000 for the validated run")
    args = ap.parse_args()
    if not math.isfinite(args.source_normalization) or args.source_normalization <= 0:
        ap.error("source normalization must be positive and finite")

    reconstructed: dict[int, int] = {}
    with args.history_comparison.open(newline="", encoding="utf-8-sig") as src:
        for row in csv.DictReader(src):
            if row["PHITS_bin"] != "":
                index = int(row["PHITS_bin"])
                reconstructed[index] = reconstructed.get(index, 0) + 1

    bins = phits_histogram(args.phits_histogram)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    mismatch_count = 0
    with args.output_csv.open("w", newline="", encoding="utf-8") as dst:
        fields = ["bin_index", "emin_MeV", "emax_MeV", "PHITS_score_per_source",
                  "PHITS_equivalent_histories", "reconstructed_histories",
                  "absolute_difference_histories", "relative_difference"]
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        for index, (low, high, score) in enumerate(bins):
            phits_count = score * args.source_normalization
            raw_count = reconstructed.get(index, 0)
            difference = raw_count - phits_count
            mismatch_count += abs(difference) > 1e-9
            writer.writerow({
                "bin_index": index, "emin_MeV": format(low, ".17g"),
                "emax_MeV": format(high, ".17g"),
                "PHITS_score_per_source": format(score, ".17g"),
                "PHITS_equivalent_histories": format(phits_count, ".17g"),
                "reconstructed_histories": raw_count,
                "absolute_difference_histories": format(abs(difference), ".17g"),
                "relative_difference": format(difference / phits_count, ".17g")
                if phits_count else "",
            })
    print(f"bins={len(bins)} mismatches={mismatch_count} "
          f"histories_in_range={sum(reconstructed.values())} "
          f"source_normalization={args.source_normalization:g}")
    return 0 if len(bins) == 150 and mismatch_count == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

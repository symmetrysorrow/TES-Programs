#!/usr/bin/env python3
"""Recheck the 1,000-history NCOL=10 omission using serialized PHITS bin edges.

Inputs are the committed derived per-history comparison, histogram-bin export,
and NCOL=10 boundary audit. This never infers history IDs from a histogram.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def bin_index(value: float, edges: list[tuple[float, float]]) -> int | None:
    for i, (low, high) in enumerate(edges):
        if value >= low and value < high:
            return i
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history_compare", type=Path)
    parser.add_argument("histogram_compare", type=Path)
    parser.add_argument("boundary_audit", type=Path)
    parser.add_argument("--source-normalization", type=float, default=1000.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    history_rows = read_csv(args.history_compare)
    histogram_rows = read_csv(args.histogram_compare)
    boundary_rows = read_csv(args.boundary_audit)
    if not history_rows or not histogram_rows or not math.isfinite(args.source_normalization) or args.source_normalization <= 0:
        parser.error("inputs must be nonempty and source normalization must be positive")

    edges = [(float(row["emin_MeV"]), float(row["emax_MeV"])) for row in histogram_rows]
    for i, (low, high) in enumerate(edges):
        if not high > low or (i and low != edges[i - 1][1]):
            parser.error("serialized PHITS bins must be increasing and contiguous")

    per_history = {int(row["NOCAS"]): row for row in history_rows}
    if len(per_history) != len(history_rows):
        parser.error("duplicate NOCAS in per-history comparison")
    corrections: dict[int, float] = collections.defaultdict(float)
    for row in boundary_rows:
        nocas = int(row["NOCAS"])
        if nocas not in per_history:
            parser.error(f"boundary audit NOCAS {nocas} absent from history comparison")
        corrections[nocas] += float(row["boundary_correction_MeV"])

    actual_counts = [0] * len(edges)
    reconstructed_counts = [0] * len(edges)
    old_raw_counts = [0] * len(edges)
    max_history_error = 0.0
    sum_abs_history_error = 0.0
    phits_values: list[float] = []
    reconstructed_values: list[float] = []
    positive_in_range = 0
    crossings: list[dict[str, object]] = []
    for nocas, row in per_history.items():
        actual = float(row["PHITS_internal_deposit_MeV"])
        reconstructed = float(row["Dumpall_reconstructed_MeV"])
        error = abs(actual - reconstructed)
        max_history_error = max(max_history_error, error)
        sum_abs_history_error += error
        phits_values.append(actual)
        reconstructed_values.append(reconstructed)
        if actual > 0.0:
            positive_in_range += 0 < actual < edges[-1][1]
        actual_bin = bin_index(actual, edges) if actual > 0.0 else None
        reconstructed_bin = bin_index(reconstructed, edges) if reconstructed > 0.0 else None
        reported_actual_bin = row["PHITS_bin"]
        reported_reconstructed_bin = row["reconstructed_bin"]
        if ("" if actual_bin is None else str(actual_bin)) != reported_actual_bin:
            parser.error(f"PHITS bin column does not match serialized edges at NOCAS {nocas}")
        if ("" if reconstructed_bin is None else str(reconstructed_bin)) != reported_reconstructed_bin:
            parser.error(f"reconstructed bin column does not match serialized edges at NOCAS {nocas}")
        if actual_bin is not None:
            actual_counts[actual_bin] += 1
        if reconstructed_bin is not None:
            reconstructed_counts[reconstructed_bin] += 1

        correction = corrections.get(nocas, 0.0)
        old_value = reconstructed - correction
        old_bin = bin_index(old_value, edges) if old_value > 0.0 else None
        if old_bin is not None:
            old_raw_counts[old_bin] += 1
        if correction and old_bin != actual_bin:
            crossings.append({
                "NOCAS": nocas,
                "omitted_NCOL10_MeV": format(correction, ".17g"),
                "old_raw_only_MeV": format(old_value, ".17g"),
                "PHITS_internal_MeV": format(actual, ".17g"),
                "old_bin": old_bin,
                "PHITS_bin": actual_bin,
            })

    phits_hist_l1 = 0.0
    reconstructed_hist_l1 = 0.0
    for i, row in enumerate(histogram_rows):
        equivalent = float(row["PHITS_score_per_source"]) * args.source_normalization
        phits_hist_l1 += abs(equivalent - actual_counts[i])
        reconstructed_hist_l1 += abs(equivalent - reconstructed_counts[i])

    old_l1 = sum(abs(old_raw_counts[i] - actual_counts[i]) for i in range(len(edges)))
    phits_total = math.fsum(phits_values)
    reconstructed_total = math.fsum(reconstructed_values)
    summary = {
        "histories": len(per_history),
        "serialized_bin_count": len(edges),
        "bin_range_MeV": [format(edges[0][0], ".17g"), format(edges[-1][1], ".17g")],
        "source_normalization": args.source_normalization,
        "PHITS_internal_total_MeV": format(phits_total, ".17g"),
        "Dumpall_reconstructed_total_MeV": format(reconstructed_total, ".17g"),
        "absolute_total_difference_MeV": format(abs(phits_total - reconstructed_total), ".17g"),
        "sum_absolute_history_difference_MeV": format(sum_abs_history_error, ".17g"),
        "maximum_history_absolute_difference_MeV": format(max_history_error, ".17g"),
        "positive_internal_histories_in_range": positive_in_range,
        "PHITS_histogram_vs_internal_history_L1": format(phits_hist_l1, ".17g"),
        "PHITS_histogram_vs_corrected_Dumpall_L1": format(reconstructed_hist_l1, ".17g"),
        "NCOL10_audited_records": len(boundary_rows),
        "NCOL10_omitted_energy_MeV": format(sum(corrections.values()), ".17g"),
        "pre_NCOL10_correction_histogram_L1": old_l1,
        "pre_correction_crossing_histories": crossings,
        "interpretation": "energy reconstruction difference; the old raw-only candidate omitted NCOL=10 EGS5 boundary loss",
    }
    payload = json.dumps(summary, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0 if phits_hist_l1 < 1e-9 and reconstructed_hist_l1 < 1e-9 and old_l1 == 6 else 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare a PHITS unit=3 deposit histogram with raw-only event ledger sums."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def read_histogram(path: Path):
    bins = []
    active = started = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.lstrip().startswith("h: n"):
            active = True
            continue
        if not active or not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        try:
            bins.append((float(fields[0]), float(fields[1]), float(fields[2])))
            started = True
        except (ValueError, IndexError):
            if started:
                break
    result = [(lo, hi, value) for lo, hi, value in bins if hi > lo]
    if not result:
        raise ValueError("no PHITS histogram data rows found")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tally", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--histories", type=int, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    totals: dict[int, float] = defaultdict(float)
    unresolved = Counter()
    with args.ledger.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            event = int(row["event_id"])
            if row["status"] != "reconstructed":
                unresolved[event] += 1
                continue
            totals[event] += float(row["E_deposit_MeV"])
    if not totals or min(totals) < 1 or max(totals) > args.histories:
        raise ValueError("ledger event IDs do not fit the requested history range")
    bins = read_histogram(args.tally)
    width = bins[0][1] - bins[0][0]
    hcounts = Counter()
    zero = under = over = 0
    values = []
    for event in range(1, args.histories + 1):
        energy = totals.get(event, 0.0)
        values.append(energy)
        if energy == 0.0:
            zero += 1
        elif energy < bins[0][0]:
            under += 1
        elif energy >= bins[-1][1]:
            over += 1
        else:
            index = math.floor((energy - bins[0][0]) / width)
            if index < 0:
                under += 1
            elif index >= len(bins):
                over += 1
            else:
                hcounts[index] += 1

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    differences = []
    with args.output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("bin", "low_MeV", "high_MeV", "PHITS_counts_unit3",
                         "raw_ledger_event_counts", "difference_counts", "relative_difference"))
        for index, (low, high, normalized_count) in enumerate(bins):
            phits_count = normalized_count * args.histories
            ledger_count = hcounts[index]
            diff = ledger_count - phits_count
            differences.append(abs(diff))
            relative = diff / phits_count if phits_count else (0.0 if ledger_count == 0 else None)
            writer.writerow((index, format(low, ".9g"), format(high, ".9g"),
                             format(phits_count, ".9g"), ledger_count, format(diff, ".9g"),
                             "" if relative is None else format(relative, ".9g")))

    mean = sum(values) / args.histories
    variance = sum((value - mean) ** 2 for value in values) / (args.histories - 1)
    peak = max(range(len(bins)), key=lambda i: hcounts[i])
    result = {
        "histories": args.histories,
        "unresolved_candidate_rows": sum(unresolved.values()),
        "events_with_unresolved_candidates": len(unresolved),
        "unresolved_event_ids": sorted(unresolved),
        "raw_ledger_sum_MeV": sum(values),
        "PHITS_unit3_binned_count_sum": sum(count * args.histories for _, _, count in bins),
        "raw_in_range_count": sum(hcounts.values()),
        "zero_histories": zero,
        "underflow_histories": under,
        "overflow_histories": over,
        "bin_width_MeV": width,
        "bin_count": len(bins),
        "bin_count_L1_difference": sum(differences),
        "max_abs_bin_count_difference": max(differences, default=0.0),
        "mean_supported_raw_deposit_MeV": mean,
        "sample_sd_MeV": math.sqrt(variance),
        "standard_error_MeV": math.sqrt(variance / args.histories),
        "peak_bin_low_MeV": bins[peak][0],
        "peak_bin_high_MeV": bins[peak][1],
        "maximum_supported_event_deposit_MeV": max(values),
    }
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare PHITS tdepstregEVENT PRE/POST trace sums with experimental ledger."""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("internal_trace_csv", type=Path,
                    help="output of parse_fortran_tally_trace.py")
    ap.add_argument("physical_ledger_csv", type=Path)
    ap.add_argument("comparison_csv", type=Path)
    ap.add_argument("summary_json", type=Path)
    ap.add_argument("--region-id", type=int, default=1000)
    ap.add_argument("--tally-unit", type=int, default=3)
    ap.add_argument("--emin", type=float, default=0.0)
    ap.add_argument("--emax", type=float, default=1.5)
    ap.add_argument("--bin-width", type=float, default=0.01)
    args = ap.parse_args()

    phits: dict[int, float] = collections.defaultdict(float)
    internal_additions = 0
    with args.internal_trace_csv.open(newline="", encoding="utf-8-sig") as src:
        for row in csv.DictReader(src):
            if int(row["tally_region_id"]) != args.region_id or int(row["tally_unit"]) != args.tally_unit:
                continue
            phits[int(row["nocas"])] += float(row["tlv_MeV"])
            internal_additions += 1

    reconstructed: dict[int, float] = collections.defaultdict(float)
    ledger_rows = 0
    unresolved_rows = 0
    with args.physical_ledger_csv.open(newline="", encoding="utf-8-sig") as src:
        for row in csv.DictReader(src):
            ledger_rows += 1
            if row["status"] != "reconstructed" or not row["E_deposit_MeV"]:
                unresolved_rows += 1
                continue
            reconstructed[int(row["event_id"])] += float(row["E_deposit_MeV"])

    histories = sorted(set(phits) | set(reconstructed))
    fields = ["NOCAS", "PHITS_internal_deposit_MeV", "Dumpall_reconstructed_MeV",
              "absolute_error_MeV", "PHITS_bin", "reconstructed_bin",
              "PHITS_boundary_distance_MeV", "reconstructed_boundary_distance_MeV"]
    args.comparison_csv.parent.mkdir(parents=True, exist_ok=True)
    max_error = 0.0
    internal_hist = collections.Counter()
    reconstructed_hist = collections.Counter()
    positive_histories = 0
    with args.comparison_csv.open("w", newline="", encoding="utf-8") as dst:
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        for history in histories:
            p = phits.get(history, 0.0)
            c = reconstructed.get(history, 0.0)
            error = abs(p - c)
            max_error = max(max_error, error)
            p_bin = c_bin = ""
            p_dist = c_dist = ""
            if p > 0.0 and args.emin <= p < args.emax:
                idx = math.floor((p - args.emin) / args.bin_width)
                p_bin = idx
                p_dist = min(p - (args.emin + idx * args.bin_width),
                             args.emin + (idx + 1) * args.bin_width - p)
                internal_hist[idx] += 1
                positive_histories += 1
            if c > 0.0 and args.emin <= c < args.emax:
                idx = math.floor((c - args.emin) / args.bin_width)
                c_bin = idx
                c_dist = min(c - (args.emin + idx * args.bin_width),
                             args.emin + (idx + 1) * args.bin_width - c)
                reconstructed_hist[idx] += 1
            writer.writerow({
                "NOCAS": history,
                "PHITS_internal_deposit_MeV": format(p, ".17g"),
                "Dumpall_reconstructed_MeV": format(c, ".17g"),
                "absolute_error_MeV": format(error, ".17g"),
                "PHITS_bin": p_bin,
                "reconstructed_bin": c_bin,
                "PHITS_boundary_distance_MeV": format(p_dist, ".17g") if p_dist != "" else "",
                "reconstructed_boundary_distance_MeV": format(c_dist, ".17g") if c_dist != "" else "",
            })

    bin_diffs = {str(i): reconstructed_hist[i] - internal_hist[i]
                 for i in set(internal_hist) | set(reconstructed_hist)
                 if reconstructed_hist[i] != internal_hist[i]}
    summary = {
        "histories": len(histories),
        "phits_tally_addition_count": internal_additions,
        "physical_ledger_row_count": ledger_rows,
        "unresolved_ledger_row_count": unresolved_rows,
        "positive_phits_histories_in_range": positive_histories,
        "positive_reconstructed_histories_in_range": sum(reconstructed_hist.values()),
        "phits_internal_sum_MeV": format(sum(phits.values()), ".17g"),
        "reconstructed_sum_MeV": format(sum(reconstructed.values()), ".17g"),
        "absolute_sum_difference_MeV": format(abs(sum(phits.values()) - sum(reconstructed.values())), ".17g"),
        "max_per_history_absolute_error_MeV": format(max_error, ".17g"),
        "per_history_rows_over_1e-12_MeV": sum(abs(phits.get(i, 0.0) - reconstructed.get(i, 0.0)) > 1e-12
                                                for i in histories),
        "histogram_L1_count_difference": sum(abs(reconstructed_hist[i] - internal_hist[i])
                                              for i in set(internal_hist) | set(reconstructed_hist)),
        "histogram_bin_count_differences": bin_diffs,
        "bin_width_MeV": args.bin_width,
        "range_MeV": [args.emin, args.emax],
    }
    args.summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))
    return 0 if unresolved_rows == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

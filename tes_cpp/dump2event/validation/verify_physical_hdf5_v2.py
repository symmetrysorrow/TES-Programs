#!/usr/bin/env python3
"""Read back HDF5 v2 and compare either internal addition rows or a PHITS histogram."""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
from pathlib import Path

import h5py
import numpy as np

from tes_cpp.event_hdf5_v2 import iter_deposits, iter_events


def xyz(text):
    if not text:
        return None
    values = tuple(float(x) for x in text.split(";"))
    return values if len(values) == 3 else None


def compare_internal(h5path: Path, reference: Path, output: Path):
    rows = {int(row["addition_id"]): (event_id, row) for event_id, row in iter_deposits(h5path)}
    if len(rows) != sum(1 for _ in iter_deposits(h5path)):
        raise ValueError("HDF5 addition_id is not unique")
    n = missing = extra = energy_bad = event_bad = particle_bad = region_bad = position_bad = 0
    maximum = max_position = 0.0
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["addition_id", "event_id", "particle_id", "region_id", "PHITS_tlv_MeV",
              "HDF5_energy_MeV", "absolute_energy_error_MeV", "kind",
              "HDF5_start_xyz_cm", "HDF5_end_xyz_cm", "HDF5_score_xyz_cm",
              "position_error_cm", "position_check"]
    with reference.open(encoding="utf-8-sig", newline="") as source, output.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        seen = set()
        for ref in csv.DictReader(source):
            key = int(ref["tally_addition_id"])
            if key in seen:
                raise ValueError(f"duplicate PHITS addition ID {key}")
            seen.add(key)
            packed = rows.get(key)
            if packed is None:
                missing += 1
                continue
            event_id, got = packed
            n += 1
            energy = float(got["energy_MeV"])
            delta = abs(float(ref["PHITS_tlv_MeV"]) - energy)
            maximum = max(maximum, delta)
            energy_bad += delta > 1e-12
            event_bad += int(ref["NOCAS"]) != event_id
            particle_bad += int(ref["particle_id"]) != int(got["particle_id"])
            region_bad += int(ref["internal_region_id"]) != int(got["region_id"])
            pe = None
            position_check = "not-applicable"
            if energy > 0 and got["kind"] == "SEGMENT":
                start_error = _maxerr(xyz(ref["internal_xyz"]), got["start_xyz_cm"])
                end_error = _maxerr(xyz(ref["internal_cxyz"]), got["end_xyz_cm"])
                pe = max(start_error, end_error)
                position_check = "segment endpoints"
            elif energy > 0 and got["kind"] == "POINT":
                internal_pos = xyz(ref["reconstructed_position"])
                pe = _maxerr(internal_pos, got["score_xyz_cm"])
                position_check = "raw child point vs internal point"
            if pe is not None:
                max_position = max(max_position, pe)
                position_bad += pe > 1e-12
            writer.writerow({"addition_id": key, "event_id": ref["NOCAS"],
                "particle_id": ref["particle_id"], "region_id": ref["internal_region_id"],
                "PHITS_tlv_MeV": ref["PHITS_tlv_MeV"], "HDF5_energy_MeV": format(energy, ".17g"),
                "absolute_energy_error_MeV": format(delta, ".17g"), "kind": got["kind"],
                "HDF5_start_xyz_cm": got["start_xyz_cm"], "HDF5_end_xyz_cm": got["end_xyz_cm"],
                "HDF5_score_xyz_cm": got["score_xyz_cm"],
                "position_error_cm": "" if pe is None else format(pe, ".17g"),
                "position_check": position_check})
    extra = len(set(rows) - seen)
    summary = {"PHITS_additions": len(seen), "HDF5_additions": len(rows), "matched": n,
               "missing": missing, "extra": extra, "energy_mismatch_gt_1e-12_MeV": energy_bad,
               "event_mismatch": event_bad, "particle_mismatch": particle_bad,
               "region_mismatch": region_bad, "position_mismatch_gt_1e-12_cm": position_bad,
               "max_energy_abs_error_MeV": maximum, "max_position_abs_error_cm": max_position,
               "passed": not any((missing, extra, energy_bad, event_bad, particle_bad, region_bad, position_bad))}
    return summary


def compare_histogram(h5path: Path, reference: Path, output: Path, normalization: int):
    event_totals = collections.defaultdict(float)
    event_counts = collections.defaultdict(int)
    for event_id, row in iter_deposits(h5path):
        if row["energy_MeV"] is not None:
            event_totals[event_id] += row["energy_MeV"]
            event_counts[event_id] += 1
    phits = []
    with reference.open(encoding="ascii", errors="replace") as stream:
        for line in stream:
            values = line.split()
            if len(values) != 4:
                continue
            try:
                low, high, score, _err = map(float, values)
            except ValueError:
                continue
            if high > low:
                phits.append((low, high, score * normalization))
    reconstructed = collections.Counter()
    for total in event_totals.values():
        if total > 0 and 0 <= total < 1.5:
            index = math.floor((total + 1e-12) / 0.01)
            reconstructed[index] += 1
    differences = []
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as target:
        fields = ["bin_index", "emin_MeV", "emax_MeV", "PHITS_equivalent_histories",
                  "HDF5_reconstructed_histories", "absolute_difference"]
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        for index, (low, high, phits_count) in enumerate(phits):
            delta = reconstructed[index] - phits_count
            differences.append(delta)
            writer.writerow({"bin_index": index, "emin_MeV": low, "emax_MeV": high,
                "PHITS_equivalent_histories": format(phits_count, ".17g"),
                "HDF5_reconstructed_histories": reconstructed[index],
                "absolute_difference": format(abs(delta), ".17g")})
    return {"histories_in_hdf5": len(event_totals), "deposits": sum(event_counts.values()),
            "phits_bins": len(phits), "histogram_L1_difference": sum(abs(x) for x in differences),
            "nonzero_bin_differences": sum(abs(x) > 1e-9 for x in differences),
            "PHITS_equivalent_histories": sum(x[2] for x in phits),
            "HDF5_in_range_positive_histories": sum(reconstructed.values()),
            "passed": len(phits) == 150 and all(abs(x) <= 1e-9 for x in differences)}


def _maxerr(a, b):
    if a is None or b is None:
        return math.inf
    return max(abs(float(x) - float(y)) for x, y in zip(a, b))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hdf5", type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--internal-additions", type=Path)
    mode.add_argument("--phits-histogram", type=Path)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--normalization", type=int, default=1000)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    # Force a full structural and offsets validation before comparisons.
    events = list(iter_events(args.hdf5))
    if args.internal_additions:
        result = compare_internal(args.hdf5, args.internal_additions, args.csv)
    else:
        result = compare_histogram(args.hdf5, args.phits_histogram, args.csv, args.normalization)
    result["hdf5_event_count"] = len(events)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

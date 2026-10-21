#!/usr/bin/env python3
"""Compare PHITS pre-histogram internal history sums with HDF5 v2 and tally bins.

The internal trace is a PHITS diagnostic oracle. It is not a raw-only input.
The report keeps incomplete HDF5 histories as lower bounds and applies the
half-open PHITS bins without any edge nudge.
"""
from __future__ import annotations

import argparse
import bisect
import csv
from collections import defaultdict
from pathlib import Path

import h5py

FIELDS = (
    "label nocas nobch rcasc particle_id ityp jcoll ncol material_id "
    "tally_region_index tally_region_id region_before region_after oldwt E EC ein ecc edep tlv "
    "tally_unit species_slot lattice_slot time_slot x_start y_start z_start x_end y_end z_end tally_before"
).split()


def read_hdf(path: Path):
    totals: dict[int, float] = defaultdict(float)
    unresolved: dict[int, int] = defaultdict(int)
    complete: dict[int, bool] = {}
    with h5py.File(path, "r") as h5:
        dep = h5["deposits"]
        for lo in range(0, len(dep["event_id"]), 100_000):
            hi = min(len(dep["event_id"]), lo + 100_000)
            ids = dep["event_id"][lo:hi]
            valid = dep["energy_valid"][lo:hi]
            energy = dep["energy_MeV"][lo:hi]
            kind = dep["kind"][lo:hi]
            for i, raw_id in enumerate(ids):
                eid = int(raw_id)
                if bool(valid[i]):
                    totals[eid] += float(energy[i])
                if int(kind[i]) == 3:
                    unresolved[eid] += 1
        for eid, is_complete in zip(h5["events/event_id"][:], h5["events/complete"][:]):
            complete[int(eid)] = bool(is_complete)
    return totals, unresolved, complete


def read_internal(path: Path):
    sums: dict[int, float] = defaultdict(float)
    counts: dict[int, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.startswith(" PRE"):
                continue
            fields = line.split()
            if len(fields) != len(FIELDS):
                raise ValueError(f"{path}:{line_no}: expected {len(FIELDS)} fields, found {len(fields)}")
            row = dict(zip(FIELDS, fields))
            if int(row["tally_region_id"]) != 1000:
                continue
            eid = int(row["nocas"])
            # Sequential IEEE binary64 accumulation, matching the PHITS trEVENT updates.
            value = float(row["tlv"])
            sums[eid] += value
            counts[eid] += 1
    return sums, counts


def read_tally(path: Path):
    rows = []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) == 4 and fields[0].endswith("E+00") or len(fields) == 4 and "E" in fields[0]:
                try:
                    low, high, val, _err = map(float, fields)
                except ValueError:
                    continue
                if high > low:
                    rows.append((low, high, val))
    if len(rows) != 150:
        raise ValueError(f"expected 150 bins from {path}, found {len(rows)}")
    return rows


def bin_index(value: float, bins):
    edges = [bins[0][0]] + [row[1] for row in bins]
    if value <= 0.0 or value < edges[0] or value >= edges[-1]:
        return None
    return bisect.bisect_right(edges, value) - 1


def write_csv(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--internal", type=Path, required=True)
    ap.add_argument("--hdf5", type=Path, required=True)
    ap.add_argument("--tally", type=Path, required=True)
    ap.add_argument("--maxcas", type=int, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    hdf, unresolved, complete = read_hdf(args.hdf5)
    internal, calls = read_internal(args.internal)
    bins = read_tally(args.tally)
    all_ids = sorted(set(hdf) | set(internal) | set(complete))
    detail = []
    internal_hist = [0] * len(bins)
    hdf_hist = [0] * len(bins)
    complete_hist = [0] * len(bins)
    incompletes = []
    for eid in all_ids:
        score = internal.get(eid, 0.0)
        subtotal = hdf.get(eid, 0.0)
        ib = bin_index(score, bins)
        hb = bin_index(subtotal, bins)
        if ib is not None:
            internal_hist[ib] += 1
        if hb is not None:
            hdf_hist[hb] += 1
        if hb is not None and complete.get(eid, False):
            complete_hist[hb] += 1
        if unresolved.get(eid, 0):
            incompletes.append((eid, score, subtotal, score - subtotal, calls.get(eid, 0), unresolved[eid], ib, hb))
        detail.append((eid, f"{score:.17g}", f"{subtotal:.17g}", f"{score-subtotal:.17g}",
                       calls.get(eid, 0), unresolved.get(eid, 0), int(complete.get(eid, False)),
                       "" if ib is None else ib, "" if hb is None else hb,
                       "" if ib is None or hb is None else int(ib == hb)))
    write_csv(args.outdir / "HISTORY_LEVEL_DISCREPANCIES.csv",
              ["NOCAS", "PHITS_internal_deposit_MeV", "HDF5_known_subtotal_MeV", "internal_minus_known_MeV",
               "internal_score_calls", "unresolved_rows", "HDF5_complete", "PHITS_internal_bin", "HDF5_known_bin",
               "same_bin"], detail)
    write_csv(args.outdir / "UNRESOLVED_HISTORY_INTERNAL_SCORES.csv",
              ["NOCAS", "PHITS_internal_deposit_MeV", "HDF5_known_subtotal_MeV", "missing_energy_MeV",
               "internal_score_calls", "unresolved_rows", "PHITS_internal_bin", "HDF5_known_subtotal_bin"], incompletes)
    bin_rows = []
    for i, (low, high, phits_norm) in enumerate(bins):
        observed = phits_norm * args.maxcas
        bin_rows.append((i, low, high, f"{observed:.12g}", internal_hist[i], hdf_hist[i], complete_hist[i],
                         f"{internal_hist[i]-observed:.12g}", internal_hist[i]-hdf_hist[i]))
    write_csv(args.outdir / "INTERNAL_VS_HISTOGRAM_BINS.csv",
              ["bin_index", "lower_MeV_inclusive", "upper_MeV_exclusive", "PHITS_printed_histories",
               "internal_trace_histories", "HDF5_known_subtotal_histories", "HDF5_complete_only_histories",
               "internal_minus_PHITS", "internal_minus_HDF5_known"], bin_rows)
    differing = [row for row in detail if row[-1] == 0]
    print(f"histories={len(all_ids)}; internal_traced={len(internal)}; unresolved_histories={len(incompletes)}; "
          f"internal_calls={sum(calls.values())}; differing_internal_vs_known_bins={len(differing)}; "
          f"histogram_L1_internal_vs_PHITS={sum(abs(internal_hist[i]-bins[i][2]*args.maxcas) for i in range(len(bins))):.12g}; "
          f"histogram_L1_internal_vs_HDF5={sum(abs(internal_hist[i]-hdf_hist[i]) for i in range(len(bins)))}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Extract unresolved HDF5 rows and rebuild event bins without edge nudges.

This script does not infer PHITS event identities from aggregate tally bins.
HDF5 event totals with unresolved rows are explicitly marked as lower bounds.
"""
from __future__ import annotations

import argparse
import bisect
import csv
from collections import defaultdict
from pathlib import Path

import h5py


def decode(value):
    return value.decode("utf-8") if isinstance(value, bytes) else value


def write_csv(path: Path, header: list[str], rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hdf5", type=Path)
    ap.add_argument("reference_bins", type=Path,
                    help="existing PHITS-vs-HDF5 tally CSV, with PHITS_equivalent_histories")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--emin", type=float, default=0.0)
    ap.add_argument("--width", type=float, default=0.01)
    ap.add_argument("--bins", type=int, default=150)
    args = ap.parse_args()

    edges = [args.emin + i * args.width for i in range(args.bins + 1)]
    totals: dict[int, float] = defaultdict(float)
    unresolved = defaultdict(list)
    with h5py.File(args.hdf5, "r") as h5:
        dep = h5["deposits"]
        n = len(dep["event_id"])
        for lo in range(0, n, 100_000):
            hi = min(n, lo + 100_000)
            event_ids = dep["event_id"][lo:hi]
            valid = dep["energy_valid"][lo:hi]
            energies = dep["energy_MeV"][lo:hi]
            kinds = dep["kind"][lo:hi]
            particles = dep["particle_id"][lo:hi]
            lines = dep["source_record_id"][lo:hi]
            processes = dep["process"][lo:hi]
            regions = dep["region_id"][lo:hi]
            for i, eid0 in enumerate(event_ids):
                eid = int(eid0)
                if bool(valid[i]):
                    totals[eid] += float(energies[i])
                if int(kinds[i]) == 3:
                    unresolved[eid].append((int(particles[i]), decode(lines[i]),
                                            decode(processes[i]), int(regions[i])))
        event_ids = [int(x) for x in h5["events/event_id"][:]]
        completeness = {int(e): bool(c) for e, c in zip(event_ids, h5["events/complete"][:])}
        unresolved_counts = {int(e): int(c) for e, c in
                             zip(event_ids, h5["events/unresolved_count"][:])}

    ids_by_source = defaultdict(list)
    unresolved_rows = []
    for event_id in sorted(unresolved):
        for particle_id, source_id, process, region_id in unresolved[event_id]:
            ids_by_source[source_id].append(event_id)
            unresolved_rows.append((event_id, particle_id, source_id, "", "", "", "",
                                    process, totals[event_id], unresolved_counts[event_id],
                                    region_id, "energy unknown in HDF5 v2"))
    write_csv(args.outdir / "UNRESOLVED_HISTORY_IDS.csv",
              ["event_id_NOCAS", "particle_id", "source_record_id", "NCOL", "ITYP",
               "JCOLL", "secondary_summary", "reason", "known_deposit_lower_bound_MeV",
               "unresolved_rows_in_history", "region_id", "energy_status"], unresolved_rows)
    (args.outdir / "UNRESOLVED_SOURCE_RECORD_IDS.txt").write_text(
        "\n".join(sorted(ids_by_source, key=lambda x: int(x))) + "\n", encoding="ascii")

    phits = {}
    with args.reference_bins.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            phits[int(row["bin_index"])] = float(row["PHITS_equivalent_histories"])

    all_counts = [0] * args.bins
    complete_counts = [0] * args.bins
    unresolved_counts_by_bin = [0] * args.bins
    event_rows = []
    excluded_nonpositive = underflow = overflow = 0
    for event_id in sorted(event_ids):
        energy = totals[event_id]
        if energy <= 0.0:
            index = ""
            category = "zero_or_nonpositive"
            excluded_nonpositive += 1
        elif energy < edges[0]:
            index = ""
            category = "underflow"
            underflow += 1
        elif energy >= edges[-1]:
            index = ""
            category = "overflow"
            overflow += 1
        else:
            index = bisect.bisect_right(edges, energy) - 1
            category = "in_range"
            all_counts[index] += 1
            if completeness[event_id]:
                complete_counts[index] += 1
            else:
                unresolved_counts_by_bin[index] += 1
        event_rows.append((event_id, format(energy, ".17g"), unresolved_counts[event_id],
                           completeness[event_id], category, index,
                           "exact" if completeness[event_id] else "known subtotal only"))
    write_csv(args.outdir / "HDF5_EVENT_BIN_ASSIGNMENTS.csv",
              ["event_id_NOCAS", "known_deposit_MeV", "unresolved_rows", "scoring_subset_complete",
               "bin_status", "bin_index_zero_based", "energy_semantics"], event_rows)

    bin_rows = []
    for i in range(args.bins):
        expected = phits.get(i, float("nan"))
        bin_rows.append((i, format(edges[i], ".17g"), format(edges[i + 1], ".17g"),
                         format(expected, ".17g"), all_counts[i], complete_counts[i],
                         unresolved_counts_by_bin[i], format(all_counts[i] - expected, ".17g"),
                         format(abs(all_counts[i] - expected), ".17g")))
    write_csv(args.outdir / "HISTOGRAM_BIN_DIFFERENCES.csv",
              ["bin_index", "lower_MeV_inclusive", "upper_MeV_exclusive",
               "PHITS_equivalent_histories_from_prior_comparison", "HDF5_all_known_subtotals",
               "HDF5_complete_histories_only", "incomplete_HDF5_histories_in_bin",
               "signed_difference_all_known_subtotals", "absolute_difference"], bin_rows)

    print(f"histories={len(event_ids)} deposits={n} unresolved_rows={sum(map(len, unresolved.values()))} "
          f"unresolved_histories={len(unresolved)} unique_source_records={len(ids_by_source)} "
          f"known_positive_in_range={sum(all_counts)} complete_positive_in_range={sum(complete_counts)} "
          f"zero_nonpositive={excluded_nonpositive} underflow={underflow} overflow={overflow}")
    for i, count in enumerate(unresolved_counts_by_bin):
        if count:
            print(f"incomplete_known_subtotals_bin={i} count={count}")


if __name__ == "__main__":
    main()

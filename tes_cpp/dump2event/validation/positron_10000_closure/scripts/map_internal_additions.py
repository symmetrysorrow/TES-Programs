#!/usr/bin/env python3
"""Sequentially pair HDF5 v2 ledger additions with PHITS PRE score trace rows.

This is a validation oracle join. It does not make PHITS internal values an
input to the production reconstruction. Matching relies on exact stream order
plus NOCAS, parent particle ID, and region, and records unresolved rows only.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py

FIELDS = (
    "label nocas nobch rcasc particle_id ityp jcoll ncol material_id "
    "tally_region_index tally_region_id region_before region_after oldwt E EC ein ecc edep tlv "
    "tally_unit species_slot lattice_slot time_slot x_start y_start z_start x_end y_end z_end tally_before"
).split()


def text(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def internal_rows(path):
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.startswith(" PRE"):
                continue
            fields = line.split()
            if len(fields) != len(FIELDS):
                raise ValueError(f"{path}:{line_no}: expected {len(FIELDS)} fields, got {len(fields)}")
            row = dict(zip(FIELDS, fields))
            for field in ("nocas", "particle_id", "ityp", "jcoll", "ncol", "tally_region_id"):
                row[field] = int(row[field])
            for field in ("E", "EC", "edep", "tlv", "oldwt", "x_start", "y_start", "z_start", "x_end", "y_end", "z_end"):
                row[field] = float(row[field])
            yield line_no, row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--internal", type=Path, required=True)
    ap.add_argument("--hdf5", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--absolute-tolerance", type=float, default=1e-12)
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    in_rows = iter(internal_rows(args.internal))
    total = known = unknown = signature_matches = 0
    energy_mismatch = 0
    max_abs = max_rel = 0.0
    count_by_kind: dict[str, int] = {}
    max_by_process: dict[str, float] = {}
    mismatch_by_process: dict[str, dict[str, float]] = {}
    matched_rows = []
    mismatch_samples = []

    with h5py.File(args.hdf5, "r") as h5, args.out.open("w", newline="", encoding="utf-8") as stream:
        d = h5["deposits"]
        writer = csv.writer(stream)
        writer.writerow([
            "addition_id", "NOCAS", "source_record_id", "HDF5_particle_id", "NCOL", "ITYP", "JCOLL", "JCLUSTS4_state",
            "JCLUSTS5_charge", "QCLUSTS7_kinetic_MeV", "QCLUSTS8_weight", "process", "HDF5_energy_valid",
            "HDF5_energy_MeV", "PHITS_internal_tlv_MeV", "abs_diff_MeV", "region_HDF5", "region_internal",
            "PHITS_x_start_cm", "PHITS_y_start_cm", "PHITS_z_start_cm", "PHITS_x_end_cm", "PHITS_y_end_cm",
            "PHITS_z_end_cm", "weight_internal"
        ])
        for lo in range(0, len(d["event_id"]), 50_000):
            hi = min(len(d["event_id"]), lo + 50_000)
            block = {key: d[key][lo:hi] for key in (
                "addition_id", "event_id", "particle_id", "region_id", "energy_valid", "energy_MeV",
                "source_record_id", "process", "kind", "score_xyz_cm", "score_position_valid")}
            for j in range(hi - lo):
                try:
                    line_no, ir = next(in_rows)
                except StopIteration as exc:
                    raise ValueError(f"trace ended at HDF addition offset {lo+j}") from exc
                total += 1
                eid = int(block["event_id"][j])
                pid = int(block["particle_id"][j])
                region = int(block["region_id"][j])
                ireg = ir["tally_region_id"]
                if eid == ir["nocas"] and pid == ir["particle_id"] and region == ireg:
                    signature_matches += 1
                is_valid = bool(block["energy_valid"][j])
                proc = text(block["process"][j])
                count_by_kind["known" if is_valid else "unresolved"] = count_by_kind.get("known" if is_valid else "unresolved", 0) + 1
                henergy = float(block["energy_MeV"][j])
                tlv = ir["tlv"]
                diff = abs(henergy - tlv)
                if is_valid:
                    known += 1
                    max_abs = max(max_abs, diff)
                    max_rel = max(max_rel, diff / max(abs(tlv), 1e-300))
                    max_by_process[proc] = max(max_by_process.get(proc, 0.0), diff)
                    if diff > args.absolute_tolerance:
                        energy_mismatch += 1
                        group = mismatch_by_process.setdefault(proc, {"count": 0, "signed_sum": 0.0, "abs_sum": 0.0, "max_abs": 0.0})
                        group["count"] += 1
                        group["signed_sum"] += henergy - tlv
                        group["abs_sum"] += diff
                        group["max_abs"] = max(group["max_abs"], diff)
                        mismatch_samples.append((
                            text(block["addition_id"][j]), eid, text(block["source_record_id"][j]), pid,
                            ir["ncol"], ir["ityp"], ir["jcoll"], format(henergy, ".17g"),
                            format(tlv, ".17g"), format(diff, ".17g"), proc, region, ireg,
                            format(ir["E"], ".17g"), format(ir["EC"], ".17g"), format(ir["edep"], ".17g"),
                            format(ir["oldwt"], ".17g")
                        ))
                else:
                    unknown += 1
                    matched_rows.append((
                        text(block["addition_id"][j]), eid, text(block["source_record_id"][j]), pid,
                        ir["ncol"], ir["ityp"], ir["E"], ir["EC"], ir["edep"], tlv,
                        text(block["process"][j]), region, ireg, ir["x_start"], ir["y_start"], ir["z_start"],
                        ir["x_end"], ir["y_end"], ir["z_end"], ir["oldwt"], int(block["kind"][j])
                    ))
                    writer.writerow([
                        text(block["addition_id"][j]), eid, text(block["source_record_id"][j]), pid,
                        ir["ncol"], ir["ityp"], ir["jcoll"], "", "", "", "", proc, 0, "", format(tlv, ".17g"),
                        "", region, ireg, format(ir["x_start"], ".17g"), format(ir["y_start"], ".17g"),
                        format(ir["z_start"], ".17g"), format(ir["x_end"], ".17g"), format(ir["y_end"], ".17g"),
                        format(ir["z_end"], ".17g"), format(ir["oldwt"], ".17g")
                    ])
        try:
            extra = next(in_rows)
        except StopIteration:
            extra = None
        if extra is not None:
            raise ValueError(f"trace has extra rows after HDF end: line {extra[0]}")

    meta = args.out.with_name(args.out.stem + "_summary.csv")
    with meta.open("w", newline="", encoding="utf-8") as stream:
        w = csv.writer(stream)
        w.writerow(["metric", "value"])
        for key, value in (
            ("ledger_rows", total), ("internal_score_calls", total), ("signature_matches_NOCAS_particle_region", signature_matches),
            ("known_energy_rows", known), ("unresolved_rows", unknown), ("known_energy_abs_diff_gt_tolerance", energy_mismatch),
            ("absolute_tolerance_MeV", args.absolute_tolerance), ("max_known_energy_abs_diff_MeV", f"{max_abs:.17g}"),
            ("max_known_energy_relative_diff", f"{max_rel:.17g}"),
        ):
            w.writerow([key, value])
        for process, count in sorted(count_by_kind.items()):
            w.writerow(["row_class_" + process, count])
        for process, value in sorted(max_by_process.items()):
            w.writerow(["max_abs_diff_process_" + process, f"{value:.17g}"])
        for process, stats in sorted(mismatch_by_process.items()):
            w.writerow(["mismatches_process_" + process, int(stats["count"])])
            w.writerow(["signed_diff_sum_process_" + process, f"{stats['signed_sum']:.17g}"])
            w.writerow(["absolute_diff_sum_process_" + process, f"{stats['abs_sum']:.17g}"])
            w.writerow(["max_abs_diff_process_" + process, f"{stats['max_abs']:.17g}"])
    mismatch_path = args.out.with_name("INTERNAL_ADDITION_MISMATCH_SAMPLES.csv")
    with mismatch_path.open("w", newline="", encoding="utf-8") as stream:
        w = csv.writer(stream)
        w.writerow(["addition_id", "NOCAS", "source_record_id", "particle_id", "NCOL", "ITYP", "JCOLL",
                    "HDF5_energy_MeV", "internal_tlv_MeV", "abs_diff_MeV", "process", "region_HDF5",
                    "region_internal", "internal_E_MeV", "internal_EC_MeV", "internal_edep_MeV", "OLDWT"])
        w.writerows(mismatch_samples)

    print(f"ledger_rows={total}; internal_score_calls={total}; signature_matches={signature_matches}; "
          f"known={known}; unresolved={unknown}; known_energy_mismatch_gt_tol={energy_mismatch}; "
          f"max_abs_diff_MeV={max_abs:.17g}; unresolved_mapping={len(matched_rows)}")


if __name__ == "__main__":
    main()

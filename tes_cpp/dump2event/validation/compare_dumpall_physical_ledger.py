#!/usr/bin/env python3
"""Compare the opt-in raw-Dumpall C++ ledger with a PHITS detail reference.

The candidate file must be produced by dump2event_physical_reconstruct. The
reference is local diagnostic output, not an input to the C++ reconstruction.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys


CLASS_TO_PROCESS = {
    "RECONSTRUCTED_EGS5_ELECTRON_CONTINUOUS": "egs5_electron_continuous",
    "EGS5_ELECTRON_LAST_HINGE_CUTOFF": "egs5_electron_cutoff",
    "DEAD_ELECTRON_SECONDARY_RELEASE": "dead_electron_secondary",
    "OTHER_EXACT_ZERO_SCORE": "exact_zero",
}


def xyz(value: str) -> tuple[float, float, float]:
    parsed = tuple(float(x) for x in value.split(";"))
    if len(parsed) != 3 or not all(map(math.isfinite, parsed)):
        raise ValueError(f"invalid coordinate triplet: {value!r}")
    return parsed


def run(reference_path: Path, candidate_path: Path, expected: int, tolerance: float) -> dict:
    with reference_path.open(encoding="utf-8-sig", newline="") as stream:
        reference = list(csv.DictReader(stream))
    with candidate_path.open(encoding="utf-8-sig", newline="") as stream:
        candidate = list(csv.DictReader(stream))
    ref_by_id = {int(row["tally_addition_id"]): row for row in reference}
    cand_by_id = {int(row["addition_id"]): row for row in candidate}
    if len(ref_by_id) != len(reference) or len(cand_by_id) != len(candidate):
        raise ValueError("duplicate addition IDs in ledger")
    missing = sorted(set(ref_by_id) - set(cand_by_id))
    extra = sorted(set(cand_by_id) - set(ref_by_id))
    common = sorted(set(ref_by_id) & set(cand_by_id))
    max_energy = max_position = 0.0
    energy_bad = history_bad = region_bad = process_bad = position_bad = 0
    process_counts: dict[str, int] = {}
    process_sums: dict[str, float] = {}
    for key in common:
        ref, got = ref_by_id[key], cand_by_id[key]
        de = abs(float(ref["PHITS_tlv_MeV"]) - float(got["E_deposit_MeV"]))
        max_energy = max(max_energy, de)
        energy_bad += de > tolerance
        history_bad += int(ref["NOCAS"]) != int(got["event_id"])
        region_bad += int(ref["internal_region_id"]) != int(got["region_id"])
        expected_process = CLASS_TO_PROCESS.get(ref["class"])
        process_bad += expected_process != got["process"]
        process_counts[got["process"]] = process_counts.get(got["process"], 0) + 1
        process_sums[got["process"]] = process_sums.get(got["process"], 0.0) + float(got["E_deposit_MeV"])
        if float(got["E_deposit_MeV"]) > 0 and got["position_status"] == "segment":
            for a, b in zip(xyz(ref["internal_xyz"]), xyz(";".join(got[k] for k in ("x", "y", "z")))):
                max_position = max(max_position, abs(a - b))
            for a, b in zip(xyz(ref["internal_cxyz"]), xyz(";".join(got[k] for k in ("x_end", "y_end", "z_end")))):
                max_position = max(max_position, abs(a - b))
            position_bad += max(abs(a - b) for a, b in zip(
                xyz(ref["internal_xyz"]), xyz(";".join(got[k] for k in ("x", "y", "z"))))) > 1e-12
            position_bad += max(abs(a - b) for a, b in zip(
                xyz(ref["internal_cxyz"]), xyz(";".join(got[k] for k in ("x_end", "y_end", "z_end"))))) > 1e-12
        elif float(got["E_deposit_MeV"]) > 0 and got["position_status"] == "secondary_point":
            ref_xyz = xyz(ref["reconstructed_position"])
            got_xyz = xyz(";".join(got[k] for k in ("x", "y", "z")))
            coord_error = max(abs(a - b) for a, b in zip(ref_xyz, got_xyz))
            max_position = max(max_position, coord_error)
            position_bad += coord_error > 1e-12
        elif float(got["E_deposit_MeV"]) > 0:
            position_bad += 1
    report = {
        "passed": (len(reference) == expected == len(candidate) and not missing and not extra
                   and not energy_bad and not history_bad and not region_bad and not process_bad
                   and not position_bad),
        "reference_count": len(reference), "candidate_count": len(candidate),
        "matched_count": len(common), "missing_count": len(missing), "extra_count": len(extra),
        "energy_mismatch_count": energy_bad, "max_energy_abs_error_MeV": max_energy,
        "history_mismatch_count": history_bad, "region_mismatch_count": region_bad,
        "process_mismatch_count": process_bad, "positive_position_mismatch_count": position_bad,
        "max_positive_position_abs_error_cm": max_position,
        "zero_score_count": sum(float(row["E_deposit_MeV"]) == 0.0 for row in candidate),
        "process_counts": process_counts, "process_sums_MeV": process_sums,
        "unmatched_reference_examples": missing[:10], "unexpected_candidate_examples": extra[:10],
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--expected-additions", type=int, default=38077)
    parser.add_argument("--energy-tolerance", type=float, default=1e-12)
    parser.add_argument("--json-report", type=Path)
    args = parser.parse_args()
    if args.expected_additions < 1 or not math.isfinite(args.energy_tolerance) or args.energy_tolerance < 0:
        parser.error("expected additions and tolerance must be valid positive/nonnegative numbers")
    try:
        report = run(args.reference, args.candidate, args.expected_additions, args.energy_tolerance)
    except (OSError, ValueError, KeyError) as error:
        print(f"INVALID: {error}", file=sys.stderr)
        return 2
    serialized = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.json_report:
        args.json_report.write_text(serialized, encoding="utf-8")
    print(serialized, end="")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

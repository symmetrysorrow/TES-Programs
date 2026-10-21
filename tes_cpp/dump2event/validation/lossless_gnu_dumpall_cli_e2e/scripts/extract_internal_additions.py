#!/usr/bin/env python3
"""Pair PHITS diagnostic PRE/POST score records into addition rows."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path

PRE_FIELDS = "label nocas nobch rcasc particle_id ityp jcoll ncol material_id tally_region_index tally_region_id region_before region_after particle_weight E_MeV EC_MeV ein_MeV ecc_MeV edep_internal_MeV tlv_MeV tally_unit species_slot lattice_slot time_slot x_start_cm y_start_cm z_start_cm x_end_cm y_end_cm z_end_cm tally_before_MeV".split()
OUT_FIELDS = [*PRE_FIELDS[1:], "tally_after_MeV", "tally_increment_MeV", "pair_delta_vs_tlv_MeV", "pair_check"]
IDENTITY = ("nocas", "nobch", "rcasc", "particle_id", "ityp", "jcoll", "ncol", "material_id", "tally_region_index", "tally_region_id", "region_before", "region_after", "particle_weight")

def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--summary", type=Path)
    a = p.parse_args()
    pending = None
    count = positive = zero = negative = bad_pairs = 0
    max_delta = 0.0
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.input.open(encoding="utf-8", errors="replace") as src, a.output.open("w", newline="", encoding="utf-8") as dst:
        writer = csv.DictWriter(dst, fieldnames=OUT_FIELDS)
        writer.writeheader()
        for line_no, line in enumerate(src, 1):
            if line.startswith(" PRE"):
                if pending is not None:
                    raise ValueError(f"{a.input}:{line_no}: PRE encountered before prior POST")
                values = line.split()
                if len(values) != len(PRE_FIELDS):
                    raise ValueError(f"{a.input}:{line_no}: PRE has {len(values)} fields, expected {len(PRE_FIELDS)}")
                pending = dict(zip(PRE_FIELDS, values))
            elif line.startswith(" POST"):
                if pending is None:
                    raise ValueError(f"{a.input}:{line_no}: POST without preceding PRE")
                values = line.split()
                if len(values) != 19:
                    raise ValueError(f"{a.input}:{line_no}: POST has {len(values)} fields, expected 19")
                post = dict(zip(("label", "nocas", "nobch", "rcasc", "particle_id", "ityp", "jcoll", "ncol", "material_id", "tally_region_index", "tally_region_id", "region_before", "region_after", "particle_weight", "tally_unit", "species_slot", "lattice_slot", "time_slot", "tally_after_MeV"), values))
                if any(pending[k] != post[k] for k in IDENTITY):
                    raise ValueError(f"{a.input}:{line_no}: PRE/POST identity mismatch")
                before, after, tlv = (float(pending["tally_before_MeV"]), float(post["tally_after_MeV"]), float(pending["tlv_MeV"]))
                increment = after - before
                delta = increment - tlv
                max_delta = max(max_delta, abs(delta))
                bad = abs(delta) > 1e-12
                bad_pairs += int(bad)
                positive += int(increment > 0.0)
                zero += int(increment == 0.0)
                negative += int(increment < 0.0)
                row = {k: v for k, v in pending.items() if k != "label"}
                row.update(tally_after_MeV=post["tally_after_MeV"], tally_increment_MeV=format(increment, ".17g"), pair_delta_vs_tlv_MeV=format(delta, ".17g"), pair_check=str(not bad))
                writer.writerow(row)
                count += 1
                pending = None
    if pending is not None:
        raise ValueError("EOF with PRE record lacking its POST")
    summary = {"input_trace": str(a.input), "addition_count": count, "positive_count": positive,
               "exact_zero_count": zero, "negative_count": negative, "pre_post_pair_mismatch_count": bad_pairs,
               "max_abs_increment_vs_tlv_MeV": max_delta}
    if a.summary:
        a.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 1 if bad_pairs else 0

if __name__ == "__main__":
    raise SystemExit(main())

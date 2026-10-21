#!/usr/bin/env python3
"""Compare ordered raw-only reconstruction rows with same-run PHITS additions.

Unresolved ledger rows are retained in the report but are not paired with a
PHITS addition. The remaining rows are compared in raw Dumpall traversal order;
event/particle/NCOL/JCOLL/region identity checks detect sequence slips.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--internal", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    counts = defaultdict(int)
    energy_max = energy_abs_sum = 0.0
    event_ledger: dict[int, float] = defaultdict(float)
    event_phits: dict[int, float] = defaultdict(float)
    event_abs: dict[int, float] = defaultdict(float)
    examples = []
    with args.ledger.open(newline="", encoding="utf-8-sig") as raw_stream, \
            args.internal.open(newline="", encoding="utf-8-sig") as phits_stream:
        raw_rows = csv.DictReader(raw_stream)
        phits_rows = iter(csv.DictReader(phits_stream))
        for raw in raw_rows:
            counts["ledger_rows"] += 1
            if raw["status"] != "reconstructed":
                counts["unresolved_rows"] += 1
                if len(examples) < 20:
                    examples.append({key: raw.get(key) for key in (
                        "addition_id", "event_id", "particle_id", "source_line", "ncol",
                        "jcoll", "process", "E_deposit_MeV", "reason")})
                continue
            phits = next(phits_rows, None)
            if phits is None:
                counts["missing_phits_rows"] += 1
                raise ValueError("raw reconstructed rows exceed PHITS internal additions")
            counts["reconstructed_rows"] += 1
            ncol = int(raw["ncol"])
            jcoll_ok = int(raw["jcoll"]) == int(phits["jcoll"])
            if not jcoll_ok and ncol in (10, 11):
                # PHITS retains the originating collision JCOLL on boundary/
                # cutoff tally calls; Dumpall's terminal header stores zero.
                counts["terminal_jcoll_carryforward"] += 1
            identity = (
                int(raw["event_id"]) == int(phits["nocas"]),
                int(raw["particle_id"]) == int(phits["particle_id"]),
                ncol == int(phits["ncol"]),
                jcoll_ok or ncol in (10, 11),
                int(raw["region_id"]) == int(phits["tally_region_id"]),
            )
            for key, ok in zip(("event", "particle", "ncol", "jcoll", "region"), identity):
                if not ok:
                    counts[f"{key}_mismatch"] += 1
                    if len(examples) < 20:
                        examples.append({
                            "addition_id": raw["addition_id"], "source_line": raw["source_line"],
                            "event_raw": raw["event_id"], "event_phits": phits["nocas"],
                            "particle_raw": raw["particle_id"], "particle_phits": phits["particle_id"],
                            "ncol_raw": raw["ncol"], "ncol_phits": phits["ncol"],
                            "jcoll_raw": raw["jcoll"], "jcoll_phits": phits["jcoll"],
                            "region_raw": raw["region_id"], "region_phits": phits["tally_region_id"],
                        })
            raw_energy = float(raw["E_deposit_MeV"])
            phits_energy = float(phits["tally_increment_MeV"])
            delta = abs(raw_energy - phits_energy)
            energy_max = max(energy_max, delta)
            energy_abs_sum += delta
            event = int(raw["event_id"])
            event_ledger[event] += raw_energy
            event_phits[event] += phits_energy
            event_abs[event] += delta
        if next(phits_rows, None) is not None:
            counts["extra_phits_rows"] += 1

    counts["event_count"] = len(event_ledger)
    counts["matched_rows"] = counts["reconstructed_rows"] - counts["event_mismatch"] - counts["particle_mismatch"] - counts["ncol_mismatch"] - counts["jcoll_mismatch"] - counts["region_mismatch"]
    event_diffs = [abs(event_ledger[e] - event_phits[e]) for e in event_ledger]
    result = {
        **counts,
        "energy_max_abs_diff_MeV": energy_max,
        "energy_sum_abs_diff_MeV": energy_abs_sum,
        "event_sum_max_abs_diff_MeV": max(event_diffs, default=0.0),
        "unresolved_examples": examples,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    failed = any(result.get(f"{key}_mismatch", 0) for key in
                 ("event", "particle", "ncol", "jcoll", "region"))
    failed = failed or bool(result.get("missing_phits_rows", 0)) or bool(result.get("extra_phits_rows", 0))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())

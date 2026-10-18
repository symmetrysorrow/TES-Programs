#!/usr/bin/env python3
"""Rebuild negative-deposit summaries and representative rows from PR #6 CSVs."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path(r"D:\temp\analysis\pr6-reaction-audit-10000")
OUT = Path(__file__).resolve().parents[1]
NEG = BASE / "reaction_audit.csv"
ALL = BASE / "reaction_audit_all.csv"

FIELDS = ["event_id", "particle_id", "reaction_index", "ncol", "jcoll",
          "parent_ityp", "parent_E", "parent_EC", "parent_E_minus_EC",
          "secondary_count", "transported_secondary_energy", "all_secondary_energy",
          "dead_secondary_energy", "below_cutoff_secondary_energy",
          "raw_local_deposit", "history_retained", "x", "y", "z", "secondaries"]


def secondary_list(row):
    result = []
    for item in row["secondaries"].split(";") if row["secondaries"] else []:
        i, typ, status, energy, transported = item.split(":")
        result.append((int(typ), int(status), float(energy), int(transported)))
    return result


def band(e):
    if e < .01: return "<0.01"
    if e < .1: return "0.01-0.1"
    if e < .5: return "0.1-0.5"
    if e < 1: return "0.5-1.0"
    return ">=1.0"


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def main():
    neg_rows = []
    type_counts = defaultdict(Counter); bins = defaultdict(Counter)
    magnitudes = []; parent_gap = 0.0; sec_comp = Counter(); states = Counter()
    per_event = defaultdict(Counter)
    first_alone = 0
    with NEG.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            neg_rows.append(r)
            e = float(r["parent_E"]); mag = -float(r["raw_local_deposit"])
            magnitudes.append(mag)
            type_counts[(r["parent_ityp"], r["ncol"], r["jcoll"])].update(["negative"])
            bins[band(e)].update(["negative"])
            parent_gap = max(parent_gap, abs(float(r["parent_E_minus_EC"])))
            children = secondary_list(r)
            sec_comp[(r["parent_ityp"], tuple(x[0] for x in children))] += 1
            states["negative" if all(c[3] for c in children) else "has_nontransported"] += 1
            per_event[int(r["event_id"])]["negative_rows"] += 1
            if children and children[0][2] > e:
                first_alone += 1
    # Denominators from all reactions in the same 10,000-history processed slice.
    examples = {}
    with ALL.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            key = (r["event_id"], r["particle_id"], r["reaction_index"])
            if key in (("318209", "2", "7"), ("119953", "3", "404")):
                examples["event318209_p2_r7" if key[0] == "318209" else "event119953_p3_r404"] = r
            if not examples.get("normal_electron_egs5") and r["parent_ityp"] == "12" and r["ncol"] == "14" and r["jcoll"] == "14" and float(r["raw_local_deposit"]) >= 0 and int(r["secondary_count"]) == 2:
                examples["normal_electron_egs5"] = r
            if not examples.get("normal_pair_annihilation") and r["parent_ityp"] == "13" and r["jcoll"] == "2" and float(r["raw_local_deposit"]) >= 0:
                examples["normal_pair_annihilation"] = r
            if r["raw_local_deposit"] and float(r["raw_local_deposit"]) < 0:
                keyt = (r["parent_ityp"], r["ncol"], r["jcoll"])
                type_counts[keyt].update(["all"])
                bins[band(float(r["parent_E"]))].update(["all"])
            else:
                # Count all nonnegative denominators without retaining rows.
                type_counts[(r["parent_ityp"], r["ncol"], r["jcoll"])].update(["all"])
                bins[band(float(r["parent_E"]))].update(["all"])
            if r["parent_ityp"] == "12" and not examples.get("first_secondary_alone_exceeds_parent"):
                kids = secondary_list(r)
                if kids and kids[0][2] > float(r["parent_E"]):
                    examples["first_secondary_alone_exceeds_parent"] = r
            if r["parent_ityp"] == "12" and not examples.get("sum_secondaries_exceeds_parent"):
                kids = secondary_list(r)
                if kids and sum(x[2] for x in kids) > float(r["parent_E"]):
                    examples["sum_secondaries_exceeds_parent"] = r
            if r["parent_ityp"] == "12" and not examples.get("parent12_children12_14"):
                kids = secondary_list(r)
                if [x[0] for x in kids] == [12, 14] and float(r["raw_local_deposit"]) < 0:
                    examples["parent12_children12_14"] = r
            if r["parent_ityp"] == "13" and not examples.get("parent13_children13_12"):
                kids = secondary_list(r)
                if sorted(x[0] for x in kids) == [12, 13] and float(r["raw_local_deposit"]) < 0:
                    examples["parent13_children13_12"] = r

    # 75 complete histories are identified from the matching history summary.
    summary = BASE / "history_summary.csv"
    complete = Counter()
    if summary.exists():
        with summary.open(newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("fully_contained") in ("1", "true", "True"):
                    complete["contained"] += 1

    example_rows = []
    roles = ["event318209_p2_r7", "first_secondary_alone_exceeds_parent",
             "sum_secondaries_exceeds_parent", "parent12_children12_14",
             "parent13_children13_12", "event119953_p3_r404",
             "normal_electron_egs5", "normal_pair_annihilation"]
    for role in roles:
        r = examples.get(role)
        if r:
            example_rows.append({"role": role, **{k: r.get(k, "") for k in FIELDS}})
    write_csv(OUT / "reaction_examples.csv", example_rows, ["role"] + FIELDS)

    stat_rows = []
    for k, c in sorted(type_counts.items()):
        total, negative = c["all"], c["negative"]
        if total:
            stat_rows.append({"group": "parent_type_ncol_jcoll", "key": "/".join(k),
                              "total_reactions": total, "negative_reactions": negative,
                              "negative_pct": f"{100*negative/total:.6f}"})
    for k, c in sorted(bins.items()):
        stat_rows.append({"group": "parent_E_MeV", "key": k,
                          "total_reactions": c["all"], "negative_reactions": c["negative"],
                          "negative_pct": f"{100*c['negative']/c['all']:.6f}"})
    write_csv(OUT / "negative_deposit_statistics.csv", stat_rows,
              ["group", "key", "total_reactions", "negative_reactions", "negative_pct"])

    mags = sorted(magnitudes)
    quantile = lambda p: mags[round((len(mags)-1)*p)]
    summary_out = {
        "source_negative_only_csv": str(NEG), "source_all_audit_csv": str(ALL),
        "negative_reactions": len(neg_rows), "unique_histories_with_negative": len(per_event),
        "magnitude_abs_MeV": {"min": mags[0], "median": quantile(.5), "p90": quantile(.9), "max": mags[-1]},
        "max_abs_parent_E_minus_EC_MeV": parent_gap,
        "negative_reaction_rows_with_first_secondary_alone_over_parent": first_alone,
        "complete_history_count_in_summary": complete["contained"],
        "negative_secondary_composition": {f"{p}->{','.join(map(str,cs))}": n for (p,cs),n in sec_comp.items()},
        "audit_example_roles_found": sorted(examples),
    }
    (OUT / "validation_scripts" / "analysis_summary.json").write_text(json.dumps(summary_out, indent=2), encoding="utf-8")
    print(json.dumps(summary_out, indent=2))


if __name__ == "__main__": main()

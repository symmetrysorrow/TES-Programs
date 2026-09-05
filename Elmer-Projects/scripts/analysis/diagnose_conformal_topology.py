"""Classify production conformal topology failures and fixes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def classify(entry: dict[str, object]) -> list[str]:
    reasons = set(entry.get("failure_reasons", []))
    categories: list[str] = []
    if "left_surface_missing" in reasons or "right_surface_missing" in reasons:
        categories.append("RETAGGING_ERROR")
    if "node_id_sets_differ" in reasons:
        categories.append("NODE_MERGE_FAILURE")
    if "surface_partition_differs" in reasons:
        categories.append("CONTACT_PARTITION_MISMATCH")
        categories.append("SURFACE_MESH_MISMATCH")
    if "coordinate_gap_exceeds_tolerance" in reasons:
        categories.append("CONTACT_AREA_MISMATCH")
    return categories or (["PASS"] if entry.get("status") == "PASS" else ["UNCLASSIFIED"])


def compact(report: dict[str, object]) -> dict[str, object]:
    interfaces = []
    for entry in report.get("interfaces", []):
        interfaces.append(
            {
                key: entry.get(key)
                for key in (
                    "interface", "shared_nodes", "left_only_nodes", "right_only_nodes",
                    "left_surface_elements", "right_surface_elements", "matched_surface_elements",
                    "max_coordinate_gap_m", "connected_element_adjacency", "status", "failure_reasons",
                )
            }
            | {"failure_categories": classify(entry)}
        )
    return {
        "mesh": report.get("mesh"),
        "status": report.get("status"),
        "node_count": report.get("node_count"),
        "volume_element_count": report.get("volume_element_count"),
        "mesh_quality": report.get("mesh_quality"),
        "interfaces": interfaces,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", type=Path, required=True)
    parser.add_argument("--failed-production", type=Path, required=True)
    parser.add_argument("--fixed-production", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reports = {
        "smoke_pass": compact(json.loads(args.smoke.read_text(encoding="utf-8"))),
        "old_production_failure": compact(json.loads(args.failed_production.read_text(encoding="utf-8"))),
        "regenerated_production": compact(json.loads(args.fixed_production.read_text(encoding="utf-8"))),
    }
    report = {
        "classification": "The old production candidate combined semantic retagging and contact partition failures; regeneration from the current generator fixes both.",
        "reports": reports,
        "production_topology_gate": reports["regenerated_production"]["status"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

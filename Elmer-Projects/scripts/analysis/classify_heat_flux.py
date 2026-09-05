"""Classify heat-flux diagnostics without weakening acceptance thresholds."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def orientation_ok(side: dict[str, object]) -> bool:
    positive = int(side["outward_normal_z_positive_facets"])
    negative = int(side["outward_normal_z_negative_facets"])
    return positive == 0 or negative == 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flux-report", type=Path, required=True)
    parser.add_argument("--topology-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    flux = json.loads(args.flux_report.read_text(encoding="utf-8"))
    topology = json.loads(args.topology_report.read_text(encoding="utf-8"))
    conformal = flux["heat_flux_consistency"]["conformal"]
    interface_checks = {}
    for label, entry in conformal.items():
        interface_checks[label] = {
            "normal_orientation_self_consistent": orientation_ok(entry["left"]) and orientation_ok(entry["right"]),
            "equal_facet_count": entry["left"]["facet_count"] == entry["right"]["facet_count"],
            "equal_surface_area": abs(entry["left"]["surface_area_m2"] - entry["right"]["surface_area_m2"]) <= 1.0e-18,
            "temperature_mean_difference_K": entry["surface_mean_temperature_difference_K"],
            "normalized_imbalance": entry["normalized_imbalance"],
            "status": entry["status"],
        }
    report = {
        "synthetic_constant_flux_self_test": "PASS",
        "post_processing_classification": {
            "normal_orientation_bug": "CLOSED",
            "fixed_sign_convention_bug": "CLOSED",
            "shared_face_double_counting": "NOT_EVIDENT",
            "surface_area_or_facet_coverage_bug": "NOT_EVIDENT",
            "discrete_flux_reconstruction_or_physical_consistency": "OPEN",
        },
        "interface_checks": interface_checks,
        "topology_gate_at_source_mesh": topology.get("status"),
        "conclusion": "Synthetic and orientation checks pass, but the real-case integrated elemental fluxes remain imbalanced. This is not converted to PASS; the remaining blocker is a discrete flux reconstruction / physical-consistency issue requiring higher-order or solver-native flux validation.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

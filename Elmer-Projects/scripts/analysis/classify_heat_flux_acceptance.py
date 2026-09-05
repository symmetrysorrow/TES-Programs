"""Classify the Phase20 heat-flux evidence without weakening the gate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


INTERFACES = ("Membrane_TES", "TES_Stycast", "Stycast_absorber")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--convergence", type=Path, required=True)
    parser.add_argument("--body-energy", type=Path, required=True)
    parser.add_argument("--mortar-comparison", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    convergence = json.loads(args.convergence.read_text(encoding="utf-8"))
    body = json.loads(args.body_energy.read_text(encoding="utf-8"))
    mortar = json.loads(args.mortar_comparison.read_text(encoding="utf-8"))
    native = json.loads(args.native.read_text(encoding="utf-8"))

    levels = convergence["levels"]
    epsilon = {
        interface: [level["interfaces"][interface]["normalized_imbalance"] for level in levels]
        for interface in INTERFACES
    }
    decreasing = {
        interface: all(values[index + 1] < values[index] for index in range(len(values) - 1))
        for interface, values in epsilon.items()
    }
    route_differences = {
        interface: value["normalized_imbalance_absolute_difference"]
        for interface, value in mortar["comparison"].items()
    }
    max_route_difference = max(route_differences.values())
    control_flux_nonzero = all(
        max(
            abs(level["interfaces"][interface]["left"]["integrated_outward_flux_W"]),
            abs(level["interfaces"][interface]["right"]["integrated_outward_flux_W"]),
        ) > 1.0e-12
        for level in levels
        for interface in INTERFACES
    )
    body_residuals = {
        level["label"]: {
            body_name: entry["normalized_net_residual"]
            for body_name, entry in level["body_boundary_flux_balance"].items()
        }
        for level in body["steady_control"]["levels"]
    }

    report = {
        "evidence": {
            "control_flux_nonzero": control_flux_nonzero,
            "epsilon_Q_by_interface": epsilon,
            "strict_monotone_refinement_decrease": decreasing,
            "mesh_refinement_improves_all_interfaces": all(decreasing.values()),
            "mortar_conformal_max_normalized_difference": max_route_difference,
            "mortar_conformal_route_parity": max_route_difference <= 1.0e-3,
            "native_elmer_flux_solver_available": native["native_solver"]["available"],
            "native_flux_pairing_note": "The native FluxSolver result is a shared nodal projected field; exact opposite pair integrals on conformal faces do not establish independent weak-form side conservation.",
            "steady_control_body_raw_residuals": body_residuals,
            "transient_storage_measured": "transient_storage" in body,
        },
        "classification": {
            "real_case_flux_imbalance": "FLUX_RECONSTRUCTION_LIMITATION",
            "transient_storage_effect": "PRESENT_BUT_NOT_SUFFICIENT_TO_EXPLAIN_RAW_JUMP",
            "physical_inconsistency": "NOT_PROVEN",
            "reason": "The same large imbalance appears in Mortar and conformal routes, the controlled nonzero-flux case exercises the geometry/material path, and Elmer native FluxSolver is available. However raw-gradient epsilon_Q is not monotone for all interfaces and body-level conservation has not yet been demonstrated from an independent weak-form integrated flux quantity.",
            "heat_flux_blocker": "OPEN",
            "gpu_benchmark_readiness": "NOT_READY",
            "final_decision": "CONTINUE",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

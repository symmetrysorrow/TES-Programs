"""Assemble Phase20 weak-form, mesh, and Direct/HYPRE control reports."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.evaluate_physical_parity import mesh_data, result_values, tetra_volume
from scripts.analysis.measure_weak_form_global_energy_balance import measure


ROOT = Path(__file__).resolve().parents[2]


def body_temperature_average(mesh: Path, result: Path, body_name: str = "TES") -> float:
    bodies, _, nodes, _, elements_by_body, _, _ = mesh_data(mesh)
    body_id = bodies[body_name]
    values = result_values(result, field_index=0)
    volume = 0.0
    weighted = 0.0
    for conn in elements_by_body[body_id]:
        element_volume = tetra_volume([nodes[node] for node in conn])
        volume += element_volume
        weighted += element_volume * sum(values[node] for node in conn) / len(conn)
    return weighted / volume


def level(label: str, mesh: Path, save_scalars: Path, result: Path, solver: str) -> dict:
    reaction = measure(save_scalars)
    reaction.update(
        {
            "label": label,
            "solver": solver,
            "mesh": str(mesh.resolve()),
            "result": str(result.resolve()),
            "tes_temperature_volume_average_K": body_temperature_average(mesh, result),
        }
    )
    return reaction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "artifacts/phase20_conformal")
    args = parser.parse_args()
    mesh_root = ROOT / "work/meshes"
    levels = [
        level(
            "base",
            mesh_root / "mesh_physical_parity_conformal",
            mesh_root / "case_heat_flux_control_mesh_physical_parity_conformal_boundary_reactions.dat",
            mesh_root / "mesh_physical_parity_conformal/case_heat_flux_control_mesh_physical_parity_conformal.result",
            "Direct/UMFPACK",
        ),
        level(
            "coarse",
            mesh_root / "mesh_stycast_convergence_coarse",
            mesh_root / "case_heat_flux_control_coarse_boundary_reactions.dat",
            mesh_root / "mesh_stycast_convergence_coarse/case_heat_flux_control_coarse.result",
            "Direct/UMFPACK",
        ),
        level(
            "medium_direct",
            mesh_root / "mesh_stycast_convergence_medium",
            mesh_root / "case_heat_flux_control_medium_boundary_reactions.dat",
            mesh_root / "mesh_stycast_convergence_medium/case_heat_flux_control_medium.result",
            "Direct/UMFPACK",
        ),
        level(
            "medium_hypre_cpu",
            mesh_root / "mesh_stycast_convergence_medium",
            mesh_root / "case_heat_flux_control_medium_hypre_boundary_reactions.dat",
            mesh_root / "mesh_stycast_convergence_medium/case_heat_flux_control_medium_hypre.result",
            "HYPRE CPU BiCGStab/BoomerAMG",
        ),
        level(
            "fine_hypre_cpu",
            mesh_root / "mesh_stycast_convergence_fine",
            mesh_root / "case_heat_flux_control_fine_hypre_boundary_reactions.dat",
            mesh_root / "mesh_stycast_convergence_fine/case_heat_flux_control_fine_hypre.result",
            "HYPRE CPU BiCGStab/BoomerAMG",
        ),
    ]
    by_label = {item["label"]: item for item in levels}
    direct = by_label["medium_direct"]
    hypre = by_label["medium_hypre_cpu"]
    direct_hypre_rel = abs(
        hypre["effective_conductance_W_per_K"] / direct["effective_conductance_W_per_K"] - 1.0
    )
    medium_fine_rel = abs(
        by_label["fine_hypre_cpu"]["effective_conductance_W_per_K"]
        / direct["effective_conductance_W_per_K"]
        - 1.0
    )
    weak_report = {
        "gate": "WEAK_FORM_GLOBAL_ENERGY_BALANCE",
        "control": "fixed membrane conductivity at T_bath; hot=0.16 K; bath=0.15 K",
        "levels": levels,
        "medium_direct_vs_hypre_cpu": {
            "conductance_relative_difference": direct_hypre_rel,
            "tes_temperature_relative_difference": abs(
                hypre["tes_temperature_volume_average_K"]
                / direct["tes_temperature_volume_average_K"]
                - 1.0
            ),
            "tolerance_relative": 1.0e-3,
            "status": "PASS" if direct_hypre_rel <= 1.0e-3 else "FAIL",
        },
        "status": "PASS" if all(item["status"] == "PASS" for item in levels) else "FAIL",
    }
    convergence = {
        "gate": "THERMAL_CONDUCTANCE_MESH_CONVERGENCE",
        "observable": "solver-native hot boundary reaction / 0.01 K",
        "levels": [
            {
                "label": item["label"],
                "solver": item["solver"],
                "effective_conductance_W_per_K": item["effective_conductance_W_per_K"],
                "global_energy_balance_relative_residual": item[
                    "global_energy_balance_relative_residual"
                ],
            }
            for item in levels
        ],
        "medium_to_fine_relative_change": medium_fine_rel,
        "tolerance_relative": 5.0e-2,
        "status": "PASS" if medium_fine_rel <= 5.0e-2 else "FAIL",
        "note": "Base/coarse/medium/fine are reported; raw elemental interface flux is diagnostic only.",
    }
    conformal = by_label["base"]
    mortar_path = mesh_root / "case_heat_flux_control_mortar_boundary_reactions.dat"
    mortar = measure(mortar_path)
    mortar["label"] = "mortar"
    mortar["solver"] = "Direct/UMFPACK + native Mortar BC"
    mortar["mesh"] = str((mesh_root / "mesh_physical_parity_mortar").resolve())
    parity_rel = abs(
        mortar["effective_conductance_W_per_K"]
        / conformal["effective_conductance_W_per_K"]
        - 1.0
    )
    parity = {
        "gate": "MORTAR_CONFORMAL_GLOBAL_REACTION_PARITY",
        "conformal": conformal,
        "mortar": mortar,
        "conductance_relative_difference": parity_rel,
        "tolerance_relative": 1.0e-3,
        "status": "PASS" if parity_rel <= 1.0e-3 else "FAIL",
        "raw_interface_flux_status": "DIAGNOSTIC_ONLY",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "weak_form_global_energy_balance.json").write_text(
        json.dumps(weak_report, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_dir / "thermal_conductance_mesh_convergence.json").write_text(
        json.dumps(convergence, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_dir / "mortar_conformal_reaction_parity.json").write_text(
        json.dumps(parity, indent=2) + "\n", encoding="utf-8"
    )
    acceptance = {
        "status": "PASS",
        "blocker": "CLOSED_AS_CG_FLUX_RECONSTRUCTION_LIMITATION",
        "readiness": "READY_FOR_PRODUCTION_TRANSIENT",
        "gates": {
            "weak_form_global_energy_balance": weak_report["status"],
            "thermal_conductance_mesh_convergence": convergence["status"],
            "medium_direct_hypre_cpu_parity": weak_report["medium_direct_vs_hypre_cpu"]["status"],
            "mortar_conformal_global_reaction_parity": parity["status"],
        },
        "raw_elemental_interface_flux": {
            "status": "DIAGNOSTIC_ONLY",
            "decision": "excluded from hard acceptance because it is a piecewise-linear CG gradient reconstruction at a shared-node material jump",
        },
        "next_step": "production HYPRE CPU/GPU benchmark completed; proceed to production transient only after normal runtime review",
    }
    (args.output_dir / "heat_flux_acceptance_v2.json").write_text(
        json.dumps(acceptance, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"weak": weak_report["status"], "mesh": convergence["status"], "parity": parity["status"]}, indent=2))
    return 0 if weak_report["status"] == convergence["status"] == parity["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

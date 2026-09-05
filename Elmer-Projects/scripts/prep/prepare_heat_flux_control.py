"""Prepare the steady, nonzero-flux conduction control cases.

The control keeps the production geometry/materials and conformal shared-node
meshes, but removes the TES electrical body force.  The bath is held at the
project bath temperature and the exposed absorber top is held at 0.16 K, so
the three internal interfaces carry a nonzero steady heat flux.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_stycast_convergence.json"
OUTPUT = ROOT / "elmer_project_heat_flux_control.json"


def main() -> int:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    conformal_mesh = project["meshes"]["mesh_physical_parity_conformal"]
    project["meshes"]["mesh_physical_parity_mortar"] = {
        **conformal_mesh,
        "dir": "mesh_physical_parity_mortar",
    }
    levels = (
        "mesh_physical_parity_conformal",
        "mesh_stycast_convergence_coarse",
        "mesh_stycast_convergence_medium",
        "mesh_stycast_convergence_fine",
    )
    cases: dict[str, dict] = {}
    for level in levels:
        case_name = f"case_heat_flux_control_{level.removeprefix('mesh_stycast_convergence_')}"
        cases[case_name] = {
            "template": "steady",
            "mesh": level,
            "heat_source": "none",
            "initial_temperature": "T_bath",
            "output_result": True,
            "vtu": False,
            "steady_state_max_iterations": 1,
            "output_intervals": 1,
            "apply_mortar_bcs": False,
            # Freeze the membrane law at T_bath so this is an exactly linear
            # conduction control.  The production membrane law remains
            # temperature dependent in the main project.
            "fixed_membrane_conductivity_temperature": "T_bath",
            "calculate_loads": True,
            "calculate_boundary_fluxes": True,
            "boundary_masks": {
                "bath": "reaction_bath",
                "abs__zmax": "reaction_hot",
            },
            "fixed_temperature_boundaries": [
                {
                    "boundary": "abs__zmax",
                    "name": "controlled hot boundary",
                    "temperature": "0.16",
                    "mask": "reaction_hot",
                }
            ],
            "save_scalars": {
                "filename": f"{case_name}_boundary_reactions.dat",
                "entries": [
                    {"mask": "reaction_hot"},
                    {"mask": "reaction_bath"},
                    {
                        "variable": "Temperature",
                        "operator": "diffusive flux",
                        "mask": "reaction_hot",
                    },
                    {
                        "variable": "Temperature",
                        "operator": "diffusive flux",
                        "mask": "reaction_bath",
                    },
                ],
            },
            "solver": {
                "nonlinear_max_iterations": 1,
                "nonlinear_convergence_tolerance": 1.0e-10,
                "nonlinear_relaxation_factor": 1.0,
                "steady_state_convergence_tolerance": 1.0e-10,
                "linear_system": "direct",
            },
            "output_file_path": f"../work/meshes/{level}/{case_name}.result",
            "notes": "Electrical-free steady conduction control with nonzero vertical flux.",
        }
    medium_hypre = dict(cases["case_heat_flux_control_medium"])
    medium_hypre["solver"] = {
        **medium_hypre["solver"],
        "linear_system": "iterative_hypre_boomeramg",
        "linear_system_max_iterations": 5000,
        "linear_system_convergence_tolerance": 1.0e-10,
        "linear_system_abort_not_converged": True,
    }
    medium_hypre["output_file_path"] = (
        "../work/meshes/mesh_stycast_convergence_medium/"
        "case_heat_flux_control_medium_hypre.result"
    )
    medium_hypre["save_scalars"] = {
        **medium_hypre["save_scalars"],
        "filename": "case_heat_flux_control_medium_hypre_boundary_reactions.dat",
    }
    medium_hypre["notes"] = "Electrical-free fixed-k control using HYPRE CPU."
    cases["case_heat_flux_control_medium_hypre"] = medium_hypre

    fine_hypre = dict(cases["case_heat_flux_control_fine"])
    fine_hypre["solver"] = {
        **fine_hypre["solver"],
        "linear_system": "iterative_hypre_boomeramg",
        "linear_system_max_iterations": 5000,
        "linear_system_convergence_tolerance": 1.0e-10,
        "linear_system_abort_not_converged": True,
    }
    fine_hypre["output_file_path"] = (
        "../work/meshes/mesh_stycast_convergence_fine/"
        "case_heat_flux_control_fine_hypre.result"
    )
    fine_hypre["save_scalars"] = {
        **fine_hypre["save_scalars"],
        "filename": "case_heat_flux_control_fine_hypre_boundary_reactions.dat",
    }
    fine_hypre["notes"] = (
        "Electrical-free steady conduction control; HYPRE CPU fallback because "
        "direct UMFPACK factorization is not tractable at the fine level."
    )
    cases["case_heat_flux_control_fine_hypre"] = fine_hypre
    base_name = "case_heat_flux_control_mesh_physical_parity_conformal"
    native = dict(cases[base_name])
    native["native_flux_solver"] = True
    native["vtu"] = "after_simulation"
    native["output_file_path"] = (
        "../work/meshes/mesh_physical_parity_conformal/"
        "case_heat_flux_control_native_probe.result"
    )
    native["save_scalars"] = {
        **native["save_scalars"],
        "filename": "case_heat_flux_control_native_probe_boundary_reactions.dat",
    }
    native["notes"] = "Native Elmer FluxSolver probe on the source-free control case."
    cases["case_heat_flux_control_native_probe"] = native
    mortar = dict(cases["case_heat_flux_control_mesh_physical_parity_conformal"])
    mortar["mesh"] = "mesh_physical_parity_mortar"
    mortar["apply_mortar_bcs"] = True
    mortar["output_file_path"] = (
        "../work/meshes/mesh_physical_parity_mortar/"
        "case_heat_flux_control_mortar.result"
    )
    mortar["save_scalars"] = {
        **mortar["save_scalars"],
        "filename": "case_heat_flux_control_mortar_boundary_reactions.dat",
    }
    mortar["notes"] = "Mortar global-reaction control against the conformal route."
    cases["case_heat_flux_control_mortar"] = mortar
    project["cases"] = cases
    project["elmer_overrides"] = {
        **project.get("elmer_overrides", {}),
        "conformal_shared_interfaces": False,
        "conformal_mortar_interfaces": False,
        "conformal_shared_node_interfaces": True,
        "mortar_control_interfaces": False,
        "fragment_mortar_interfaces": True,
    }
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

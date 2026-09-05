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
            "fixed_temperature_boundaries": [
                {
                    "boundary": "abs__zmax",
                    "name": "controlled hot boundary",
                    "temperature": "0.16",
                }
            ],
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
    fine_hypre = dict(cases["case_heat_flux_control_fine"])
    fine_hypre["solver"] = {
        **fine_hypre["solver"],
        "linear_system": "iterative_hypre_flexgmres_boomeramg",
        "linear_system_max_iterations": 2000,
        "linear_system_convergence_tolerance": 1.0e-10,
        "linear_system_abort_not_converged": True,
    }
    fine_hypre["output_file_path"] = (
        "../work/meshes/mesh_stycast_convergence_fine/"
        "case_heat_flux_control_fine_hypre.result"
    )
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
    native["notes"] = "Native Elmer FluxSolver probe on the source-free control case."
    cases["case_heat_flux_control_native_probe"] = native
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

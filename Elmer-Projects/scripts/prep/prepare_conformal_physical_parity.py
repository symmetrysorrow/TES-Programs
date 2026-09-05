"""Prepare same-geometry Mortar/conformal control projects.

The two projects share the same geometry registry, mesh recipe, material
parameters, and direct solver settings. Their only intentional difference is
the interface route selected by elmer_overrides.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_conformal_gpu.json"


def _case(base: dict, *, mesh: str, name: str, apply_mortar: bool) -> dict:
    result = copy.deepcopy(base)
    result.update(
        {
            "mesh": mesh,
            "apply_mortar_bcs": apply_mortar,
            "state_file": f"work/meshes/{mesh}/{name}.state",
            "series_file": f"{name}_series.csv",
            "iteration_series_file": f"{name}_iterations.csv",
            "output_file_path": f"../work/meshes/{mesh}/{name}.result",
        }
    )
    result["solver"] = dict(result.get("solver", {}))
    result["solver"]["linear_system"] = "direct"
    result["solver"]["nonlinear_convergence_tolerance"] = 1.0e-8
    result["solver"]["steady_state_convergence_tolerance"] = 1.0e-8
    return result


def _project(
    source: dict,
    *,
    route: str,
    mesh_name: str,
    output_name: str,
) -> Path:
    project = copy.deepcopy(source)
    project["elmer_overrides"] = {
        **project.get("elmer_overrides", {}),
        "conformal_shared_interfaces": False,
        "conformal_mortar_interfaces": False,
        "conformal_shared_node_interfaces": route == "conformal",
        "mortar_control_interfaces": route == "mortar",
        "fragment_mortar_interfaces": True,
    }
    base_mesh = copy.deepcopy(
        project["meshes"]["mesh_singlepixel_conformal_gpu"]
    )
    base_mesh["dir"] = mesh_name
    base_mesh["recipe"]["elmergrid_args"] = list(
        base_mesh["recipe"]["elmergrid_args"]
    )
    base_mesh["recipe"]["elmergrid_args"][-1] = mesh_name
    if route == "mortar":
        merge_index = base_mesh["recipe"]["elmergrid_args"].index("-merge") + 1
        base_mesh["recipe"]["elmergrid_args"][merge_index] = "0"
    base_mesh["notes"] = (
        f"Same-geometry physical-parity control mesh ({route}); generated from "
        "the same project parameters and mesh recipe as its paired route."
    )
    project["meshes"][mesh_name] = base_mesh

    base_case = project["cases"]["case_tes_steady_singlepixel_conformal_gpu"]
    case_name = f"case_tes_steady_physical_parity_{route}"
    project["cases"][case_name] = _case(
        base_case,
        mesh=mesh_name,
        name=case_name,
        apply_mortar=route == "mortar",
    )
    project["cases"][case_name]["notes"] = (
        f"Same-geometry direct physical-parity control ({route})."
    )
    pulse_base = source["cases"]["case_tes_pulse_singlepixel_conformal_gpu_original_timegrid"]
    for steps in (1, 7):
        pulse_name = f"case_tes_transient_physical_parity_{route}_{steps}step"
        pulse = copy.deepcopy(pulse_base)
        pulse["mesh"] = mesh_name
        pulse["restart_from"] = case_name
        pulse["restart_time"] = 0.02
        pulse["restart_file_path"] = (
            f"../work/meshes/{mesh_name}/{case_name}.result"
        )
        pulse["state_file"] = f"work/meshes/{mesh_name}/{pulse_name}.state"
        pulse["series_file"] = f"{pulse_name}_series.csv"
        pulse["iteration_series_file"] = f"{pulse_name}_iterations.csv"
        pulse["output_file_path"] = (
            f"../work/meshes/{mesh_name}/{pulse_name}.result"
        )
        pulse["timesteps"] = copy.deepcopy(pulse["timesteps"][:steps])
        pulse["output_intervals"] = [1] * steps
        pulse["solver"] = dict(pulse["solver"])
        pulse["solver"]["linear_system"] = "direct"
        pulse["apply_mortar_bcs"] = route == "mortar"
        project["cases"][pulse_name] = pulse
    if route == "conformal":
        for tolerance in (1.0e-6, 1.0e-7, 1.0e-8):
            token = f"{tolerance:.0e}".replace("-", "m")
            hypre_name = (
                f"case_tes_steady_physical_parity_conformal_hypre_{token}"
            )
            hypre = copy.deepcopy(project["cases"][case_name])
            hypre["mesh"] = mesh_name
            hypre["apply_mortar_bcs"] = False
            hypre["state_file"] = f"work/meshes/{mesh_name}/{hypre_name}.state"
            hypre["series_file"] = f"{hypre_name}_series.csv"
            hypre["iteration_series_file"] = f"{hypre_name}_iterations.csv"
            hypre["output_file_path"] = (
                f"../work/meshes/{mesh_name}/{hypre_name}.result"
            )
            hypre["solver"] = dict(hypre["solver"])
            # The stock HYPRE GPU build does not contain the custom HeatSolve
            # inner-circuit hook.  Use the portable external circuit UDF for
            # HYPRE CPU/GPU observables, while direct transient cases continue
            # to use the fully coupled custom HeatSolve route.
            hypre["heat_source"] = "circuit_parallel"
            hypre["parallel_circuit_iterations"] = 1
            hypre["solver"].update(
                {
                    "linear_system": "iterative_hypre_flexgmres_boomeramg",
                    "linear_system_max_iterations": 1000,
                    "linear_system_convergence_tolerance": tolerance,
                    "linear_system_abort_not_converged": True,
                    "omit_hypre_gpu_when_false": True,
                }
            )
            project["cases"][hypre_name] = hypre
            if tolerance == 1.0e-7:
                gpu_name = hypre_name + "_gpu"
                gpu = copy.deepcopy(hypre)
                gpu["solver"] = dict(gpu["solver"])
                gpu["solver"]["linear_system"] = (
                    "iterative_hypre_flexgmres_boomeramg_gpu"
                )
                gpu["state_file"] = f"work/meshes/{mesh_name}/{gpu_name}.state"
                gpu["series_file"] = f"{gpu_name}_series.csv"
                gpu["iteration_series_file"] = f"{gpu_name}_iterations.csv"
                gpu["output_file_path"] = (
                    f"../work/meshes/{mesh_name}/{gpu_name}.result"
                )
                project["cases"][gpu_name] = gpu
    project["cases"] = {
        key: value
        for key, value in project["cases"].items()
        if key == case_name
        or key.startswith(f"case_tes_transient_physical_parity_{route}_")
        or key.startswith("case_tes_steady_physical_parity_conformal_hypre_")
    }
    output = ROOT / output_name
    output.write_text(json.dumps(project, indent=2) + chr(10), encoding="utf-8")
    return output


def main() -> int:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    outputs = [
        _project(
            source,
            route="conformal",
            mesh_name="mesh_physical_parity_conformal",
            output_name="elmer_project_physical_parity_conformal.json",
        ),
        _project(
            source,
            route="mortar",
            mesh_name="mesh_physical_parity_mortar",
            output_name="elmer_project_physical_parity_mortar.json",
        ),
    ]
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

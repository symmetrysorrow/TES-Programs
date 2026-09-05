"""Prepare a strict common-window Mortar/conformal replacement benchmark.

The two formulations are kept in separate JSON files because their global
interface overrides differ.  They share the physical-parity geometry,
parameters, pulse, and one/ seven-step time grid.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MORTAR_SOURCE = ROOT / "elmer_project_physical_parity_mortar.json"
CONFORMAL_SOURCE = ROOT / "elmer_project_physical_parity_conformal.json"
MORTAR_OUTPUT = ROOT / "elmer_project_phase21_mortar_replacement.json"
CONFORMAL_OUTPUT = ROOT / "elmer_project_phase21_conformal_replacement.json"


def rename_case(case: dict, *, name: str, mesh: str, restart_from: str | None) -> dict:
    result = copy.deepcopy(case)
    result["mesh"] = mesh
    result["restart_from"] = restart_from
    result["state_file"] = f"work/meshes/{mesh}/{name}.state"
    result["series_file"] = f"{name}_series.csv"
    result["iteration_series_file"] = f"{name}_iterations.csv"
    result["output_file_path"] = f"../work/meshes/{mesh}/{name}.result"
    result["vtu"] = False
    result["output_result"] = True
    result["max_output_level"] = 6
    if restart_from:
        result["restart_file_path"] = f"../work/meshes/{mesh}/{restart_from}.result"
    else:
        result.pop("restart_file_path", None)
    return result


def hypre_solver(base: dict, *, gpu: bool) -> dict:
    result = copy.deepcopy(base)
    result.update(
        {
            "linear_system": (
                "iterative_hypre_flexgmres_boomeramg_gpu"
                if gpu
                else "iterative_hypre_flexgmres_boomeramg"
            ),
            "linear_system_max_iterations": 5000,
            "linear_system_convergence_tolerance": 1.0e-8,
            "linear_system_abort_not_converged": True,
            "linear_system_residual_output": 6,
            "omit_hypre_gpu_when_false": True,
        }
    )
    return result


def main() -> int:
    mortar = json.loads(MORTAR_SOURCE.read_text(encoding="utf-8"))
    conformal = json.loads(CONFORMAL_SOURCE.read_text(encoding="utf-8"))
    mortar_steady = "case_phase21_mortar_steady"
    mortar_cases = {
        mortar_steady: rename_case(
            mortar["cases"]["case_tes_steady_physical_parity_mortar"],
            name=mortar_steady,
            mesh="mesh_physical_parity_mortar",
            restart_from=None,
        ),
    }
    for window in ("1step", "7step"):
        name = f"case_phase21_mortar_{window}"
        mortar_cases[name] = rename_case(
            mortar["cases"][f"case_tes_transient_physical_parity_mortar_{window}"],
            name=name,
            mesh="mesh_physical_parity_mortar",
            restart_from=mortar_steady,
        )
    mortar["cases"] = mortar_cases
    mortar["performance_metadata"] = {
        "phase": "21",
        "role": "validated CPU Mortar replacement reference",
        "common_windows": ["steady", "1step", "7step"],
    }
    MORTAR_OUTPUT.write_text(json.dumps(mortar, indent=2) + "\n", encoding="utf-8")

    conformal_cases: dict[str, dict] = {}
    for backend, gpu in (("cpu", False), ("gpu", True)):
        steady_name = f"case_phase21_conformal_{backend}_steady"
        base_steady = conformal["cases"]["case_tes_steady_physical_parity_conformal"]
        steady = rename_case(
            base_steady,
            name=steady_name,
            mesh="mesh_physical_parity_conformal",
            restart_from=None,
        )
        steady["solver"] = hypre_solver(steady["solver"], gpu=gpu)
        conformal_cases[steady_name] = steady
        for window in ("1step", "7step"):
            name = f"case_phase21_conformal_{backend}_{window}"
            transient = rename_case(
                conformal["cases"][f"case_tes_transient_physical_parity_conformal_{window}"],
                name=name,
                mesh="mesh_physical_parity_conformal",
                restart_from=steady_name,
            )
            transient["solver"] = hypre_solver(transient["solver"], gpu=gpu)
            conformal_cases[name] = transient
    conformal["cases"] = conformal_cases
    conformal["performance_metadata"] = {
        "phase": "21",
        "role": "conformal shared-node CPU/GPU replacement candidate",
        "common_windows": ["steady", "1step", "7step"],
        "mortar_project": str(MORTAR_OUTPUT.name),
    }
    CONFORMAL_OUTPUT.write_text(json.dumps(conformal, indent=2) + "\n", encoding="utf-8")
    print(MORTAR_OUTPUT)
    print(CONFORMAL_OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

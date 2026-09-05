"""Prepare identical production-size HYPRE CPU/GPU benchmark cases."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
OUTPUT = ROOT / "elmer_project_phase20_production_hypre_benchmark.json"
MESH = "mesh_singlepixel_conformal_gpu_fine"


def make_case(source: dict, *, name: str, gpu: bool) -> dict:
    case = copy.deepcopy(source)
    case.update(
        {
            "mesh": MESH,
            "heat_source": "circuit_parallel",
            "parallel_circuit_iterations": 1,
            "apply_mortar_bcs": False,
            "output_result": True,
            "vtu": False,
            "state_file": f"work/meshes/{MESH}/{name}.state",
            "series_file": f"{name}_series.csv",
            "iteration_series_file": f"{name}_iterations.csv",
            "output_file_path": f"../work/meshes/{MESH}/{name}.result",
        }
    )
    case["solver"] = {
        **case["solver"],
        "linear_system": (
            "iterative_hypre_flexgmres_boomeramg_gpu"
            if gpu
            else "iterative_hypre_flexgmres_boomeramg"
        ),
        "linear_system_max_iterations": 2000,
        "linear_system_convergence_tolerance": 1.0e-8,
        "linear_system_abort_not_converged": True,
        "omit_hypre_gpu_when_false": True,
    }
    return case


def main() -> int:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    cpu_name = "case_phase20_production_hypre_cpu"
    gpu_name = "case_phase20_production_hypre_gpu"
    project["cases"] = {
        cpu_name: make_case(
            project["cases"]["case_tes_steady_singlepixel_conformal_hypre_cpu"],
            name=cpu_name,
            gpu=False,
        ),
        gpu_name: make_case(
            project["cases"]["case_tes_steady_singlepixel_conformal_hypre_gpu"],
            name=gpu_name,
            gpu=True,
        ),
    }
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

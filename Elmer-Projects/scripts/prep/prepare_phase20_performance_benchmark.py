"""Prepare Phase20 performance benchmark and transient-prefix projects.

The performance cases keep the validated conformal physics and differ only in
mesh size, HYPRE backend, repetition tag, or transient prefix length.  The
HYPRE residual verbosity is deliberately raised so the native setup timer and
iteration count are present in the solver log.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
MESH_SOURCE = ROOT / "elmer_project_stycast_convergence.json"
OUTPUT = ROOT / "elmer_project_phase20_performance_benchmark.json"
TRANSIENT_OUTPUT = ROOT / "elmer_project_phase20_performance_transient.json"

MESHES = {
    "small": "mesh_physical_parity_conformal",
    "production": "mesh_singlepixel_conformal_gpu_fine",
    "medium": "mesh_stycast_convergence_medium",
    "fine": "mesh_stycast_convergence_fine",
}


def mesh_registry(source: dict) -> dict:
    convergence = json.loads(MESH_SOURCE.read_text(encoding="utf-8"))
    meshes = copy.deepcopy(source["meshes"])
    for name in (
        "mesh_physical_parity_conformal",
        "mesh_stycast_convergence_medium",
        "mesh_stycast_convergence_fine",
    ):
        meshes[name] = copy.deepcopy(convergence["meshes"][name])
    return meshes


def solver_config(base: dict, *, gpu: bool) -> dict:
    result = copy.deepcopy(base)
    result.update(
        {
            "linear_system": (
                "iterative_hypre_flexgmres_boomeramg_gpu"
                if gpu
                else "iterative_hypre_flexgmres_boomeramg"
            ),
            "linear_system_max_iterations": 2000,
            "linear_system_convergence_tolerance": 1.0e-8,
            "linear_system_abort_not_converged": True,
            # HYPRE verbosity >= 6 emits setup time, iterations, and residual.
            "linear_system_residual_output": 6,
            "omit_hypre_gpu_when_false": True,
        }
    )
    return result


def steady_case(base: dict, *, name: str, mesh: str, gpu: bool) -> dict:
    case = copy.deepcopy(base)
    case.update(
        {
            "template": "steady",
            "mesh": mesh,
            "heat_source": "circuit_parallel",
            "parallel_circuit_iterations": 1,
            "apply_mortar_bcs": False,
            "max_output_level": 6,
            "output_result": True,
            "vtu": False,
            "steady_state_max_iterations": 1,
            "state_file": f"work/meshes/{mesh}/{name}.state",
            "series_file": f"{name}_series.csv",
            "iteration_series_file": f"{name}_iterations.csv",
            "output_file_path": f"../work/meshes/{mesh}/{name}.result",
            "solver": solver_config(case["solver"], gpu=gpu),
            "performance_role": "steady_size_crossover_or_repeat",
        }
    )
    return case


def main() -> int:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    production = source["cases"]["case_tes_steady_singlepixel_conformal_gpu_fine"]
    project = copy.deepcopy(source)
    project["meshes"] = mesh_registry(source)
    cases: dict[str, dict] = {}
    for size, mesh in MESHES.items():
        for backend, gpu in (("cpu", False), ("gpu", True)):
            for repeat in range(1, 4):
                name = f"case_phase20_perf_{size}_{backend}_r{repeat}"
                cases[name] = steady_case(
                    production, name=name, mesh=mesh, gpu=gpu
                )
    project["cases"] = cases
    project["performance_metadata"] = {
        "mesh_points": MESHES,
        "repeats": 3,
        "backends": ["cpu", "gpu"],
        "tolerance": 1.0e-8,
        "purpose": "wall-time breakdown, repeated timing, and GPU crossover",
    }
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")

    pulse_source = source["cases"]["case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step"]
    expanded_timesteps = []
    for dt, count in pulse_source["timesteps"]:
        expanded_timesteps.extend([[dt, 1] for _ in range(int(count))])
    transient = copy.deepcopy(source)
    transient["meshes"] = mesh_registry(source)
    transient_cases: dict[str, dict] = {}
    existing_restart = {
        "cpu": "../work/meshes/mesh_singlepixel_conformal_gpu_fine/case_phase20_production_hypre_cpu.result",
        "gpu": "../work/meshes/mesh_singlepixel_conformal_gpu_fine/case_phase20_production_hypre_gpu.result",
    }
    for backend, gpu in (("cpu", False), ("gpu", True)):
        for prefix, step_count in (("7step", 7), ("50step", 50)):
            name = f"case_phase20_perf_transient_{backend}_{prefix}"
            case = copy.deepcopy(pulse_source)
            case.update(
                {
                    "template": "pulse",
                    "mesh": MESHES["production"],
                    "heat_source": "circuit_parallel",
                    "parallel_circuit_iterations": 1,
                    "apply_mortar_bcs": False,
                    "max_output_level": 6,
                    "restart_from": None,
                    "restart_file_base": "external_phase20_steady",
                    "restart_file_path": existing_restart[backend],
                    "state_file": f"work/meshes/{MESHES['production']}/{name}.state",
                    "series_file": f"{name}_series.csv",
                    "iteration_series_file": f"{name}_iterations.csv",
                    "output_file_path": f"../work/meshes/{MESHES['production']}/{name}.result",
                    # The source time grid is grouped (7 groups, 176 actual
                    # steps).  Performance prefixes are deliberately actual
                    # timestep counts, so expand it before truncating.
                    "timesteps": copy.deepcopy(expanded_timesteps[:step_count]),
                    "output_intervals": [1] * step_count,
                    "solver": solver_config(pulse_source["solver"], gpu=gpu),
                    "performance_role": "transient_prefix",
                }
            )
            transient_cases[name] = case
        reuse_name = f"case_phase20_perf_transient_{backend}_7step_reuse_probe"
        reuse = copy.deepcopy(transient_cases[f"case_phase20_perf_transient_{backend}_7step"])
        reuse["series_file"] = f"{reuse_name}_series.csv"
        reuse["iteration_series_file"] = f"{reuse_name}_iterations.csv"
        reuse["state_file"] = f"work/meshes/{MESHES['production']}/{reuse_name}.state"
        reuse["output_file_path"] = f"../work/meshes/{MESHES['production']}/{reuse_name}.result"
        reuse["solver"] = copy.deepcopy(reuse["solver"])
        reuse["solver"]["linear_system_refactorize"] = False
        reuse["performance_role"] = "transient_unsafe_reuse_probe"
        transient_cases[reuse_name] = reuse
    transient["cases"] = transient_cases
    transient["performance_metadata"] = {
        "mesh": MESHES["production"],
        "prefixes": [7, 50],
        "backends": ["cpu", "gpu"],
        "tolerance": 1.0e-8,
        "restart_inputs": existing_restart,
        "purpose": "production transient CPU/GPU prefix timing and setup reuse study",
    }
    TRANSIENT_OUTPUT.write_text(
        json.dumps(transient, indent=2) + "\n", encoding="utf-8"
    )
    print(OUTPUT)
    print(TRANSIENT_OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

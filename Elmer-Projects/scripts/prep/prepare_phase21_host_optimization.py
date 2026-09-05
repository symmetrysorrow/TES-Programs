"""Prepare Phase21 host-overhead, fine-transient, and I/O benchmarks.

The cases keep the Phase20 conformal/shared-node physics and HYPRE settings.
They differ only in the host-side UDF cache revision, output switches, mesh
size, and backend.  The script deliberately does not add an unconditional
HYPRE hierarchy reuse option: Phase20 already rejected that experiment on
convergence grounds.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
PHASE20 = ROOT / "elmer_project_phase20_performance_benchmark.json"
OUTPUT = ROOT / "elmer_project_phase21_host_optimization.json"

PRODUCTION_MESH = "mesh_singlepixel_conformal_gpu_fine"
FINE_MESH = "mesh_stycast_convergence_fine"


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
            "linear_system_residual_output": 6,
            "omit_hypre_gpu_when_false": True,
        }
    )
    return result


def expanded_grid(pulse_source: dict, count: int) -> list[list[object]]:
    expanded: list[list[object]] = []
    for dt, stages in pulse_source["timesteps"]:
        expanded.extend([[dt, 1] for _ in range(int(stages))])
    return copy.deepcopy(expanded[:count])


def transient_case(
    pulse_source: dict,
    *,
    name: str,
    mesh: str,
    gpu: bool,
    steps: int,
    restart_from: str | None = None,
    external_restart: str | None = None,
    output_result: bool = True,
    vtu: object = False,
    write_series: bool = True,
    write_iteration_series: bool = True,
    role: str,
) -> dict:
    case = copy.deepcopy(pulse_source)
    case.update(
        {
            "template": "pulse",
            "mesh": mesh,
            "heat_source": "circuit_parallel",
            "parallel_circuit_iterations": 1,
            "apply_mortar_bcs": False,
            "max_output_level": 6,
            "restart_from": restart_from,
            "restart_file_path": (
                external_restart
                if external_restart
                else (
                    f"../work/meshes/{mesh}/{restart_from}.result"
                    if restart_from
                    else None
                )
            ),
            "restart_file_base": (
                "external_phase21_steady" if external_restart else restart_from
            ),
            "state_file": f"work/meshes/{mesh}/{name}.state",
            "series_file": f"{name}_series.csv",
            "iteration_series_file": f"{name}_iterations.csv",
            "output_file_path": f"../work/meshes/{mesh}/{name}.result",
            "timesteps": expanded_grid(pulse_source, steps),
            "output_intervals": [1] * steps,
            "output_result": output_result,
            "vtu": vtu,
            "write_series": write_series,
            "write_iteration_series": write_iteration_series,
            "solver": solver_config(pulse_source["solver"], gpu=gpu),
            "performance_role": role,
        }
    )
    if restart_from is not None:
        case.pop("preexisting_restart", None)
    return case


def main() -> int:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    phase20 = json.loads(PHASE20.read_text(encoding="utf-8"))
    pulse_source = source["cases"]["case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step"]
    production_restart = {
        "cpu": "../work/meshes/mesh_singlepixel_conformal_gpu_fine/case_phase20_production_hypre_cpu.result",
        "gpu": "../work/meshes/mesh_singlepixel_conformal_gpu_fine/case_phase20_production_hypre_gpu.result",
    }

    project = copy.deepcopy(source)
    project["meshes"] = copy.deepcopy(phase20["meshes"])
    project["elmer_overrides"] = {
        **project.get("elmer_overrides", {}),
        "conformal_shared_interfaces": False,
        "conformal_mortar_interfaces": False,
        "conformal_shared_node_interfaces": True,
        "mortar_control_interfaces": False,
        "fragment_mortar_interfaces": True,
    }
    cases: dict[str, dict] = {}

    # Re-run the production prefix after rebuilding the cached UDFs.  The
    # Phase20 steady results are the fixed restart input, so this isolates the
    # host-side change from steady-state initialization.
    for backend, gpu in (("cpu", False), ("gpu", True)):
        for prefix, steps in (("7step", 7), ("50step", 50)):
            name = f"case_phase21_host_transient_{backend}_{prefix}"
            cases[name] = transient_case(
                pulse_source,
                name=name,
                mesh=PRODUCTION_MESH,
                gpu=gpu,
                steps=steps,
                external_restart=production_restart[backend],
                role="host_cache_optimized_transient",
            )

    # Fine transient crossover: construct a matching fine-mesh steady restart
    # in this same project, then use short 7/20-step windows.
    for backend, gpu in (("cpu", False), ("gpu", True)):
        steady_name = f"case_phase21_fine_steady_{backend}"
        steady_base = phase20["cases"][f"case_phase20_perf_fine_{backend}_r1"]
        steady = copy.deepcopy(steady_base)
        steady.update(
            {
                "mesh": FINE_MESH,
                "state_file": f"work/meshes/{FINE_MESH}/{steady_name}.state",
                "series_file": f"{steady_name}_series.csv",
                "iteration_series_file": f"{steady_name}_iterations.csv",
                "output_file_path": f"../work/meshes/{FINE_MESH}/{steady_name}.result",
                "performance_role": "fine_transient_restart",
            }
        )
        cases[steady_name] = steady
        for prefix, steps in (("7step", 7), ("20step", 20)):
            name = f"case_phase21_fine_transient_{backend}_{prefix}"
            cases[name] = transient_case(
                pulse_source,
                name=name,
                mesh=FINE_MESH,
                gpu=gpu,
                steps=steps,
                restart_from=steady_name,
                role="fine_transient_crossover",
            )

    # I/O matrix on the production 7-step window.  ``full_io`` turns VTU on;
    # the other rows remove exactly one output class at a time.  Results are
    # intentionally retained in the full and no-series rows, so the report
    # cannot mistake an output-disabled run for a production replacement.
    io_modes = {
        "full_io": {"output_result": True, "vtu": "after_timestep", "write_series": True, "write_iteration_series": True},
        "no_vtu": {"output_result": True, "vtu": False, "write_series": True, "write_iteration_series": True},
        "no_result": {"output_result": False, "vtu": False, "write_series": True, "write_iteration_series": True},
        "no_iteration_csv": {"output_result": True, "vtu": False, "write_series": True, "write_iteration_series": False},
        "no_series_csv": {"output_result": True, "vtu": False, "write_series": False, "write_iteration_series": False},
    }
    for backend, gpu in (("cpu", False), ("gpu", True)):
        for mode, switches in io_modes.items():
            name = f"case_phase21_io_{backend}_{mode}"
            cases[name] = transient_case(
                pulse_source,
                name=name,
                mesh=PRODUCTION_MESH,
                gpu=gpu,
                steps=7,
                external_restart=production_restart[backend],
                role="io_overhead_matrix",
                **switches,
            )

    project["cases"] = cases
    project["performance_metadata"] = {
        "phase": "21",
        "baseline": "Phase20 artifacts and fixed Phase20 steady restarts",
        "backends": ["cpu", "gpu"],
        "transient_prefixes": [7, 20, 50],
        "io_modes": sorted(io_modes),
        "physics_policy": "no model simplification; same conformal/shared-node pulse and HYPRE tolerance",
        "reuse_policy": "unconditional Linear System Refactorize=False remains rejected",
    }
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

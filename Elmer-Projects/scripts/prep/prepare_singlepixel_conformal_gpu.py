"""Create an isolated all-tetra, conformal single-pixel GPU project."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_comsol_timegrid.json"
OUTPUT = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
MESH = "mesh_singlepixel_conformal_gpu"
STEADY = "case_tes_steady_singlepixel_conformal_gpu"
PULSE = "case_tes_pulse_singlepixel_conformal_gpu_original_timegrid"
SMOKE = f"{PULSE}_amgx_smoke_7step"
FINE_MESH = "mesh_singlepixel_conformal_gpu_fine"
FINE_STEADY = "case_tes_steady_singlepixel_conformal_gpu_fine"
FINE_PULSE = "case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step"
REFERENCE_STEADY = "case_tes_steady_singlepixel_conformal_mortar_reference"
HYPRE_CPU_STEADY = "case_tes_steady_singlepixel_conformal_hypre_cpu"
HYPRE_GPU_STEADY = "case_tes_steady_singlepixel_conformal_hypre_gpu"

PULSE_PREFIX = [
    ["18[us]", 1], ["1[us]", 2], ["1[ns]", 1], ["10[ns]", 10],
    ["100[ns]", 9], ["1[us]", 9], ["0.625[us]", 144],
]


def timestep_seconds(token: str) -> float:
    value, unit = token[:-1].split("[")
    return float(value) * {"s": 1.0, "ms": 1e-3, "us": 1e-6, "ns": 1e-9}[unit]


def hybrid_timesteps(original: list[list[object]]) -> list[list[object]]:
    early = copy.deepcopy(PULSE_PREFIX)
    early_end = sum(timestep_seconds(str(token)) * int(count) for token, count in early)
    cumulative = 0.0
    for index, (token, count) in enumerate(original):
        group_end = cumulative + timestep_seconds(str(token)) * int(count)
        if group_end > early_end:
            result = early
            bridge = group_end - early_end
            if bridge > 1e-18:
                result.append([f"{bridge:.17g}[s]", 1])
            result.extend(copy.deepcopy(original[index + 1 :]))
            return result
        cumulative = group_end
    raise ValueError("original time grid ends before hybrid rise grid")


def truncate_timesteps(timesteps: list[list[object]], max_steps: int) -> list[list[object]]:
    remaining = max_steps
    result: list[list[object]] = []
    for token, count in timesteps:
        take = min(int(count), remaining)
        if take:
            result.append([token, take])
            remaining -= take
        if remaining == 0:
            return result
    raise ValueError(f"time grid contains fewer than {max_steps} steps")


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    project.setdefault("elmer_overrides", {})[
        "conformal_shared_interfaces"
    ] = False
    project["elmer_overrides"]["conformal_mortar_interfaces"] = False
    project["elmer_overrides"]["conformal_shared_node_interfaces"] = True
    project["elmer_overrides"]["fragment_mortar_interfaces"] = True

    mesh = copy.deepcopy(project["meshes"]["mesh_refined_3x"])
    mesh["dir"] = MESH
    mesh["recipe"]["elmergrid_args"] = [
        "14", "2", "gmsh/project.msh", "-merge", "1e-10", "-out", MESH,
    ]
    mesh["notes"] = (
        "Independent all-tetra 3x shared-node experiment.  The three contact "
        "footprints are imprinted and face meshes are paired before ElmerGrid "
        "node merging; this route must pass the post-conversion connectivity gate."
    )
    project["meshes"][MESH] = mesh
    fine_mesh = copy.deepcopy(mesh)
    fine_mesh["dir"] = FINE_MESH
    fine_mesh["recipe"]["mesh_overrides"] = {
        "mesh_min": 4.5e-5,
        "mesh_max": 9.0e-5,
        "mesh_min_mode": "fixed",
        "mesh_max_mode": "fixed",
    }
    fine_mesh["recipe"]["elmergrid_args"][-1] = FINE_MESH
    fine_mesh["notes"] = (
        "Fine shared-node contact mesh targeting production-v2-like spatial "
        "resolution; generated independently from the Mortar production path."
    )
    project["meshes"][FINE_MESH] = fine_mesh

    steady = {
        "template": "steady",
        "mesh": MESH,
        "heat_source": "circuit_inner",
        "initial_temperature": "T_0",
        "output_result": True,
        "vtu": False,
        "steady_state_max_iterations": 1,
        "output_intervals": 1,
        "apply_mortar_bcs": False,
        "solver": {
            "nonlinear_max_iterations": 120,
            "nonlinear_convergence_tolerance": 1e-8,
            "nonlinear_relaxation_factor": 1.0,
            "steady_state_convergence_tolerance": 1e-8,
            # The stock Windows Elmer 26.1 binary used for this validation
            # does not ship MUMPS; use the direct UMFPACK backend for the
            # CPU parity gate and keep the solver class (direct) explicit.
            "linear_system": "direct",
        },
        "state_file": f"work/meshes/{MESH}/{STEADY}.state",
        "series_file": f"{STEADY}_series.csv",
        "iteration_series_file": f"{STEADY}_iterations.csv",
        "output_file_path": f"../work/meshes/{MESH}/{STEADY}.result",
    }
    project["cases"][STEADY] = steady

    reference = copy.deepcopy(steady)
    reference.update(
        {
            "mesh": "mesh_refined_3x",
            "apply_mortar_bcs": True,
            "state_file": f"work/meshes/mesh_refined_3x/{REFERENCE_STEADY}.state",
            "series_file": f"{REFERENCE_STEADY}_series.csv",
            "iteration_series_file": f"{REFERENCE_STEADY}_iterations.csv",
            "output_file_path": f"../work/meshes/mesh_refined_3x/{REFERENCE_STEADY}.result",
            "solver": dict(steady["solver"]),
        }
    )
    reference["solver"]["nonlinear_convergence_tolerance"] = 1.0e-6
    reference["solver"]["steady_state_convergence_tolerance"] = 1.0e-6
    project["cases"][REFERENCE_STEADY] = reference

    for case_name, linear_system in (
        (HYPRE_CPU_STEADY, "iterative_hypre_flexgmres_boomeramg"),
        (HYPRE_GPU_STEADY, "iterative_hypre_flexgmres_boomeramg_gpu"),
    ):
        iterative = copy.deepcopy(steady)
        iterative["solver"] = dict(steady["solver"])
        iterative["solver"].update(
            {
                "linear_system": linear_system,
                "linear_system_max_iterations": 1000,
                # This primal mesh reaches 1.94e-8 in the stock HYPRE
                # FlexGMRES/AMG build at 1000 iterations; keep the first
                # conformal gate honest and attainable, then tighten in a
                # separate convergence study.
                "linear_system_convergence_tolerance": 1.0e-7,
                "linear_system_abort_not_converged": True,
                # The stock Windows Elmer build does not register the
                # HYPRE GPU keyword.  Omitting an explicit false value keeps
                # the CPU case portable; the GPU case still emits True.
                "omit_hypre_gpu_when_false": True,
            }
        )
        iterative.update(
            {
                "state_file": f"work/meshes/{MESH}/{case_name}.state",
                "series_file": f"{case_name}_series.csv",
                "iteration_series_file": f"{case_name}_iterations.csv",
                "output_file_path": f"../work/meshes/{MESH}/{case_name}.result",
            }
        )
        project["cases"][case_name] = iterative

    pulse = copy.deepcopy(
        project["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]
    )
    pulse.update(
        {
            "mesh": MESH,
            "restart_from": STEADY,
            "restart_time": 0.02,
            "restart_file_path": f"../work/meshes/{MESH}/{STEADY}.result",
            "state_file": f"work/meshes/{MESH}/{STEADY}.state",
            "series_file": f"{PULSE}_series.csv",
            "iteration_series_file": f"{PULSE}_iterations.csv",
            "output_file_path": f"../work/meshes/{MESH}/{PULSE}.result",
            "vtu": False,
            "apply_mortar_bcs": False,
            "comparison_time_grid": {
                "source_case": "case_tes_mpi_comsol_grid_full_uniform_continuous",
                "source_project": SOURCE.name,
                "description": "Conformal GPU mesh with the original COMSOL time grid.",
            },
        }
    )
    pulse["solver"] = dict(pulse["solver"])
    pulse["solver"]["linear_system"] = "direct"
    project["cases"][PULSE] = pulse

    smoke = copy.deepcopy(pulse)
    smoke["timesteps"] = truncate_timesteps(smoke["timesteps"], 7)
    smoke["output_intervals"] = [999999] * len(smoke["timesteps"])
    smoke["output_intervals"][-1] = 1
    smoke["series_file"] = f"{SMOKE}_series.csv"
    smoke["iteration_series_file"] = f"{SMOKE}_iterations.csv"
    smoke["output_file_path"] = f"../work/meshes/{MESH}/{SMOKE}.result"
    project["cases"][SMOKE] = smoke

    fine_steady = copy.deepcopy(steady)
    fine_steady.update(
        {
            "mesh": FINE_MESH,
            "state_file": f"work/meshes/{FINE_MESH}/{FINE_STEADY}.state",
            "series_file": f"{FINE_STEADY}_series.csv",
            "iteration_series_file": f"{FINE_STEADY}_iterations.csv",
            "output_file_path": f"../work/meshes/{FINE_MESH}/{FINE_STEADY}.result",
        }
    )
    project["cases"][FINE_STEADY] = fine_steady

    fine_pulse = copy.deepcopy(pulse)
    fine_pulse.update(
        {
            "mesh": FINE_MESH,
            "restart_from": FINE_STEADY,
            "restart_file_path": f"../work/meshes/{FINE_MESH}/{FINE_STEADY}.result",
            "state_file": f"work/meshes/{FINE_MESH}/{FINE_STEADY}.state",
            "series_file": f"{FINE_PULSE}_series.csv",
            "iteration_series_file": f"{FINE_PULSE}_iterations.csv",
            "output_file_path": f"../work/meshes/{FINE_MESH}/{FINE_PULSE}.result",
        }
    )
    fine_pulse["timesteps"] = truncate_timesteps(
        hybrid_timesteps(fine_pulse["timesteps"]), 177
    )
    fine_pulse["output_intervals"] = [999999] * len(fine_pulse["timesteps"])
    fine_pulse["output_intervals"][-1] = 1
    project["cases"][FINE_PULSE] = fine_pulse

    # Elmer prepends the mesh-name component to result paths.  The existing
    # production workflow uses this empty traversal anchor as well.
    (ROOT / MESH).mkdir(exist_ok=True)
    (ROOT / FINE_MESH).mkdir(exist_ok=True)
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()

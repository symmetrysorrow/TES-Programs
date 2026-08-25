"""Diagnostic: isolate global_mesh_size_um's effect on the decay tail.

Identical to production-v2 except global_mesh_size_um is reverted to the
legacy 50 um. Also drops the validated-but-unnecessary-for-this-check 0.625 us
mid section, jumping straight from the pulse prefix into the same 10/100 us
/ 1 ms tail schedule that already matched COMSOL closely on mesh_refined_3x
(case_tes_mpi_comsol_grid_full_uniform_continuous). This roughly halves the
step count (129 vs 273) for a faster A/B read.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"

MESH_NAME = "mesh_singlepixel_prod_v2_g50"
MESH_FILE = "project_hybrid_prism_prod_v2_g50.msh"
STEADY_NAME = "case_tes_steady_singlepixel_prod_v2_g50"
PULSE_NAME = "case_tes_pulse_singlepixel_prod_v2_g50_fast"
SOURCE_STEADY = "case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight"
SOURCE_PULSE = "case_p19_pulse_phase23_tight"

PULSE_PREFIX = [["18[us]", 1], ["1[us]", 2], ["1[ns]", 1], ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9]]
# No 0.625 us mid section this time: go straight into the proven fast tail.
TAIL = [["10[us]", 9], ["100[us]", 9], ["1[ms]", 79]]


def main() -> None:
    project = json.loads(PROJECT.read_text(encoding="utf-8"))

    gen_command = (
        "python generate_hybrid_prism_geometry.py elmer_project_comsol_timegrid.json "
        f"--output gmsh/{MESH_FILE} "
        "--global-mesh-size 50e-6 "
        "--stack-local-size 12.5e-6 --stack-local-half-width 400e-6 "
        "--absorber-local-size 25e-6 --absorber-local-radius 25e-6 "
        "--disable-mesh-size-extend-from-boundary --mesh-algorithm 5 "
        "--stycast-layers 32 --tes-layers 1 --sinx-layers 1 "
        "--sio2-1-layers 2 --si-1-layers 2 --si-2-layers 2"
    )
    grid_command = f"ElmerGrid 14 2 gmsh/{MESH_FILE} -merge 1e-10 -out {MESH_NAME}"
    project["meshes"][MESH_NAME] = {
        "geometry": "single_pixel",
        "dir": MESH_NAME,
        "recipe": {
            "generator": "generate_hybrid_prism_geometry.py",
            "commands": [gen_command, grid_command],
        },
        "notes": (
            "Diagnostic mesh: same as mesh_singlepixel_prod_v2 but "
            "global_mesh_size_um reverted 100->50 to isolate its effect on the "
            "20-80ms decay-tail mismatch against COMSOL."
        ),
    }

    steady = copy.deepcopy(project["cases"][SOURCE_STEADY])
    steady.update({
        "mesh": MESH_NAME,
        "state_file": f"{MESH_NAME}/prod_v2_g50_steady.state",
        "output_file_path": f"../work/meshes/{MESH_NAME}/{STEADY_NAME}.result",
        "series_file": f"{STEADY_NAME}_series.csv",
        "iteration_series_file": f"{STEADY_NAME}_iterations.csv",
    })
    project["cases"][STEADY_NAME] = steady

    pulse = copy.deepcopy(project["cases"][SOURCE_PULSE])
    timesteps = [*PULSE_PREFIX, *TAIL]
    pulse.update({
        "mesh": MESH_NAME,
        "restart_from": STEADY_NAME,
        "restart_file_path": f"../work/meshes/{MESH_NAME}/{STEADY_NAME}.result",
        "state_file": f"{MESH_NAME}/prod_v2_g50_steady.state",
        "output_file_path": f"../work/meshes/{MESH_NAME}/{PULSE_NAME}.result",
        "series_file": f"{PULSE_NAME}_series.csv",
        "iteration_series_file": f"{PULSE_NAME}_iterations.csv",
        "timesteps": timesteps,
        "output_intervals": [1] * len(timesteps),
    })
    pulse.pop("solver_comment", None)
    pulse.pop("comparison_time_grid", None)
    project["cases"][PULSE_NAME] = pulse

    PROJECT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"added mesh={MESH_NAME!r} steady={STEADY_NAME!r} pulse={PULSE_NAME!r}")
    print(f"pulse timesteps ({len(timesteps)} entries):", timesteps)


if __name__ == "__main__":
    main()

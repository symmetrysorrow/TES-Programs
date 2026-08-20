"""One-off: add the validated production-v2 mesh, steady case, and full-decay
pulse case to elmer_project_hybrid_prism_phase23_pulse_tight.json.

Settings come from the singlepixel_resolution_optimization.json sensitivity
search (global_mesh_size_um=100, stack_local_size_um=12.5,
stack_local_half_width_um=400, absorber_local_size_um=25,
absorber_local_radius_um=25, stycast_layers=32, tes/sinx/sio2_1... layers).
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"

MESH_NAME = "mesh_singlepixel_prod_v2"
MESH_FILE = "project_hybrid_prism_prod_v2.msh"
STEADY_NAME = "case_tes_steady_singlepixel_prod_v2"
PULSE_NAME = "case_tes_pulse_singlepixel_prod_v2_full"
SOURCE_STEADY = "case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight"
SOURCE_PULSE = "case_p19_pulse_phase23_tight"

# PULSE_PREFIX from run_singlepixel_resolution_pilot.py, ending 10.001 us after the pulse.
PULSE_PREFIX = [["18[us]", 1], ["1[us]", 2], ["1[ns]", 1], ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9]]
# Validated 0.625 us step through the 100 us comparison window (10.001 -> 100.001 us).
FINE_STEP = ["0.625[us]", 144]
# Coarsening tail mirroring case_tes_pulse_20ms_3x_refined's established post-pulse schedule,
# reaching ~80 ms after the pulse (full return toward baseline).
TAIL = [["10[us]", 9], ["100[us]", 9], ["1[ms]", 79]]


def main() -> None:
    project = json.loads(PROJECT.read_text(encoding="utf-8"))

    gen_command = (
        "python generate_hybrid_prism_geometry.py elmer_project_comsol_timegrid.json "
        f"--output gmsh/{MESH_FILE} "
        "--global-mesh-size 100e-6 "
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
            "Production-v2 mesh: settings from the singlepixel_resolution_optimization.json "
            "sensitivity search (global=100um, stack_local=12.5um/half-width 400um, "
            "absorber_local=25um/radius 25um, stycast=32 layers, tes/sinx/sio2_1/si_1/si_2 "
            "as validated). Uses --mesh-algorithm 5 to match what was actually validated."
        ),
    }

    steady = copy.deepcopy(project["cases"][SOURCE_STEADY])
    steady.update({
        "mesh": MESH_NAME,
        "state_file": f"{MESH_NAME}/prod_v2_steady.state",
        # Elmer's default "Output File" resolution writes to a repo-root-relative
        # <meshname>/ directory, ignoring the "work/meshes" Mesh DB prefix. This
        # explicit override routes it back to where the mesh actually lives,
        # matching the pattern run_singlepixel_resolution_pilot.py already uses.
        "output_file_path": f"../work/meshes/{MESH_NAME}/{STEADY_NAME}.result",
        "series_file": f"{STEADY_NAME}_series.csv",
        "iteration_series_file": f"{STEADY_NAME}_iterations.csv",
    })
    project["cases"][STEADY_NAME] = steady

    pulse = copy.deepcopy(project["cases"][SOURCE_PULSE])
    timesteps = [*PULSE_PREFIX, FINE_STEP, *TAIL]
    pulse.update({
        "mesh": MESH_NAME,
        "restart_from": STEADY_NAME,
        "restart_file_path": f"../work/meshes/{MESH_NAME}/{STEADY_NAME}.result",
        "state_file": f"{MESH_NAME}/prod_v2_steady.state",
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

"""Prepare an 8-layer Stycast short-pulse validation against COMSOL.

The established Phase23 one-layer hybrid result is retained as the baseline.
This diagnostic changes only the number of prism elements through the 20 um
Stycast thickness and stops after pulse +225 us, which is sufficient to cover
the COMSOL and Elmer 10/90-percent rise crossings.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"
OUTPUT = ROOT / "elmer_project_stycast_z8_validation.json"

MESH = "mesh_hybrid_phase19_stycast_z8"
STEADY = "case_stycast_z8_steady_tight"
PULSE = "case_stycast_z8_pulse_225us_tight"
SOURCE_STEADY = "case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight"
SOURCE_PULSE = "case_p19_pulse_phase23_tight"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    project["meshes"][MESH] = {
        "geometry": "single_pixel",
        "dir": MESH,
        "recipe": {
            "generator": "generate_hybrid_prism_geometry.py",
            "commands": [
                "python generate_hybrid_prism_geometry.py "
                "elmer_project_comsol_timegrid.json "
                "--output gmsh/project_hybrid_phase19_stycast_z8.msh "
                "--stack-local-size 16.6666666666667e-6 "
                "--absorber-local-size 35e-6 --absorber-local-radius 50e-6 "
                "--disable-mesh-size-extend-from-boundary --stycast-layers 8",
                "ElmerGrid 14 2 gmsh/project_hybrid_phase19_stycast_z8.msh "
                f"-merge 1e-10 -out work/meshes/{MESH}",
            ],
        },
        "notes": (
            "Phase19 geometry and in-plane fields with eight equal prism "
            "layers through the 20 um Stycast thickness."
        ),
    }

    steady = copy.deepcopy(project["cases"][SOURCE_STEADY])
    steady.update({
        "mesh": MESH,
        "state_file": f"work/meshes/{MESH}/stycast_z8_steady.state",
        "output_file_path": f"../work/meshes/{MESH}/{STEADY}.result",
        "series_file": f"{STEADY}_series.csv",
        "iteration_series_file": f"{STEADY}_iterations.csv",
    })
    project["cases"][STEADY] = steady

    pulse = copy.deepcopy(project["cases"][SOURCE_PULSE])
    pulse.update({
        "mesh": MESH,
        "restart_from": STEADY,
        "restart_file_path": f"../work/meshes/{MESH}/{STEADY}.result",
        "state_file": f"work/meshes/{MESH}/stycast_z8_steady.state",
        "series_file": f"{PULSE}_series.csv",
        "iteration_series_file": f"{PULSE}_iterations.csv",
        "timesteps": [
            ["18[us]", 1], ["1[us]", 2], ["1[ns]", 1],
            ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9],
            ["10[us]", 9], ["5[us]", 25],
        ],
        "output_intervals": [1] * 8,
        "comparison_time_grid": {
            "mode": "Stycast through-thickness resolution diagnostic",
            "reference_case": SOURCE_PULSE,
            "only_changed_factor": "Stycast prism layers 1 -> 8",
            "end_after_pulse": "225[us]",
        },
    })
    project["cases"][PULSE] = pulse

    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()

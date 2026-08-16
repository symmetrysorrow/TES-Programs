"""Prepare a 32-layer Stycast spatial-resolution check at fixed 1.25 us dt."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"
OUTPUT = ROOT / "elmer_project_stycast_z32_validation.json"

MESH = "mesh_hybrid_phase19_stycast_z32"
STEADY = "case_stycast_z32_steady_tight"
PULSE = "case_stycast_z32_pulse_105us_fine1p25us_tight"
PULSE_FINE = "case_stycast_z32_pulse_105us_fine0p625us_tight"
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
                "--output gmsh/project_hybrid_phase19_stycast_z32.msh "
                "--stack-local-size 16.6666666666667e-6 "
                "--absorber-local-size 35e-6 --absorber-local-radius 50e-6 "
                "--disable-mesh-size-extend-from-boundary --stycast-layers 32",
                "ElmerGrid 14 2 gmsh/project_hybrid_phase19_stycast_z32.msh "
                f"-merge 1e-10 -out work/meshes/{MESH}",
            ],
        },
        "notes": "Phase19 geometry with thirty-two 0.625 um Stycast prism layers.",
    }

    steady = copy.deepcopy(project["cases"][SOURCE_STEADY])
    steady.update({
        "mesh": MESH,
        "state_file": f"work/meshes/{MESH}/stycast_z32_steady.state",
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
        "state_file": f"work/meshes/{MESH}/stycast_z32_steady.state",
        "series_file": f"{PULSE}_series.csv",
        "iteration_series_file": f"{PULSE}_iterations.csv",
        "timesteps": [
            ["18[us]", 1], ["1[us]", 2], ["1[ns]", 1],
            ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9], ["1.25[us]", 76],
        ],
        "output_intervals": [1] * 7,
        "comparison_time_grid": {
            "mode": "Stycast z32 spatial-resolution diagnostic",
            "reference_case": "case_stycast_z16_pulse_105us_fine1p25us_tight",
            "only_changed_factor": "Stycast prism layers 16 -> 32",
            "end_after_pulse": "105[us]",
        },
    })
    project["cases"][PULSE] = pulse

    pulse_fine = copy.deepcopy(pulse)
    pulse_fine.update({
        "series_file": f"{PULSE_FINE}_series.csv",
        "iteration_series_file": f"{PULSE_FINE}_iterations.csv",
        "timesteps": [
            ["18[us]", 1], ["1[us]", 2], ["1[ns]", 1],
            ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9], ["0.625[us]", 152],
        ],
        "comparison_time_grid": {
            "mode": "Stycast z32 time-convergence diagnostic",
            "reference_case": PULSE,
            "only_changed_factor": "post-pulse timestep 1.25 us -> 0.625 us",
            "end_after_pulse": "105[us]",
        },
    })
    project["cases"][PULSE_FINE] = pulse_fine

    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()

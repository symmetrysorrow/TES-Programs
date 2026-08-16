"""Prepare Stycast z16 spatial- and temporal-resolution diagnostics."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"
OUTPUT = ROOT / "elmer_project_stycast_z16_validation.json"

MESH = "mesh_hybrid_phase19_stycast_z16"
STEADY = "case_stycast_z16_steady_tight"
PULSE_COARSE = "case_stycast_z16_pulse_105us_tight"
PULSE_FINE = "case_stycast_z16_pulse_105us_fine5us_tight"
PULSE_FINER = "case_stycast_z16_pulse_105us_fine2p5us_tight"
PULSE_FINER_LONG = "case_stycast_z16_pulse_225us_fine2p5us_tight"
PULSE_FINEST = "case_stycast_z16_pulse_105us_fine1p25us_tight"
PULSE_FINEST_LONG = "case_stycast_z16_pulse_225us_fine1p25us_tight"
SOURCE_STEADY = "case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight"
SOURCE_PULSE = "case_p19_pulse_phase23_tight"


def pulse_case(project: dict, name: str, timesteps: list[list[object]]) -> dict:
    pulse = copy.deepcopy(project["cases"][SOURCE_PULSE])
    pulse.update({
        "mesh": MESH,
        "restart_from": STEADY,
        "restart_file_path": f"../work/meshes/{MESH}/{STEADY}.result",
        "state_file": f"work/meshes/{MESH}/stycast_z16_steady.state",
        "series_file": f"{name}_series.csv",
        "iteration_series_file": f"{name}_iterations.csv",
        "timesteps": timesteps,
        "output_intervals": [1] * len(timesteps),
    })
    return pulse


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
                "--output gmsh/project_hybrid_phase19_stycast_z16.msh "
                "--stack-local-size 16.6666666666667e-6 "
                "--absorber-local-size 35e-6 --absorber-local-radius 50e-6 "
                "--disable-mesh-size-extend-from-boundary --stycast-layers 16",
                "ElmerGrid 14 2 gmsh/project_hybrid_phase19_stycast_z16.msh "
                f"-merge 1e-10 -out work/meshes/{MESH}",
            ],
        },
        "notes": "Phase19 geometry with sixteen 1.25 um Stycast prism layers.",
    }

    steady = copy.deepcopy(project["cases"][SOURCE_STEADY])
    steady.update({
        "mesh": MESH,
        "state_file": f"work/meshes/{MESH}/stycast_z16_steady.state",
        "output_file_path": f"../work/meshes/{MESH}/{STEADY}.result",
        "series_file": f"{STEADY}_series.csv",
        "iteration_series_file": f"{STEADY}_iterations.csv",
    })
    project["cases"][STEADY] = steady

    common_early = [
        ["18[us]", 1], ["1[us]", 2], ["1[ns]", 1],
        ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9],
    ]
    coarse = pulse_case(
        project,
        PULSE_COARSE,
        common_early + [["10[us]", 9], ["5[us]", 1]],
    )
    coarse["comparison_time_grid"] = {
        "mode": "Stycast z-resolution diagnostic",
        "reference_case": SOURCE_PULSE,
        "only_changed_factor": "Stycast prism layers 1 -> 16",
        "end_after_pulse": "105[us]",
    }
    project["cases"][PULSE_COARSE] = coarse

    fine = pulse_case(project, PULSE_FINE, common_early + [["5[us]", 19]])
    fine["comparison_time_grid"] = {
        "mode": "Stycast z16 time-resolution diagnostic",
        "reference_case": PULSE_COARSE,
        "only_changed_factor": "post-pulse 10--100 us timestep 10 us -> 5 us",
        "end_after_pulse": "105[us]",
    }
    project["cases"][PULSE_FINE] = fine

    finer = pulse_case(project, PULSE_FINER, common_early + [["2.5[us]", 38]])
    finer["comparison_time_grid"] = {
        "mode": "Stycast z16 time-resolution diagnostic",
        "reference_case": PULSE_FINE,
        "only_changed_factor": "post-pulse 5 us timestep -> 2.5 us",
        "end_after_pulse": "105[us]",
    }
    project["cases"][PULSE_FINER] = finer

    finer_long = pulse_case(project, PULSE_FINER_LONG, common_early + [["2.5[us]", 86]])
    finer_long["comparison_time_grid"] = {
        "mode": "Stycast z16 10--90-percent fall-time diagnostic",
        "reference_case": PULSE_FINER,
        "only_changed_factor": "2.5 us tail extended from 105 us to 225 us",
        "end_after_pulse": "225[us]",
    }
    project["cases"][PULSE_FINER_LONG] = finer_long

    finest = pulse_case(project, PULSE_FINEST, common_early + [["1.25[us]", 76]])
    finest["comparison_time_grid"] = {
        "mode": "Stycast z16 time-convergence diagnostic",
        "reference_case": PULSE_FINER,
        "only_changed_factor": "post-pulse 2.5 us timestep -> 1.25 us",
        "end_after_pulse": "105[us]",
    }
    project["cases"][PULSE_FINEST] = finest

    finest_long = pulse_case(project, PULSE_FINEST_LONG, common_early + [["1.25[us]", 172]])
    finest_long["comparison_time_grid"] = {
        "mode": "Stycast z16 10--90-percent fall-time diagnostic",
        "reference_case": PULSE_FINEST,
        "only_changed_factor": "1.25 us tail extended from 105 us to 225 us",
        "end_after_pulse": "225[us]",
    }
    project["cases"][PULSE_FINEST_LONG] = finest_long

    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()

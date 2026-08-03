"""Prepare Phase26 with Stycast/absorber contact area matched to all-tet."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase26_stycast_area_tight.json"
MESH = "mesh_hybrid_abs_tet_layers_prism_phase26_sty493_abs35r50_noextend"
SOURCE_STEADY = "case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight"
SOURCE_PULSE = "case_p19_pulse_phase23_tight"
STEADY = "case_tes_steady_hybrid_prism_phase26_sty493_tight"
PULSE = "case_p19_pulse_phase26_sty493_tight"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    project["meshes"][MESH] = {
        "geometry": "single_pixel",
        "dir": MESH,
        "recipe": {
            "generator": "generate_hybrid_prism_geometry.py",
            "commands": [
                "python generate_hybrid_prism_geometry.py elmer_project_comsol_timegrid.json "
                "--output gmsh/project_hybrid_prism_phase26_sty493_abs35r50_noextend.msh "
                "--stack-local-size 16.6666666666667e-6 "
                "--absorber-local-size 35e-6 --absorber-local-radius 50e-6 "
                "--stycast-diameter 493.465731183667e-6 "
                "--disable-mesh-size-extend-from-boundary",
                "ElmerGrid 14 2 gmsh/project_hybrid_prism_phase26_sty493_abs35r50_noextend.msh "
                f"-merge 1e-10 -out {MESH}",
            ],
        },
        "notes": "Phase26: Stycast diameter sets the Stycast/absorber face area to the all-tet value.",
    }
    steady = copy.deepcopy(project["cases"][SOURCE_STEADY])
    steady.update({
        "mesh": MESH,
        "state_file": f"{MESH}/phase26_steady.state",
        "series_file": f"{STEADY}_series.csv",
        "iteration_series_file": f"{STEADY}_iterations.csv",
    })
    project["cases"][STEADY] = steady
    pulse = copy.deepcopy(project["cases"][SOURCE_PULSE])
    pulse.update({
        "mesh": MESH,
        "restart_from": STEADY,
        "restart_file_path": f"{STEADY}.result",
        "state_file": f"{MESH}/phase26_steady.state",
        "series_file": f"{PULSE}_series.csv",
        "iteration_series_file": f"{PULSE}_iterations.csv",
    })
    project["cases"][PULSE] = pulse
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

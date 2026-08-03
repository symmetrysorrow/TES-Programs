"""Prepare the Phase24 stack-resolution probe from the Phase23 tight project.

Only the central prism-stack/contact Box field changes: 16.7 um (Phase19) to
14.2857 um.  The 35 um / 50 um absorber-centre field, 5 us pulse grid and
strict nonlinear coupling are retained.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase24_stack14_tight.json"

MESH = "mesh_hybrid_abs_tet_layers_prism_stack14_abs35r50_noextend"
STEADY = "case_tes_steady_hybrid_prism_stack14_abs35r50_noextend_phase24_tight"
PULSE = "case_p19_pulse_phase24_stack14_tight"
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
                "python generate_hybrid_prism_geometry.py elmer_project_comsol_timegrid.json "
                "--output gmsh/project_hybrid_prism_stack14_abs35r50_noextend.msh "
                "--stack-local-size 14.2857142857143e-6 "
                "--absorber-local-size 35e-6 --absorber-local-radius 50e-6 "
                "--disable-mesh-size-extend-from-boundary",
                "ElmerGrid 14 2 gmsh/project_hybrid_prism_stack14_abs35r50_noextend.msh "
                f"-merge 1e-10 -out {MESH}",
            ],
        },
        "notes": (
            "Phase24: Phase19 absorber-local field fixed; central prism stack/contact "
            "Box field only is refined from 16.7 um to 14.2857 um."
        ),
    }

    steady = copy.deepcopy(project["cases"][SOURCE_STEADY])
    steady.update({
        "mesh": MESH,
        "state_file": f"{MESH}/phase24_steady.state",
        "series_file": f"{STEADY}_series.csv",
        "iteration_series_file": f"{STEADY}_iterations.csv",
    })
    project["cases"][STEADY] = steady

    pulse = copy.deepcopy(project["cases"][SOURCE_PULSE])
    pulse.update({
        "mesh": MESH,
        "restart_from": STEADY,
        # The solver resolves an unqualified restart name in the mesh DB
        # directory, where run.py stores the dependency result.
        "restart_file_path": f"{STEADY}.result",
        "state_file": f"{MESH}/phase24_steady.state",
        "series_file": f"{PULSE}_series.csv",
        "iteration_series_file": f"{PULSE}_iterations.csv",
    })
    project["cases"][PULSE] = pulse

    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    print(f"Cases: {STEADY}, {PULSE}")


if __name__ == "__main__":
    main()

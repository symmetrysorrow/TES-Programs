"""Create a short post-pulse transient case from the refined no-mortar restart."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "artifacts/comparison/nomortar_refinement_probe/steady_project.json"
OUTDIR = ROOT / "artifacts/comparison/nomortar_refinement_probe"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-name", default="case_tes_pulse_singlepixel_conformal_gpu_refine20_short")
    parser.add_argument("--vtu", action="store_true")
    parser.add_argument("--points", type=int, default=2, choices=(2, 10))
    args = parser.parse_args()
    project = json.loads(BASE.read_text(encoding="utf-8"))
    mesh_name = "mesh_singlepixel_conformal_gpu_refine20"
    steady_name = "case_tes_steady_singlepixel_conformal_gpu_refine20"
    mesh_dir = f"work/meshes/{mesh_name}"
    pulse = {
        "template": "pulse",
        "mesh": mesh_name,
        "restart_file_base": "external_refine20_steady",
        "restart_file_path": f"../{mesh_dir}/{steady_name}.result",
        "restart_time": 0.02,
        "series_file": f"{args.case_name}_series.csv",
        "iteration_series_file": f"{args.case_name}_iterations.csv",
        "initial_temperature": "T_0",
        "lumped_mass": True,
        "bdf_order": 1,
        "vtu": "after_timestep" if args.vtu else False,
        # The steady restart is already at 20 ms.  Aligning the pulse with
        # that restart removes the inert pre-event interval while preserving
        # the post-event response in elapsed time.
        "timesteps": [["1[ns]", 1], ["100[ns]", args.points - 1]],
        "output_intervals": [1, 1],
        "pulse": {
            "energy": "1332[keV]",
            "start": "20[ms]",
            "duration": "1[ns]",
            "sigma": "50[um]",
            "center": "auto",
        },
        "steady_state_max_iterations": 1,
        "solver": {
            "nonlinear_max_iterations": 25,
            "nonlinear_convergence_tolerance": 3e-7,
            "nonlinear_relaxation_factor": 1.0,
            "steady_state_convergence_tolerance": 1e-9,
            "linear_system": "mumps",
        },
        "heat_source": "circuit_inner",
        "state_file": f"{mesh_dir}/{steady_name}.state",
        "output_result": True,
        "restart_from": None,
        "output_file_path": f"../{mesh_dir}/{args.case_name}.result",
        "apply_mortar_bcs": False,
    }
    project["cases"] = {args.case_name: pulse}
    OUTDIR.mkdir(parents=True, exist_ok=True)
    output = OUTDIR / "transient_project.json"
    output.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"project={output}")
    print(f"case={args.case_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

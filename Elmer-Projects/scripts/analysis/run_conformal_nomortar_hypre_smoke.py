"""Prepare a short backend run on the conformal no-mortar mesh.

The case deliberately reuses the completed conformal shared-node steady
restart and changes only the linear backend.  It is a solver-isolation test:
it must not be interpreted as a production COMSOL qualification.
"""
from __future__ import annotations

import copy
import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
OUT_DIR = ROOT / "artifacts" / "conformal_nomortar_hypre_smoke"
PROJECT = OUT_DIR / "project.json"
SOURCE_CASE = "case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step"
STAGES = [
    ["18[us]", 1],
    ["1[us]", 2],
    ["1[ns]", 1],
    ["10[ns]", 10],
    ["100[ns]", 9],
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("hypre", "mumps"), default="hypre")
    parser.add_argument("--tolerance", type=float, default=6e-7)
    parser.add_argument("--max-iterations", type=int, default=4000)
    args = parser.parse_args()
    token = f"{args.tolerance:.0e}".replace("+", "")
    case_name = (
        "case_tes_pulse_singlepixel_conformal_gpu_fine_"
        f"{args.backend}_nomortar_"
        f"smoke_5stage_tol{token}"
    )
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    source = project["cases"][SOURCE_CASE]
    case = copy.deepcopy(source)
    case.update(
        {
            "restart_from": None,
            "preexisting_restart": True,
            "restart_file_base": "case_tes_steady_singlepixel_conformal_gpu_fine",
            "state_file": (
                "work/meshes/mesh_singlepixel_conformal_gpu_fine/"
                "case_tes_steady_singlepixel_conformal_gpu_fine.state"
            ),
            "timesteps": STAGES,
            "output_intervals": [999999, 999999, 999999, 999999, 1],
            "series_file": f"{case_name}_series.csv",
            "iteration_series_file": f"{case_name}_iterations.csv",
            "output_file_path": (
                "../work/meshes/mesh_singlepixel_conformal_gpu_fine/"
                f"{case_name}.result"
            ),
            "apply_mortar_bcs": False,
            "solver_comment": (
                f"Conformal shared-node no-mortar {args.backend} isolation; "
                f"0.9 us window; tolerance={args.tolerance:g}"
            ),
        }
    )
    case["solver"] = copy.deepcopy(case["solver"])
    case["solver"].update(
        {
            "linear_system": (
                "iterative_hypre_flexgmres_boomeramg"
                if args.backend == "hypre" else "mumps"
            ),
            "nonlinear_max_iterations": 120,
            "nonlinear_convergence_tolerance": 1e-8,
            "steady_state_convergence_tolerance": 1e-9,
        }
    )
    if args.backend == "hypre":
        case["solver"].update(
            {
                "linear_system_max_iterations": args.max_iterations,
                "linear_system_convergence_tolerance": args.tolerance,
            }
        )
    project["cases"][case_name] = case
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT.write_text(
        json.dumps(project, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(PROJECT)
    print(case_name)


if __name__ == "__main__":
    main()

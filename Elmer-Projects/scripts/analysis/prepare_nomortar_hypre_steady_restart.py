"""Prepare a HYPRE steady solve initialized from the converged MUMPS result.

This is an isolated solver check for the large no-mortar refine20 mesh.  The
restart result and TES circuit state are copied by the caller into distinct
files so the reference MUMPS artifacts remain untouched.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "artifacts/comparison/nomortar_refinement_probe/steady_project.json"
SOURCE = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
OUTDIR = ROOT / "artifacts/comparison/nomortar_refinement_probe"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-name",
        default="case_tes_steady_singlepixel_conformal_gpu_refine20_hypre_restart",
    )
    parser.add_argument(
        "--linear-system",
        choices=("iterative_hypre_boomeramg", "iterative_hypre_flexgmres_boomeramg"),
        default="iterative_hypre_flexgmres_boomeramg",
    )
    parser.add_argument("--linear-tolerance", type=float, default=1e-10)
    parser.add_argument("--linear-max-iterations", type=int, default=4000)
    parser.add_argument("--nonlinear-tolerance", type=float, default=1e-7)
    parser.add_argument("--nonlinear-max-iterations", type=int, default=10)
    args = parser.parse_args()

    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    # The refine20 mesh is an isolated generated artifact and is registered in
    # the refinement probe project, not in the checked-in source project.
    registry = OUTDIR / "project.json"
    if registry.is_file():
        project["meshes"] = json.loads(registry.read_text(encoding="utf-8"))["meshes"]
    source_name = "case_tes_steady_singlepixel_conformal_gpu_fine"
    source = copy.deepcopy(project["cases"].get(source_name))
    if source is None:
        raise SystemExit(f"missing source case {source_name!r} in {BASE}")

    mesh = "mesh_singlepixel_conformal_gpu_refine20"
    mumps_case = "case_tes_steady_singlepixel_conformal_gpu_refine20"
    state = f"work/meshes/{mesh}/{args.case_name}.state"
    result = f"../work/meshes/{mesh}/{args.case_name}.result"
    # Restart paths are resolved relative to the mesh directory by Elmer.
    restart = f"../work/meshes/{mesh}/{mumps_case}.result"
    solver = {
        **source["solver"],
        "linear_system": args.linear_system,
        "linear_system_convergence_tolerance": args.linear_tolerance,
        "linear_system_max_iterations": args.linear_max_iterations,
        "nonlinear_convergence_tolerance": args.nonlinear_tolerance,
        "nonlinear_max_iterations": args.nonlinear_max_iterations,
        "steady_state_convergence_tolerance": args.nonlinear_tolerance,
    }
    case = {
        **source,
        "mesh": mesh,
        "restart_from": None,
        "restart_file_base": mumps_case,
        "restart_file_path": restart,
        "restart_position": 0,
        "state_file": state,
        "series_file": f"{args.case_name}_series.csv",
        "iteration_series_file": f"{args.case_name}_iterations.csv",
        "output_file_path": result,
        "vtu": False,
        "steady_state_max_iterations": 1,
        "solver": solver,
        "solver_comment": (
            "HYPRE steady solve initialized from converged MUMPS refine20 "
            f"restart; backend={args.linear_system}; tol={args.linear_tolerance:g}"
        ),
    }
    project["cases"] = {args.case_name: case}
    OUTDIR.mkdir(parents=True, exist_ok=True)
    output = OUTDIR / "hypre_restart_steady_project.json"
    output.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"project={output}")
    print(f"case={args.case_name}")
    print(f"restart={restart}")
    print(
        "run_command="
        f"python run.py {args.case_name} --project {output} "
        '--elmer-solver "D:/Github/TES-Programs/tools/elmer-hypre/install-stage11/bin/ElmerSolver.exe" '
        '--runtime-bin "D:/Github/TES-Programs/tools/elmer-hypre/install-stage11/bin" '
        '--toolchain-bin "C:/msys64/ucrt64/bin"'
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

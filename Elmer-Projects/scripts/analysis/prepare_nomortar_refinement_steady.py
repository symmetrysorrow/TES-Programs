"""Create an isolated steady-state case for the refined no-mortar mesh.

The mesh is intentionally large, so this helper keeps the existing project
and result names untouched.  It is paired with ``run.py`` and is also useful
for producing a VTU-enabled variant later when interface field checks are
needed.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
OUTDIR = ROOT / "artifacts/comparison/nomortar_refinement_probe"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-name", default="mesh_singlepixel_conformal_gpu_refine20")
    parser.add_argument("--case-name", default="case_tes_steady_singlepixel_conformal_gpu_refine20")
    parser.add_argument("--vtu", action="store_true", help="request a final VTU field output")
    parser.add_argument("--steady-tol", type=float, default=1e-7)
    parser.add_argument("--nonlinear-tol", type=float, default=1e-7)
    parser.add_argument("--max-iterations", type=int, default=40)
    parser.add_argument(
        "--linear-system",
        choices=("mumps", "iterative_hypre_boomeramg", "iterative_hypre_flexgmres_boomeramg"),
        default="mumps",
    )
    parser.add_argument("--linear-tolerance", type=float, default=1e-10)
    parser.add_argument("--linear-max-iterations", type=int, default=4000)
    args = parser.parse_args()

    mesh_project = OUTDIR / "project.json"
    source_project = mesh_project if mesh_project.exists() else SOURCE
    project = json.loads(source_project.read_text(encoding="utf-8"))
    source_case = copy.deepcopy(project["cases"]["case_tes_steady_singlepixel_conformal_gpu_fine"])
    solver = {
        **source_case["solver"],
        "steady_state_convergence_tolerance": args.steady_tol,
        "nonlinear_convergence_tolerance": args.nonlinear_tol,
        "nonlinear_max_iterations": args.max_iterations,
        "linear_system": args.linear_system,
    }
    if args.linear_system != "mumps":
        solver.update(
            {
                "linear_system_convergence_tolerance": args.linear_tolerance,
                "linear_system_max_iterations": args.linear_max_iterations,
            }
        )
    source_case.update(
        {
            "mesh": args.mesh_name,
            "series_file": f"{args.case_name}_series.csv",
            "iteration_series_file": f"{args.case_name}_iterations.csv",
            "state_file": f"work/meshes/{args.mesh_name}/{args.case_name}.state",
            "output_file_path": f"../work/meshes/{args.mesh_name}/{args.case_name}.result",
            "vtu": "after_simulation" if args.vtu else False,
            "steady_state_max_iterations": 1,
            "solver": solver,
        }
    )
    project["cases"] = {args.case_name: source_case}
    OUTDIR.mkdir(parents=True, exist_ok=True)
    output = OUTDIR / "steady_project.json"
    output.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"project={output}")
    print(f"case={args.case_name}")
    print(
        "run_command="
        f"python run.py {args.case_name} --project {output} --skip-sync "
        '--elmer-solver "D:/Github/TES-Programs/tools/elmer-hypre/install-stage11/bin/ElmerSolver.exe" '
        '--runtime-bin "D:/Github/TES-Programs/tools/elmer-hypre/install-stage11/bin" '
        '--toolchain-bin "C:/msys64/ucrt64/bin"'
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

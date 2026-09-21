"""Prepare an isolated MUMPS steady solve on the validated production-v2 mesh.

Only the mortar flag is changed from the historical production-v2 steady case.
This separates the hybrid-prism mesh effect from the mortar effect without
touching the existing production result.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_prod_v2_original_timegrid.json"
OUTDIR = ROOT / "artifacts/comparison/prod_v2_nomortar_probe"
CASE = "case_tes_steady_prod_v2_nomortar"
MESH = "mesh_singlepixel_prod_v2"


def main() -> int:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_case = project["cases"]["case_tes_steady_singlepixel_prod_v2_original_timegrid_hybrid"]
    case = copy.deepcopy(source_case)
    case.update(
        {
            "mesh": MESH,
            "apply_mortar_bcs": False,
            "state_file": f"work/meshes/{MESH}/{CASE}.state",
            "series_file": f"{CASE}_series.csv",
            "iteration_series_file": f"{CASE}_iterations.csv",
            "output_file_path": f"../work/meshes/{MESH}/{CASE}.result",
            "solver_comment": "Production-v2 hybrid-prism mesh, no-mortar MUMPS isolation",
        }
    )
    case["solver"] = {
        **case["solver"],
        "linear_system": "mumps",
        "nonlinear_max_iterations": 120,
        "nonlinear_convergence_tolerance": 1e-8,
        "steady_state_convergence_tolerance": 1e-8,
    }
    project["cases"] = {CASE: case}
    OUTDIR.mkdir(parents=True, exist_ok=True)
    output = OUTDIR / "project.json"
    output.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"project={output}")
    print(f"case={CASE}")
    print(
        "run_command="
        f"python run.py {CASE} --project {output} "
        '--elmer-solver "D:/Github/TES-Programs/tools/elmer-hypre/install-stage11/bin/ElmerSolver.exe" '
        '--runtime-bin "D:/Github/TES-Programs/tools/elmer-hypre/install-stage11/bin" '
        '--toolchain-bin "C:/msys64/ucrt64/bin"'
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

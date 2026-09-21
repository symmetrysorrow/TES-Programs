"""Prepare a current native-HYPRE steady solve on the legacy production mesh."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_prod_v2_original_timegrid.json"
OUTDIR = ROOT / "artifacts/comparison/prod_v2_hypre_steady_probe"
CASE = "case_tes_steady_prod_v2_hypre"
MESH = "mesh_singlepixel_prod_v2"


def main() -> int:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    source = copy.deepcopy(
        project["cases"]["case_tes_steady_singlepixel_prod_v2_original_timegrid_hybrid"]
    )
    source.update(
        {
            "mesh": MESH,
            "apply_mortar_bcs": True,
            "state_file": f"work/meshes/{MESH}/{CASE}.state",
            "series_file": f"{CASE}_series.csv",
            "iteration_series_file": f"{CASE}_iterations.csv",
            "output_file_path": f"../work/meshes/{MESH}/{CASE}.result",
            "solver_comment": "Legacy production-v2 hybrid-prism + mortar native HYPRE steady probe",
        }
    )
    source["solver"] = {
        **source["solver"],
        "linear_system": "iterative_hypre_flexgmres_boomeramg",
        "linear_system_max_iterations": 4000,
        "linear_system_convergence_tolerance": 1e-10,
        "nonlinear_max_iterations": 120,
        "nonlinear_convergence_tolerance": 1e-8,
        "steady_state_convergence_tolerance": 1e-8,
    }
    project["cases"] = {CASE: source}
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

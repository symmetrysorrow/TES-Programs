"""Prepare the validated COMSOL-qualification profile.

The current Phase24 hybrid mesh/HYPRE case is a performance experiment, not a
COMSOL reference: its pre-pulse operating point is different.  This profile
reuses the already validated production-v2 spatial mesh and direct MUMPS
backend, preserving the fine rise-time grid and 100-us tail.  HYPRE cases are
left untouched and can still be used for speed/memory experiments.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_prod_v2_original_timegrid.json"
OUTPUT = ROOT / "elmer_project_phase24_comsol_qualification.json"
SOURCE_CASE = "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid"
CASE = "case_phase24_comsol_qualification_100us"


def main() -> int:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    candidate = copy.deepcopy(project["cases"][SOURCE_CASE])
    candidate.update({
        "series_file": f"{CASE}_series.csv",
        "iteration_series_file": f"{CASE}_iterations.csv",
        "output_file_path": f"../work/meshes/mesh_singlepixel_prod_v2/{CASE}.result",
        "solver_comment": (
            "COMSOL qualification profile: production-v2 mesh, direct MUMPS, "
            "validated fine rise grid and 100-us tail"
        ),
        "comparison_time_grid": {
            "mode": "COMSOL qualification",
            "reference": "docs/Single-Pixel.txt",
            "purpose": "absolute baseline plus baseline-corrected 0..100 us waveform",
        },
        "phase24_smoke": {
            "purpose": "COMSOL qualification, not HYPRE performance",
            "backend": "Direct MUMPS",
            "no_matrix_dump": True,
            "no_vtu": True,
        },
    })
    candidate["solver"] = copy.deepcopy(candidate["solver"])
    candidate["solver"]["linear_system"] = "mumps"
    candidate["phase24_vector_assembly"] = False
    candidate["phase24_static_matrix_reuse"] = False
    candidate["phase24_operator_lagging"] = "disabled"
    candidate["phase24_preconditioner_lagging"] = "disabled"
    project["cases"] = {CASE: candidate}
    OUTPUT.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"prepared {OUTPUT}")
    print(f"case={CASE}")
    print("backend=direct MUMPS")
    print("mesh=mesh_singlepixel_prod_v2")
    print("qualification=existing result: artifacts/comparison/prod_v2_mumps_vs_comsol")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

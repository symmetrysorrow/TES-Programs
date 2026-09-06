"""Create the small Phase24 production-path before/after smoke project."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hypre_gpu_phase19.json"
OUTPUT = ROOT / "elmer_project_phase24_production.json"
BASE_CASE = "case_p19_hypre_flexgmres_boomeramg_cpu_time5us_smoke_7step"
CASE = "case_phase24_hypre_cpu_smoke_7step"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    base = project["cases"][BASE_CASE]
    candidate = copy.deepcopy(base)
    candidate["restart_from"] = None
    candidate["preexisting_restart"] = True
    candidate["series_file"] = f"{CASE}_series.csv"
    candidate["iteration_series_file"] = f"{CASE}_iterations.csv"
    candidate["output_file_path"] = (
        f"../work/meshes/{candidate['mesh']}/{CASE}.result"
    )
    candidate["phase24_smoke"] = {
        "purpose": "same production HYPRE CPU smoke with Phase24 opt-in",
        "baseline_case": BASE_CASE,
        "path": "scalar 3D P1 tetra, static/dynamic split, cached CSR insertion",
    }
    candidate["phase24_vector_assembly"] = True
    candidate["phase24_wall_profiling"] = True
    candidate["solver"] = dict(candidate["solver"])
    candidate["solver"]["matrix_dump_prefix"] = CASE
    project["cases"][CASE] = candidate
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

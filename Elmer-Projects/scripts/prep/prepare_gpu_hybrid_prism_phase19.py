"""Create a GPU/AMGX copy of the Phase19 hybrid-prism 5-us benchmark."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase19_time5us.json"
OUTPUT = ROOT / "elmer_project_gpu_hybrid_prism_phase19.json"
SOURCE_CASE = "case_p19_pulse_time5us"
GPU_CASE = "case_p19_gpu_amgx_phase19_time5us"
SMOKE_CASE = f"{GPU_CASE}_smoke_7step"


def truncate(groups: list[list[object]], n: int) -> list[list[object]]:
    out: list[list[object]] = []
    left = n
    for token, count in groups:
        take = min(left, int(count))
        if take:
            out.append([token, take])
            left -= take
        if not left:
            return out
    raise ValueError("source time grid is shorter than requested smoke run")


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    source = copy.deepcopy(project["cases"][SOURCE_CASE])
    source.update(
        {
            # Reuse the already validated Phase19 steady result.  This avoids
            # rerunning a steady solve with the WSL build just to benchmark
            # the transient AMGX path.
            "restart_from": None,
            "preexisting_restart": True,
            "restart_file_base": "case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight",
            "restart_file_path": "../work/meshes/mesh_hybrid_abs_tet_layers_prism_stack17_abs35r50_noextend/case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight.result",
            "state_file": "work/meshes/mesh_hybrid_abs_tet_layers_prism_stack17_abs35r50_noextend/phase19_steady.state",
            "series_file": f"{GPU_CASE}_series.csv",
            "iteration_series_file": f"{GPU_CASE}_iterations.csv",
            "output_file_path": f"../work/meshes/{source['mesh']}/{GPU_CASE}.result",
            "comparison_time_grid": {
                "mode": "Phase19 hybrid-prism GPU AMGX benchmark",
                "purpose": "same mesh, restart, source, circuit, and 5-us rise grid as CPU Phase19",
                "reference_case": SOURCE_CASE,
                "linear_solver": "AMGX on RTX 3060 Ti",
            },
        }
    )
    project["cases"][GPU_CASE] = source
    smoke = copy.deepcopy(source)
    smoke["timesteps"] = truncate(source["timesteps"], 7)
    smoke["output_intervals"] = [999999] * len(smoke["timesteps"])
    smoke["output_intervals"][-1] = 1
    smoke["series_file"] = f"{SMOKE_CASE}_series.csv"
    smoke["iteration_series_file"] = f"{SMOKE_CASE}_iterations.csv"
    smoke["output_file_path"] = f"../work/meshes/{source['mesh']}/{SMOKE_CASE}.result"
    project["cases"][SMOKE_CASE] = smoke
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()

"""Create the Phase22 RMS-equivalent uniform-sphere source probe.

The Phase19 mesh, 5-us time grid, restart, circuit and total deposited energy
are inherited unchanged.  Only the spatial source profile changes from the
50-us Gaussian to a uniform sphere with equal <r^2> (sqrt(5)*sigma).
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase19_time5us.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase22_heatshape.json"
BASE_CASE = "case_p19_pulse_time5us"
CASE = "case_p19_pulse_rms_sphere_time5us"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = copy.deepcopy(project["cases"][BASE_CASE])
    case["series_file"] = f"{CASE}_series.csv"
    case["iteration_series_file"] = f"{CASE}_iterations.csv"
    case["pulse"].update(
        {
            "shape": "uniform_sphere",
            "radius": "111.8033988749895[um]",
        }
    )
    case["comparison_time_grid"] = {
        "mode": "Phase22 source-shape sensitivity",
        "reference_case": BASE_CASE,
        "only_changed_factor": "Gaussian sigma=50 um -> RMS-equivalent uniform sphere radius=sqrt(5)*50 um",
        "time_step": "5[us]",
    }
    project["cases"][CASE] = case
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

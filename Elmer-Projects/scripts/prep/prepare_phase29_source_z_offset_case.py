"""Create a +Z 20-um heat-source offset probe from the Phase28 C2 case."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase28_c2_pulse.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase29_source_z20um.json"
BASE_CASE = "case_p19_pulse_phase28_c2_01ps"
CASE = "case_p19_pulse_phase29_c2_zplus20um"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = copy.deepcopy(project["cases"][BASE_CASE])
    case["series_file"] = f"{CASE}_series.csv"
    case["iteration_series_file"] = f"{CASE}_iterations.csv"
    # The mesh-specific absorber centroid is z=562.16 um.  Retain automatic
    # x/y centering and move only the Gaussian centre along +Z by 20 um.
    case["pulse"]["center"] = {"x": "auto", "y": "auto", "z": "582.16[um]"}
    case["comparison_time_grid"] = {
        **case["comparison_time_grid"],
        "mode": "Phase29 source-position sensitivity",
        "reference_case": BASE_CASE,
        "only_changed_factor": "Gaussian source centre z: absorber centroid + 20[um]",
    }
    project["cases"][CASE] = case
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

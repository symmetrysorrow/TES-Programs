"""Create a Phase29 rerun with sparse transient VTU temperature output."""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase29_source_z20um.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase30_abs_temperature.json"
BASE = "case_p19_pulse_phase29_c2_zplus20um"
CASE = "case_p19_pulse_phase30_abs_temperature"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = copy.deepcopy(project["cases"][BASE])
    case["series_file"] = f"{CASE}_series.csv"
    case["iteration_series_file"] = f"{CASE}_iterations.csv"
    # Keep every edge and early transient output; save every 10th 5-us step
    # thereafter, which is sufficient for the ~345-us peak trajectory.
    case["vtu"] = {"exec_intervals": [1, 1, 1, 1, 1, 1, 1, 1, 1, 10]}
    case["comparison_time_grid"] = {
        **case["comparison_time_grid"],
        "mode": "Phase30 absorber temperature extraction",
        "reference_case": BASE,
        "result_output": "VTU after selected timesteps; absorber volume average",
    }
    project["cases"][CASE] = case
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

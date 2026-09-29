"""Split the 1-ns heat-deposition plateau into 10-ps steps through 1 us."""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase31_abs_temperature_fine1us.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase32_abs_temperature_splitpulse.json"
BASE = "case_p19_pulse_phase31_abs_temperature_fine1us"
CASE = "case_p19_pulse_phase32_abs_temperature_splitpulse"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = copy.deepcopy(project["cases"][BASE])
    case["series_file"] = f"{CASE}_series.csv"
    case["iteration_series_file"] = f"{CASE}_iterations.csv"
    # The 999-ps plateau between the C2 edges is now 10 ps x 100.  The
    # post-pulse 1-us interval remains 1 ns x 100.
    case["timesteps"] = [
        ["18[us]", 1], ["1999.9995[ns]", 1], ["0.1[ps]", 10],
        ["10[ps]", 100], ["0.1[ps]", 10], ["1[ns]", 100],
    ]
    case["output_intervals"] = [1] * len(case["timesteps"])
    case["vtu"] = {"exec_intervals": [1, 1, 1, 10, 1, 10]}
    case["comparison_time_grid"] = {
        **case["comparison_time_grid"],
        "mode": "Phase32 heat-deposition plateau refinement",
        "reference_case": BASE,
        "pulse_plateau_step": "10[ps]",
    }
    project["cases"][CASE] = case
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

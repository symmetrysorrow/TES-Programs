"""Create a fine-step absorber-temperature probe through 1 us after the pulse."""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase30_abs_temperature.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase31_abs_temperature_fine1us.json"
BASE = "case_p19_pulse_phase30_abs_temperature"
CASE = "case_p19_pulse_phase31_abs_temperature_fine1us"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = copy.deepcopy(project["cases"][BASE])
    case["series_file"] = f"{CASE}_series.csv"
    case["iteration_series_file"] = f"{CASE}_iterations.csv"
    # Resolve both 1-ps source edges as before, then use 1-ns steps for the
    # first microsecond after the nominal pulse.  The final time is t0+1 us.
    case["timesteps"] = [
        ["18[us]", 1], ["1999.9995[ns]", 1], ["0.1[ps]", 10],
        ["999[ps]", 1], ["0.1[ps]", 10], ["1[ns]", 100],
    ]
    case["output_intervals"] = [1] * len(case["timesteps"])
    case["vtu"] = {"exec_intervals": [1, 1, 1, 1, 1, 1]}
    case["comparison_time_grid"] = {
        **case["comparison_time_grid"],
        "mode": "Phase31 absorber-temperature time-step refinement",
        "reference_case": BASE,
        "fine_step": "1[ns] through pulse + 1[us]",
    }
    project["cases"][CASE] = case
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

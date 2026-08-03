"""Create the Phase28 COMSOL-style C2 pulse-edge comparison case.

The Phase19 hybrid mesh and Phase23 nonlinear tolerance are retained.  Only
the temporal source profile and its immediate time grid are changed: COMSOL's
1-ps transition zone with two continuous derivatives is sampled at 0.1 ps.
The run stops 500 us after the nominal pulse end, comfortably after the
approximately 345-us current peak observed in Phase19.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase28_c2_pulse.json"
BASE_CASE = "case_p19_pulse_phase23_tight"
CASE = "case_p19_pulse_phase28_c2_01ps"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = copy.deepcopy(project["cases"][BASE_CASE])
    case["series_file"] = f"{CASE}_series.csv"
    case["iteration_series_file"] = f"{CASE}_iterations.csv"
    case["pulse"]["transition_zone"] = "1[ps]"
    # Starting from the existing 20-ms restart: land at the leading-edge
    # transition boundary, resolve both 1-ps zones in ten steps, then stop
    # 500 us after the nominal 1-ns pulse ends.
    case["timesteps"] = [
        ["18[us]", 1],
        ["1999.9995[ns]", 1],
        ["0.1[ps]", 10],
        ["999[ps]", 1],
        ["0.1[ps]", 10],
        ["10[ns]", 10],
        ["100[ns]", 9],
        ["1[us]", 9],
        ["10[us]", 9],
        ["5[us]", 80],
    ]
    case["output_intervals"] = [1] * len(case["timesteps"])
    case["comparison_time_grid"] = {
        "mode": "Phase28 COMSOL C2 temporal-pulse comparison",
        "reference_case": BASE_CASE,
        "only_changed_factor": "rectangular temporal pulse -> 1 ps COMSOL flc2hs C2 transition zone",
        "edge_time_step": "0.1[ps]",
        "end_time": "pulse end + 500[us]",
        "peak_coverage": "Phase19 peak was approximately pulse + 345[us]",
    }
    project["cases"][CASE] = case
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

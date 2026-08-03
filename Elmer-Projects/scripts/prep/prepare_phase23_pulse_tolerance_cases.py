"""Create paired tight nonlinear-coupling probes on the Phase19 5-us grid."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase19_time5us.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"


def make_case(project: dict, source: str, target: str) -> None:
    case = copy.deepcopy(project["cases"][source])
    case["series_file"] = f"{target}_series.csv"
    case["iteration_series_file"] = f"{target}_iterations.csv"
    case["solver"]["nonlinear_convergence_tolerance"] = 1e-8
    case["solver"]["nonlinear_max_iterations"] = 120
    case["comparison_time_grid"] = {
        "mode": "Phase23 nonlinear coupling tolerance probe",
        "reference_case": source,
        "only_changed_factor": "nonlinear_convergence_tolerance 3e-7 -> 1e-8; max_iterations 25 -> 120",
        "time_step": "5[us]",
    }
    project["cases"][target] = case


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    make_case(
        project,
        "case_tes_pulse_3x_phase19_time5us",
        "case_tes_pulse_3x_phase23_tight",
    )
    make_case(project, "case_p19_pulse_time5us", "case_p19_pulse_phase23_tight")
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

"""Prepare a Phase25 pulse case with only the Stycast/absorber mortar reversed."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"
OUTPUT = ROOT / "elmer_project_hybrid_prism_phase25_mortar_orientation.json"
SOURCE_CASE = "case_p19_pulse_phase23_tight"
CASE = "case_p19_pulse_phase25_abs_slave_tight"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = copy.deepcopy(project["cases"][SOURCE_CASE])
    case.update({
        "reverse_stycast_abs_mortar": True,
        "series_file": f"{CASE}_series.csv",
        "iteration_series_file": f"{CASE}_iterations.csv",
    })
    project["cases"][CASE] = case
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()

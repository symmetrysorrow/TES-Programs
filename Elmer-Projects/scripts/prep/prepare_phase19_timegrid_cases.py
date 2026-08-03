"""Add paired post-pulse time-grid probes to the Phase19 project.

The probes leave mesh, restart state, source position, and circuit settings
unchanged.  They replace only the final 6 x 100-us stage (20.120001--
20.720001 ms in the legacy definition) with a uniform probe grid, ending at the
same 20.620001-ms comparison boundary used by the Phase19 report.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / "elmer_project_hybrid_prism.json"

REFERENCE = "case_tes_pulse_3x_from_20ms_phase13_tight_1rank_step_commit"
CANDIDATE = "case_p19_pulse"


def make_probe(project: dict, source_name: str, probe_name: str, step_us: int) -> dict:
    case = copy.deepcopy(project["cases"][source_name])
    # The first seven stages resolve the 1-ns heat input and the initial
    # 100 us response.  The final stage is the only 100-us sampling interval.
    assert case["timesteps"][-1] == ["100[us]", 6]
    assert len(case["timesteps"]) == len(case["output_intervals"])
    case["timesteps"][-1] = [f"{step_us}[us]", 500 // step_us]
    case["output_intervals"][-1] = 1
    case["series_file"] = f"{probe_name}_series.csv"
    case["iteration_series_file"] = f"{probe_name}_iterations.csv"
    # Elmer prepends the Mesh DB directory to this path.  Its mesh directory
    # name is relative to the repository root, so step out once before
    # returning to work/meshes.
    case["restart_file_path"] = (
        f"../work/meshes/{case['mesh']}/{case['restart_from']}.result"
    )
    case["comparison_time_grid"] = {
        "mode": f"paired {step_us}-us post-pulse time-grid probe",
        "purpose": "separate 100-us peak-sampling/timestep effects from mesh effects",
        "same_mesh_restart_pulse_and_circuit_as": source_name,
        "post_pulse_step": f"{step_us}[us]",
        "comparison_end": "20.620001[ms]",
    }
    project["cases"][probe_name] = case
    return case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--post-pulse-step-us", type=int, default=10)
    args = parser.parse_args()
    if args.post_pulse_step_us <= 0 or 500 % args.post_pulse_step_us:
        raise ValueError("--post-pulse-step-us must be a positive divisor of 500")
    step_us = args.post_pulse_step_us
    output = args.output or ROOT / f"elmer_project_hybrid_prism_phase19_time{step_us}us.json"

    project = json.loads(args.source.read_text(encoding="utf-8"))
    for name in (REFERENCE, CANDIDATE):
        if name not in project["cases"]:
            raise KeyError(f"missing source case: {name}")
    make_probe(project, REFERENCE, f"case_tes_pulse_3x_phase19_time{step_us}us", step_us)
    make_probe(project, CANDIDATE, f"case_p19_pulse_time{step_us}us", step_us)
    output.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

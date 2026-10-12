#!/usr/bin/env python3
"""Compare the time-integration and solver controls in the two run SIFs."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
single_path = ROOT.parent / "wsl-8min-replay-20261013" / "replay.sif"
dual_path = ROOT / "inputs" / "dualtes_wsl_transient_20261013_newton.sif"
keys = [
    "Timestepping Method", "BDF Order", "Phase24 BDF2 Predictor",
    "Phase24 BDF1 Step Start Times(2)", "Timestep Sizes(26)",
    "Timestep Intervals(26)", "Output Intervals(26)",
    "Nonlinear System Max Iterations", "Nonlinear System Convergence Tolerance",
    "Nonlinear System Relaxation Factor", "Linear System Use Hypre",
    "Linear System Solver", "Linear System Iterative Method",
    "Linear System Preconditioning", "Linear System Max Iterations",
    "Linear System Convergence Tolerance", "BoomerAMG Relax Type",
    "Phase24 Prism Lumped Preconditioner", "Phase24 Precond Matrix Reuse",
    "Phase24 Coupled Temperature Tolerance", "Phase24 Vector Assembly",
    "TES Inner Circuit Newton", "TES Inner Circuit Predict From Temperature",
    "TES Inner Circuit Step Commit", "Apply Mortar BCs", "Restart Time",
    "Steady State Convergence Tolerance",
]


def read_sif(path):
    text = path.read_text(encoding="utf-8")
    result = {}
    for key in keys:
        match = re.search(r"^\s*\"?" + re.escape(key) + r"\"?\s*=\s*(.*?)\s*$", text, flags=re.M | re.I)
        result[key] = match.group(1) if match else None
    cache = re.search(r"^\s*Phase24 Dynamic K Node Cache\s*=\s*(.*?)\s*$", text, flags=re.M | re.I)
    result["Phase24 Dynamic K Node Cache"] = cache.group(1) if cache else "not specified; solver default OFF"
    return result


single = read_sif(single_path)
dual = read_sif(dual_path)
same = {key: single[key] == dual[key] for key in keys}
for key in ("Nonlinear System Convergence Tolerance", "Nonlinear System Relaxation Factor", "Linear System Convergence Tolerance", "Phase24 Coupled Temperature Tolerance", "Restart Time", "Steady State Convergence Tolerance"):
    try:
        same[key] = float(single[key].split()[-1]) == float(dual[key].split()[-1])
    except (AttributeError, ValueError, IndexError):
        pass
active_keys = [key for key in keys if key != "Steady State Convergence Tolerance"]
report = {
    "singlepixel_sif": str(single_path.relative_to(ROOT.parent.parent)),
    "singlepixel_sif_sha256": __import__("hashlib").sha256(single_path.read_bytes()).hexdigest(),
    "dualtes_sif": str(dual_path.relative_to(ROOT.parent.parent)),
    "dualtes_sif_sha256": __import__("hashlib").sha256(dual_path.read_bytes()).hexdigest(),
    "singlepixel": single,
    "dualtes": dual,
    "exactly_equal_solver_controls": same,
    "all_listed_controls_match": all(same.values()),
    "all_active_transient_controls_match": all(same[key] for key in active_keys),
    "transient_inactive_steady_tolerance": {
        "singlepixel": single["Steady State Convergence Tolerance"],
        "dualtes": dual["Steady State Convergence Tolerance"],
        "applies_to_transient_steps": False,
    },
    "model_specific_exceptions": ["Mesh, body IDs, two TES circuits, circuit constants, pulse/source, and input-specific restart result differ by model."],
}
(ROOT / "solver_settings_comparison.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"all_active_transient_controls_match": report["all_active_transient_controls_match"], "all_listed_controls_match": report["all_listed_controls_match"], "differences": [k for k, v in same.items() if not v]}, indent=2))

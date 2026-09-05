"""Audit the reduced TES operating point, Jacobian, and stability gate.

This diagnostic does not fit any noise parameter and does not add an
experimental residual source.  It compares the analytic linearization used by
``MakeNoise`` with a central finite-difference Jacobian of its canonical
nonlinear RHS.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))
import PoST_Simulation as simulation  # noqa: E402


def run(parameters: dict):
    point = simulation.tes_operating_point(parameters)
    analytic = simulation.tes_linearized_time_matrix(parameters)
    numeric = simulation.numerical_jacobian(
        lambda state: simulation.tes_nonlinear_rhs(state, parameters),
        point["state"],
    )
    difference = analytic - numeric
    scale = np.maximum(np.abs(numeric), np.finfo(float).tiny)
    absolute_error = np.abs(difference)
    relative_error = absolute_error / scale
    stability = simulation.diagnose_linear_stability(analytic)
    return {
        "operating_current_A": point["current_A"],
        "operating_temperature_K": point["tes_temperature_K"],
        "bath_power_W": point["bath_power_W"],
        "matrix_convention": "M(omega) = -A + i*omega*I; dx/dt = A*x",
        "jacobian_max_absolute_error_per_s": float(np.max(np.abs(difference))),
        "jacobian_max_relative_error": float(np.max(np.abs(difference) / scale)),
        "jacobian_element_errors": [
            [
                {
                    "absolute_per_s": float(absolute_error[row, column]),
                    "relative": float(relative_error[row, column]),
                }
                for column in range(analytic.shape[1])
            ]
            for row in range(analytic.shape[0])
        ],
        "eigenvalues_per_s": [
            {"real": float(value.real), "imag": float(value.imag)}
            for value in stability["eigenvalues_per_s"]
        ],
        "max_real_part_per_s": stability["max_real_part_per_s"],
        "unstable_mode": stability["unstable_mode"],
        "pole_frequency_scale_hz": [
            float(value) for value in stability["pole_frequency_scale_hz"]
        ],
        "noise_sources": {
            "tes_johnson": "ASD sqrt(4*k_B*T*R*(1+2*beta)*(1+M^2)); injected as one correlated voltage source",
            "load_johnson": "ASD sqrt(4*k_B*T_bath*R_l); independent voltage source",
            "thermal_links": "ASD sqrt(4*k_B*T_c^2*G*F); shared link source enters connected nodes with opposite signs",
            "aggregation": "independent columns add in PSD; CH0/CH1 cross PSD retains shared-source complex transfer",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=SIMULATION_ROOT / "input.json")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = run(json.loads(args.input.read_text(encoding="utf-8")))
    print(json.dumps(result, indent=2, allow_nan=False))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

"""Compute the corrected Gate3 one-shot thermal image from a native capture."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.support.analyze_phase24_saddle_sensitivity import (
    indexed_vector,
    load_matrix,
    load_result_field,
    tes_weights,
)


def circuit_at(temperature: float) -> dict[str, float]:
    ibias, rsh, r0, rmin = 715e-6, 3.9e-3, 15.527e-3, 1e-6
    alpha, beta, i0, t0 = 256.46, 5.03, 143.537344932311e-6, 168.57e-3
    a = r0 * (1.0 + alpha * (temperature - t0) / t0 - beta)
    b = r0 * beta / i0
    current = (np.sqrt((rsh + a) ** 2 + 4.0 * b * ibias * rsh) - (rsh + a)) / (2.0 * b)
    current = float(np.clip(current, 0.0, ibias))
    resistance = max(a + b * abs(current), rmin)
    if resistance == rmin:
        current = ibias * rsh / (rsh + resistance)
    return {"temperature_K": temperature, "raw_current_A": current,
            "resistance_ohm": resistance, "raw_power_W": current * current * resistance}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--restart-result", type=Path, required=True)
    parser.add_argument("--mesh-elements", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metadata = json.loads((args.capture / "metadata.json").read_text(encoding="utf-8"))
    total = int(metadata["runtime"]["total_rows"])
    primal = int(metadata["runtime"]["primal_rows"])
    matrix = load_matrix(args.capture / "full_A_before.dat", (total, total)).tocsc()
    rhs = indexed_vector(args.capture / "full_b_before.dat", total)
    restart = indexed_vector(args.capture / "full_x_before.dat", total)
    field, permutation = load_result_field(args.restart_result)
    weights = tes_weights(args.mesh_elements, permutation, 101)
    direct = np.asarray(splu(matrix).solve(rhs), dtype=np.float64)
    restart_t = float(weights.dot(restart[:primal]))
    direct_t = float(weights.dot(direct[:primal]))
    residual = np.asarray(matrix @ direct - rhs).ravel()
    circuit = circuit_at(direct_t)
    saved_t = restart_t
    saved_circuit = circuit_at(saved_t)
    result = {
        "method": "scipy sparse LU of corrected full_A_before; no nonlinear feedback",
        "restart_tes_temperature_K": saved_t,
        "thermal_only_tes_temperature_K": direct_t,
        "delta_T_mK": (direct_t - saved_t) * 1.0e3,
        "saved_circuit": saved_circuit,
        "field_re_evaluated_raw_circuit": circuit,
        "delta_I_uA": (circuit["raw_current_A"] - saved_circuit["raw_current_A"]) * 1.0e6,
        "full_residual": {"l2": float(np.linalg.norm(residual)), "max_abs": float(np.max(np.abs(residual)))},
        "dimensions": {"total_rows": total, "primal_rows": primal, "constraint_rows": total - primal},
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"delta_T_mK": result["delta_T_mK"], "delta_I_uA": result["delta_I_uA"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

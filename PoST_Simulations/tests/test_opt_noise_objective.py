from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))

import Opt_noise as optimizer  # noqa: E402


def test_absolute_level_residual_is_log10_ratio() -> None:
    residual = optimizer.absolute_level_residual(2.0e-12, 1.0e-12)
    assert residual == pytest.approx(np.log10(2.0))


def test_absorber_heat_capacity_prior_is_factor_two_about_material_value() -> None:
    nominal = optimizer.C_ABS_MATERIAL_J_PER_K
    assert optimizer.C_ABS_FIT_MIN_J_PER_K == pytest.approx(0.5 * nominal)
    assert optimizer.C_ABS_FIT_MAX_J_PER_K == pytest.approx(2.0 * nominal)


def test_boundary_diagnostics_uses_log_position_for_log_bound() -> None:
    bound = optimizer.Bound(1.0e-10, 1.0e-8, logarithmic=True)
    candidate = {"L": 1.0e-9}
    result = optimizer.parameter_boundary_diagnostics(
        candidate,
        ["L"],
        {"L": bound},
    )["L"]
    assert result["position_fraction"] == pytest.approx(0.5)
    assert result["nearest_bound"] == "lower"
    assert result["within_1pct_of_bound"] is False

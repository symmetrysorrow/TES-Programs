from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_c2_physical_replacement_diagnostic as diag


def test_magnicon_bessel_magnitude_is_unity_at_dc():
    value = diag.magnicon_bessel_magnitude(
        np.array([0.0]),
        10_000.0,
        1.0 / np.sqrt(3.0),
    )
    assert value[0] == pytest.approx(1.0)


def test_lead_lag_magnitude_limits():
    zero = 2_000.0
    pole = 200_000.0
    low = diag.lead_lag_magnitude(
        np.array([0.0]),
        zero,
        pole,
    )[0]
    high = diag.lead_lag_magnitude(
        np.array([1.0e12]),
        zero,
        pole,
    )[0]
    assert low == pytest.approx(1.0)
    assert high == pytest.approx(pole / zero, rel=1e-6)


def test_lead_lag_rejects_nonlead_ordering():
    with pytest.raises(ValueError):
        diag.lead_lag_magnitude(
            np.array([1_000.0]),
            10_000.0,
            5_000.0,
        )


def test_equivalent_c2_zero_uses_reference_scale():
    zero = diag.equivalent_c2_zero_hz(100.0, 40_000.0)
    assert zero == pytest.approx(4_000.0)
    assert diag.equivalent_c2_zero_hz(0.0, 40_000.0) is None


def _fit(score, rms):
    return {
        "joint_shape_score": float(score),
        "joint_continuum_rms_dB": float(rms),
    }


def test_replacement_summary_accepts_near_c2_candidate():
    screen = {
        "max_score_ratio_to_c2_reference": 1.5,
        "max_rms_delta_from_c2_reference_dB": 0.15,
        "max_score_ratio_to_c2_zero_for_material_improvement": 0.8,
        "min_rms_improvement_vs_c2_zero_dB": 0.05,
    }
    result = diag.replacement_summary(
        _fit(0.011, 0.55),
        _fit(0.08, 2.9),
        _fit(0.010, 0.50),
        screen,
    )
    assert result["material_improvement_vs_c2_zero"]
    assert result["within_c2_reference_screen"]
    assert result["classification"] == (
        "candidate_replaces_c2_within_screen"
    )


def test_replacement_summary_distinguishes_partial_improvement():
    screen = {
        "max_score_ratio_to_c2_reference": 1.5,
        "max_rms_delta_from_c2_reference_dB": 0.15,
        "max_score_ratio_to_c2_zero_for_material_improvement": 0.8,
        "min_rms_improvement_vs_c2_zero_dB": 0.05,
    }
    result = diag.replacement_summary(
        _fit(0.03, 1.0),
        _fit(0.08, 2.9),
        _fit(0.01, 0.50),
        screen,
    )
    assert result["material_improvement_vs_c2_zero"]
    assert not result["within_c2_reference_screen"]
    assert result["classification"] == (
        "candidate_materially_improves_c2_zero_but_not_c2_reference"
    )

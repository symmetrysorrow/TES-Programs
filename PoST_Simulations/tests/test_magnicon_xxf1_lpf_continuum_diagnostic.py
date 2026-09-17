from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import Opt_noise as opt
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as diag


def test_second_order_bessel_canonical_matches_analog_response_mag():
    frequency = np.geomspace(10.0, 200_000.0, 300)
    coordinates = diag.second_order_bessel_canonical(10_000.0, "mag")

    canonical = diag.effective.effective_transfer_magnitude(
        frequency,
        coordinates["pole_Hz"],
        coordinates["pole_Q"],
        0,
        {},
        40_000.0,
    )
    direct = opt.hardware_filter_magnitude(
        frequency,
        cutoff_hz=10_000.0,
        order=2,
        norm="mag",
    )

    np.testing.assert_allclose(canonical, direct, rtol=1e-11, atol=1e-13)
    assert coordinates["pole_Q"] == pytest.approx(1.0 / np.sqrt(3.0))
    assert coordinates["dc_gain"] == pytest.approx(1.0)


def test_second_order_bessel_canonical_matches_analog_response_phase():
    frequency = np.geomspace(10.0, 200_000.0, 300)
    coordinates = diag.second_order_bessel_canonical(10_000.0, "phase")

    canonical = diag.effective.effective_transfer_magnitude(
        frequency,
        coordinates["pole_Hz"],
        coordinates["pole_Q"],
        0,
        {},
        40_000.0,
    )
    direct = opt.hardware_filter_magnitude(
        frequency,
        cutoff_hz=10_000.0,
        order=2,
        norm="phase",
    )

    np.testing.assert_allclose(canonical, direct, rtol=1e-11, atol=1e-13)


def test_equivalent_cutoff_round_trip():
    coordinates = diag.second_order_bessel_canonical(10_000.0, "mag")
    assert diag.equivalent_cutoff_from_pole(
        coordinates["pole_Hz"], coordinates
    ) == pytest.approx(10_000.0)

    shifted = diag.second_order_bessel_canonical(10_250.0, "mag")
    assert diag.equivalent_cutoff_from_pole(
        shifted["pole_Hz"], coordinates
    ) == pytest.approx(10_250.0)


def _row(score, full_rms, low_rms, high_rms):
    def metrics(value):
        return {"residual_metrics": {"rms_residual_dB": float(value)}}

    return {
        "shape_score": float(score),
        "continuum_metrics_full": metrics(full_rms),
        "continuum_metrics_1_40k": metrics(low_rms),
        "continuum_metrics_40_200k": metrics(high_rms),
    }


def test_comparison_reports_region_deltas():
    reference = _row(1.0, 0.50, 0.10, 0.90)
    candidate = _row(1.2, 0.53, 0.12, 0.91)
    result = diag._comparison(candidate, reference)

    assert result["shape_score_ratio_to_free_pole"] == pytest.approx(1.2)
    assert result["continuum_rms_delta_to_free_pole_dB"] == pytest.approx(0.03)
    assert result["continuum_1_40k_rms_delta_to_free_pole_dB"] == pytest.approx(0.02)
    assert result["continuum_40_200k_rms_delta_to_free_pole_dB"] == pytest.approx(0.01)


def test_direct_pass_requires_score_and_rms_screens():
    reference = _row(1.0, 0.50, 0.10, 0.90)
    screen = {
        "known_filter_max_rms_degradation_dB": 0.05,
        "known_filter_max_score_ratio_to_free_pole": 1.25,
    }
    assert diag._direct_pass(
        _row(1.20, 0.54, 0.12, 0.91), reference, screen
    )
    assert not diag._direct_pass(
        _row(1.40, 0.54, 0.12, 0.91), reference, screen
    )
    assert not diag._direct_pass(
        _row(1.20, 0.57, 0.12, 0.91), reference, screen
    )

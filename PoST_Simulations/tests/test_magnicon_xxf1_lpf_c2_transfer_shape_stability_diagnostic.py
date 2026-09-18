from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_lpf_c2_transfer_shape_stability_diagnostic as diag


def test_readout_components_decompose_total():
    frequency = np.geomspace(1_000.0, 200_000.0, 200)
    readout = {
        "pole_Hz": 10_200.0,
        "pole_Q": 1.0 / np.sqrt(3.0),
        "c2": 20.0,
        "c4": 0.0,
    }
    components = diag.readout_components(
        frequency,
        readout,
        scale_hz=40_000.0,
        reference_hz=1_000.0,
    )
    np.testing.assert_allclose(
        components["total"],
        components["denominator_only"] * components["numerator_only"],
        rtol=1e-12,
        atol=1e-12,
    )


def test_difference_db_identity_zero():
    values = np.array([1.0, 2.0, 3.0])
    np.testing.assert_allclose(
        diag.difference_db(values, values),
        0.0,
        atol=0.0,
    )


def test_band_metrics():
    frequency = np.array([1_000.0, 2_000.0, 3_000.0, 10_000.0])
    values = np.array([0.0, 1.0, -1.0, 2.0])
    rows = diag.band_metrics(
        frequency,
        values,
        [
            {"name": "low", "min": 1_000.0, "max": 3_000.0},
            {"name": "high", "min": 10_000.0, "max": 10_000.0},
        ],
    )
    assert rows["low"]["mean_difference_dB"] == pytest.approx(0.0)
    assert rows["low"]["rms_difference_dB"] == pytest.approx(
        np.sqrt(2.0 / 3.0)
    )
    assert rows["low"]["max_abs_difference_dB"] == pytest.approx(1.0)
    assert rows["high"]["rms_difference_dB"] == pytest.approx(2.0)


def _bands(value):
    return {
        "a": {
            "rms_difference_dB": float(value),
            "max_abs_difference_dB": float(value),
        }
    }


def _screen():
    return {
        "max_full_band_rms_difference_dB": 0.5,
        "max_any_band_rms_difference_dB": 0.75,
        "max_abs_difference_dB": 1.5,
    }


def test_classify_direct_shape_repeatable():
    result = diag.classify(
        {
            "rms_difference_dB": 0.2,
            "max_abs_difference_dB": 0.4,
        },
        _bands(0.3),
        _screen(),
        cross_application_passes=True,
    )
    assert result["direct_transfer_shape_passes_screen"]
    assert result["classification"] == "direct_readout_transfer_shape_repeatable"


def test_classify_detector_refit_absorbs_difference():
    result = diag.classify(
        {
            "rms_difference_dB": 0.8,
            "max_abs_difference_dB": 2.0,
        },
        _bands(1.0),
        _screen(),
        cross_application_passes=True,
    )
    assert not result["direct_transfer_shape_passes_screen"]
    assert result["classification"] == (
        "direct_readout_transfer_differs_but_detector_refit_absorbs_difference"
    )


def test_classify_not_repeatable_without_cross_application():
    result = diag.classify(
        {
            "rms_difference_dB": 0.8,
            "max_abs_difference_dB": 2.0,
        },
        _bands(1.0),
        _screen(),
        cross_application_passes=False,
    )
    assert result["classification"] == "direct_readout_transfer_not_repeatable"

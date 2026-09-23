from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_filter_order_thermal_backbone_diagnostic as orderdiag
from subScript import magnicon_xxf1_filter_state_thermal_backbone_diagnostic as diag


def test_order_zero_is_exact_unity_for_filter_off():
    frequency = np.geomspace(1.0, 2_000_000.0, 400)
    actual = orderdiag.analog_bessel_magnitude(
        frequency,
        order=0,
        cutoff_hz=10_000.0,
        norm="phase",
    )
    np.testing.assert_array_equal(actual, np.ones_like(frequency))


def test_order_zero_normalized_response_is_unity():
    frequency = np.asarray([2_000.0, 10_000.0, 100_000.0])
    actual = orderdiag.normalized_magnicon_magnitude(
        frequency,
        order=0,
        cutoff_hz=9_750.0,
        norm="phase",
    )
    np.testing.assert_array_equal(actual, np.ones_like(frequency))


def test_negative_order_remains_invalid():
    with pytest.raises(ValueError):
        orderdiag.analog_bessel_magnitude(
            np.asarray([10_000.0]),
            order=-1,
            cutoff_hz=10_000.0,
            norm="phase",
        )


def test_compare_states_marks_material_off_improvement():
    off = {"joint_shape_score": 0.4, "joint_continuum_rms_dB": 0.55}
    on = {"joint_shape_score": 1.0, "joint_continuum_rms_dB": 1.10}
    cfg = {
        "materiality_screen": {
            "max_off_over_on_score_ratio": 0.8,
            "min_off_rms_improvement_dB": 0.20,
            "target_good_fit_rms_dB": 0.65,
        }
    }
    result = diag.compare_states(off, on, cfg)
    assert result["filter_off_materially_better"] is True
    assert result["filter_off_reaches_target_good_fit"] is True
    assert result["classification"] == (
        "filter_off_reaches_good_fit_and_beats_documented_on"
    )


def test_config_preserves_documented_on_and_explicit_off_state():
    path = (
        ROOT
        / "config"
        / "magnicon_xxf1_filter_state_thermal_backbone_diagnostic_config.json"
    )
    cfg = json.loads(path.read_text(encoding="utf-8"))
    assert cfg["states"] == ["off", "on"]
    assert cfg["magnicon_filter"]["documented_order"] == 2
    assert cfg["magnicon_filter"]["cutoff_Hz"]["min"] == pytest.approx(9750.0)
    assert cfg["magnicon_filter"]["cutoff_Hz"]["max"] == pytest.approx(10250.0)
    assert "unity" in cfg["magnicon_filter"]["off_semantics"].lower()

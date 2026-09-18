from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as diag


def test_equivalent_first_order_zero_hz():
    assert diag.equivalent_first_order_zero_hz(
        16.0, 40_000.0
    ) == pytest.approx(10_000.0)
    assert diag.equivalent_first_order_zero_hz(
        20.0, 40_000.0
    ) == pytest.approx(40_000.0 / (20.0**0.5))
    assert diag.equivalent_first_order_zero_hz(0.0, 40_000.0) is None


def test_equivalent_first_order_zero_rejects_negative_c2():
    with pytest.raises(ValueError):
        diag.equivalent_first_order_zero_hz(-1.0, 40_000.0)


def test_symmetric_fractional_difference():
    assert diag.symmetric_fractional_difference(
        9_000.0, 9_000.0
    ) == pytest.approx(0.0)
    assert diag.symmetric_fractional_difference(
        9_000.0, 11_000.0
    ) == pytest.approx(0.2)
    assert diag.symmetric_fractional_difference(0.0, 1.0) is None


def test_cross_apply_pass_uses_rms_and_score():
    screen = {
        "max_cross_applied_rms_degradation_dB": 0.05,
        "max_cross_applied_shape_score_ratio": 1.25,
    }
    assert diag._cross_apply_pass(
        {
            "continuum_rms_delta_dB": 0.04,
            "shape_score_ratio": 1.20,
        },
        screen,
    )
    assert not diag._cross_apply_pass(
        {
            "continuum_rms_delta_dB": 0.06,
            "shape_score_ratio": 1.20,
        },
        screen,
    )
    assert not diag._cross_apply_pass(
        {
            "continuum_rms_delta_dB": 0.04,
            "shape_score_ratio": 1.30,
        },
        screen,
    )


def _summary(zero_hz, required=True):
    return {
        "c2_required_by_nested_screen": bool(required),
        "equivalent_first_order_zero_Hz": float(zero_hz),
    }


def _comparison(rms=0.02, score=1.05):
    return {
        "continuum_rms_delta_dB": float(rms),
        "shape_score_ratio": float(score),
    }


def _screen():
    return {
        "max_cross_applied_rms_degradation_dB": 0.05,
        "max_cross_applied_shape_score_ratio": 1.25,
        "max_equivalent_zero_fractional_difference": 0.20,
    }


def test_classify_repeatable_c2_candidate():
    result = diag._classify(
        reference_summary=_summary(9_000.0),
        repeat_summary=_summary(9_500.0),
        cross_reference_from_repeat=_comparison(),
        cross_repeat_from_reference=_comparison(),
        screen=_screen(),
    )
    assert result["c2_required_on_both_days"]
    assert result["equivalent_zero_repeatable_within_screen"]
    assert result["bidirectional_cross_application_passes"]
    assert result["classification"] == (
        "c2_residual_repeatable_across_adjacent_days_"
        "supports_stable_readout_shape_candidate"
    )


def test_classify_required_but_not_cross_day_repeatable():
    result = diag._classify(
        reference_summary=_summary(9_000.0),
        repeat_summary=_summary(9_500.0),
        cross_reference_from_repeat=_comparison(rms=0.08),
        cross_repeat_from_reference=_comparison(),
        screen=_screen(),
    )
    assert result["c2_required_on_both_days"]
    assert not result["bidirectional_cross_application_passes"]
    assert result["classification"] == "c2_required_but_not_cross_day_repeatable"


def test_classify_requirement_not_reproduced():
    result = diag._classify(
        reference_summary=_summary(9_000.0, required=False),
        repeat_summary=_summary(9_500.0, required=True),
        cross_reference_from_repeat=_comparison(),
        cross_repeat_from_reference=_comparison(),
        screen=_screen(),
    )
    assert not result["c2_required_on_both_days"]
    assert result["classification"] == "c2_requirement_not_reproduced_on_both_days"

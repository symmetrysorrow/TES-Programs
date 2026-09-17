from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_lpf_c2_necessity_diagnostic as diag


def _metrics(value):
    return {"residual_metrics": {"rms_residual_dB": float(value)}}


def _row(score, full_rms, low_rms, high_rms):
    return {
        "shape_score": float(score),
        "continuum_metrics_full": _metrics(full_rms),
        "continuum_metrics_1_40k": _metrics(low_rms),
        "continuum_metrics_40_200k": _metrics(high_rms),
    }


def _screen():
    return {
        "max_c2_free_shape_score_ratio_to_c2_zero": 0.80,
        "min_c2_rms_improvement_dB": 0.05,
        "max_no_c2_rms_degradation_to_current_free_pole_dB": 0.05,
        "max_no_c2_shape_score_ratio_to_current_free_pole": 1.25,
        "max_c2_free_rms_degradation_to_current_free_pole_dB": 0.05,
        "max_c2_free_shape_score_ratio_to_current_free_pole": 1.25,
        "cutoff_bound_relative_tolerance": 0.0001,
    }


def test_nested_gain_reports_improvements_with_expected_sign():
    zero = _row(1.0, 0.60, 0.20, 1.00)
    free = _row(0.70, 0.52, 0.15, 0.94)

    gain = diag._nested_gain(free, zero)

    assert gain["shape_score_ratio_c2_free_to_c2_zero"] == pytest.approx(0.70)
    assert gain["continuum_rms_improvement_dB"] == pytest.approx(0.08)
    assert gain["continuum_1_40k_rms_improvement_dB"] == pytest.approx(0.05)
    assert gain["continuum_40_200k_rms_improvement_dB"] == pytest.approx(0.06)


def test_c2_material_requires_both_score_and_rms_thresholds():
    screen = _screen()
    assert diag._c2_material(
        {
            "shape_score_ratio_c2_free_to_c2_zero": 0.75,
            "continuum_rms_improvement_dB": 0.06,
        },
        screen,
    )
    assert not diag._c2_material(
        {
            "shape_score_ratio_c2_free_to_c2_zero": 0.85,
            "continuum_rms_improvement_dB": 0.06,
        },
        screen,
    )
    assert not diag._c2_material(
        {
            "shape_score_ratio_c2_free_to_c2_zero": 0.75,
            "continuum_rms_improvement_dB": 0.04,
        },
        screen,
    )


def test_passes_current_uses_matching_prefix_thresholds():
    screen = _screen()
    reference = _row(1.0, 0.50, 0.10, 0.90)

    assert diag._passes_current(
        _row(1.20, 0.54, 0.12, 0.91),
        reference,
        screen,
        "no_c2",
    )
    assert not diag._passes_current(
        _row(1.30, 0.54, 0.12, 0.91),
        reference,
        screen,
        "no_c2",
    )
    assert not diag._passes_current(
        _row(1.20, 0.56, 0.12, 0.91),
        reference,
        screen,
        "c2_free",
    )


def _branch(*, material, zero_pass, free_pass):
    return {
        "norm": "phase",
        "profiled_c2_material_improvement": bool(material),
        "profiled_c2_zero_passes_current_screen": bool(zero_pass),
        "profiled_c2_free_passes_current_screen": bool(free_pass),
        "profiled_c2_zero_cutoff_Hz": 10_100.0,
        "profiled_c2_free_cutoff_Hz": 10_200.0,
        "profiled_c2_zero_cutoff_boundary_state": {
            "interior": True,
            "at_lower": False,
            "at_upper": False,
        },
        "profiled_c2_free_cutoff_boundary_state": {
            "interior": True,
            "at_lower": False,
            "at_upper": False,
        },
        "profiled_c2_free": {
            "readout": {"c2": 12.0},
            "readout_boundary_hits": {
                "c2": {"at_lower": False, "at_upper": False}
            },
        },
        "profiled_c2_gain": {
            "shape_score_ratio_c2_free_to_c2_zero": 0.95,
            "continuum_rms_improvement_dB": 0.01,
            "continuum_1_40k_rms_improvement_dB": 0.01,
            "continuum_40_200k_rms_improvement_dB": 0.00,
        },
    }


def test_classification_says_c2_not_required_when_zero_passes_and_gain_small():
    result = diag._classify_primary(
        _branch(material=False, zero_pass=True, free_pass=True),
        _screen(),
    )
    assert result["classification"] == (
        "c2_not_materially_required_under_documented_lpf"
    )
    assert result["profiled_c2_boundary_hit"] is False


def test_classification_says_c2_required_when_only_material_free_branch_passes():
    result = diag._classify_primary(
        _branch(material=True, zero_pass=False, free_pass=True),
        _screen(),
    )
    assert result["classification"] == (
        "c2_required_to_meet_current_fit_screen_under_documented_lpf"
    )


def test_classification_preserves_material_but_not_required_for_screen_case():
    result = diag._classify_primary(
        _branch(material=True, zero_pass=True, free_pass=True),
        _screen(),
    )
    assert result["classification"] == (
        "c2_materially_improves_fit_but_is_not_required_for_current_screen"
    )

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import c2_bound_profile_continuum_diagnostic as diag


def _fit_row(c2, rms, low, high):
    return {
        "shape_score": rms**2,
        "readout": {
            "pole_Hz": 11800.0,
            "pole_Q": 0.66,
            "c2": float(c2),
            "c4": 34.23,
        },
        "profiled_white_asd_A_rtHz": 6.0e-11,
        "continuum_metrics_full": {
            "residual_metrics": {"rms_residual_dB": float(rms)}
        },
        "continuum_metrics_1_40k": {
            "residual_metrics": {"rms_residual_dB": float(low)}
        },
        "continuum_metrics_40_200k": {
            "residual_metrics": {"rms_residual_dB": float(high)}
        },
        "detector_candidate": {
            "alpha": 160.0,
            "beta": 0.02,
            "C_tes": 7.0e-13,
            "L": 1.1e-10,
            "T_bath": 0.217,
        },
    }


def _result(c2, rms=0.5, low=0.14, high=0.9):
    return {
        "fits": {
            "pole_section_plus_c2": _fit_row(
                c2, rms, low, high
            ),
            "full_order2": _fit_row(
                c2 * 0.8, rms * 0.8, low, high
            ),
        },
        "selection": {
            "selected_family": "full_order2"
        },
    }


def _config():
    return {
        "interpretation": {
            "interior_fraction_threshold": 0.9,
            "material_total_rms_improvement_dB": 0.05,
            "negligible_last_step_rms_improvement_dB": 0.01,
        }
    }


def test_profile_row_uses_fixed_pole_section_family():
    row = diag._profile_row(600.0, _result(420.0))
    assert row["profiled_family"] == "pole_section_plus_c2"
    assert row["c2_value"] == 420.0
    assert row["c2_fraction_of_upper_bound"] == 0.7
    assert not row["c2_at_upper_bound"]


def test_classify_finite_interior_c2():
    rows = [
        {
            "c2_at_upper_bound": True,
            "c2_fraction_of_upper_bound": 1.0,
            "continuum_rms_dB": 0.55,
        },
        {
            "c2_at_upper_bound": False,
            "c2_fraction_of_upper_bound": 0.72,
            "continuum_rms_dB": 0.48,
        },
        {
            "c2_at_upper_bound": False,
            "c2_fraction_of_upper_bound": 0.43,
            "continuum_rms_dB": 0.48,
        },
    ]
    result = diag.classify_profile(rows, _config())
    assert result["classification"] == "finite_interior_c2_supported"


def test_classify_c2_tracks_bound_when_gain_continues():
    rows = [
        {
            "c2_at_upper_bound": True,
            "c2_fraction_of_upper_bound": 1.0,
            "continuum_rms_dB": 0.55,
        },
        {
            "c2_at_upper_bound": True,
            "c2_fraction_of_upper_bound": 1.0,
            "continuum_rms_dB": 0.48,
        },
        {
            "c2_at_upper_bound": True,
            "c2_fraction_of_upper_bound": 1.0,
            "continuum_rms_dB": 0.45,
        },
    ]
    result = diag.classify_profile(rows, _config())
    assert (
        result["classification"]
        == "c2_tracks_bound_missing_shape_freedom_candidate"
    )


def test_classify_bound_limited_but_saturated():
    rows = [
        {
            "c2_at_upper_bound": True,
            "c2_fraction_of_upper_bound": 1.0,
            "continuum_rms_dB": 0.55,
        },
        {
            "c2_at_upper_bound": True,
            "c2_fraction_of_upper_bound": 1.0,
            "continuum_rms_dB": 0.49,
        },
        {
            "c2_at_upper_bound": True,
            "c2_fraction_of_upper_bound": 1.0,
            "continuum_rms_dB": 0.485,
        },
    ]
    result = diag.classify_profile(rows, _config())
    assert (
        result["classification"]
        == "c2_bound_limited_but_fit_saturated"
    )

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_lpf_c4zero_cutoff_profile_diagnostic as diag


def _metrics(rms):
    return {"residual_metrics": {"rms_residual_dB": float(rms)}}


def _row(score, full, low=0.1, high=0.9):
    return {
        "shape_score": float(score),
        "continuum_metrics_full": _metrics(full),
        "continuum_metrics_1_40k": _metrics(low),
        "continuum_metrics_40_200k": _metrics(high),
    }


def test_cutoff_boundary_state_distinguishes_interior_and_edges():
    interior = diag._cutoff_boundary_state(
        10_000.0, 9_750.0, 10_250.0, 1.0e-4
    )
    assert interior["interior"]
    assert not interior["at_lower"]
    assert not interior["at_upper"]

    lower = diag._cutoff_boundary_state(
        9_750.0, 9_750.0, 10_250.0, 1.0e-4
    )
    assert lower["at_lower"]
    assert not lower["interior"]

    upper = diag._cutoff_boundary_state(
        10_250.0, 9_750.0, 10_250.0, 1.0e-4
    )
    assert upper["at_upper"]
    assert not upper["interior"]


def test_comparison_to_reports_full_and_region_deltas():
    reference = _row(1.0, 0.50, 0.10, 0.90)
    candidate = _row(1.1, 0.52, 0.13, 0.91)
    result = diag._comparison_to(candidate, reference)

    assert result["shape_score_ratio"] == pytest.approx(1.1)
    assert result["continuum_rms_delta_dB"] == pytest.approx(0.02)
    assert result["continuum_1_40k_rms_delta_dB"] == pytest.approx(0.03)
    assert result["continuum_40_200k_rms_delta_dB"] == pytest.approx(0.01)


def test_classify_supported_when_best_is_interior_and_passes_both_screens():
    branches = {
        "mag": {
            "norm": "mag",
            "profiled": {
                **_row(1.05, 0.52),
                "readout": {"c2": 8.0},
            },
            "fitted_cutoff_Hz": 10_040.0,
            "cutoff_boundary_state": {"interior": True},
            "passes_current_free_pole_screen": True,
            "passes_free_pole_c4zero_screen": True,
        },
        "phase": {
            "norm": "phase",
            "profiled": {
                **_row(1.10, 0.53),
                "readout": {"c2": 15.0},
            },
            "fitted_cutoff_Hz": 10_100.0,
            "cutoff_boundary_state": {"interior": True},
            "passes_current_free_pole_screen": True,
            "passes_free_pole_c4zero_screen": True,
        },
    }
    result = diag._classify(branches, {})
    assert (
        result["classification"]
        == "documented_lpf_c4zero_within_tolerance_replaces_free_pole_supported"
    )
    assert result["best_normalization"] == "mag"


def test_classify_flags_manual_tolerance_boundary():
    branches = {
        "mag": {
            "norm": "mag",
            "profiled": {
                **_row(1.05, 0.52),
                "readout": {"c2": 8.0},
            },
            "fitted_cutoff_Hz": 9_750.0,
            "cutoff_boundary_state": {
                "interior": False,
                "at_lower": True,
                "at_upper": False,
            },
            "passes_current_free_pole_screen": True,
            "passes_free_pole_c4zero_screen": True,
        }
    }
    result = diag._classify(branches, {})
    assert (
        result["classification"]
        == "documented_lpf_c4zero_numerically_close_but_cutoff_hits_manual_tolerance"
    )


def test_reference_screens_apply_requested_thresholds():
    current = _row(1.0, 0.50)
    zero = _row(1.1, 0.51)
    candidate = _row(1.20, 0.54)
    screen = {
        "max_rms_degradation_to_current_free_pole_dB": 0.05,
        "max_shape_score_ratio_to_current_free_pole": 1.25,
        "max_rms_degradation_to_free_pole_c4zero_dB": 0.05,
    }
    assert diag._passes_current_reference(candidate, current, screen)
    assert diag._passes_c4zero_reference(candidate, zero, screen)

    too_bad = _row(1.20, 0.57)
    assert not diag._passes_current_reference(too_bad, current, screen)

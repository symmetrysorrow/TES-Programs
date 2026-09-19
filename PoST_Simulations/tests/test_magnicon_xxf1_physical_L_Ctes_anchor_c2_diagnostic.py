from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_physical_L_Ctes_anchor_c2_diagnostic as diag


def test_multiplier_key():
    assert diag._multiplier_key(1.0) == "1x"
    assert diag._multiplier_key(0.25) == "0.25x"
    assert diag._multiplier_key(2.0) == "2x"


def test_problem_with_fixed_detector_values_removes_parameters():
    problem = {
        "baseline_detector": {
            "alpha": 100.0,
            "C_tes": 5.0e-13,
            "L": 5.0e-9,
        },
        "detector_names": ("alpha", "C_tes", "L"),
        "detector_bounds": {
            "alpha": (10.0, 200.0),
            "C_tes": (1.0e-13, 5.0e-12),
            "L": (1.0e-10, 1.23e-8),
        },
    }
    result = diag._problem_with_fixed_detector_values(
        problem,
        {
            "L": 1.0e-10,
            "C_tes": 8.0e-13,
        },
    )
    assert result["baseline_detector"]["L"] == pytest.approx(1.0e-10)
    assert result["baseline_detector"]["C_tes"] == pytest.approx(8.0e-13)
    assert result["detector_names"] == ("alpha",)
    assert problem["detector_names"] == ("alpha", "C_tes", "L")


def test_problem_with_fixed_detector_values_rejects_out_of_bounds():
    problem = {
        "baseline_detector": {"C_tes": 5.0e-13},
        "detector_names": ("C_tes",),
        "detector_bounds": {"C_tes": (1.0e-13, 5.0e-12)},
    }
    with pytest.raises(ValueError):
        diag._problem_with_fixed_detector_values(
            problem,
            {"C_tes": 1.0e-14},
        )


def test_equivalent_zero_hz():
    assert diag._equivalent_zero_hz(4.0, 40_000.0) == pytest.approx(
        20_000.0
    )
    assert diag._equivalent_zero_hz(0.0, 40_000.0) is None


def test_transfer_passes():
    metrics = {
        "full_1_200k": {
            "rms_difference_dB": 0.2,
            "max_abs_difference_dB": 0.8,
        }
    }
    screen = {
        "max_full_band_rms_difference_dB": 0.5,
        "max_abs_difference_dB": 1.5,
    }
    assert diag._transfer_passes(metrics, screen)


def test_transfer_rejects_max_abs():
    metrics = {
        "full_1_200k": {
            "rms_difference_dB": 0.2,
            "max_abs_difference_dB": 1.8,
        }
    }
    screen = {
        "max_full_band_rms_difference_dB": 0.5,
        "max_abs_difference_dB": 1.5,
    }
    assert not diag._transfer_passes(metrics, screen)


def test_classify_both_not_material():
    assert diag._classify(False, False, False) == (
        "physical_L_Ctes_anchors_remove_material_c2_need_on_both_days"
    )


def test_classify_repeatable_material():
    assert diag._classify(True, True, True) == (
        "physical_L_Ctes_anchors_leave_repeatable_material_c2_residual"
    )


def test_classify_nonrepeatable_material():
    assert diag._classify(True, True, False) == (
        "physical_L_Ctes_anchors_leave_nonrepeatable_material_c2_residual"
    )


def test_classify_day_dependent():
    assert diag._classify(True, False, True) == (
        "physical_L_Ctes_anchors_leave_day_dependent_c2_requirement"
    )

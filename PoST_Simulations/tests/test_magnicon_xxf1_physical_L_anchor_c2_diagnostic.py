from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_physical_L_anchor_c2_diagnostic as diag


def test_multiplier_key():
    assert diag._multiplier_key(1.0) == "1x"
    assert diag._multiplier_key(3.0) == "3x"
    assert diag._multiplier_key(0.5) == "0.5x"


def test_problem_with_fixed_L_removes_L_from_nuisance():
    problem = {
        "baseline_detector": {
            "alpha": 100.0,
            "L": 5.0e-9,
        },
        "detector_names": ("alpha", "L"),
        "detector_bounds": {
            "alpha": (10.0, 200.0),
            "L": (1.0e-10, 1.23e-8),
        },
    }
    result = diag._problem_with_fixed_L(problem, 3.0e-10)
    assert result["baseline_detector"]["L"] == pytest.approx(3.0e-10)
    assert result["detector_names"] == ("alpha",)
    assert problem["baseline_detector"]["L"] == pytest.approx(5.0e-9)


def test_problem_with_fixed_L_rejects_out_of_bounds():
    problem = {
        "baseline_detector": {"L": 5.0e-9},
        "detector_names": ("L",),
        "detector_bounds": {"L": (1.0e-10, 1.23e-8)},
    }
    with pytest.raises(ValueError):
        diag._problem_with_fixed_L(problem, 1.0e-11)


def test_transfer_repeatability_passes():
    metrics = {
        "full_1_200k": {
            "rms_difference_dB": 0.3,
            "max_abs_difference_dB": 1.0,
        }
    }
    screen = {
        "max_full_band_rms_difference_dB": 0.5,
        "max_abs_difference_dB": 1.5,
    }
    assert diag.transfer_repeatability_passes(metrics, screen)


def test_transfer_repeatability_rejects_rms():
    metrics = {
        "full_1_200k": {
            "rms_difference_dB": 0.6,
            "max_abs_difference_dB": 1.0,
        }
    }
    screen = {
        "max_full_band_rms_difference_dB": 0.5,
        "max_abs_difference_dB": 1.5,
    }
    assert not diag.transfer_repeatability_passes(metrics, screen)


def test_classify_both_c2_not_material():
    assert diag.classify_primary(False, False, False) == (
        "physical_L_anchor_removes_material_c2_need_on_both_days"
    )


def test_classify_repeatable_material_c2():
    assert diag.classify_primary(True, True, True) == (
        "physical_L_anchor_leaves_repeatable_material_c2_residual"
    )


def test_classify_nonrepeatable_material_c2():
    assert diag.classify_primary(True, True, False) == (
        "physical_L_anchor_leaves_nonrepeatable_material_c2_residual"
    )


def test_classify_day_dependent_requirement():
    assert diag.classify_primary(True, False, True) == (
        "physical_L_anchor_leaves_day_dependent_c2_requirement"
    )

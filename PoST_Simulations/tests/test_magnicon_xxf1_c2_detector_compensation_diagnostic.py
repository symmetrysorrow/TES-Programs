from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as diag


def test_recovery_summary_identifies_top_single():
    summary = diag.recovery_summary(
        c2_free_rms=0.50,
        frozen_c2_zero_rms=1.20,
        all_detector_c2_zero_rms=0.70,
        single_rows={
            "alpha": 0.80,
            "beta": 1.05,
            "L": 0.90,
        },
    )
    assert summary["c2_removal_penalty_with_detector_frozen_dB"] == pytest.approx(
        0.70
    )
    assert summary["all_detector_rms_recovery_dB"] == pytest.approx(0.50)
    assert summary["remaining_rms_penalty_after_all_detector_refit_dB"] == pytest.approx(
        0.20
    )
    assert summary["top_single_parameter"] == "alpha"
    assert summary["top_single"]["fraction_of_all_detector_recovery"] == pytest.approx(
        0.80
    )


def test_classify_strong_single_parameter():
    summary = diag.recovery_summary(
        c2_free_rms=0.50,
        frozen_c2_zero_rms=1.20,
        all_detector_c2_zero_rms=0.70,
        single_rows={"alpha": 0.80, "beta": 1.10},
    )
    screen = {
        "min_material_all_detector_rms_recovery_dB": 0.05,
        "single_parameter_strong_fraction_of_all_recovery": 0.70,
        "single_parameter_mixed_fraction_of_all_recovery": 0.40,
    }
    result = diag.classify_recovery(summary, screen)
    assert result["all_detector_recovery_is_material"]
    assert result["top_single_parameter"] == "alpha"
    assert result["classification"] == (
        "one_parameter_recovers_most_detector_compensation"
    )


def test_classify_mixed_compensation():
    summary = diag.recovery_summary(
        c2_free_rms=0.50,
        frozen_c2_zero_rms=1.20,
        all_detector_c2_zero_rms=0.70,
        single_rows={"alpha": 0.92, "beta": 0.95},
    )
    screen = {
        "min_material_all_detector_rms_recovery_dB": 0.05,
        "single_parameter_strong_fraction_of_all_recovery": 0.70,
        "single_parameter_mixed_fraction_of_all_recovery": 0.40,
    }
    result = diag.classify_recovery(summary, screen)
    assert result["classification"] == "mixed_compensation_with_leading_parameter"


def test_classify_distributed_compensation():
    summary = diag.recovery_summary(
        c2_free_rms=0.50,
        frozen_c2_zero_rms=1.20,
        all_detector_c2_zero_rms=0.70,
        single_rows={"alpha": 1.05, "beta": 1.10},
    )
    screen = {
        "min_material_all_detector_rms_recovery_dB": 0.05,
        "single_parameter_strong_fraction_of_all_recovery": 0.70,
        "single_parameter_mixed_fraction_of_all_recovery": 0.40,
    }
    result = diag.classify_recovery(summary, screen)
    assert result["classification"] == (
        "distributed_or_correlated_detector_compensation"
    )


def test_classify_nonmaterial_detector_recovery():
    summary = diag.recovery_summary(
        c2_free_rms=0.50,
        frozen_c2_zero_rms=0.58,
        all_detector_c2_zero_rms=0.55,
        single_rows={"alpha": 0.56},
    )
    screen = {
        "min_material_all_detector_rms_recovery_dB": 0.05,
        "single_parameter_strong_fraction_of_all_recovery": 0.70,
        "single_parameter_mixed_fraction_of_all_recovery": 0.40,
    }
    result = diag.classify_recovery(summary, screen)
    assert not result["all_detector_recovery_is_material"]
    assert result["classification"] == (
        "detector_nuisance_does_not_materially_compensate_c2"
    )


def test_parameter_shift_summary():
    rows = diag.parameter_shift_summary(
        {"alpha": 100.0, "L": 1.0e-9},
        {"alpha": 120.0, "L": 2.0e-9},
        ("alpha", "L"),
        {
            "alpha": (10.0, 200.0),
            "L": (1.0e-10, 1.0e-8),
        },
        {
            "alpha": {"at_lower": False, "at_upper": False},
            "L": {"at_lower": False, "at_upper": False},
        },
    )
    assert rows["alpha"]["delta"] == pytest.approx(20.0)
    assert rows["alpha"]["ratio_c2_zero_over_c2_free"] == pytest.approx(1.2)
    assert rows["L"]["ratio_c2_zero_over_c2_free"] == pytest.approx(2.0)

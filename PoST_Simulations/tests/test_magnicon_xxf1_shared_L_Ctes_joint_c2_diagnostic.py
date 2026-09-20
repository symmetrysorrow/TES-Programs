from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_shared_L_Ctes_joint_c2_diagnostic as diag


def test_encoded_boundary_linear():
    row = diag._encoded_boundary(1.0, 0.0, 1.0, log=False)
    assert row["at_upper"]
    assert not row["at_lower"]


def test_encoded_boundary_log():
    row = diag._encoded_boundary(
        1.0e-10,
        1.0e-10,
        1.0e-9,
        log=True,
    )
    assert row["at_lower"]
    assert not row["at_upper"]


def test_combined_rms():
    ref = {"continuum_rms_dB": 3.0}
    rep = {"continuum_rms_dB": 4.0}
    assert diag._combined_rms(ref, rep) == pytest.approx(
        np.sqrt(12.5)
    )


def test_nested_summary_material():
    zero = {
        "joint_shape_score": 1.0,
        "joint_continuum_rms_dB": 1.0,
        "solution": {"shared_c2_boundary_hit": None},
    }
    free = {
        "joint_shape_score": 0.5,
        "joint_continuum_rms_dB": 0.8,
        "solution": {
            "shared_c2_boundary_hit": {
                "at_lower": False,
                "at_upper": False,
            }
        },
    }
    screen = {
        "max_shared_c2_free_joint_score_ratio_to_c2_zero": 0.8,
        "min_joint_rms_improvement_dB": 0.05,
    }
    result = diag._nested_summary(zero, free, screen)
    assert result["shared_c2_material_improvement"]
    assert result["classification"] == (
        "shared_c2_material_after_shared_L_Ctes_profile"
    )


def test_nested_summary_not_material():
    zero = {
        "joint_shape_score": 1.0,
        "joint_continuum_rms_dB": 1.0,
        "solution": {"shared_c2_boundary_hit": None},
    }
    free = {
        "joint_shape_score": 0.95,
        "joint_continuum_rms_dB": 0.98,
        "solution": {
            "shared_c2_boundary_hit": {
                "at_lower": False,
                "at_upper": False,
            }
        },
    }
    screen = {
        "max_shared_c2_free_joint_score_ratio_to_c2_zero": 0.8,
        "min_joint_rms_improvement_dB": 0.05,
    }
    result = diag._nested_summary(zero, free, screen)
    assert not result["shared_c2_material_improvement"]
    assert result["classification"] == (
        "shared_L_Ctes_profile_removes_material_c2_need"
    )


def test_nested_summary_upper_bound_limited():
    zero = {
        "joint_shape_score": 1.0,
        "joint_continuum_rms_dB": 2.0,
        "solution": {"shared_c2_boundary_hit": None},
    }
    free = {
        "joint_shape_score": 0.2,
        "joint_continuum_rms_dB": 0.5,
        "solution": {
            "shared_c2_boundary_hit": {
                "at_lower": False,
                "at_upper": True,
            }
        },
    }
    screen = {
        "max_shared_c2_free_joint_score_ratio_to_c2_zero": 0.8,
        "min_joint_rms_improvement_dB": 0.05,
    }
    result = diag._nested_summary(zero, free, screen)
    assert result["classification"] == (
        "shared_c2_material_but_upper_bound_limited_after_"
        "shared_L_Ctes_profile"
    )


def test_decode_shared_parameters_and_c2():
    spec = {
        "names": (
            "shared_log10_C_tes",
            "shared_log10_L",
            "reference_alpha",
            "reference_beta",
            "reference_T_bath",
            "repeat_alpha",
            "repeat_beta",
            "repeat_T_bath",
            "reference_log10_white_scale",
            "repeat_log10_white_scale",
            "shared_c2",
        ),
        "day_names": ("alpha", "beta", "T_bath"),
        "c2_free": True,
    }
    problem = {
        "baseline_detector": {
            "alpha": 100.0,
            "beta": 0.1,
            "T_bath": 0.215,
            "C_tes": 5.0e-13,
            "L": 5.0e-9,
            "R": 0.02,
        }
    }
    vector = np.array(
        [
            np.log10(4.0e-13),
            np.log10(3.0e-10),
            150.0,
            0.2,
            0.216,
            160.0,
            0.3,
            0.217,
            np.log10(1.2),
            np.log10(0.8),
            42.0,
        ]
    )
    result = diag._decode(
        vector,
        spec=spec,
        reference_problem=problem,
        repeat_problem=problem,
        readout_reference={
            "pole_Hz": 10000.0,
            "pole_Q": 0.577,
            "c2": 0.0,
            "c4": 0.0,
        },
    )
    assert result["shared_C_tes"] == pytest.approx(4.0e-13)
    assert result["shared_L"] == pytest.approx(3.0e-10)
    assert result["readout"]["c2"] == pytest.approx(42.0)
    assert result["detectors"]["reference"]["alpha"] == pytest.approx(
        150.0
    )
    assert result["detectors"]["repeat"]["alpha"] == pytest.approx(
        160.0
    )
    assert result["white_scales"]["reference"] == pytest.approx(1.2)
    assert result["white_scales"]["repeat"] == pytest.approx(0.8)

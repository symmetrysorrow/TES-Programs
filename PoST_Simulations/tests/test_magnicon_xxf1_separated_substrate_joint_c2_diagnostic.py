from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from lib.tes_noise_model import THERMAL_EXTENSION_TES_BATH_SERIES
from subScript import magnicon_xxf1_separated_substrate_joint_c2_diagnostic as diag


def test_split_series_conductance_preserves_equivalent():
    g_eq = 4.0e-8
    for ratio in (0.03, 0.3, 1.0, 3.0, 30.0):
        g1, g2 = diag.split_series_conductance(g_eq, ratio)
        assert g1 / g2 == pytest.approx(ratio)
        assert diag.series_equivalent(g1, g2) == pytest.approx(g_eq)


def test_split_series_conductance_rejects_nonpositive():
    with pytest.raises(ValueError):
        diag.split_series_conductance(0.0, 1.0)
    with pytest.raises(ValueError):
        diag.split_series_conductance(1.0, 0.0)


def test_decode_inserts_substrate_node_and_preserves_bath_G():
    spec = {
        "names": (
            "shared_log10_C_tes",
            "shared_log10_L",
            "shared_log10_C_substrate",
            "shared_log10_G_ratio",
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
            "C_tes": 8.0e-13,
            "L": 1.0e-10,
            "G_tes-bath": 4.0e-8,
            "thermal_link_model": "stycast_node",
        }
    }
    vector = np.array(
        [
            np.log10(8.0e-13),
            np.log10(3.0e-10),
            np.log10(7.0e-12),
            np.log10(2.0),
            150.0,
            0.2,
            0.216,
            160.0,
            0.3,
            0.217,
            np.log10(1.2),
            np.log10(0.8),
            12.0,
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

    detector = result["detectors"]["reference"]
    assert detector["thermal_extension"] == THERMAL_EXTENSION_TES_BATH_SERIES
    assert detector["C_substrate"] == pytest.approx(7.0e-12)
    assert diag.series_equivalent(
        detector["G_tes-substrate"],
        detector["G_substrate-bath"],
    ) == pytest.approx(detector["G_tes-bath"])
    assert (
        detector["G_tes-substrate"] / detector["G_substrate-bath"]
        == pytest.approx(2.0)
    )
    assert result["readout"]["c2"] == pytest.approx(12.0)


def _branch(score, rms, *, c2=0.0, c2_upper=False, boundary=False):
    hit = {
        "value": 1.0,
        "lower": 0.0,
        "upper": 2.0,
        "at_lower": bool(boundary),
        "at_upper": False,
    }
    return {
        "joint_shape_score": float(score),
        "joint_continuum_rms_dB": float(rms),
        "solution": {
            "shared_c2_boundary_hit": (
                {
                    "value": c2,
                    "lower": 0.0,
                    "upper": 1000.0,
                    "at_lower": False,
                    "at_upper": bool(c2_upper),
                }
                if c2 != 0.0
                else None
            ),
            "shared_physical_boundary_hits": {
                "C_tes": dict(hit),
                "L": {**hit, "at_lower": False},
                "C_substrate": {**hit, "at_lower": False},
                "G_ratio": {**hit, "at_lower": False},
            },
        },
    }


def test_nested_summary_reports_removed_c2_need():
    zero = _branch(1.0, 2.0)
    free = _branch(0.95, 1.98, c2=10.0)
    screen = {
        "max_shared_c2_free_joint_score_ratio_to_c2_zero": 0.8,
        "min_joint_rms_improvement_dB": 0.05,
    }
    result = diag._nested_summary(zero, free, screen)
    assert not result["shared_c2_material_improvement"]
    assert result["classification"] == (
        "separated_substrate_node_removes_material_c2_need"
    )


def test_nested_summary_marks_profile_boundary():
    zero = _branch(1.0, 3.0, boundary=True)
    free = _branch(0.1, 0.5, c2=20.0)
    screen = {
        "max_shared_c2_free_joint_score_ratio_to_c2_zero": 0.8,
        "min_joint_rms_improvement_dB": 0.05,
    }
    result = diag._nested_summary(zero, free, screen)
    assert result["shared_c2_material_improvement"]
    assert result["shared_physical_profile_boundary_active"]
    assert result["classification"] == (
        "shared_c2_material_but_separated_substrate_profile_"
        "boundary_active"
    )

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_zero_vs_bypass_thermal_backbone_diagnostic as diag


def test_direct_comparison_prefers_pure_zero_when_materially_better():
    pure = {"joint_shape_score": 0.010, "joint_continuum_rms_dB": 0.60}
    bypass = {"joint_shape_score": 0.020, "joint_continuum_rms_dB": 0.90}
    screen = {
        "direct_pair_max_score_ratio": 0.80,
        "direct_pair_min_rms_advantage_dB": 0.15,
    }
    result = diag.compare_pure_zero_and_bypass(
        pure, bypass, screen=screen
    )
    assert result["pure_zero_materially_better"] is True
    assert result["bypass_noise_materially_better"] is False
    assert result["classification"] == (
        "pure_zero_materially_better_than_bypass_noise"
    )


def test_direct_comparison_prefers_bypass_when_materially_better():
    pure = {"joint_shape_score": 0.025, "joint_continuum_rms_dB": 1.00}
    bypass = {"joint_shape_score": 0.010, "joint_continuum_rms_dB": 0.70}
    screen = {
        "direct_pair_max_score_ratio": 0.80,
        "direct_pair_min_rms_advantage_dB": 0.15,
    }
    result = diag.compare_pure_zero_and_bypass(
        pure, bypass, screen=screen
    )
    assert result["pure_zero_materially_better"] is False
    assert result["bypass_noise_materially_better"] is True
    assert result["classification"] == (
        "bypass_noise_materially_better_than_pure_zero"
    )


def test_replacement_metrics_can_match_c2_and_improve_zero():
    candidate = {
        "joint_shape_score": 0.012,
        "joint_continuum_rms_dB": 0.70,
    }
    zero = {
        "joint_shape_score": 0.030,
        "joint_continuum_rms_dB": 1.30,
    }
    c2 = {
        "joint_shape_score": 0.010,
        "joint_continuum_rms_dB": 0.62,
    }
    screen = {
        "max_score_ratio_to_c2_reference": 1.50,
        "max_rms_delta_from_c2_reference_dB": 0.15,
        "max_score_ratio_to_c2_zero_for_material_improvement": 0.80,
        "min_rms_improvement_vs_c2_zero_dB": 0.05,
    }
    result = diag.replacement_metrics(
        candidate, zero=zero, c2=c2, screen=screen
    )
    assert result["within_c2_reference_screen"] is True
    assert result["material_improvement_vs_c2_zero"] is True


def test_config_keeps_bypass_location_after_magnicon_filter():
    path = (
        ROOT
        / "config"
        / "magnicon_xxf1_zero_vs_bypass_thermal_backbone_diagnostic_config.json"
    )
    cfg = json.loads(path.read_text(encoding="utf-8"))
    semantics = cfg["bypass_electronics_noise"]["semantics"].lower()
    assert "bypasses the magnicon 10 khz lpf" in semantics
    assert cfg["base_thermal_backbone_config"] == (
        "magnicon_xxf1_thermal_backbone_pure_zero_diagnostic_config.json"
    )
    bounds = cfg["bypass_electronics_noise"][
        "shared_asd_scale_to_geometric_mean_tracked_white"
    ]
    assert bounds["min"] == pytest.approx(0.01)
    assert bounds["max"] == pytest.approx(20.0)

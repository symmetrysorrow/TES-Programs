from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_high_L_thermal_backbone_c2_diagnostic as diag


def test_parse_elmer_nh_expression():
    assert diag._parse_nh_expression("12.3[nH]") == pytest.approx(12.3)
    assert diag._parse_nh_expression(" 1e2 [nH] ") == pytest.approx(100.0)


def test_fixed_L_profile_is_extremely_narrow():
    base = {
        "shared_physical_profile": {
            "L_H": {"min": 1e-10, "max": 1.3e-9}
        }
    }
    fixed = diag._base_config_at_fixed_L(
        base, 12.3e-9, 1e-6
    )
    bounds = fixed["shared_physical_profile"]["L_H"]
    assert bounds["nominal"] == pytest.approx(12.3e-9)
    assert bounds["min"] == pytest.approx(12.3e-9 * (1 - 1e-6))
    assert bounds["max"] == pytest.approx(12.3e-9 * (1 + 1e-6))
    assert base["shared_physical_profile"]["L_H"]["max"] == pytest.approx(1.3e-9)


def test_nested_metrics_reports_material_c2():
    zero = {
        "joint_shape_score": 0.030,
        "joint_continuum_rms_dB": 1.4,
        "solution": {"shared_c2": 0.0},
    }
    c2 = {
        "joint_shape_score": 0.010,
        "joint_continuum_rms_dB": 0.8,
        "solution": {
            "shared_c2": 16.0,
            "readout_boundary_hit": {"at_upper": False},
        },
    }
    screen = {
        "max_c2_over_zero_shape_score_ratio": 0.8,
        "min_joint_rms_improvement_dB": 0.05,
        "target_good_fit_rms_dB": 0.65,
        "_scale_hz": 40000.0,
    }
    result = diag._nested_metrics(zero, c2, screen)
    assert result["c2_material_improvement"] is True
    assert result["equivalent_first_order_zero_Hz"] == pytest.approx(10000.0)
    assert result["classification"] == "fixed_L_still_requires_material_c2"


def test_config_contains_elmer_and_extended_L_grid():
    path = (
        ROOT
        / "config"
        / "magnicon_xxf1_high_L_thermal_backbone_c2_diagnostic_config.json"
    )
    cfg = json.loads(path.read_text(encoding="utf-8"))
    grid = cfg["L_grid_nH"]
    assert 12.3 in grid
    assert max(grid) >= 100.0
    assert cfg["elmer_reference"]["expected_parameter_expression"] == "12.3[nH]"

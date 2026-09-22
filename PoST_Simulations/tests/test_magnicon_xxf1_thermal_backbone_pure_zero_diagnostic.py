from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_thermal_backbone_pure_zero_diagnostic as diag


def test_pb_absorber_nominal_0p5mm_values():
    row = diag.pb_absorber_values(
        0.5e-3,
        length_m=20e-3,
        width_m=1e-3,
        rho_kg_per_m3=9860.0,
        cp_J_per_kgK=3.26e-5,
        k_W_per_mK=1.68e-2,
    )
    assert row["C_abs_J_per_K"] == pytest.approx(3.21436e-9)
    assert row["G_abs_abs_W_per_K"] == pytest.approx(4.2e-7)


def test_pb_absorber_thickness_scales_C_and_G_linearly():
    low = diag.pb_absorber_values(
        0.4e-3,
        length_m=20e-3,
        width_m=1e-3,
        rho_kg_per_m3=9860.0,
        cp_J_per_kgK=3.26e-5,
        k_W_per_mK=1.68e-2,
    )
    high = diag.pb_absorber_values(
        0.6e-3,
        length_m=20e-3,
        width_m=1e-3,
        rho_kg_per_m3=9860.0,
        cp_J_per_kgK=3.26e-5,
        k_W_per_mK=1.68e-2,
    )
    assert high["C_abs_J_per_K"] / low["C_abs_J_per_K"] == pytest.approx(1.5)
    assert high["G_abs_abs_W_per_K"] / low["G_abs_abs_W_per_K"] == pytest.approx(1.5)


def test_pure_zero_c2_roundtrip():
    scale = 40_000.0
    zero = 560.0
    c2 = diag.pure_zero_to_c2(zero, scale)
    assert c2 == pytest.approx((40_000.0 / 560.0) ** 2)
    assert diag.c2_to_pure_zero(c2, scale) == pytest.approx(zero)


def test_config_has_downward_tbath_and_measured_pb_thickness():
    config_path = (
        ROOT
        / "config"
        / "magnicon_xxf1_thermal_backbone_pure_zero_diagnostic_config.json"
    )
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    tbath = cfg["day_specific_bounds"]["T_bath_K"]
    thickness = cfg["shared_thermal_backbone_profile"]["Pb_absorber"][
        "thickness_m"
    ]
    assert tbath["min"] < 0.215
    assert tbath["max"] > 0.215
    assert thickness["min"] == pytest.approx(0.4e-3)
    assert thickness["nominal"] == pytest.approx(0.5e-3)
    assert thickness["max"] == pytest.approx(0.6e-3)

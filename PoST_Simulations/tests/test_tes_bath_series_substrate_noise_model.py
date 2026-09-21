"""Regression tests for the diagnostic TES--substrate--bath series node."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
POST = ROOT / "PoST_Simulations"
if str(POST) not in sys.path:
    sys.path.insert(0, str(POST))

from lib.tes_noise_model import (  # noqa: E402
    THERMAL_EXTENSION_TES_BATH_SERIES,
    linearized_matrix,
    noise_components,
    operating_point,
    source_matrix,
)


def _parameters():
    g_eq = 4.0e-8
    ratio = 2.0
    g_tes_substrate = g_eq * (1.0 + ratio)
    g_substrate_bath = g_eq * (1.0 + ratio) / ratio
    return {
        "C_abs": 4.5e-9,
        "C_tes": 8.0e-13,
        "C_stycast": 1.1e-11,
        "C_substrate": 7.0e-12,
        "G_abs-abs": 5.9e-7,
        "G_abs-tes": 2.6e-8,
        "G_tes-stycast": 3.0e-8,
        "G_stycast-abs": 3.0e-8,
        "G_tes-bath": g_eq,
        "G_tes-substrate": g_tes_substrate,
        "G_substrate-bath": g_substrate_bath,
        "R": 0.0175,
        "R_l": 0.004,
        "T_c": 0.25,
        "T_bath": 0.215,
        "alpha": 0.5,
        "beta": 1.0,
        "L": 5.0e-9,
        "n": 3.0,
        "excess_johnson_M": 0.0,
        "thermal_link_model": "stycast_node",
        "thermal_extension": THERMAL_EXTENSION_TES_BATH_SERIES,
    }


def test_substrate_series_node_preserves_dc_bath_conductance():
    p = _parameters()
    point = operating_point(p)
    expected = p["G_tes-bath"]
    assert point["G_substrate_equivalent_W_per_K"] == pytest.approx(
        expected
    )
    assert point["G_substrate_equivalent_over_G_tes_bath"] == pytest.approx(
        1.0
    )


def test_substrate_series_node_has_expected_state_and_source_counts():
    p = _parameters()
    assert linearized_matrix(p, 1_000.0).shape == (9, 9)
    assert source_matrix(p).shape == (9, 12)


def test_substrate_series_node_noise_components_are_finite():
    p = _parameters()
    point = operating_point(p)
    assert point["valid"]
    assert point["stable"]

    result = noise_components(
        p,
        np.array([1_000.0, 10_000.0, 100_000.0]),
    )
    assert result["thermal_link_model"] == "stycast_node"
    assert result["thermal_extension"] == THERMAL_EXTENSION_TES_BATH_SERIES
    assert len(result["source_names"]) == 12
    assert result["components_ch0"].shape == (3, 12)
    assert "TES_substrate_TFN" in result["source_class_indices"]
    assert "substrate_bath_TFN" in result["source_class_indices"]
    assert np.all(np.isfinite(result["total_ch0"]))
    assert np.all(result["total_ch0"] > 0.0)


def test_substrate_series_node_rejects_nonpositive_capacity():
    p = _parameters()
    p["C_substrate"] = 0.0
    point = operating_point(p)
    assert not point["valid"]
    assert not point["stable"]

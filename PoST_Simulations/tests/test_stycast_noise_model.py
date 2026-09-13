"""Regression tests for the optional seven-state Stycast noise model."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
POST = ROOT / "PoST_Simulations"
if str(POST) not in sys.path:
    sys.path.insert(0, str(POST))

from lib.tes_noise_model import (  # noqa: E402
    linearized_matrix,
    noise_components,
    operating_point,
    source_matrix,
)


def _parameters():
    return {
        "C_abs": 4.5e-9,
        "C_tes": 8.0e-13,
        "C_stycast": 1.1e-11,
        "G_abs-abs": 5.9e-7,
        "G_abs-tes": 2.6e-8,
        "G_tes-stycast": 3.0e-8,
        "G_stycast-abs": 3.0e-8,
        "G_tes-bath": 4.0e-8,
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
    }


def test_stycast_model_has_expected_state_and_source_counts():
    parameters = _parameters()
    assert linearized_matrix(parameters, 1_000.0).shape == (7, 7)
    assert source_matrix(parameters).shape == (7, 10)


def test_stycast_noise_components_are_finite_and_use_dynamic_sources():
    parameters = _parameters()
    point = operating_point(parameters)
    assert point["valid"]
    assert point["stable"]

    result = noise_components(
        parameters,
        np.array([1_000.0, 10_000.0, 100_000.0]),
    )
    assert result["thermal_link_model"] == "stycast_node"
    assert len(result["source_names"]) == 10
    assert result["components_ch0"].shape == (3, 10)
    assert np.all(np.isfinite(result["total_ch0"]))
    assert np.all(result["total_ch0"] > 0.0)

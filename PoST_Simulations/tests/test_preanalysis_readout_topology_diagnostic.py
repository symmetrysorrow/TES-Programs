from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import preanalysis_readout_topology_diagnostic as diag  # noqa: E402


def test_general_second_order_unity_when_pole_and_zero_match():
    frequency = np.geomspace(100.0, 500_000.0, 401)
    response = diag.second_order_magnitude(
        frequency,
        pole_hz=12_000.0,
        pole_q=0.25,
        zero_hz=12_000.0,
        zero_q=0.25,
    )
    np.testing.assert_allclose(response, np.ones_like(response), rtol=1e-13)


def test_overdamped_second_order_matches_equivalent_real_roots():
    frequency = np.geomspace(100.0, 300_000.0, 501)
    center = 18_000.0
    q = 0.25
    roots = diag.equivalent_real_roots(center, q)
    assert roots is not None
    second_order_denominator = diag.second_order_magnitude(
        frequency,
        pole_hz=center,
        pole_q=q,
        zero_hz=center,
        zero_q=1.0,
    )
    numerator = np.sqrt(
        (1.0 - (frequency / center) ** 2) ** 2
        + (frequency / center) ** 2
    )
    expected_denominator = (
        1.0
        / diag.real_pole_zero_magnitude(
            frequency,
            poles_hz=roots,
            zeros_hz=[],
        )
    )
    reconstructed = numerator / expected_denominator
    np.testing.assert_allclose(
        second_order_denominator,
        reconstructed,
        rtol=1e-12,
        atol=1e-12,
    )


def test_real_2p2z_unity_when_corners_match():
    frequency = np.geomspace(100.0, 500_000.0, 401)
    response = diag.real_pole_zero_magnitude(
        frequency,
        poles_hz=[8_000.0, 90_000.0],
        zeros_hz=[8_000.0, 90_000.0],
    )
    np.testing.assert_allclose(response, np.ones_like(response), rtol=1e-13)


def test_root_regime_classification():
    assert diag.root_regime(0.2) == "overdamped_real_roots"
    assert diag.root_regime(0.5) == "critical"
    assert diag.root_regime(2.0) == "complex_conjugate"


def test_decode_real_family_sorts_poles_and_zeros():
    vector = np.log10(
        np.array([80_000.0, 10_000.0, 100_000.0, 20_000.0])
    )
    params = diag.decode_family("real_2p2z_4p", vector)
    assert params["poles_Hz"] == pytest.approx([10_000.0, 80_000.0])
    assert params["zeros_Hz"] == pytest.approx([20_000.0, 100_000.0])


def test_pre_analysis_model_applies_transfer_before_alias_and_not_to_white(
    monkeypatch,
):
    monkeypatch.setattr(
        diag.opt,
        "hardware_filter_magnitude",
        lambda frequency, **kwargs: np.ones_like(
            np.asarray(frequency, dtype=float)
        ),
    )
    context = {
        "frequency_Hz": np.array([1_000.0, 20_000.0]),
        "alias_frequency_Hz": np.array([99_000.0, 80_000.0]),
        "main_intrinsic_asd": np.array([2.0, 2.0]),
        "alias_intrinsic_asd": np.array([1.0, 1.0]),
        "same_bin": np.array([False, False]),
        "post_filter_white_asd_A_rtHz": 3.0,
    }
    params = {
        "poles_Hz": [12_000.0, 150_000.0],
        "zeros_Hz": [20_000.0, 80_000.0],
    }
    result = diag.pre_analysis_model(
        context,
        "real_2p2z_4p",
        params,
    )
    main_h = diag.real_pole_zero_magnitude(
        context["frequency_Hz"],
        params["poles_Hz"],
        params["zeros_Hz"],
    )
    alias_h = diag.real_pole_zero_magnitude(
        context["alias_frequency_Hz"],
        params["poles_Hz"],
        params["zeros_Hz"],
    )
    absolute = np.sqrt(
        (context["main_intrinsic_asd"] * main_h) ** 2
        + (context["alias_intrinsic_asd"] * alias_h) ** 2
        + 3.0**2
    )
    expected = absolute / absolute[0]
    np.testing.assert_allclose(result, expected)


def test_family_bounds_match_core_parameter_counts():
    for family in diag.CORE_FAMILIES:
        names, bounds = diag.family_bounds(
            family,
            1_000.0,
            300_000.0,
            0.55,
            0.10,
            20.0,
        )
        assert len(names) == 4
        assert len(bounds) == 4

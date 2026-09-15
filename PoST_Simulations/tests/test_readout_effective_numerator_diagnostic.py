from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_effective_numerator_diagnostic as diag  # noqa: E402
from subScript import preanalysis_hybrid_pole_profile_diagnostic as hybrid  # noqa: E402


def test_polynomial_coefficients_match_direct_factorization():
    latent = {"u": -1.7, "v": 2.3, "w": 0.8}
    coeff = diag.polynomial_coefficients(3, latent)
    x = np.linspace(0.0, 5.0, 100)
    direct = (
        ((1.0 + latent["u"] * x**2) ** 2 + latent["v"] ** 2 * x**2)
        * (1.0 + latent["w"] ** 2 * x**2)
    )
    polynomial = (
        1.0
        + coeff["c2"] * x**2
        + coeff["c4"] * x**4
        + coeff["c6"] * x**6
    )
    np.testing.assert_allclose(direct, polynomial, rtol=1e-13, atol=1e-13)


def test_nested_orders_are_positive_for_all_frequencies():
    frequency = np.geomspace(1.0, 1.0e7, 500)
    examples = {
        0: {},
        1: {"v": 4.0},
        2: {"u": -3.0, "v": 2.0},
        3: {"u": -3.0, "v": 2.0, "w": 5.0},
    }
    for order, latent in examples.items():
        values = diag.effective_numerator_squared(
            frequency,
            order,
            latent,
            40000.0,
        )
        assert np.all(np.isfinite(values))
        assert np.all(values > 0.0)


def test_hybrid_to_effective_order3_is_algebraically_exact():
    parameters = {
        "pole_Hz": 11820.989566184204,
        "pole_Q": 0.6569188322650389,
        "zero_Hz": 29821.758720326587,
        "zero_Q": 0.14308634813973375,
        "leadlag_zero_Hz": 81127.75255825827,
    }
    frequency = np.geomspace(1000.0, 500000.0, 1000)
    mapped = diag.hybrid_to_effective(parameters, 40000.0)
    effective = diag.effective_transfer_magnitude(
        frequency,
        mapped["pole_Hz"],
        mapped["pole_Q"],
        3,
        mapped["latent"],
        40000.0,
    )
    expected = hybrid.hybrid_transfer_magnitude(
        frequency,
        parameters,
        None,
    )
    np.testing.assert_allclose(
        effective,
        expected,
        rtol=5e-13,
        atol=1e-13,
    )


def test_vector_spec_matches_nested_parameter_count():
    optimizer = {
        "pole_min_Hz": 1000.0,
        "pole_max_Hz": 300000.0,
        "pole_Q_min": 0.1,
        "pole_Q_max": 20.0,
        "u_min": -30.0,
        "u_max": 30.0,
        "v_min": 0.001,
        "v_max": 30.0,
        "w_min": 0.001,
        "w_max": 30.0,
    }
    for order in range(4):
        names, bounds = diag.vector_spec(order, optimizer)
        assert len(names) == 2 + order
        assert len(bounds) == 2 + order


def test_warm_vector_order3_roundtrips_hybrid_parameters():
    parameters = {
        "pole_Hz": 13000.0,
        "pole_Q": 0.7,
        "zero_Hz": 70000.0,
        "zero_Q": 1.2,
        "leadlag_zero_Hz": 5000.0,
    }
    vector = diag.warm_vector(3, parameters, 40000.0)
    decoded = diag.decode_vector(3, vector)
    mapped = diag.hybrid_to_effective(parameters, 40000.0)
    assert decoded["pole_Hz"] == pytest.approx(mapped["pole_Hz"])
    assert decoded["pole_Q"] == pytest.approx(mapped["pole_Q"])
    assert decoded["latent"]["u"] == pytest.approx(mapped["latent"]["u"])
    assert decoded["latent"]["v"] == pytest.approx(mapped["latent"]["v"])
    assert decoded["latent"]["w"] == pytest.approx(mapped["latent"]["w"])


def test_annotate_complexity_uses_score_and_rms_screen():
    rows = [
        {
            "order": 0,
            "shape_score": 2.0,
            "residual_metrics": {"rms_residual_dB": 0.13},
        },
        {
            "order": 1,
            "shape_score": 7.0,
            "residual_metrics": {"rms_residual_dB": 0.12},
        },
        {
            "order": 2,
            "shape_score": 2.0,
            "residual_metrics": {"rms_residual_dB": 0.20},
        },
    ]
    result = diag.annotate_complexity(
        rows,
        free_score=1.0,
        free_rms=0.10,
        screen={
            "near_free_score_ratio": 5.0,
            "near_free_rms_delta_dB": 0.05,
        },
    )
    assert result[0]["near_free_hybrid"] is True
    assert result[1]["near_free_hybrid"] is False
    assert result[2]["near_free_hybrid"] is False


def test_default_git_tracked_inputs_parse():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    manifest_path = diag.ident.resolve_config_path(
        config["manifest"],
        diag.DEFAULT_CONFIG,
    )
    reference_path = diag.ident.resolve_config_path(
        config["reference_transfer"],
        diag.DEFAULT_CONFIG,
    )
    snapshot_path = diag.ident.resolve_config_path(
        config["lowmid_free_snapshot"],
        diag.DEFAULT_CONFIG,
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))

    assert manifest["cases"]
    assert reference["best_profile_row"]["parameters"]["pole_Hz"] > 0.0
    assert snapshot["transfer"]["parameters"]["pole_Hz"] > 0.0
    assert config["effective_numerator"]["orders"] == [0, 1, 2, 3]
    assert (
        config["effective_numerator"]["reference_scale_Hz"]
        == 40000.0
    )


def test_fit_and_holdout_regions_are_disjoint():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    frequency = np.geomspace(1000.0, 200000.0, 601)
    fit = config["fit_region_Hz"]
    hold = config["holdout_region_Hz"]
    fit_mask = (
        (frequency >= float(fit["min"]))
        & (frequency < float(fit["max_exclusive"]))
    )
    hold_mask = (
        (frequency >= float(hold["min"]))
        & (frequency <= float(hold["max"]))
    )
    assert not np.any(fit_mask & hold_mask)


def test_effective_numerator_families_are_exactly_nested():
    frequency = np.geomspace(1000.0, 500000.0, 300)
    scale = 40000.0

    order0 = diag.effective_numerator_squared(
        frequency, 0, {}, scale
    )
    order1_zero = diag.effective_numerator_squared(
        frequency, 1, {"v": 0.0}, scale
    )
    np.testing.assert_allclose(order0, order1_zero, rtol=0.0, atol=0.0)

    latent1 = {"v": 3.2}
    order1 = diag.effective_numerator_squared(
        frequency, 1, latent1, scale
    )
    order2_u_zero = diag.effective_numerator_squared(
        frequency, 2, {"u": 0.0, "v": 3.2}, scale
    )
    np.testing.assert_allclose(order1, order2_u_zero, rtol=0.0, atol=0.0)

    latent2 = {"u": -1.1, "v": 2.4}
    order2 = diag.effective_numerator_squared(
        frequency, 2, latent2, scale
    )
    order3_w_zero = diag.effective_numerator_squared(
        frequency,
        3,
        {"u": -1.1, "v": 2.4, "w": 0.0},
        scale,
    )
    np.testing.assert_allclose(order2, order3_w_zero, rtol=0.0, atol=0.0)

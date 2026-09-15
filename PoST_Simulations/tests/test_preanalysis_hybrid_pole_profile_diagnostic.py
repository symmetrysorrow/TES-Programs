from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import preanalysis_hybrid_pole_profile_diagnostic as diag  # noqa: E402


def test_leadlag_infinite_pole_is_exact_numerator_limit():
    frequency = np.geomspace(100.0, 500_000.0, 401)
    zero_hz = 5_000.0
    infinite = diag.leadlag_magnitude(
        frequency,
        zero_hz,
        None,
    )
    expected = np.sqrt(1.0 + (frequency / zero_hz) ** 2)
    np.testing.assert_allclose(infinite, expected, rtol=1e-14)


def test_large_finite_pole_converges_to_infinite_limit():
    frequency = np.geomspace(100.0, 500_000.0, 401)
    infinite = diag.leadlag_magnitude(
        frequency,
        5_000.0,
        None,
    )
    finite = diag.leadlag_magnitude(
        frequency,
        5_000.0,
        1.0e12,
    )
    np.testing.assert_allclose(finite, infinite, rtol=2e-13)


def test_component_fraction_sample_sums_to_one():
    frequency = np.geomspace(1_000.0, 200_000.0, 601)
    components = {
        "main_asd": np.full_like(frequency, 2.0),
        "alias_asd": np.full_like(frequency, 3.0),
        "white_asd": np.full_like(frequency, 4.0),
    }
    rows = diag.component_fraction_sample(frequency, components)
    assert rows
    for row in rows:
        total = (
            row["main_psd_fraction"]
            + row["first_alias_psd_fraction"]
            + row["post_filter_white_psd_fraction"]
        )
        assert total == pytest.approx(1.0, abs=1e-14)
        assert row["alias_over_main_asd"] == pytest.approx(1.5)


def test_profile_trend_flags_pole_to_infinity():
    rows = [
        {
            "fixed_leadlag_pole_Hz": 150_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.010,
        },
        {
            "fixed_leadlag_pole_Hz": 300_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.007,
        },
        {
            "fixed_leadlag_pole_Hz": 1_000_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.005,
        },
        {
            "fixed_leadlag_pole_Hz": 5_000_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.004,
        },
        {
            "fixed_leadlag_pole_Hz": None,
            "leadlag_pole_is_infinite": True,
            "shape_score": 0.00399,
        },
    ]
    trend = diag.profile_trend(rows)
    assert trend["score_nonincreasing_from_300k_to_5M"] is True
    assert trend["infinite_is_best_profile_point"] is True
    assert trend["interpretation_flags"][
        "supports_pole_to_infinity_nonidentifiability"
    ] is True
    assert trend["interpretation_flags"][
        "supports_finite_high_frequency_corner"
    ] is False


def test_profile_trend_flags_finite_corner():
    rows = [
        {
            "fixed_leadlag_pole_Hz": 150_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.010,
        },
        {
            "fixed_leadlag_pole_Hz": 300_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.006,
        },
        {
            "fixed_leadlag_pole_Hz": 500_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.004,
        },
        {
            "fixed_leadlag_pole_Hz": 1_000_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.006,
        },
        {
            "fixed_leadlag_pole_Hz": 5_000_000.0,
            "leadlag_pole_is_infinite": False,
            "shape_score": 0.008,
        },
        {
            "fixed_leadlag_pole_Hz": None,
            "leadlag_pole_is_infinite": True,
            "shape_score": 0.009,
        },
    ]
    trend = diag.profile_trend(rows)
    assert trend["best_finite_pole_Hz"] == pytest.approx(500_000.0)
    assert trend["interpretation_flags"][
        "supports_finite_high_frequency_corner"
    ] is True
    assert trend["interpretation_flags"][
        "supports_pole_to_infinity_nonidentifiability"
    ] is False


def test_parse_pole_grid():
    assert diag.parse_pole_grid(
        "150000, 300000,1000000"
    ) == pytest.approx((150000.0, 300000.0, 1000000.0))
    with pytest.raises(ValueError):
        diag.parse_pole_grid("150000,-1")


def test_fit_fixed_pole_keeps_exact_warm_start_candidate(monkeypatch):
    from types import SimpleNamespace

    target = np.asarray([10.0, 1.0, 20.0, 2.0, 5.0])
    warm = {
        "pole_Hz": 10.0,
        "pole_Q": 1.0,
        "zero_Hz": 20.0,
        "zero_Q": 2.0,
        "leadlag_zero_Hz": 5.0,
    }
    bad_vector = np.log10(
        np.asarray([100.0, 3.0, 100.0, 3.0, 100.0])
    )

    monkeypatch.setattr(
        diag,
        "pre_analysis_model",
        lambda context, parameters, fixed_pole_hz: np.asarray(
            [
                parameters["pole_Hz"],
                parameters["pole_Q"],
                parameters["zero_Hz"],
                parameters["zero_Q"],
                parameters["leadlag_zero_Hz"],
            ],
            dtype=float,
        ),
    )
    monkeypatch.setattr(
        diag.opt,
        "fit_score",
        lambda model, target_value, frequency, args: float(
            np.sum(
                np.log10(
                    np.asarray(model) / np.asarray(target_value)
                )
                ** 2
            )
        ),
    )
    monkeypatch.setattr(
        diag.opt,
        "weighted_residual_vector",
        lambda model, target_value, frequency, args: np.log10(
            np.asarray(model) / np.asarray(target_value)
        ),
    )
    monkeypatch.setattr(
        diag,
        "differential_evolution",
        lambda *args, **kwargs: SimpleNamespace(
            x=bad_vector.copy(),
            success=True,
            nfev=1,
        ),
    )
    monkeypatch.setattr(
        diag,
        "least_squares",
        lambda fun, x0, **kwargs: SimpleNamespace(
            x=np.asarray(x0, dtype=float).copy(),
            success=True,
            nfev=1,
        ),
    )
    monkeypatch.setattr(
        diag.base,
        "residual_db_metrics",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr(
        diag.base,
        "band_summary",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr(
        diag.topology,
        "named_boundary_positions",
        lambda *args, **kwargs: {},
    )

    result = diag.fit_fixed_pole(
        None,
        {},
        target,
        np.arange(1.0, 6.0),
        SimpleNamespace(),
        center_min_hz=1.0,
        center_max_hz=300_000.0,
        general_q_min=0.1,
        q_max=20.0,
        seed=1,
        de_maxiter=1,
        rms_screen_db=1.0,
        max_screen_db=3.0,
        initial_parameters=warm,
    )

    assert result["warm_start_used"] is True
    assert result["best_candidate_source"] == "warm_start_exact"
    assert result["shape_score"] == pytest.approx(0.0, abs=1e-15)
    assert result["_parameters_raw"] == pytest.approx(warm)

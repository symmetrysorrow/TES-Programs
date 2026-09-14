from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import transfer_model_falsification_diagnostic as diag  # noqa: E402


def args():
    return SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )


def test_pole_zero_transfer_is_normalized_at_reference():
    frequency = np.array([1_000.0, 10_000.0, 100_000.0])
    response = diag.pole_zero_transfer(
        frequency,
        poles_hz=(70_000.0, 180_000.0),
        zeros_hz=(8_000.0,),
        reference_hz=1_000.0,
    )
    assert response[0] == pytest.approx(1.0, abs=1e-14)
    assert np.all(response > 0.0)


def test_all_tested_families_are_proper():
    for _name, n_poles, n_zeros in diag.FAMILIES:
        assert n_poles >= n_zeros


def test_fit_family_recovers_exact_synthetic_transfer(monkeypatch):
    frequency = np.geomspace(1_000.0, 200_000.0, 121)
    exact_pole = 70_000.0
    exact_zero = 8_000.0
    required = diag.pole_zero_transfer(
        frequency,
        poles_hz=(exact_pole,),
        zeros_hz=(exact_zero,),
    )
    model = np.ones_like(frequency)
    target = required.copy()
    exact_x = np.log10([exact_pole, exact_zero])

    def fake_de(*_args, **_kwargs):
        return SimpleNamespace(
            x=exact_x.copy(),
            success=True,
            nfev=5,
        )

    def fake_ls(*_args, **_kwargs):
        return SimpleNamespace(
            x=exact_x.copy(),
            success=True,
            nfev=2,
        )

    monkeypatch.setattr(diag, "differential_evolution", fake_de)
    monkeypatch.setattr(diag, "least_squares", fake_ls)

    row = diag.fit_family(
        "1p1z",
        1,
        1,
        frequency,
        required,
        model,
        target,
        args(),
        corner_min_hz=300.0,
        corner_max_hz=2_000_000.0,
        seed=1,
        de_maxiter=5,
        rms_threshold_db=1.0,
        max_threshold_db=3.0,
    )

    assert row["poles_Hz"][0] == pytest.approx(exact_pole)
    assert row["zeros_Hz"][0] == pytest.approx(exact_zero)
    assert row["rms_residual_dB"] == pytest.approx(0.0, abs=1e-12)
    assert row["max_abs_residual_dB"] == pytest.approx(0.0, abs=1e-12)
    assert row["corrected_shape_score"] == pytest.approx(0.0, abs=1e-12)
    assert row["passes_screen_tolerance"] is True


def test_unity_family_reports_required_transfer_mismatch():
    frequency = np.geomspace(1_000.0, 200_000.0, 81)
    required = np.ones_like(frequency)
    required[frequency > 40_000.0] = 0.5
    model = np.ones_like(frequency)
    target = required.copy()

    row = diag.fit_family(
        "unity",
        0,
        0,
        frequency,
        required,
        model,
        target,
        args(),
        corner_min_hz=300.0,
        corner_max_hz=2_000_000.0,
        seed=1,
        de_maxiter=5,
        rms_threshold_db=1.0,
        max_threshold_db=3.0,
    )

    assert row["n_parameters"] == 0
    assert row["rms_residual_dB"] > 1.0
    assert row["passes_screen_tolerance"] is False

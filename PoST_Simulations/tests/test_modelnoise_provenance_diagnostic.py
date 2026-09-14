from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import modelnoise_provenance_diagnostic as diag  # noqa: E402


def test_shape_difference_ignores_absolute_scale():
    frequency = np.arange(50_001, dtype=float) * 5.0
    fresh = 1.0 / np.sqrt(1.0 + (frequency / 30_000.0) ** 2)
    stored = 23.0 * fresh
    result = diag.shape_difference_analysis(
        stored,
        fresh,
        frequency,
        smooth_width_hz=1_000.0,
    )
    assert result["total"]["rms_dB"] == pytest.approx(
        0.0, abs=1e-12
    )
    assert result["broad_component"]["rms_dB"] == pytest.approx(
        0.0, abs=1e-12
    )
    assert result["narrow_component"]["rms_dB"] == pytest.approx(
        0.0, abs=1e-12
    )


def test_shape_difference_separates_narrow_excursion():
    frequency = np.arange(50_001, dtype=float) * 5.0
    fresh = np.ones_like(frequency)
    broad_db = 0.25 * np.log10(
        np.maximum(frequency, 1_000.0) / 1_000.0
    )
    stored = fresh * 10.0 ** (broad_db / 20.0)
    spike_index = int(np.argmin(np.abs(frequency - 70_000.0)))
    stored[spike_index] *= 10.0 ** (6.0 / 20.0)

    result = diag.shape_difference_analysis(
        stored,
        fresh,
        frequency,
        smooth_width_hz=1_000.0,
    )
    assert result["broad_component"]["rms_dB"] > 0.05
    assert result["narrow_component"]["max_abs_dB"] > 5.0
    excursions = result["narrow_component"][
        "top_separated_excursions"
    ]
    assert excursions
    assert excursions[0]["frequency_Hz"] == pytest.approx(70_000.0)


def test_optimizer_stored_target_matches_historical_frequency_semantics():
    stored = np.linspace(1.0, 2.0, 11)
    rate = 100.0
    fit_frequency = np.array([10.0, 25.0, 40.0])
    target = diag.optimizer_stored_target(
        stored,
        rate,
        fit_frequency,
    )
    historical_frequency = (
        np.arange(len(stored), dtype=float)
        * (rate / 2.0)
        / len(stored)
    )
    expected = np.interp(
        fit_frequency,
        historical_frequency,
        diag.opt.normalize_at(historical_frequency, stored),
    )
    np.testing.assert_allclose(target, expected)


def test_model_target_metrics_zero_for_exact_shape(monkeypatch):
    frequency = np.geomspace(1_000.0, 200_000.0, 101)
    model = np.ones_like(frequency)
    target = np.ones_like(frequency)
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.30,
    )
    result = diag.model_target_metrics(
        model,
        target,
        frequency,
        args,
    )
    assert result["shape_score"] == pytest.approx(0.0, abs=1e-14)
    for row in result["bands"].values():
        assert row["mean_model_over_target_dB"] == pytest.approx(
            0.0, abs=1e-14
        )

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import empirical_digital_transfer_diagnostic as diag  # noqa: E402


def test_bessel_magnitude_has_unity_dc_and_rolls_off():
    frequency = np.array([0.0, 1_000.0, 10_000.0, 40_000.0])
    response = diag.bessel_magnitude(
        frequency,
        100_000.0,
        10_000.0,
        order=2,
        passes=2,
    )
    assert response[0] == pytest.approx(1.0, abs=1e-14)
    assert np.all(response > 0.0)
    assert response[-1] < response[2] < response[1]


def test_residual_metrics_are_zero_for_exact_prediction():
    frequency = np.geomspace(1_000.0, 200_000.0, 301)
    empirical = 1.0 / np.sqrt(1.0 + (frequency / 20_000.0) ** 2)
    metrics = diag.residual_metrics(
        empirical,
        empirical.copy(),
        frequency,
    )
    assert metrics["absolute_rms_residual_dB"] == pytest.approx(
        0.0, abs=1e-14
    )
    assert metrics[
        "shape_normalized_rms_residual_dB"
    ] == pytest.approx(0.0, abs=1e-14)
    for row in metrics["bands"].values():
        assert row["rms_residual_dB"] == pytest.approx(0.0, abs=1e-14)


def test_stored_shape_comparison_ignores_absolute_calibration():
    frequency = np.linspace(0.0, 250_000.0, 50_001)
    recomputed = 1.0 / (
        1.0 + (frequency / 30_000.0) ** 2
    )
    stored = 17.0 * recomputed
    metrics = diag.stored_shape_metrics(
        stored,
        recomputed,
        frequency,
    )
    assert metrics["rms_shape_difference_dB"] == pytest.approx(
        0.0, abs=1e-12
    )
    assert metrics["max_abs_shape_difference_dB"] == pytest.approx(
        0.0, abs=1e-12
    )


def test_paired_finite_record_transfer_returns_matching_grids():
    sample = 512
    rate = 100_000.0
    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    input_asd = np.ones_like(frequency)
    pre, post, transfer = diag.paired_finite_record_transfer(
        input_asd,
        sample,
        rate,
        cutoff_hz=10_000.0,
        records=4,
        seed=123,
    )
    assert pre.shape == frequency.shape
    assert post.shape == frequency.shape
    assert transfer.shape == frequency.shape
    assert np.all(np.isfinite(transfer))
    assert np.all(transfer >= 0.0)

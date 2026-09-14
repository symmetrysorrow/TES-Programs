from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import measurement_chain_convention_diagnostic as diag  # noqa: E402


def test_digital_baseline_matches_production_analysis_response():
    frequency = np.geomspace(1_000.0, 200_000.0, 101)
    rate = 500_000.0
    actual = diag.digital_bessel_magnitude(
        frequency,
        rate,
        10_000.0,
        2,
        2,
    )
    expected = diag.opt.analysis_filter_magnitude(
        frequency,
        rate,
        cutoff_hz=10_000.0,
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-15)


def test_grid_has_648_cases_and_one_production_baseline():
    configs = []
    for analog_cutoff in diag.ANALOG_CUTOFFS_HZ:
        for analog_order in diag.ANALOG_ORDERS:
            for analog_norm in diag.ANALOG_NORMS:
                for digital_cutoff in diag.DIGITAL_CUTOFFS_HZ:
                    for digital_order in diag.DIGITAL_ORDERS:
                        for digital_passes in diag.DIGITAL_PASSES:
                            configs.append(
                                {
                                    "analog_cutoff_Hz": float(analog_cutoff),
                                    "analog_order": int(analog_order),
                                    "analog_norm": str(analog_norm),
                                    "digital_cutoff_Hz": float(digital_cutoff),
                                    "digital_order": int(digital_order),
                                    "digital_passes": int(digital_passes),
                                }
                            )
    assert len(configs) == 648
    assert sum(c == diag.BASELINE_CHAIN for c in configs) == 1


def test_chain_model_uses_only_measurement_chain_context(monkeypatch):
    frequency = np.array([1_000.0, 10_000.0, 100_000.0])
    context = {
        "frequency_Hz": frequency,
        "alias_frequency_Hz": 500_000.0 - frequency,
        "main_intrinsic_asd": np.array([1.0, 2.0, 3.0]),
        "alias_intrinsic_asd": np.array([0.2, 0.3, 0.4]),
        "same_bin": np.array([False, False, False]),
        "post_filter_white_asd_A_rtHz": 0.5,
        "rate_Hz": 500_000.0,
    }

    def fake_hardware(frequency_hz, **_kwargs):
        return np.ones_like(np.asarray(frequency_hz, dtype=float))

    monkeypatch.setattr(diag.opt, "hardware_filter_magnitude", fake_hardware)
    monkeypatch.setattr(
        diag,
        "digital_bessel_magnitude",
        lambda frequency_hz, *_args, **_kwargs: np.ones_like(
            np.asarray(frequency_hz, dtype=float)
        ),
    )

    model = diag.chain_model(context, diag.BASELINE_CHAIN)
    raw = np.sqrt(
        context["main_intrinsic_asd"] ** 2
        + context["alias_intrinsic_asd"] ** 2
        + context["post_filter_white_asd_A_rtHz"] ** 2
    )
    expected = raw / raw[0]
    np.testing.assert_allclose(model, expected)


def test_digital_bessel_rejects_nonphysical_settings():
    frequency = np.array([1_000.0, 10_000.0])
    with pytest.raises(ValueError):
        diag.digital_bessel_magnitude(
            frequency, 500_000.0, 260_000.0, 2, 2
        )
    with pytest.raises(ValueError):
        diag.digital_bessel_magnitude(
            frequency, 500_000.0, 10_000.0, 0, 2
        )

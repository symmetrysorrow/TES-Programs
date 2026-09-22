from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_filter_order_thermal_backbone_diagnostic as orderdiag
from subScript import magnicon_xxf1_free_Q_thermal_backbone_diagnostic as diag


def test_documented_bessel_Q_matches_phase_normalized_second_order():
    frequency = np.geomspace(100.0, 200_000.0, 400)
    q_bessel = 1.0 / np.sqrt(3.0)
    expected = orderdiag.analog_bessel_magnitude(
        frequency,
        order=2,
        cutoff_hz=10_000.0,
        norm="phase",
    )
    actual = diag.second_order_magnitude(
        frequency,
        10_000.0,
        q_bessel,
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-13)


def test_minus3db_frequency_is_consistent_with_transfer():
    f0 = 10_000.0
    for q in (0.3, 1.0 / np.sqrt(3.0), 0.9, 1.5):
        f3 = diag.minus3db_frequency(f0, q)
        magnitude = diag.second_order_magnitude(
            np.asarray([f3]),
            f0,
            q,
        )[0]
        assert magnitude == pytest.approx(1.0 / np.sqrt(2.0), rel=1e-10)


def test_free_Q_keeps_two_pole_high_frequency_slope():
    q_values = (0.25, 1.0 / np.sqrt(3.0), 1.5)
    for q in q_values:
        magnitude = diag.second_order_magnitude(
            np.asarray([100_000.0, 1_000_000.0]),
            10_000.0,
            q,
        )
        decade_db = 20.0 * np.log10(magnitude[1] / magnitude[0])
        assert decade_db == pytest.approx(-40.0, abs=0.5)


def test_config_profiles_Q_but_preserves_documented_order_two():
    path = (
        ROOT
        / "config"
        / "magnicon_xxf1_free_Q_thermal_backbone_diagnostic_config.json"
    )
    cfg = json.loads(path.read_text(encoding="utf-8"))
    filt = cfg["magnicon_filter"]
    assert filt["documented_order"] == 2
    assert filt["documented_bessel_Q"] == pytest.approx(1.0 / np.sqrt(3.0))
    assert filt["free_Q"]["min"] < filt["documented_bessel_Q"]
    assert filt["free_Q"]["max"] > filt["documented_bessel_Q"]
    assert filt["natural_frequency_Hz"]["min"] == pytest.approx(9750.0)
    assert filt["natural_frequency_Hz"]["max"] == pytest.approx(10250.0)

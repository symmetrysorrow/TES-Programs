from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_filter_order_thermal_backbone_diagnostic as diag
from subScript import magnicon_xxf1_c2_physical_replacement_diagnostic as replacement
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1


def test_first_order_phase_bessel_is_single_pole_at_cutoff():
    cutoff = 10_000.0
    frequency = np.asarray([cutoff])
    magnitude = diag.analog_bessel_magnitude(
        frequency,
        order=1,
        cutoff_hz=cutoff,
        norm="phase",
    )[0]
    assert magnitude == pytest.approx(1.0 / np.sqrt(2.0), rel=1e-10)


def test_second_order_phase_bessel_matches_existing_canonical_response():
    cutoff = 10_000.0
    frequency = np.geomspace(100.0, 200_000.0, 300)
    canonical = xxf1.second_order_bessel_canonical(cutoff, "phase")
    expected = replacement.magnicon_bessel_magnitude(
        frequency,
        canonical["pole_Hz"],
        canonical["pole_Q"],
    )
    actual = diag.analog_bessel_magnitude(
        frequency,
        order=2,
        cutoff_hz=cutoff,
        norm="phase",
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-13)


def test_first_over_second_gains_one_asymptotic_pole():
    frequency = np.asarray([100_000.0, 1_000_000.0])
    first = diag.normalized_magnicon_magnitude(
        frequency,
        order=1,
        cutoff_hz=10_000.0,
        norm="phase",
    )
    second = diag.normalized_magnicon_magnitude(
        frequency,
        order=2,
        cutoff_hz=10_000.0,
        norm="phase",
    )
    ratio_db = 20.0 * np.log10(first / second)
    assert ratio_db[1] - ratio_db[0] == pytest.approx(20.0, abs=0.5)


def test_config_keeps_manual_second_order_as_documented_branch():
    path = (
        ROOT
        / "config"
        / "magnicon_xxf1_filter_order_thermal_backbone_diagnostic_config.json"
    )
    cfg = json.loads(path.read_text(encoding="utf-8"))
    assert cfg["orders"] == [1, 2]
    assert cfg["magnicon_filter"]["documented_order"] == 2
    assert cfg["magnicon_filter"]["cutoff_Hz"]["min"] == pytest.approx(9750.0)
    assert cfg["magnicon_filter"]["cutoff_Hz"]["max"] == pytest.approx(10250.0)
    assert "counterfactual" in cfg["magnicon_filter"]["semantics"].lower()

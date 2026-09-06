"""Numerical checks for the shared TES measurement-chain helper."""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))

from PoST_Simulations.subScript.noise_measurement_model import (  # noqa: E402
    ANALYSIS_BESSEL_CUTOFF_HZ,
    HARDWARE_BESSEL_CUTOFF_HZ,
    analysis_filter_magnitude,
    finite_record_post_analysis_asd,
    fold_hardware_asd,
    hardware_filter_magnitude,
)
from PoST_Simulations.lib import general  # noqa: E402


def test_hardware_and_analysis_cutoffs_are_distinct():
    assert HARDWARE_BESSEL_CUTOFF_HZ == 100_000.0
    assert ANALYSIS_BESSEL_CUTOFF_HZ == 10_000.0


def test_analysis_filter_is_forward_backward_bessel():
    frequency = np.array([1_000.0, 10_000.0, 50_000.0])
    measured = analysis_filter_magnitude(
        frequency,
        500_000.0,
        cutoff_hz=10_000.0,
    )
    expected = general.BesselMagnitudeResponse(
        frequency,
        500_000.0,
        10_000.0,
        passes=2,
    )
    np.testing.assert_allclose(measured, expected)
    assert measured[0] > measured[1] > measured[2]


def test_nyquist_is_not_double_counted_in_alias_fold():
    rate = 500_000.0
    frequency = np.array([rate / 2.0])
    folded = fold_hardware_asd(
        frequency,
        np.ones(1),
        np.ones(1),
        rate,
        cutoff_hz=100_000.0,
        order=4,
    )
    expected = general.AnalogBesselMagnitudeResponse(
        frequency,
        100_000.0,
        order=4,
    )
    np.testing.assert_allclose(folded, expected)


def test_phase_and_mag_hardware_normalizations_are_distinct_at_100_khz():
    frequency = np.array([100_000.0])
    phase = hardware_filter_magnitude(
        frequency,
        cutoff_hz=100_000.0,
        order=4,
        norm="phase",
    )[0]
    mag = hardware_filter_magnitude(
        frequency,
        cutoff_hz=100_000.0,
        order=4,
        norm="mag",
    )[0]
    assert phase < mag
    np.testing.assert_allclose(mag, 1.0 / np.sqrt(2.0), rtol=2e-3, atol=0.0)


def test_hardware_bypass_is_unity_and_keeps_alias_psd_fold():
    frequency = np.array([100_000.0, 250_000.0])
    response = hardware_filter_magnitude(
        frequency,
        cutoff_hz=100_000.0,
        order=4,
        norm="phase",
        bypass=True,
    )
    np.testing.assert_allclose(response, np.ones_like(frequency))

    folded = fold_hardware_asd(
        frequency,
        np.array([3.0, 3.0]),
        np.array([4.0, 4.0]),
        500_000.0,
        cutoff_hz=100_000.0,
        order=4,
        bypass=True,
    )
    np.testing.assert_allclose(folded, np.array([5.0, 3.0]))


def test_finite_record_path_applies_software_lowpass():
    sample = 4096
    rate = 100_000.0
    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    input_asd = np.ones_like(frequency)
    estimated = finite_record_post_analysis_asd(
        input_asd,
        sample,
        rate,
        analysis_cutoff_hz=10_000.0,
        records=64,
        seed=1234,
    )
    reference = estimated[np.argmin(np.abs(frequency - 1_000.0))]
    at_2k = estimated[np.argmin(np.abs(frequency - 2_000.0))] / reference
    at_30k = estimated[np.argmin(np.abs(frequency - 30_000.0))] / reference
    assert 0.7 < at_2k < 1.3
    assert at_30k < 0.2

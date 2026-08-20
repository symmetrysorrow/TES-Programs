import sys
from pathlib import Path

import numpy as np
from scipy import signal


SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))

from lib import general  # noqa: E402


def test_bessel_magnitude_response_matches_freqz():
    rate = 500_000.0
    cutoff = 10_000.0
    frequency = np.array([0.0, 1_000.0, 30_000.0, 100_000.0])

    b, a = signal.bessel(2, cutoff / (rate / 2), "low")
    _, expected = signal.freqz(
        b,
        a,
        worN=2 * np.pi * frequency / rate,
    )

    actual = general.BesselMagnitudeResponse(
        frequency,
        rate,
        cutoff,
        passes=2,
    )

    np.testing.assert_allclose(actual, np.abs(expected) ** 2)


def test_bessel_magnitude_response_rejects_out_of_band_frequency():
    with np.testing.assert_raises(ValueError):
        general.BesselMagnitudeResponse([250_001.0], 500_000.0, 10_000.0)


def test_analog_bessel_response_stays_finite_at_digital_nyquist():
    frequency = np.array([100_000.0, 250_000.0])

    response = general.AnalogBesselMagnitudeResponse(
        frequency,
        100_000.0,
    )

    np.testing.assert_allclose(response[0], 1 / np.sqrt(3), rtol=1e-12)
    assert response[1] > 0.1


def test_sim965_four_pole_bessel_normalization():
    # SIM965 manual table 1.1: for the 24 dB/octave (four-pole) Bessel
    # setting, the actual -3 dB frequency is 0.6604 times the display.
    response = general.AnalogBesselMagnitudeResponse(
        [0.6604 * 100_000.0],
        100_000.0,
        order=4,
    )

    np.testing.assert_allclose(response[0], 1 / np.sqrt(2), rtol=2e-4)

import numpy as np

from Analyze_Experimental_Data.tes_analysis.noise_utils import (
    one_sided_asd_from_power,
    voltage_asd_to_pA,
)


def test_white_noise_asd_is_recovered_from_power_average():
    sample = 512
    rate = 2000.0
    rng = np.random.default_rng(1234)
    sigma = 3.0e-3
    window = np.hanning(sample)
    window_power_gain = np.sqrt(np.mean(window**2))
    power = np.zeros(sample // 2 + 1)

    for _ in range(128):
        spectrum = np.fft.rfft(rng.normal(scale=sigma, size=sample) * window)
        power += np.abs(spectrum) ** 2

    asd = one_sided_asd_from_power(
        power / 128,
        sample,
        rate,
        window_power_gain,
    )
    expected = sigma * np.sqrt(2.0 / rate)
    np.testing.assert_allclose(np.median(asd[2:-2]), expected, rtol=0.08)


def test_one_sided_asd_keeps_dc_and_nyquist_without_sqrt2_gain():
    sample = 8
    rate = 8.0
    expected_scale = 1.0 / sample
    asd = one_sided_asd_from_power(
        np.ones(sample // 2 + 1),
        sample,
        rate,
    )

    np.testing.assert_allclose(asd[[0, -1]], expected_scale)
    np.testing.assert_allclose(asd[1:-1], expected_scale * np.sqrt(2.0))


def test_power_average_is_not_amplitude_average():
    # Two record magnitudes 1 and 3 produce sqrt(mean(power))=sqrt(5), not 2.
    mean_power = np.array([5.0, 5.0, 5.0])
    asd = one_sided_asd_from_power(mean_power, sample=4, rate=4.0)
    np.testing.assert_allclose(asd[1], np.sqrt(10.0) / 4.0)
    assert not np.isclose(asd[1], 2.0 * np.sqrt(2.0) / 4.0)


def test_voltage_asd_to_pa_uses_eta_in_uA_per_V():
    np.testing.assert_allclose(
        voltage_asd_to_pA(np.array([2.0e-6]), eta_uA_per_V=3.5),
        np.array([7.0]),
    )

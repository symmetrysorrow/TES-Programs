import json
import sys
from pathlib import Path

import numpy as np
import pytest


SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))

import PoST_Simulation as simulation  # noqa: E402


def test_tes_johnson_m_zero_is_standard_expression():
    expected = np.sqrt(4 * simulation.k_b * 0.1 * 0.01 * (1 + 2 * 1.6))
    np.testing.assert_allclose(
        simulation.tes_johnson_voltage_asd(0.1, 0.01, 1.6, 0.0), expected
    )


def test_tes_johnson_m_one_increases_asd_by_sqrt_two():
    standard = simulation.tes_johnson_voltage_asd(0.1, 0.01, 1.6, 0.0)
    excess = simulation.tes_johnson_voltage_asd(0.1, 0.01, 1.6, 1.0)
    np.testing.assert_allclose(excess / standard, np.sqrt(2.0))


def test_excess_johnson_validation_and_default():
    assert simulation.resolve_excess_johnson_M({}) == 0.0
    for value in (-1, np.inf, "not-a-number"):
        with pytest.raises(ValueError):
            simulation.resolve_excess_johnson_M({"excess_johnson_M": value})

    for args in ((0, 1, 0, 0), (1, 0, 0, 0), (1, 1, -0.6, 0), (1, 1, 0, -1)):
        with pytest.raises(ValueError):
            simulation.tes_johnson_voltage_asd(*args)


def test_constant_johnson_hook_preserves_scalar_and_frequency_independent_values():
    parameters = {"T_c": 0.1, "R": 0.01, "beta": 1.6, "excess_johnson_M": 0.25}
    scalar = simulation.tes_johnson_voltage_asd(
        parameters["T_c"], parameters["R"], parameters["beta"], 0.25
    )
    assert simulation.resolve_tes_johnson_model(parameters) == "constant_M"
    assert simulation.resolve_tes_johnson_voltage_asd(parameters) == scalar
    np.testing.assert_array_equal(
        simulation.resolve_tes_johnson_voltage_asd(
            parameters, np.array([0.0, 1.0, 1000.0])
        ),
        np.full(3, scalar),
    )


def test_resistance_fluctuation_validation_rejects_nonphysical_parameters():
    for updates in (
        {"tes_resistance_fluctuation_model": "lorentzian", "resistance_fluctuation_M0": -1.0, "resistance_fluctuation_tau_s": 1e-5},
        {"tes_resistance_fluctuation_model": "lorentzian", "resistance_fluctuation_M0": 1.0, "resistance_fluctuation_tau_s": 0.0},
    ):
        with pytest.raises(ValueError):
            simulation.resolve_resistance_fluctuation_model(updates)


def test_asd_to_time_and_back_preserves_one_sided_asd():
    sample = 1024
    rate = 1024.0
    frequency = np.fft.rfftfreq(sample, d=1 / rate)
    input_asd = 2.0e-9 * (1.0 + 0.2 * frequency / frequency[-1])

    noise = simulation.generate_noise_from_asd(
        input_asd,
        sample,
        rate,
        rng=np.random.default_rng(1234),
    )
    recovered_asd = simulation.asd_from_rfft(
        np.fft.rfft(noise),
        sample,
        rate,
    )

    # The generator sets the requested Fourier magnitudes exactly; the
    # random phases only change the time-domain realization.
    np.testing.assert_allclose(recovered_asd[1:-1], input_asd[1:-1], rtol=1e-12)


def test_finite_record_spectrum_uses_power_average_and_asd_normalization():
    sample = 512
    rate = 512.0
    frequency = np.fft.rfftfreq(sample, d=1 / rate)
    input_asd = np.full(frequency.shape, 2.0e-9)

    estimated_asd = simulation.finite_record_simulation_spectrum(
        input_asd,
        sample,
        rate,
        cutoff=100.0,
        records=128,
        seed=1234,
    )

    # Frequencies well inside the low-pass passband should retain the input
    # ASD.  A power average has finite-record scatter, so use a robust median.
    passband = estimated_asd[2:12] / input_asd[2:12]
    assert 0.9 < np.median(passband) < 1.1


def test_show_noise_spectrum_converts_current_asd_to_microamps_without_eta(
    tmp_path, monkeypatch
):
    sample = 32
    rate = 32.0
    frequency = np.fft.rfftfreq(sample, d=1 / rate)
    input_path = tmp_path / "input.json"
    input_path.write_text(
        json.dumps({"samples": sample, "rate": rate, "cutoff": 10.0}),
        encoding="utf-8",
    )

    monkeypatch.setattr(simulation, "output", str(tmp_path))
    monkeypatch.setattr(
        simulation,
        "_simulation_pre_analysis_asd",
        lambda freq, _rate, _para: np.ones_like(freq),
    )
    monkeypatch.setattr(
        simulation,
        "finite_record_simulation_spectrum",
        lambda *_args, **_kwargs: np.full(frequency.shape, 3.0),
    )

    simulation._ShowNoiseSpectrum()
    saved_asd = np.loadtxt(tmp_path / "noise_total-bessel100k.dat")

    np.testing.assert_allclose(saved_asd, 3.0e6)

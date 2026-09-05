import json
import sys
from pathlib import Path

import numpy as np
import pytest


SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))

import PoST_Simulation as simulation  # noqa: E402
from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    estimate_one_sided_asd,
)


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


def test_experimental_and_simulation_estimators_have_absolute_and_shape_parity():
    sample = 1024
    rate = 1024.0
    frequency = np.fft.rfftfreq(sample, d=1 / rate)
    input_asd = 2.0e-9 * (1.0 + 0.35 * frequency / frequency[-1])
    records = [
        simulation.generate_noise_from_asd(
            input_asd,
            sample,
            rate,
            rng=np.random.default_rng(seed),
        )
        for seed in range(512)
    ]

    experimental_style, count = estimate_one_sided_asd(
        records,
        sample,
        rate,
        cutoff=256.0,
        remove_mean=True,
    )
    simulation_style = simulation.finite_record_simulation_spectrum(
        input_asd,
        sample,
        rate,
        cutoff=256.0,
        records=512,
        seed=9876,
    )
    assert count == 512
    passband = slice(4, 32)
    # Both absolute level and frequency-dependent shape must agree; the
    # tolerances cover finite-record scatter, not a normalization convention.
    np.testing.assert_allclose(
        np.median(experimental_style[passband] / input_asd[passband]),
        1.0,
        rtol=0.04,
    )
    np.testing.assert_allclose(
        np.median(simulation_style[passband] / input_asd[passband]),
        1.0,
        rtol=0.04,
    )
    shape_band = slice(8, 64)
    exp_shape = experimental_style[shape_band] / experimental_style[32]
    sim_shape = simulation_style[shape_band] / simulation_style[32]
    assert np.median(np.abs(np.log(exp_shape / sim_shape))) < 0.03
    np.testing.assert_allclose(
        np.median(experimental_style[shape_band] / simulation_style[shape_band]),
        1.0,
        rtol=0.03,
    )


def test_reduced_nonlinear_rhs_jacobian_matches_make_noise_linearization():
    parameters = json.loads(
        (SIMULATION_ROOT / "input.json").read_text(encoding="utf-8")
    )
    operating_point = simulation.tes_operating_point(parameters)
    analytic = simulation.tes_linearized_time_matrix(parameters)
    numerical = simulation.numerical_jacobian(
        lambda state: simulation.tes_nonlinear_rhs(state, parameters),
        operating_point["state"],
    )
    np.testing.assert_allclose(analytic, numerical, rtol=2.0e-5, atol=1.0e-5)


def test_stability_diagnostic_reports_poles_and_strict_gate():
    parameters = json.loads(
        (SIMULATION_ROOT / "input.json").read_text(encoding="utf-8")
    )
    stable = simulation.diagnose_linear_stability(
        simulation.tes_linearized_time_matrix(parameters)
    )
    assert stable["unstable_mode"] is False
    assert stable["max_real_part_per_s"] < 0.0
    assert len(stable["eigenvalues_per_s"]) == 5
    assert np.all(stable["pole_frequency_scale_hz"] >= 0.0)

    unstable_parameters = dict(parameters)
    unstable_parameters.update({"alpha": 1000.0, "L": 1.0e-6})
    unstable = simulation.diagnose_linear_stability(
        simulation.tes_linearized_time_matrix(unstable_parameters)
    )
    assert unstable["unstable_mode"] is True
    assert unstable["max_real_part_per_s"] > 0.0


def test_make_noise_strict_stability_gate_runs_before_frequency_solve(
    tmp_path, monkeypatch
):
    parameters = json.loads(
        (SIMULATION_ROOT / "input.json").read_text(encoding="utf-8")
    )
    parameters.update(
        {
            "alpha": 1000.0,
            "L": 1.0e-6,
            "samples": 64,
            "rate": 100_000.0,
            "stability_mode": "strict",
        }
    )
    (tmp_path / "input.json").write_text(
        json.dumps(parameters), encoding="utf-8"
    )
    monkeypatch.setattr(simulation, "output", str(tmp_path))
    with pytest.raises(ValueError, match="unstable TES operating point"):
        simulation.MakeNoise()


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

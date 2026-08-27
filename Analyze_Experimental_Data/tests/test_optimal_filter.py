"""Regression tests for the pulse-height estimation chain.

These cover the pieces that directly set the energy resolution: the baseline
definition, the peak-search window, the optimal-filter normalization, and the
baseline (temperature) calibration.
"""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from Analyze_Experimental_Data.tes_analysis import analysis_utils as general


RATE = 1.0e5
SAMPLE = 1024
PRESAMPLE = 200


def make_config(**analysis):
    config = {
        "Readout": {"Rate": RATE, "Sample": SAMPLE, "PreSample": PRESAMPLE},
        "Analysis": {
            "CutoffFrequency": 1.0e4,
            "PeakAveragePreSample": 5,
            "PeakAveragePostSample": 5,
            "RiseHighRatio": 0.9,
            "RiseLowRatio": 0.2,
            "DecayHighRatio": 0.9,
            "DecayLowRatio": 0.1,
        },
    }
    config["Analysis"].update(analysis)
    return config


def make_pulse(amplitude=1.0):
    """Single-exponential pulse starting at PRESAMPLE."""
    pulse = np.zeros(SAMPLE)
    n = np.arange(SAMPLE - PRESAMPLE)
    pulse[PRESAMPLE:] = amplitude * (1.0 - np.exp(-n / 8.0)) * np.exp(-n / 120.0)
    return pulse


def make_asd(length):
    """Pink-ish one-sided ASD, strictly positive."""
    fq = np.arange(length, dtype=float)
    return 1.0e-3 / np.sqrt(1.0 + fq) + 1.0e-5


# ---------------------------------------------------------------------------
# Baseline definition (item 6)
# ---------------------------------------------------------------------------

def test_baseline_window_defaults_to_full_pretrigger():
    assert general.BaselineWindow(make_config()) == (0, PRESAMPLE)


def test_baseline_window_follows_getpara_base_x_base_w():
    # getpara: data[presamples - base_x : presamples - base_x + base_w]
    config = make_config(BaseStart=150, BaseWidth=50)
    assert general.BaselineWindow(config) == (PRESAMPLE - 150, PRESAMPLE - 150 + 50)


def test_compute_baseline_uses_configured_window_only():
    pulse = np.zeros(SAMPLE)
    pulse[:100] = 5.0  # only inside the default window
    assert general.ComputeBaseline(pulse, make_config()) == pytest.approx(2.5)
    # A window that starts after the contaminated region sees nothing.
    config = general.NormalizeConfig(make_config(BaseStart=100.0, BaseWidth=100.0))
    assert general.ComputeBaseline(pulse, config) == pytest.approx(0.0)


def test_normalize_config_casts_new_baseline_fields_to_int():
    config = general.NormalizeConfig(make_config(BaseStart=1000.0, BaseWidth=500.0))
    assert isinstance(config["Analysis"]["BaseStart"], int)
    assert isinstance(config["Analysis"]["BaseWidth"], int)


# ---------------------------------------------------------------------------
# Peak search window
# ---------------------------------------------------------------------------

def test_peak_search_window_is_whole_record_without_config():
    assert general.PeakSearchWindow(make_config(), SAMPLE) == (0, SAMPLE)


def test_peak_search_window_ignores_pretrigger_spike():
    pulse = make_pulse()
    pulse[10] = 10.0  # pre-trigger spike, much larger than the real pulse

    unrestricted, _ = general.PeakHeight(pulse, make_config())
    restricted, index = general.PeakHeight(
        pulse, make_config(PeakSearchSample=SAMPLE)
    )

    assert index >= PRESAMPLE
    assert restricted < unrestricted
    assert restricted == pytest.approx(np.max(pulse[PRESAMPLE:]), rel=0.2)


# ---------------------------------------------------------------------------
# Optimal filter template (items 3, 4, 5)
# ---------------------------------------------------------------------------

def test_method_list_exposes_all_comparison_variants():
    assert general.CURRENT_METHOD in general.OPTIMAL_FILTER_METHODS
    assert general.CURRENT_NO_TEMPLATE_BESSEL_METHOD in general.OPTIMAL_FILTER_METHODS
    assert general.PSD_OPTIMAL_METHOD in general.OPTIMAL_FILTER_METHODS
    assert general.LEGACY_FFT_METHOD in general.OPTIMAL_FILTER_METHODS


def test_template_bessel_variant_differs_from_plain_variant():
    config = make_config()
    avg = make_pulse()
    asd = make_asd(SAMPLE // 2 + 1)

    with_bessel = general.OptimalFilterTemplate(
        asd, avg, config, method=general.CURRENT_METHOD, plot=False
    )
    without = general.OptimalFilterTemplate(
        asd, avg, config,
        method=general.CURRENT_NO_TEMPLATE_BESSEL_METHOD, plot=False,
    )

    assert with_bessel.shape == without.shape == (SAMPLE,)
    assert not np.allclose(with_bessel, without)


def test_psd_optimal_estimator_is_normalized_to_unit_amplitude():
    config = make_config()
    avg = make_pulse()
    asd = make_asd(SAMPLE // 2 + 1)

    filt = general.OptimalFilterTemplate(
        asd, avg, config, method=general.PSD_OPTIMAL_METHOD, plot=False
    )
    assert np.sum(avg * filt) == pytest.approx(1.0, rel=1e-9)


def test_psd_optimal_estimator_is_linear_in_pulse_amplitude():
    config = make_config()
    avg = make_pulse()
    asd = make_asd(SAMPLE // 2 + 1)

    filt = general.OptimalFilterTemplate(
        asd, avg, config, method=general.PSD_OPTIMAL_METHOD, plot=False,
        AmplitudeScale=3.0,
    )
    assert np.sum(avg * filt) == pytest.approx(3.0, rel=1e-9)
    assert np.sum(2.0 * avg * filt) == pytest.approx(6.0, rel=1e-9)


def test_psd_optimal_template_is_invariant_to_noise_normalization():
    # Only the *shape* of the PSD may matter; a global unit/eta factor must not
    # change the estimator.  This is what makes the ASD/PSD ambiguity harmless
    # for the normalized method.
    config = make_config()
    avg = make_pulse()
    asd = make_asd(SAMPLE // 2 + 1)

    a = general.OptimalFilterTemplate(
        asd, avg, config, method=general.PSD_OPTIMAL_METHOD, plot=False
    )
    b = general.OptimalFilterTemplate(
        asd * 137.0, avg, config, method=general.PSD_OPTIMAL_METHOD, plot=False
    )
    np.testing.assert_allclose(a, b, rtol=1e-9, atol=0)


def test_psd_optimal_beats_asd_weighting_on_simulated_noise():
    """S*/PSD must not be worse than the 1/ASD weighting currently used."""
    config = make_config()
    avg = make_pulse()
    length = SAMPLE // 2 + 1
    asd = make_asd(length)

    psd_filt = general.OptimalFilterTemplate(
        asd, avg, config, method=general.PSD_OPTIMAL_METHOD, plot=False
    )
    asd_filt = general.OptimalFilterTemplate(
        asd, avg, config,
        method=general.CURRENT_NO_TEMPLATE_BESSEL_METHOD, plot=False,
    )
    # Put both filters on the same signal gain so the noise widths compare.
    asd_filt = asd_filt / np.sum(avg * asd_filt)

    rng = np.random.default_rng(7)
    psd_out = []
    asd_out = []
    for _ in range(400):
        phase = np.exp(1j * rng.uniform(0, 2 * np.pi, length))
        spectrum = asd * phase
        spectrum[0] = spectrum[0].real
        spectrum[-1] = spectrum[-1].real
        noise = np.fft.irfft(spectrum, n=SAMPLE)
        psd_out.append(np.sum(noise * psd_filt))
        asd_out.append(np.sum(noise * asd_filt))

    assert np.std(psd_out) <= np.std(asd_out)


def test_unknown_method_is_rejected():
    with pytest.raises(ValueError):
        general.OptimalFilterTemplate(
            make_asd(SAMPLE // 2 + 1), make_pulse(), make_config(),
            method="nope", plot=False,
        )


# ---------------------------------------------------------------------------
# Baseline (temperature) calibration (item 1)
# ---------------------------------------------------------------------------

def _calib_frame():
    base = np.linspace(-0.1, 0.1, 200)
    # PeakOpt sits exactly on a line in Base, so a correct calibration flattens
    # it perfectly.  Peak is deliberately unrelated to Base, so dividing the
    # wrong column by the fitted gain would leave an obvious Base dependence.
    peak_opt = 2.0 + 5.0 * base
    peak = np.full(base.size, 10.0)
    return pd.DataFrame({"key": np.arange(base.size), "Base": base,
                         "Peak": peak, "PeakOpt": peak_opt})


def test_tempcalib_calibrates_the_column_it_fitted():
    df = general.TempCalib(_calib_frame(), plot=False)

    # Base vs PeakOpt was fitted, so PeakOpt/f(Base) must be exactly flat.
    np.testing.assert_allclose(
        df["PeakOptTemp"], df["PeakOpt"].mean(), rtol=1e-9
    )
    # Had the code divided Peak (not PeakOpt) by the fitted gain, the result
    # would still slope with Base.  Confirm that alternative is different.
    wrong = df["Peak"] / (2.0 + 5.0 * df["Base"]) * df["PeakOpt"].mean()
    assert not np.allclose(df["PeakOptTemp"], wrong)


def test_tempcalib_honours_value_and_result_keys():
    df = general.TempCalib(
        _calib_frame(), ValueKey="Peak", ResultKey="PeakTemp", plot=False
    )
    assert "PeakTemp" in df.columns
    assert "PeakOptTemp" not in df.columns
    assert df["PeakTemp"].mean() == pytest.approx(df["Peak"].mean(), rel=1e-3)


def test_tempcalib_can_be_applied_to_every_estimator_column():
    frame = _calib_frame()
    frame["PeakOptPSD"] = frame["PeakOpt"] * 3.0
    for column in ("Peak", "PeakOpt", "PeakOptPSD"):
        frame = general.TempCalib(
            frame, ValueKey=column, ResultKey=f"{column}Temp", plot=False
        )
    assert {"PeakTemp", "PeakOptTemp", "PeakOptPSDTemp"} <= set(frame.columns)
    np.testing.assert_allclose(
        frame["PeakOptPSDTemp"], 3.0 * frame["PeakOptTemp"], rtol=1e-9
    )


# ---------------------------------------------------------------------------
# Resolution comparison (item 7)
# ---------------------------------------------------------------------------

def test_resolution_summary_reports_fwhm_and_ratio():
    rng = np.random.default_rng(11)
    df = pd.DataFrame({
        "key": np.arange(4000),
        "PeakOpt": rng.normal(loc=1.0, scale=0.01, size=4000),
        "PeakOptPSD": rng.normal(loc=1.0, scale=0.005, size=4000),
    })

    summary = general.ResolutionSummary(df, columns=["PeakOpt", "PeakOptPSD"])

    assert list(summary["column"]) == ["PeakOpt", "PeakOptPSD"]
    expected = 0.01 * 2 * np.sqrt(2 * np.log(2))
    assert summary.loc[0, "FWHM"] == pytest.approx(expected, rel=0.1)
    assert summary.loc[1, "FWHM"] < summary.loc[0, "FWHM"]
    np.testing.assert_allclose(
        summary["FWHM/mean"], summary["FWHM"] / summary["fit_mean"].abs(), rtol=1e-12
    )


def test_resolution_summary_survives_a_column_that_cannot_be_fitted():
    df = pd.DataFrame({"key": [0, 1, 2], "PeakOpt": [1.0, 1.0, 1.0]})
    summary = general.ResolutionSummary(df, columns=["PeakOpt"])
    assert len(summary) == 1
    assert summary.loc[0, "mean"] == pytest.approx(1.0)

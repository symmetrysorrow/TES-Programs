"""Contract tests for the exploratory high-frequency parameter search."""
from pathlib import Path


def _source_text() -> str:
    source = (
        Path(__file__).parents[1]
        / "subScript"
        / "explore_high_frequency_noise_parameters.py"
    )
    return source.read_text(encoding="utf-8")


def test_search_is_separate_from_frozen_target_inference():
    text = _source_text()
    assert "exploratory_noise_guided_parameter_search" in text
    assert '"strict_target_parameter_estimate_allowed": False' in text
    assert '"separate_from_frozen_ensemble": True' in text
    assert '"noise_used_for_ranking": True' in text
    assert '"noise_residual_fit": False' in text


def test_search_reuses_original_stage_a_ranges_and_pulse_gate():
    text = _source_text()
    assert "0.5" in text and "2.0" in text
    assert "scenario_slow_pole_time_constant_range_s" in text
    assert "linear_modes" in text
    assert "operating_point" in text


def test_search_uses_correct_measurement_chain_and_1_to_10khz_score():
    text = _source_text()
    assert "expected_post_analysis_asd" in text
    assert "finite_record_post_analysis_asd" in text
    assert "HARDWARE_BESSEL_CUTOFF_HZ" in text
    assert '"selection_band_Hz": [1000.0, 10000.0]' in text
    assert '"low_frequency_excluded_from_score": True' in text


def test_main_plot_compares_only_finite_record_curves_at_nyquist():
    text = _source_text()
    assert "frozen_finite_asd = finite_record_post_analysis_asd" in text
    assert "Frozen best finite" in text
    assert "Frozen best expected" not in text
    assert "plt.xlim(rate / sample, rate / 2.0)" in text
    assert "same finite-record estimator as the experiment" in text

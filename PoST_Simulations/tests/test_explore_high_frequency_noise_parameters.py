"""Contract tests for the exploratory multi-band parameter search."""
from pathlib import Path

from PoST_Simulations.subScript.explore_high_frequency_noise_parameters import (
    SCORE_FIELDS,
    pareto_front,
)


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


def test_search_has_three_independent_frequency_objectives():
    text = _source_text()
    assert '"mid": (1_000.0, 10_000.0)' in text
    assert '"high": (10_000.0, 100_000.0)' in text
    assert '"all": (1_000.0, 100_000.0)' in text
    assert "rms_log_ratio_1_10_kHz" in text
    assert "rms_log_ratio_10_100_kHz" in text
    assert "rms_log_ratio_1_100_kHz" in text
    assert '"low_frequency_below_1kHz_excluded_from_score": True' in text


def test_hardware_is_fixed_and_not_searched():
    text = _source_text()
    assert "HARDWARE_BESSEL_ORDER = 4" in text
    assert 'params["hardware_bessel_order"] = HARDWARE_BESSEL_ORDER' in text
    assert "HARDWARE_BESSEL_CUTOFF_HZ" in text
    assert '"hardware_parameters_searched": False' in text


def test_finite_pool_covers_each_objective_and_uses_same_estimator():
    text = _source_text()
    assert "deterministic_candidate_pool" in text
    assert "for field in SCORE_FIELDS.values()" in text
    assert "finite_record_post_analysis_asd" in text
    assert "finite_best_by_objective" in text
    assert "finite_record_common_seed" in text


def test_pareto_front_keeps_only_non_dominated_rows():
    rows = [
        {"trial_id": "a", SCORE_FIELDS["mid"]: 0.10, SCORE_FIELDS["high"]: 0.40},
        {"trial_id": "b", SCORE_FIELDS["mid"]: 0.20, SCORE_FIELDS["high"]: 0.20},
        {"trial_id": "c", SCORE_FIELDS["mid"]: 0.30, SCORE_FIELDS["high"]: 0.30},
        {"trial_id": "d", SCORE_FIELDS["mid"]: 0.40, SCORE_FIELDS["high"]: 0.10},
    ]
    front = pareto_front(rows)
    assert [row["trial_id"] for row in front] == ["a", "b", "d"]


def test_plots_include_ratio_and_pareto_diagnostics_and_nyquist_limit():
    text = _source_text()
    assert "high_frequency_parameter_search_ratio.png" in text
    assert "high_frequency_parameter_search_pareto.png" in text
    assert "plt.xlim(rate / sample, rate / 2.0)" in text
    assert "Exploratory high-band best" in text
    assert "Exploratory balanced best" in text

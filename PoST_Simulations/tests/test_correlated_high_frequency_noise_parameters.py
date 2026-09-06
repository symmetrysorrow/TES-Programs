"""Contract tests for the correlated multivariate TES noise search."""
import ast
from pathlib import Path


def _text() -> str:
    path = (
        Path(__file__).parents[1]
        / "subScript"
        / "explore_high_frequency_noise_parameters_correlated.py"
    )
    return path.read_text(encoding="utf-8")


def test_correlated_script_is_valid_python():
    ast.parse(_text())


def test_correlated_search_keeps_all_13_unfixed_parameters():
    text = _text()
    assert "SEARCH_PARAMETERS" in text
    assert "searched_parameter_count" in text
    assert "specs_from_envelope" in text
    assert "T_bath_K" in text


def test_covariance_adaptation_learns_parameter_combinations():
    text = _text()
    assert "learned_covariance" in text
    assert "rank_weights" in text
    assert "covariance_refine" in text
    assert "multivariate_normal" in text
    assert "diverse_anchors" in text
    assert "covariance_correlations" in text


def test_full_residual_curve_is_locally_optimized():
    text = _text()
    assert "deterministic_residual" in text
    assert "residual_from_unit" in text
    assert "optimize.least_squares" in text
    assert 'loss="soft_l1"' in text
    assert "least_squares_plan" in text
    assert "least_squares_pareto_seeds" in text


def test_expected_interaction_pairs_are_mapped():
    text = _text()
    for pair in (
        '("L", "R_l")',
        '("L", "beta")',
        '("L", "alpha")',
        '("alpha", "C_tes")',
        '("C_tes", "G_tes-bath")',
    ):
        assert pair in text
    assert "interaction_scan" in text
    assert "correlated_high_frequency_parameter_interactions.png" in text


def test_confirmed_filter_chain_is_not_searched():
    text = _text()
    assert "HARDWARE_BESSEL_ORDER" in text
    assert '"hardware_bessel_cutoff_Hz": 100_000.0' in text
    assert '"analysis_bessel_cutoff_Hz": cutoff' in text
    assert "realize_finite" in text


def test_correlated_search_reports_physical_combination_coordinates():
    text = _text()
    assert "derived_coordinates" in text
    assert '"electrical_pole_Hz"' in text
    assert '"loop_gain"' in text
    assert '"G_eff_W_per_K"' in text


def test_correlated_outputs_include_main_ratio_pareto_interactions_and_json():
    text = _text()
    assert "correlated_high_frequency_parameter_search.png" in text
    assert "correlated_high_frequency_parameter_search_ratio.png" in text
    assert "correlated_high_frequency_parameter_pareto.png" in text
    assert "correlated_high_frequency_parameter_interactions.png" in text
    assert "correlated_high_frequency_parameter_search.json" in text


def test_exploratory_firewall_is_preserved():
    text = _text()
    assert '"strict_target_parameter_estimate_allowed": False' in text
    assert '"noise_residual_fit": False' in text
    assert '"additive_noise_parameter_fit": False' in text
    assert '"simulation_amplitude_rescale": False' in text

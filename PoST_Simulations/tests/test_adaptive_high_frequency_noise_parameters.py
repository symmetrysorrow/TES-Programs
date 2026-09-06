"""Contract tests for the adaptive all-parameter TES noise search."""
from pathlib import Path


def _text() -> str:
    path = (
        Path(__file__).parents[1]
        / "subScript"
        / "explore_high_frequency_noise_parameters_adaptive.py"
    )
    return path.read_text(encoding="utf-8")


def test_adaptive_search_covers_all_unfixed_tes_noise_model_parameters():
    text = _text()
    block = text.split("SEARCH_PARAMETERS = (", 1)[1].split(")", 1)[0]
    for name in (
        "T_c", "R", "R_l", "alpha", "beta", "L", "n", "C_tes",
        "C_abs", "G_tes-bath", "G_abs-tes", "G_abs-abs",
        "excess_johnson_M",
    ):
        assert f'"{name}"' in block
    assert '"T_bath"' not in block
    assert '"R_SH"' not in block


def test_confirmed_hardware_and_analysis_chain_stays_fixed():
    text = _text()
    assert "HARDWARE_BESSEL_ORDER" in text
    assert '"hardware_bessel_cutoff_Hz": 100_000.0' in text
    assert '"analysis_bessel_cutoff_Hz": cutoff' in text
    assert "realize_finite" in text


def test_r_sh_is_provenance_only_and_correlated_with_r():
    text = _text()
    assert "correlated_r_sh" in text
    assert "provenance-only; tes_noise_model does not consume R_SH" in text


def test_search_is_more_than_plain_random_sampling():
    text = _text()
    assert "latin_hypercube" in text
    assert "elite_refinement" in text
    assert "expand_bounds" in text
    assert "polish" in text
    assert "sensitivity_scan" in text
    assert "pareto_front" in text


def test_l_search_range_covers_measurement_band_electrical_poles():
    text = _text()
    assert "300_000.0" in text
    assert "1_000.0" in text
    assert "electrical-pole-informed exploratory range" in text


def test_adaptive_outputs_include_ratio_sensitivity_and_pareto():
    text = _text()
    assert "adaptive_high_frequency_parameter_search_ratio.png" in text
    assert "adaptive_high_frequency_parameter_sensitivity.png" in text
    assert "adaptive_high_frequency_parameter_pareto.png" in text
    assert "adaptive_high_frequency_parameter_search.json" in text


def test_search_preserves_exploratory_firewall():
    text = _text()
    assert '"strict_target_parameter_estimate_allowed": False' in text
    assert '"noise_residual_fit": False' in text
    assert '"simulation_amplitude_rescale": False' in text

"""Contract tests for lower-envelope hardware-convention refinement."""
import ast
from pathlib import Path


def _text() -> str:
    path = (
        Path(__file__).parents[1]
        / "subScript"
        / "refine_hardware_bessel_lower_envelope.py"
    )
    return path.read_text(encoding="utf-8")


def test_lower_envelope_refinement_is_valid_python():
    ast.parse(_text())


def test_raw_experiment_is_preserved_and_lower_target_is_explicit():
    text = _text()
    assert "asymmetric_lower_baseline" in text
    assert "lower_envelope_targets" in text
    assert '"raw_experiment_preserved_for_reporting": True' in text
    assert "Lower-baseline optimization target" in text
    assert "Simulation / raw experiment" in text
    assert "Simulation / optimization target" in text


def test_high_frequency_bump_has_reduced_leverage_not_manual_masking():
    text = _text()
    assert "asymmetric_least_squares_lower_baseline_in_log_frequency_log_ASD" in text
    assert "lower-asymmetry" in text
    assert "lower-smoothness" in text
    assert "lower-blend-end-hz" in text
    assert "spsolve" in text
    assert "soft_l1" in text
    assert "exclude" not in text.lower()


def test_all_three_hardware_conventions_are_supported():
    text = _text()
    assert 'HARDWARE_MODES = ("phase", "mag", "bypass")' in text
    assert '"hardware_bypass_is_diagnostic_only": True' in text
    assert '"hardware_bessel_order": HARDWARE_BESSEL_ORDER' in text
    assert '"hardware_bessel_cutoff_Hz": HARDWARE_BESSEL_CUTOFF_HZ' in text


def test_confirmed_software_analysis_chain_is_unchanged():
    text = _text()
    assert "finite_record_post_analysis_asd" in text
    assert "analysis_cutoff_hz=cutoff" in text
    assert '"analysis_bessel_cutoff_Hz": cutoff' in text
    assert "records = args.finite_records if args.finite_records > 0 else len(exp_paths)" in text


def test_refinement_uses_correlated_search_seeds_and_all_tes_parameters():
    text = _text()
    assert "correlated_high_frequency_parameter_search.json" in text
    assert "SEARCH_PARAMETERS" in text
    assert "pareto_front" in text
    assert "best_deterministic" in text
    assert "parameters_from_unit" in text
    assert "unit_from_parameters" in text


def test_no_additive_noise_or_amplitude_rescale_is_introduced():
    text = _text()
    assert '"noise_residual_fit": False' in text
    assert '"additive_noise_parameter_fit": False' in text
    assert '"simulation_amplitude_rescale": False' in text


def test_expected_outputs_are_named():
    text = _text()
    for name in (
        "lower_envelope_hardware_refinement.json",
        "lower_envelope_hardware_refinement.png",
        "lower_envelope_hardware_refinement_raw_ratio.png",
        "lower_envelope_hardware_refinement_target_ratio.png",
        "lower_envelope_hardware_refinement_target.png",
    ):
        assert name in text

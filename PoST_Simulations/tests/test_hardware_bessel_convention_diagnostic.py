"""Contract tests for the hardware-Bessel convention diagnostic."""
import ast
from pathlib import Path


def _text() -> str:
    path = (
        Path(__file__).parents[1]
        / "subScript"
        / "compare_hardware_bessel_conventions.py"
    )
    return path.read_text(encoding="utf-8")


def test_diagnostic_script_is_valid_python():
    ast.parse(_text())


def test_diagnostic_compares_phase_bypass_and_mag_only():
    text = _text()
    assert '"key": "phase"' in text
    assert '"key": "bypass"' in text
    assert '"key": "mag"' in text
    assert '"norm": "phase"' in text
    assert '"norm": "mag"' in text
    assert '"bypass": True' in text


def test_balanced_parameters_are_held_fixed_without_refit():
    text = _text()
    assert "load_balanced_parameters" in text
    assert '"best_finite"' in text
    assert '"best_deterministic"' in text
    assert '"parameters_refit_for_this_diagnostic": False' in text


def test_every_case_uses_same_finite_record_estimator_and_seed():
    text = _text()
    assert "hardware_sampled_asd" in text
    assert "finite_record_post_analysis_asd" in text
    assert "analysis_cutoff_hz=analysis_cutoff" in text
    assert "seed=int(args.finite_seed)" in text
    assert '"same_finite_record_seed_for_all_cases": True' in text


def test_bypass_is_explicitly_diagnostic_not_physical_replacement():
    text = _text()
    assert "bypass is diagnostic only" in text.lower()
    assert '"bypass_is_diagnostic_only": True' in text


def test_outputs_show_post_analysis_ratio_preanalysis_and_response():
    text = _text()
    assert "hardware_bessel_convention_comparison.png" in text
    assert "hardware_bessel_convention_ratio.png" in text
    assert "hardware_bessel_convention_pre_analysis.png" in text
    assert "hardware_bessel_convention_response.png" in text
    assert "hardware_bessel_convention_diagnostic.json" in text

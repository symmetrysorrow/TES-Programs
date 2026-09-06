"""Contract and small numerical tests for intermediate-target refinement."""
import ast
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))

from PoST_Simulations.subScript.refine_hardware_bessel_intermediate_target import (  # noqa: E402
    BAND_WEIGHTS,
    DEFAULT_LAMBDAS,
    geometric_intermediate,
    scan_grids,
    symmetric_log_trend,
)


def _text() -> str:
    path = (
        Path(__file__).parents[1]
        / "subScript"
        / "refine_hardware_bessel_intermediate_target.py"
    )
    return path.read_text(encoding="utf-8")


def test_intermediate_script_is_valid_python():
    ast.parse(_text())


def test_geometric_intermediate_has_expected_endpoints_and_midpoint():
    lower = np.array([1.0, 4.0, 9.0])
    upper = np.array([4.0, 16.0, 36.0])
    np.testing.assert_allclose(geometric_intermediate(lower, upper, 0.0), lower)
    np.testing.assert_allclose(geometric_intermediate(lower, upper, 1.0), upper)
    np.testing.assert_allclose(
        geometric_intermediate(lower, upper, 0.5),
        np.sqrt(lower * upper),
    )


def test_symmetric_log_trend_reduces_a_narrow_spike():
    x = np.linspace(0.0, 1.0, 101)
    y = -3.0 * x
    y[50] += 1.5
    smooth = symmetric_log_trend(y, smoothness=500.0)
    assert smooth[50] < y[50]
    assert np.isfinite(smooth).all()


def test_scan_explicitly_resolves_80_to_120_khz_edge_band():
    grids = scan_grids(120_000.0)
    assert np.isclose(grids["edge"][0], 80_000.0)
    assert np.isclose(grids["edge"][-1], 120_000.0)
    assert BAND_WEIGHTS["edge"] > BAND_WEIGHTS["mid"]


def test_default_lambda_scan_brackets_midpoint():
    assert 0.5 in DEFAULT_LAMBDAS
    assert min(DEFAULT_LAMBDAS) < 0.5 < max(DEFAULT_LAMBDAS)


def test_selection_is_fixed_to_midpoint_not_raw_convex_feature():
    text = _text()
    assert '"selection_lambda": 0.5' in text
    assert "midpoint_validation_scores" in text
    assert "best_lambda_candidates" in text
    assert "raw_trend_smoothness" in text
    assert "geometric_log_ASD_between_lower_and_smoothed_raw_trend" in text


def test_previous_lower_solution_is_used_as_an_endpoint_seed():
    text = _text()
    assert "previous_lower_" in text
    assert "lower_envelope_hardware_refinement.json" in text
    assert "seed_pool_for_mode" in text


def test_bypass_remains_diagnostic_and_physical_selection_excludes_it():
    text = _text()
    assert '"hardware_bypass_is_diagnostic_only": True' in text
    assert 'if m != "bypass"' in text
    assert "best_physical_hardware_mode" in text
    assert "best_including_bypass_diagnostic" in text


def test_output_includes_lambda_and_midpoint_diagnostics():
    text = _text()
    for name in (
        "intermediate_target_hardware_refinement.json",
        "intermediate_target_hardware_refinement.png",
        "intermediate_target_hardware_refinement_raw_ratio.png",
        "intermediate_target_hardware_refinement_midpoint_ratio.png",
        "intermediate_target_hardware_refinement_lambda_scan.png",
        "intermediate_target_hardware_refinement_target_family.png",
    ):
        assert name in text


def test_exploratory_firewall_is_preserved():
    text = _text()
    assert '"strict_target_parameter_estimate_allowed": False' in text
    assert '"noise_residual_fit": False' in text
    assert '"additive_noise_parameter_fit": False' in text
    assert '"simulation_amplitude_rescale": False' in text

"""Regression checks for the two-stage conditional proxy workflow."""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"
STAGE_A = ROOT / "PoST_Simulations/subScript/build_proxy_envelope.py"
STAGE_B = ROOT / "PoST_Simulations/subScript/run_noise_blind_ensemble.py"
COMPARE = ROOT / "PoST_Simulations/subScript/compare_proxy_ensemble_to_experiment.py"
PULSE_AUDIT = ROOT / "PoST_Simulations/subScript/target_pulse_pole_audit.py"


def _load(name: str) -> dict:
    return json.loads((CASE / name).read_text(encoding="utf-8"))


def test_stage_a_has_no_experimental_spectrum_reader_or_comparison_import():
    source = STAGE_A.read_text(encoding="utf-8")
    assert "comparison_summary" not in source
    assert "CH0_noise" not in source
    tree = ast.parse(source)
    imported = [node.module or "" for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
    assert not any("noise" in module.lower() for module in imported)


def test_stage_b_consumes_frozen_stage_a_outputs_only():
    source = STAGE_B.read_text(encoding="utf-8")
    assert "--generic-input" not in source
    assert "proxy_parameter_envelope.json" in source
    assert "proxy_scenarios.json" in source


def test_every_parameter_has_provenance_and_strict_flags():
    envelope = _load("proxy_parameter_envelope.json")
    required = {"T_c", "T_bath", "R", "R_SH", "R_l", "alpha", "beta", "L", "n", "C_tes", "C_abs", "G_tes-bath", "G_abs-tes", "G_abs-abs", "K"}
    assert required <= set(envelope["parameters"])
    for row in envelope["parameters"].values():
        assert row["source_class"]
        assert "derivation" in row and "confidence" in row and "correlations" in row
        assert row["allowed_for_strict_target"] is False


def test_proxy_and_generic_reference_are_not_target_truth():
    scenarios = _load("proxy_scenarios.json")
    assert scenarios["freeze_status"] == "frozen"
    assert all(not row["strict_target_allowed"] for row in scenarios["scenarios"])
    assert all("simulation_reference_only" in row["source_class_by_parameter"].values() for row in scenarios["scenarios"])


def test_fixed_seed_reproducibility_and_stage_a_manifest():
    manifest = _load("noise_blind_sweep_manifest.json")
    scenarios = _load("proxy_scenarios.json")
    assert manifest["freeze_status"] == "frozen"
    assert manifest["experimental_spectrum_read"] is False
    assert manifest["parameter_ranges_frozen_before_comparison"] is True
    assert scenarios["seed"] == manifest["seed"]
    assert scenarios["sample_count"] == manifest["sample_count"]


def test_stability_exclusions_have_reason_schema_and_sources_are_baseline_only():
    envelope = _load("proxy_noise_envelope.json")
    assert set(envelope["excluded_scenario_schema"]) == {"scenario_id", "reason", "excluded_before_noise"}
    assert all(row["reason"] for row in envelope["excluded_scenarios"])
    assert envelope["empirical_white_noise_A_rtHz"] == 0.0
    assert envelope["empirical_readout_floor_A_rtHz"] == 0.0
    assert envelope["resistance_fluctuation_model"] == "none"
    assert envelope["target_hanging_model"] == "disabled"


def test_comparison_cannot_construct_or_adjust_parameters():
    source = COMPARE.read_text(encoding="utf-8")
    assert "build_proxy_envelope" not in source
    assert "range_adjustment_performed" in source
    summary = _load("conditional_comparison_summary.json")
    assert summary["parameter_generation_called"] is False
    assert summary["range_adjustment_performed"] is False


def test_pre_post_semantics_and_strict_conclusion_are_preserved():
    summary = _load("conditional_comparison_summary.json")
    assert summary["strict_target_conclusion"].startswith("C")
    assert len(summary["rows"]) == len(summary["post_analysis"]["rows"]) == 7
    assert "Bessel" in summary["post_analysis"]["filter"]
    assert _load("provenance.json")["existing_data_search"]["conclusion"] == "C_frozen"


def test_post_factor_matches_production_bessel_helper_numerically():
    sys.path.insert(0, str(ROOT / "PoST_Simulations/subScript"))
    sys.path.insert(0, str(ROOT / "PoST_Simulations"))
    import compare_proxy_ensemble_to_experiment as comparison
    from lib import general

    expected = general.BesselMagnitudeResponse(comparison.FREQUENCIES, 500000.0, 10000.0, passes=2)
    np.testing.assert_allclose(comparison._post_factor(), expected, rtol=0.0, atol=1e-14)


def test_proxy_uses_eight_independent_sources_and_production_flink():
    envelope = _load("proxy_parameter_envelope.json")
    params = envelope["sensitivity_reference"]
    sys.path.insert(0, str(ROOT / "PoST_Simulations/subScript"))
    import proxy_physics

    components, metadata = proxy_physics.noise_components(params, np.array([10.0, 1000.0]))
    assert components.shape == (2, 8)
    assert metadata["F_LINK"] == 0.9
    assert len(metadata["source_names"]) == 8


def test_pulse_audit_is_noise_blind_and_stationarity_result_is_conservative():
    source = PULSE_AUDIT.read_text(encoding="utf-8")
    assert "CH0_noise" not in source and "CH1_noise" not in source
    audit = _load("target_pulse_pole_audit.json")
    assert audit["noise_records_read"] is False
    stationarity = _load("low_frequency_stationarity_audit.json")
    assert stationarity["classification"] in {"stationary_detector_like", "stationary_channel_specific", "inconclusive"}

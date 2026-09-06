"""Structural and scientific-firewall checks for pulse-contamination v3."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUB = ROOT / "PoST_Simulations/subScript"
CASE = ROOT / "PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"


def load(name):
    return json.loads((CASE / name).read_text(encoding="utf-8"))


def test_v3_required_artifacts_exist_and_input_is_unresolved():
    required = {
        "pulse_template_library_v3.json", "pulse_detection_decision_policy_v3.json",
        "pulse_only_block_bootstrap_null.json", "negative_polarity_null_control.json",
        "pulse_detection_dual_null_comparison.json", "noise_record_pulse_classification_v3.json",
        "pulse_detector_injection_recovery_v3.json", "pulse_partitioned_noise_spectra_v3.json",
        "pulse_selection_bias_audit.json", "pulse_mask_robustness.json",
        "record_length_scaling_v3.json", "detrend_sensitivity_v3.json",
        "auxiliary_ch0_ch1_cross_spectrum.json", "clean_v3_simulation_comparison.json",
        "pulse_contamination_v3_summary.json", "pulse_contamination_v3_summary.md",
    }
    assert required <= {p.name for p in CASE.iterdir()}
    target = load("input.json")
    assert target["T_c"] is None and target["R"] is None
    assert target["post_filter_white_asd_A_rtHz"] == 0.0
    assert target["readout_white_asd_A_rtHz"] == 0.0


def test_v3_policy_has_resolution_compatible_categories():
    policy = load("pulse_detection_decision_policy_v3.json")
    assert policy["target_primary_fwer"] == 1e-3
    assert policy["one_e_minus_4_category"] == "removed"
    assert "definite" in policy["operational_boundaries"]
    classification = load("noise_record_pulse_classification_v3.json")
    classes = {row["CH0"]["classification"] for row in classification["records"].values()}
    assert not any("1e-4" in value for value in classes)


def test_v3_raw_injection_and_null_construction_firewalls():
    source = (SUB / "pulse_contamination_v3.py").read_text(encoding="utf-8")
    assert "inject_waveform(raw, waveform" in source
    assert "production_preprocessing_passes_after_injection" in source
    assert "double_filtering" in source
    assert "np.resize" not in source
    assert "negative_template_is_physically_impossible" in source
    injection = load("pulse_detector_injection_recovery_v3.json")
    assert injection["raw_injection"] is True
    assert injection["production_preprocessing_passes_after_injection"] == 1
    assert injection["double_filtering"] is False


def test_v3_templates_and_spectra_semantics_are_explicit():
    library = load("pulse_template_library_v3.json")
    for channel in ("CH0", "CH1"):
        for template in library["channels"][channel]["templates"].values():
            assert template["sample_step"] == 1
            assert template["fullrate"] is True
            assert "raw_values" in template and "processed_values" in template
    spectra = load("pulse_partitioned_noise_spectra_v3.json")
    for subset in spectra["subsets"].values():
        assert "pre_analysis" in subset and "post_analysis" in subset
    assert "no Bessel" in spectra["semantics"]["pre_analysis"]
    assert "Bessel" in spectra["semantics"]["post_analysis"]


def test_v3_recovery_gate_blocks_pc3_and_record_length_promotion():
    recovery = load("pulse_detector_injection_recovery_v3.json")
    summary = load("pulse_contamination_v3_summary.json")
    scaling = load("record_length_scaling_v3.json")
    assert recovery["criteria_pass"] is False
    assert summary["pc_classification"] == "PC4"
    assert summary["stationary_physical_source_investigation_allowed"] is False
    assert scaling["status"] == "blocked_recovery_criteria_failed"


def test_v3_auxiliary_coherence_is_not_production_accepted():
    cross = load("auxiliary_ch0_ch1_cross_spectrum.json")
    assert cross["production_accepted"] is False
    assert cross["pairing"] == "exact event key"
    assert cross["groups"]["all"]["bootstrap"]["replicates"] == 200

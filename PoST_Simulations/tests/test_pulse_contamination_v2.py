"""Regression checks for the record-wise pulse-contamination audit."""

from __future__ import annotations

import ast
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUB = ROOT / "PoST_Simulations/subScript"
CASE = ROOT / "PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"


def _load(name):
    return json.loads((CASE / name).read_text(encoding="utf-8"))


def test_v2_has_required_artifacts_and_strict_input_is_unchanged():
    required = {
        "recordwise_pulse_detection_thresholds.json",
        "pulse_detection_familywise_null.json",
        "pulse_tail_template_library_v2.json",
        "noise_record_pulse_classification_v2.json",
        "pulse_detector_injection_recovery.json",
        "pulse_partitioned_noise_spectra_v2.json",
        "low_frequency_shape_diagnostic.json",
        "record_length_scaling.json",
        "detrend_sensitivity.json",
        "slow_transient_time_domain_audit.json",
        "low_frequency_record_dominance.json",
        "ch1_acceptance_provenance.json",
        "clean_v2_simulation_comparison.json",
        "experimental_pulse_decomposition_reconstruction.json",
        "pulse_contamination_v2_summary.md",
    }
    assert required <= {p.name for p in CASE.iterdir()}
    target = _load("input.json")
    assert target["T_c"] is None and target["R"] is None
    assert target["post_filter_white_asd_A_rtHz"] == 0.0
    assert target["readout_white_asd_A_rtHz"] == 0.0


def test_v2_calibration_is_full_record_familywise_and_pulse_only():
    calibration = _load("recordwise_pulse_detection_thresholds.json")
    assert calibration["all_lags"] is True
    assert calibration["all_templates"] is True
    assert calibration["simulation_spectrum_read"] is False
    assert calibration["noise_spectrum_read"] is False
    assert calibration["primary_fpr"] == 1e-3
    assert calibration["record_duration_s"] == 0.2
    for channel in ("CH0", "CH1"):
        data = calibration["channels"][channel]
        assert data["null_sample_count"] >= 8
        assert set(data["null_score_samples"]) == {"full_pulse", "post_peak_tail", "slow_tail"}


def test_v2_templates_are_full_rate_and_injection_has_tail_cases():
    library = _load("pulse_tail_template_library_v2.json")
    for channel in ("CH0", "CH1"):
        templates = library["channels"][channel]["templates"]
        assert all(row["sample_step"] == 1 and row["fullrate"] for row in templates.values())
        assert library["channels"][channel]["tail_valid_length_ms"] >= 50.0
    injection = _load("pulse_detector_injection_recovery.json")
    assert set(injection["tail_efficiency_by_age_ms"]) >= {"0", "10", "50", "100", "150"}
    assert injection["simulation_spectrum_used"] is False


def test_v2_clean_spectrum_has_one_estimator_and_no_subtraction():
    spectra = _load("pulse_partitioned_noise_spectra_v2.json")
    assert "mean removal" in spectra["estimator"]
    assert "Hann" in spectra["estimator"]
    assert "power average" in spectra["estimator"]
    assert spectra["subsets"]["clean_v2"]["record_count"] > 0
    reconstruction = _load("experimental_pulse_decomposition_reconstruction.json")
    assert reconstruction["pulse_subtraction_from_clean_v2"] is False


def test_v2_script_does_not_fit_noise_residual_or_add_floor():
    source = (SUB / "pulse_contamination_v2.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert "least_squares" not in source
    assert "curve_fit" not in source
    assert "empirical_white" not in source
    assert "readout_white" not in source
    assert ast.dump(tree)

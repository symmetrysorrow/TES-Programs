"""Structural and artifact checks for the pulse-contamination phase."""

from __future__ import annotations

import ast
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUB = ROOT / "PoST_Simulations/subScript"
CASE = ROOT / "PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"


def _load(name):
    return json.loads((CASE / name).read_text(encoding="utf-8"))


def _source(name):
    return (SUB / name).read_text(encoding="utf-8")


def test_template_and_threshold_stages_do_not_read_noise_or_simulation_spectra():
    template = _source("build_pulse_template_library.py")
    threshold = _source("build_pulse_detection_thresholds.py")
    assert "CH0_noise" not in template and "CH1_noise" not in template
    assert "simulation_spectrum" not in template
    assert "comparison_summary" not in template
    assert "proxy_noise_envelope" not in threshold
    assert "proxy_noise_envelope" not in threshold


def test_detection_is_fixed_fpr_and_event_key_pairing_is_explicit():
    source = _source("classify_noise_record_pulses.py")
    assert "record_key" in source
    assert "accepted_noise_indices" in source
    thresholds = _load("pulse_detection_thresholds.json")
    assert thresholds["primary_fpr"] == 1e-3
    assert thresholds["fprs_reported"] == [1e-2, 1e-3, 1e-4]


def test_independent_acceptance_and_same_estimator_are_recorded():
    classification = _load("noise_record_pulse_classification.json")
    assert set(classification["accepted_counts"]) == {"CH0", "CH1"}
    partition = _load("pulse_partitioned_noise_spectra.json")
    assert "Hann" in partition["estimator"]
    assert "mean removal" in partition["estimator"]
    assert set(partition["subsets"]["CH0"]) >= {"all_accepted", "pulse_free", "definitely_contaminated", "ambiguous"}


def test_primary_clean_is_classified_records_without_subtraction():
    source = _source("partition_noise_spectra.py")
    assert "pulse_free" in source and "estimate_one_sided_asd" in source
    assert "filtfilt" not in source
    partition = _load("pulse_partitioned_noise_spectra.json")
    assert partition["subsets"]["CH0"]["pulse_free"]["record_count"] > 0


def test_false_positive_controls_and_channel_coincidence_are_present():
    false_positive = _load("pulse_false_positive_audit.json")
    assert set(false_positive["controls_included"]) >= {"pulse pretrigger baseline", "time-reversed template", "sign-inverted template"}
    coincidence = _load("pulse_channel_coincidence.json")
    assert "exact event key" in coincidence["pairing"]


def test_pulse_only_prediction_uses_event_statistics_without_spectral_fit():
    prediction = _load("pulse_contamination_psd_prediction.json")
    assert prediction["stationary_residual_used"] is False
    assert prediction["amplitude_fit_to_spectrum"] is False
    assert prediction["fixed_seed"] == 20260908
    assert _load("full_acquisition_reproduction.json")["pc_classification"] == "PC3"


def test_strict_target_and_zero_empirical_floors_remain():
    case_input = _load("input.json")
    assert case_input["T_c"] is None and case_input["R"] is None
    envelope = _load("proxy_noise_envelope.json")
    assert envelope["empirical_white_noise_A_rtHz"] == 0.0
    assert envelope["empirical_readout_floor_A_rtHz"] == 0.0
    assert _load("phase_next_final_classification.json")["pc_classification"] == "PC3"

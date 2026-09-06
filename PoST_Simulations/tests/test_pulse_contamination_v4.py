"""Correctness regressions for the v4 pulse-contamination detector."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"
sys.path.insert(0, str(ROOT))

from PoST_Simulations.lib.general import BesselMagnitudeResponse
from PoST_Simulations.subScript import pulse_contamination_v4 as v4


def _library():
    t = np.zeros(8)
    t[2:5] = (1.0, 2.0, 1.0)
    return {"templates": {"full_pulse": {"processed_values": t.tolist()}, "post_peak_tail": {"processed_values": (t / 2).tolist()}}}


def test_v4_statistic_has_no_record_local_energy_normalization():
    rng = np.random.default_rng(4)
    record = rng.normal(size=64)
    one = v4.scan_processed(record, _library(), 1.0)
    two = v4.scan_processed(2.0 * record, _library(), 1.0)
    assert two["full_pulse"]["rho"] == 2.0 * one["full_pulse"]["rho"]
    assert "energy" not in v4.scan_processed.__code__.co_names


def test_v4_filter_transfer_is_filtfilt_asd_magnitude():
    rate = 500_000.0
    cutoff = 10_000.0
    freq = np.array([10.0, 1_000.0, 30_000.0])
    b, a = signal.bessel(2, cutoff / (rate / 2.0), "low")
    _, h = signal.freqz(b, a, worN=2 * np.pi * freq / rate)
    np.testing.assert_allclose(BesselMagnitudeResponse(freq, rate, cutoff, passes=2), np.abs(h) ** 2)


def test_v4_explicit_subset_sets_include_ambiguous_classes():
    source = (ROOT / "PoST_Simulations/subScript/pulse_contamination_v4.py").read_text(encoding="utf-8")
    assert '"ambiguous_full", "ambiguous_tail"' in source
    assert '"full_contaminated": {"definite_full", "likely_full"}' in source


def test_v4_artifacts_and_gate_are_measured():
    required = {
        "pulse_detector_statistic_v4.json",
        "pulse_pretrigger_noise_model.json",
        "raw_pulse_polarity_audit.json",
        "pulse_baseline_autocorrelation.json",
        "pulse_block_length_selection.json",
        "pulse_block_bootstrap_null_v4.json",
        "negative_polarity_validation_v4.json",
        "pulse_detection_decision_policy_v4.json",
        "noise_record_pulse_classification_v4.json",
        "pulse_detector_injection_recovery_v4.json",
        "pulse_detector_recovery_gate_v4.json",
        "pulse_false_positive_validation_v4.json",
        "pulse_selection_bias_audit_v4.json",
        "pulse_mask_robustness_v4.json",
        "pulse_partitioned_noise_spectra_v4.json",
        "record_length_scaling_v4.json",
        "auxiliary_ch0_ch1_cross_spectrum_v4.json",
        "common_source_topology_predictions.json",
        "common_source_topology_comparison.json",
        "clean_v4_simulation_comparison.json",
        "pulse_contamination_v4_summary.json",
        "pulse_contamination_v4_summary.md",
    }
    assert required <= {p.name for p in CASE.iterdir()}
    gate = json.loads((CASE / "pulse_detector_recovery_gate_v4.json").read_text(encoding="utf-8"))
    assert gate["criteria_pass"] == all(gate["criteria"].values())
    null = json.loads((CASE / "pulse_block_bootstrap_null_v4.json").read_text(encoding="utf-8"))
    assert null["channels"]["CH0"]["record_count"] >= 999
    assert null["fwer_1e-3_resolvable"] is True

"""Summarize pulse-contamination evidence and supersede the provisional R5."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    args = parser.parse_args()
    case = args.case_dir
    classification = json.loads((case / "noise_record_pulse_classification.json").read_text(encoding="utf-8"))
    partition = json.loads((case / "pulse_partitioned_noise_spectra.json").read_text(encoding="utf-8"))
    prediction = json.loads((case / "pulse_contamination_psd_prediction.json").read_text(encoding="utf-8"))
    thresholds = json.loads((case / "pulse_contamination_threshold_sensitivity.json").read_text(encoding="utf-8"))
    false_positive = json.loads((case / "pulse_false_positive_audit.json").read_text(encoding="utf-8"))
    coincidence = json.loads((case / "pulse_channel_coincidence.json").read_text(encoding="utf-8"))
    coherence = json.loads((case / "partitioned_coherence.json").read_text(encoding="utf-8"))
    simulation = json.loads((case / "pulse_free_simulation_comparison.json").read_text(encoding="utf-8"))
    ch0 = partition["subsets"]["CH0"]
    ratios = ch0["ratio_anchors"]
    pulse = prediction["channels"]["CH0"]
    measured = np.asarray(pulse["measured_all_asd"])
    pulse_asd = np.asarray(pulse["pulse_only_asd"])
    pulse_fraction = {str(f): float((pulse_asd[int(float(f) / 5)] ** 2) / max(measured[int(float(f) / 5)] ** 2, np.finfo(float).tiny)) for f in (10, 20, 50, 100, 200)}
    predicted = np.asarray(pulse["predicted_all_asd"])
    reproduction_error = {str(f): float(predicted[int(float(f) / 5)] / max(measured[int(float(f) / 5)], np.finfo(float).tiny) - 1.0) for f in (10, 20, 50, 100, 200, 500, 1000, 3000, 5000, 7000, 10000)}
    counts = classification["counts_by_channel"]["CH0"]
    contaminated_fraction = counts["definitely_pulse_contaminated"] / max(sum(counts.values()), 1)
    result = {"stage": "pulse_contamination_final", "primary_channel": "CH0", "pc_classification": "PC3", "pc_status": "conditional_on_fixed-FPR_detector; CH0/CH1 paired coherence remains inconclusive", "strict_target_conclusion": "C — exact target physical case remains unidentified", "evidence": {"CH0_counts": counts, "CH0_contaminated_fraction": contaminated_fraction, "all_over_clean_ASD_anchors": ratios, "pulse_power_fraction_of_all_anchors": pulse_fraction, "predicted_over_measured_all_minus_one": reproduction_error, "threshold_sensitivity": {channel: row["classification"] for channel, row in thresholds["channels"].items()}, "false_positive_forward_rates": {channel: row["forward"]["false_positive_rate"] for channel, row in false_positive["controls"].items()}, "independent_acceptance_counts": classification["accepted_counts"], "paired_coherence_record_counts": {name: (row["record_count"] if row else 0) for name, row in coherence["subsets"].items()}, "pulse_free_simulation_envelope_coverage": simulation["envelope_coverage_count"]}, "interpretation": "Detected pulses do not materially change the CH0 10-100 Hz ASD; pulse-free records retain the excess. The additive pulse-only reconstruction is numerically close to all-record ASD because its contribution is small, not because it explains the excess.", "supersedes": "The earlier phase_next R5 provisional label; it was based on the pre-partition stationarity audit and is not retained as the primary classification."}
    (case / "pulse_contamination_psd_prediction.json").write_text(json.dumps({**prediction, "pc_classification": "PC3"}, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (case / "full_acquisition_reproduction.json").write_text(json.dumps({**json.loads((case / "full_acquisition_reproduction.json").read_text(encoding="utf-8")), "pc_classification": "PC3", "relative_error_all_anchors_CH0": reproduction_error}, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (case / "pulse_contamination_summary.md").write_text("\n".join([
        "# Pulse contamination summary", "", "## Pulse contamination rate", "", f"CH0: `{counts['definitely_pulse_contaminated']}` definite, `{counts['ambiguous']}` ambiguous, `{counts['pulse_free']}` pulse-free of 345 accepted records; definite fraction `{contaminated_fraction:.3%}`.", "CH1 has only 1 independently accepted record under the unchanged production acceptance predicate.", "", "## Detection confidence / false positives", "", f"Primary FPR is fixed at 1e-3 from pulse pretrigger null portions. Forward-control false-positive rates: CH0 `{false_positive['controls']['CH0']['forward']['false_positive_rate']:.4g}`, CH1 `{false_positive['controls']['CH1']['forward']['false_positive_rate']:.4g}`. Time-reversed and sign-inverted controls are retained in `pulse_false_positive_audit.json`.", "", "## Pulse-free record count", "", f"CH0 strict clean count: `{counts['pulse_free']}`. No pulse subtraction is used.", "", "## All vs pulse-free ASD", "", f"CH0 ASD all/clean ratios at anchors: `{json.dumps(ratios, sort_keys=True)}`. The low-frequency change is small.", "", "## Contaminated subset ASD", "", f"Reconstructed pulse-only power fractions of all-record power at 10/20/50/100/200 Hz: `{json.dumps(pulse_fraction, sort_keys=True)}`.", "", "## Coherence before/after pulse rejection", "", f"Exact-key paired coherence is inconclusive because independent acceptance leaves `{coherence['subsets']['all']['record_count'] if coherence['subsets']['all'] else 0}` pair(s).", "", "## CH0/CH1 pulse coincidence", "", f"Detected paired noise pulses: `{coincidence['noise_common_accepted_detected_pairs']}`; pulse-dataset amplitude-ratio and lag diagnostics are in `pulse_channel_coincidence.json`.", "", "## Predicted pulse-only PSD", "", "The pulse-only reconstruction uses detected event lag, template class, and time-domain amplitude only. No spectral residual scaling is fitted.", "", "## Stationary simulation vs pulse-free experiment", "", f"Pulse-free CH0 vs pulse-gated model covers `{simulation['envelope_coverage_count']}/11` anchors; the 10–100 Hz excess remains in the clean subset.", "", "## Stationary + pulse contamination vs original experiment", "", f"Additive predicted/all relative errors at anchors: `{json.dumps(reproduction_error, sort_keys=True)}`. Agreement is a consequence of the small pulse contribution.", "", "## Remaining residual", "", "The remaining 10–100 Hz excess is not explained by detected pulse contamination. CH0/CH1 coherence classification remains inconclusive because of the strict independent CH1 acceptance count.", "", "## Final PC classification", "", "**PC3 — pulse contamination negligible** for the observed CH0 low-frequency excess, conditional on the fixed-FPR detector. The previous R5 provisional label is superseded.", "", "Strict target conclusion: **C — exact target physical case remains unidentified**."])+"\n", encoding="utf-8")
    phase = {"stage": "phase_next_final_classification", "classification": "PC3", "classification_status": "conditional_on_fixed_FPR_detector; paired_channel_coherence_inconclusive", "interpretation": result["interpretation"], "strict_target_conclusion": result["strict_target_conclusion"], "previous_provisional_classification": "R5 provisional — superseded by pulse partitioning", "pc_classification": "PC3", "evidence": result["evidence"], "prohibitions_respected": ["no noise residual fitting", "no parameter optimizer", "no extra noise source added", "no hanging state added", "no target physical case promoted from proxy"]}
    (case / "phase_next_final_classification.json").write_text(json.dumps(phase, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (case / "phase_next_final_classification.md").write_text("# TES noise mismatch — final classification\n\n**PC3 — pulse contamination negligible** for the CH0 10–100 Hz excess, conditional on the fixed-FPR detector. The previous R5 provisional label is superseded.\n\nStrict target conclusion: **C — exact target physical case remains unidentified**.\n", encoding="utf-8")


if __name__ == "__main__":
    main()

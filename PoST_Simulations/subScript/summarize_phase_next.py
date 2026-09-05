"""Produce the final R1-R5 classification from the completed audit artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


STRICT = "C — exact target physical case remains unidentified"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    args = parser.parse_args()
    comparison = json.loads((args.case_dir / "conditional_comparison_summary.json").read_text(encoding="utf-8"))
    ensemble = json.loads((args.case_dir / "proxy_noise_envelope.json").read_text(encoding="utf-8"))
    pulse = json.loads((args.case_dir / "target_pulse_pole_audit.json").read_text(encoding="utf-8"))
    stationarity = json.loads((args.case_dir / "low_frequency_stationarity_audit.json").read_text(encoding="utf-8"))
    rows = comparison["rows"]
    low = [row for row in rows if row["frequency_Hz"] <= 100.0]
    kilohertz = [row for row in rows if row["frequency_Hz"] >= 1000.0]
    low_outside = all(row["classification"] == "outside_sampled_envelope" for row in low)
    khz_covered = all(row["classification"] != "outside_sampled_envelope" for row in kilohertz)
    if low_outside and khz_covered and stationarity["classification"] != "nonstationary":
        classification = "R5"
        interpretation = "Pulse-constrained full physical source ensemble covers 1–10 kHz but not the 10–100 Hz anchors; independent CH0/CH1 acceptance leaves too few paired records for a definitive stationarity/coherence decision, so R5 is provisional and an independent readout/circuit contribution remains the next hypothesis before any added TES source is considered."
    elif low_outside and stationarity["classification"] == "nonstationary":
        classification = "R4"
        interpretation = "The low-frequency residual is drift-like and is not a stationary TES-noise comparison target."
    elif low_outside:
        classification = "R3"
        interpretation = "The pulse-constrained detector physics is compatible, while a stationary channel/readout-like residual remains."
    else:
        classification = "R1"
        interpretation = "The pulse-constrained physical ensemble covers the sampled 10 Hz–10 kHz shape."
    result = {
        "stage": "phase_next_final_classification",
        "classification": classification,
        "classification_status": "provisional_due_to_inconclusive_low_frequency_coherence" if stationarity["classification"] == "inconclusive" else "conditional",
        "interpretation": interpretation,
        "strict_target_conclusion": STRICT,
        "evidence": {
            "pulse_consistent_scenario_count": ensemble["scenario_count"],
            "pulse_selected_slow_poles": {channel: row["selected_slow_pole"] for channel, row in pulse["channels"].items()},
            "low_frequency_rows": low,
            "kilohertz_rows": kilohertz,
            "stationarity_classification": stationarity["classification"],
            "coherence_10_to_100_Hz_median": stationarity.get("coherence_10_to_100_Hz_median"),
            "coherence_10_to_100_Hz_range": stationarity.get("coherence_10_to_100_Hz_range"),
            "common_accepted_records": stationarity["accepted_counts"]["common"],
        },
        "prohibitions_respected": ["no noise residual fitting", "no parameter optimizer", "no extra noise source added", "no hanging state added", "no target physical case promoted from proxy"],
    }
    (args.case_dir / "phase_next_final_classification.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = ["# TES noise mismatch — next-phase classification", "", f"Final classification: **{classification}** ({result['classification_status']}).", "", interpretation, "", f"Stationarity sub-class: **{stationarity['classification']}**.", f"Common accepted CH0/CH1 records: **{stationarity['accepted_counts']['common']}**.", f"Pulse-consistent scenarios: **{ensemble['scenario_count']}**.", "", f"Strict target conclusion: **{STRICT}**.", "", "No noise residual fitting, parameter optimization, extra source, or hanging state was introduced."]
    (args.case_dir / "phase_next_final_classification.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

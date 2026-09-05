"""Classify accepted noise records using fixed pulse-only template thresholds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Analyze_Experimental_Data.tes_analysis.noise_utils import accepted_noise_indices
from pulse_contamination_common import noise_paths, read_record, record_key, scan_templates


RATE_HZ = 500000.0
SAMPLES = 100000
CUTOFF_HZ = 10000.0


def _range_ok(values: np.ndarray) -> bool:
    return float(np.max(values) - np.min(values)) <= 0.04


def _accepted(paths: list[Path]) -> set[str]:
    def records():
        for path in paths:
            try:
                yield read_record(path)
            except ValueError:
                yield np.empty(0, dtype=float)
    indices = accepted_noise_indices(records(), SAMPLES, RATE_HZ, cutoff=CUTOFF_HZ, remove_mean=True, accept_raw=_range_ok, accept_processed=_range_ok)
    return {record_key(paths[index]) for index in indices}


def _threshold_for(row: dict, template: str, fpr: float) -> float:
    return float(row["thresholds_by_template"][template][f"fpr_{fpr:g}"])


def _classify(best: dict, thresholds: dict) -> str:
    template = best["template"]
    if template is None:
        return "pulse_free"
    score = best["score"]
    if score >= _threshold_for(thresholds, template, 1e-4):
        return "definitely_pulse_contaminated"
    if score >= _threshold_for(thresholds, template, 1e-3):
        return "likely_pulse_contaminated"
    if score >= _threshold_for(thresholds, template, 1e-2):
        return "ambiguous"
    return "pulse_free"


def _scan_channel(target_root: Path, channel: str, accepted: set[str], library: dict, thresholds: dict) -> dict[str, dict]:
    result = {}
    for path in noise_paths(target_root, channel):
        key = record_key(path)
        if key not in accepted:
            continue
        try:
            values = read_record(path)
        except ValueError:
            continue
        best = scan_templates(values - np.mean(values), library["templates"])
        classification = _classify(best, thresholds)
        morphology = "full_pulse" if best["template"] == "full_pulse" else ("tail_only" if best["template"] in {"post_peak_tail", "slow_tail"} else "possible_pulse")
        result[key] = {"event_key": key, "channel": channel, "classification": classification, "morphology": morphology if classification != "pulse_free" else "clean", "best_template": best["template"], "best_lag_samples": best["lag_samples"], "best_lag_s": best["lag_s"], "matched_score": best["score"], "estimated_amplitude_normalized": best["amplitude_normalized"]}
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    library = json.loads(args.library.read_text(encoding="utf-8"))
    threshold_file = json.loads(args.thresholds.read_text(encoding="utf-8"))
    paths = {channel: noise_paths(args.target_root, channel) for channel in ("CH0", "CH1")}
    accepted = {channel: _accepted(paths[channel]) for channel in paths}
    scans = {channel: _scan_channel(args.target_root, channel, accepted[channel], library["channels"][channel], threshold_file["channels"][channel]) for channel in paths}
    keys = sorted(set().union(*[set(paths[channel] and accepted[channel]) for channel in paths]), key=lambda value: int(value))
    records = {}
    for key in keys:
        rows = {channel: scans[channel].get(key, {"event_key": key, "channel": channel, "classification": "not_accepted", "morphology": "not_accepted"}) for channel in paths}
        records[key] = {"event_key": key, "CH0": rows["CH0"], "CH1": rows["CH1"], "paired_classification": "paired_pulse" if rows["CH0"]["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"} and rows["CH1"]["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"} else ("CH0_only_pulse" if rows["CH0"]["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"} else ("CH1_only_pulse" if rows["CH1"]["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"} else "neither"))}
    def counts(channel):
        values = [row[channel]["classification"] for row in records.values() if row[channel]["classification"] != "not_accepted"]
        return {name: values.count(name) for name in ("pulse_free", "ambiguous", "likely_pulse_contaminated", "definitely_pulse_contaminated")}
    result = {"stage": "noise_record_pulse_classification", "target_root": args.target_root.as_posix(), "record_duration_s": SAMPLES / RATE_HZ, "threshold_source": args.thresholds.as_posix(), "threshold_selected_before_noise_psd": True, "accepted_event_key_sets": {channel: sorted(accepted[channel], key=int) for channel in accepted}, "accepted_counts": {channel: len(accepted[channel]) for channel in accepted}, "records": records, "counts_by_channel": {channel: counts(channel) for channel in paths}, "pair_categories": {name: sum(row["paired_classification"] == name for row in records.values()) for name in ("CH0_only_pulse", "CH1_only_pulse", "paired_pulse", "neither")}, "contamination_rate_basis": "time-domain detections per accepted-record duration; no spectrum or residual used", "provenance": "classification uses exact event keys and independent production acceptance for CH0 and CH1"}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "noise_record_pulse_classification.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

"""Compare CH0/CH1 pulse arrival and amplitude ratios by exact event key."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pulse_contamination_common import noise_paths, pulse_paths, read_record, record_key, normalized_pulse


def _quantiles(values):
    return {name: float(np.quantile(values, q)) for name, q in (("q05", 0.05), ("q50", 0.5), ("q95", 0.95))} if values else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--classification", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    classified = json.loads(args.classification.read_text(encoding="utf-8"))
    pulse_amplitudes = {channel: {} for channel in ("CH0", "CH1")}
    for channel in pulse_amplitudes:
        for path in pulse_paths(args.target_root, channel):
            try:
                pulse_amplitudes[channel][record_key(path)] = normalized_pulse(read_record(path))["amplitude_raw"]
            except ValueError:
                continue
    pulse_pair_keys = sorted(set(pulse_amplitudes["CH0"]) & set(pulse_amplitudes["CH1"]), key=int)
    ratios = np.asarray([pulse_amplitudes["CH0"][key] / max(abs(pulse_amplitudes["CH1"][key]), np.finfo(float).tiny) for key in pulse_pair_keys])
    pair_lags, pair_ratios = [], []
    for key, row in classified["records"].items():
        ch0, ch1 = row["CH0"], row["CH1"]
        accepted = ch0["classification"] != "not_accepted" and ch1["classification"] != "not_accepted"
        detected = ch0["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"} and ch1["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"}
        if accepted and detected:
            pair_lags.append(float(ch0["best_lag_s"] - ch1["best_lag_s"]))
            pair_ratios.append(float(ch0["estimated_amplitude_normalized"] / max(abs(ch1["estimated_amplitude_normalized"]), np.finfo(float).tiny)))
    result = {"stage": "pulse_channel_coincidence", "pairing": "exact event key; array position is not used", "pulse_dataset_amplitude_ratio_CH0_over_CH1": _quantiles(ratios.tolist()), "noise_common_accepted_detected_pairs": len(pair_lags), "noise_detection_lag_delta_t_CH0_minus_CH1_s": _quantiles(pair_lags), "noise_detection_amplitude_ratio_CH0_over_CH1": _quantiles(pair_ratios), "interpretation": "Coincidence confidence is inconclusive when the independently accepted paired sample is too small; no noise PSD is used."}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pulse_channel_coincidence.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

"""Predict pulse-contamination PSD from detected event templates and amplitudes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Analyze_Experimental_Data.tes_analysis.noise_utils import estimate_one_sided_asd
from pulse_contamination_common import noise_paths, read_record, record_key


RATE_HZ = 500000.0
SAMPLES = 100000
ANCHORS = (10, 20, 50, 100, 200, 500, 1000, 3000, 5000, 7000, 10000)


def _pulse_record(row: dict, template_library: dict) -> np.ndarray:
    output = np.zeros(SAMPLES, dtype=float)
    if row["classification"] not in {"definitely_pulse_contaminated", "likely_pulse_contaminated"}:
        return output
    name = row["best_template"]
    if not name:
        return output
    template = np.asarray(template_library["templates"][name]["values"], dtype=float)
    step = int(template_library["templates"][name].get("sample_step", 1))
    lag = int(row["best_lag_samples"])
    amplitude = float(row["estimated_amplitude_normalized"])
    if step == 1:
        left, right = max(0, lag), min(SAMPLES, lag + len(template))
        if right > left:
            output[left:right] = amplitude * template[left - lag:right - lag]
    else:
        indices = lag + step * np.arange(len(template))
        valid = (indices >= 0) & (indices < SAMPLES)
        output[indices[valid].astype(int)] = amplitude * template[valid]
    return output


def _asd(records):
    return estimate_one_sided_asd(records, SAMPLES, RATE_HZ, cutoff=0.0, remove_mean=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--classification", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--partitioned-spectra", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    classified = json.loads(args.classification.read_text(encoding="utf-8"))
    library = json.loads(args.library.read_text(encoding="utf-8"))
    partitioned = json.loads(args.partitioned_spectra.read_text(encoding="utf-8"))
    result = {"stage": "pulse_contamination_psd_prediction", "method": "method_A_detected_event_templates_and_measured_amplitudes", "fixed_seed": 20260908, "seed_role": "deterministic event ordering provenance; method A has no random sampling", "stationary_residual_used": False, "amplitude_fit_to_spectrum": False, "channels": {}}
    for channel in ("CH0", "CH1"):
        rows = classified["records"]
        keys = [key for key, row in rows.items() if row[channel]["classification"] != "not_accepted"]
        paths = {record_key(path): path for path in noise_paths(args.target_root, channel)}
        channel_rows = [rows[key][channel] for key in keys]
        pulse_record_count = sum(row["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"} for row in channel_rows)
        pulse_asd, pulse_accepted = _asd((_pulse_record(row, library["channels"][channel]) for row in channel_rows))
        frequencies = np.arange(pulse_asd.size) * RATE_HZ / SAMPLES
        measured = partitioned["subsets"][channel]["all_accepted"]["asd"]
        clean = partitioned["subsets"][channel]["pulse_free"]["asd"]
        if measured is None or clean is None:
            predicted = None
        else:
            predicted = np.sqrt(np.asarray(clean, dtype=float) ** 2 + pulse_asd**2)
        result["channels"][channel] = {"accepted_records": len(keys), "records_with_detected_pulses": pulse_record_count, "pulse_only_estimator_accepted_count": pulse_accepted, "frequencies_Hz": frequencies.tolist(), "pulse_only_asd": pulse_asd.tolist(), "predicted_all_asd": predicted.tolist() if predicted is not None else None, "measured_all_asd": measured, "clean_asd": clean, "anchors": {"pulse_only_asd": {str(anchor): float(pulse_asd[int(np.argmin(np.abs(frequencies - anchor)))]) for anchor in ANCHORS}, "predicted_all_asd": ({str(anchor): float(predicted[int(np.argmin(np.abs(frequencies - anchor)))]) for anchor in ANCHORS} if predicted is not None else None), "measured_all_asd": ({str(anchor): float(np.asarray(measured)[int(np.argmin(np.abs(frequencies - anchor)))]) for anchor in ANCHORS} if measured is not None else None)}}
    result["interpretation"] = "Predicted all-record ASD is sqrt(clean PSD + independently reconstructed pulse-only PSD); no residual or spectral amplitude scaling was fitted."
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pulse_contamination_psd_prediction.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    reproduction = {"stage": "full_acquisition_reproduction", "model": "pulse-free experimental ASD plus measured pulse-contamination PSD", "channels": {channel: {"predicted_all_asd": row["anchors"]["predicted_all_asd"], "measured_all_asd": row["anchors"]["measured_all_asd"]} for channel, row in result["channels"].items()}, "strict_target_conclusion": "C — exact target physical case remains unidentified"}
    (args.output_dir / "full_acquisition_reproduction.json").write_text(json.dumps(reproduction, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

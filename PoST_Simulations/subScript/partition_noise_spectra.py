"""Compute all/clean/contaminated/ambiguous ASD with one canonical estimator."""

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
BANDS = {"5-30_Hz": (5, 30), "30-150_Hz": (30, 150), "150-500_Hz": (150, 500), "0.5-1_kHz": (500, 1000), "1-10_kHz": (1000, 10000)}


def _asd(paths: list[Path], keys: set[str]) -> tuple[np.ndarray, np.ndarray, int]:
    selected = [path for path in paths if record_key(path) in keys]
    values, accepted = estimate_one_sided_asd((read_record(path) for path in selected), SAMPLES, RATE_HZ, cutoff=0.0, remove_mean=True)
    return np.arange(values.size) * RATE_HZ / SAMPLES, values, accepted


def _points(frequencies: np.ndarray, values: np.ndarray) -> dict[str, float]:
    return {str(anchor): float(values[int(np.argmin(np.abs(frequencies - anchor)))]) for anchor in ANCHORS}


def _band_summary(frequencies: np.ndarray, values: np.ndarray) -> dict:
    return {name: {"frequency_range_Hz": [left, right], "median_asd": float(np.median(values[(frequencies >= left) & (frequencies <= right)])), "mean_asd": float(np.mean(values[(frequencies >= left) & (frequencies <= right)]))} for name, (left, right) in BANDS.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--classification", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    classification = json.loads(args.classification.read_text(encoding="utf-8"))
    result = {"stage": "pulse_partitioned_noise_spectra", "target_root": args.target_root.as_posix(), "estimator": "per-record mean removal, Hann window, power average, one-sided ASD; no inverse filter", "rate_Hz": RATE_HZ, "samples": SAMPLES, "record_duration_s": SAMPLES / RATE_HZ, "anchors_Hz": list(ANCHORS), "subsets": {}}
    for channel in ("CH0", "CH1"):
        rows = classification["records"]
        keys = {key for key, row in rows.items() if row[channel]["classification"] != "not_accepted"}
        clean = {key for key in keys if rows[key][channel]["classification"] == "pulse_free"}
        contaminated = {key for key in keys if rows[key][channel]["classification"] == "definitely_pulse_contaminated"}
        ambiguous = {key for key in keys if rows[key][channel]["classification"] == "ambiguous"}
        subsets = {"all_accepted": keys, "pulse_free": clean, "definitely_contaminated": contaminated, "ambiguous": ambiguous}
        output = {}
        paths = noise_paths(args.target_root, channel)
        for name, subset in subsets.items():
            if not subset:
                output[name] = {"record_count": 0, "frequencies_Hz": None, "asd": None, "anchors_asd": None, "bands": None}
                continue
            frequencies, asd, accepted = _asd(paths, subset)
            output[name] = {"record_count": len(subset), "estimator_accepted_count": accepted, "frequencies_Hz": frequencies.tolist(), "asd": asd.tolist(), "anchors_asd": _points(frequencies, asd), "bands": _band_summary(frequencies, asd)}
        all_asd, clean_asd = output["all_accepted"]["asd"], output["pulse_free"]["asd"]
        if all_asd is not None and clean_asd is not None:
            all_array, clean_array = np.asarray(all_asd), np.asarray(clean_asd)
            output["ratios"] = {"ASD_all_over_clean": (all_array / np.maximum(clean_array, np.finfo(float).tiny)).tolist(), "ASD_contaminated_over_clean": ((np.asarray(output["definitely_contaminated"]["asd"]) / np.maximum(clean_array, np.finfo(float).tiny)).tolist() if output["definitely_contaminated"]["asd"] is not None else None), "PSD_all_over_clean": ((all_array**2 / np.maximum(clean_array**2, np.finfo(float).tiny)).tolist())}
            output["ratio_anchors"] = {name: float((all_array / np.maximum(clean_array, np.finfo(float).tiny))[int(np.argmin(np.abs(frequencies - float(name))))]) for name in map(str, ANCHORS)}
        result["subsets"][channel] = output
    result["provenance"] = "All subsets use the same canonical estimator; primary clean is the strict pulse_free classification and no pulse subtraction is applied."
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pulse_partitioned_noise_spectra.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

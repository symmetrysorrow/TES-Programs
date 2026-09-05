"""Frequency-resolved CH0/CH1 coherence for all pulse partitions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pulse_contamination_common import noise_paths, read_record, record_key


RATE_HZ = 500000.0
SAMPLES = 100000
SUBSETS = ("all", "pulse_free", "pulse_contaminated", "ambiguous")
POINTS = (10, 20, 50, 100)


def _calculate(pairs: list[tuple[Path, Path]]) -> dict | None:
    if not pairs:
        return None
    window = np.hanning(SAMPLES)
    p00 = p11 = p01 = None
    for path0, path1 in pairs:
        x0, x1 = read_record(path0), read_record(path1)
        f0 = np.fft.rfft((x0 - np.mean(x0)) * window)
        f1 = np.fft.rfft((x1 - np.mean(x1)) * window)
        if p00 is None:
            p00, p11, p01 = np.zeros_like(np.abs(f0) ** 2), np.zeros_like(np.abs(f1) ** 2), np.zeros_like(f0 * np.conjugate(f1))
        p00 += np.abs(f0) ** 2
        p11 += np.abs(f1) ** 2
        p01 += f0 * np.conjugate(f1)
    p00 /= len(pairs)
    p11 /= len(pairs)
    p01 /= len(pairs)
    frequency = np.arange(p00.size) * RATE_HZ / SAMPLES
    coherence = np.abs(p01) ** 2 / np.maximum(p00 * p11, np.finfo(float).tiny)
    phase = np.angle(p01)
    return {"record_count": len(pairs), "frequencies_Hz": frequency.tolist(), "coherence": coherence.tolist(), "cross_phase_rad": phase.tolist(), "S00": p00.tolist(), "S11": p11.tolist(), "S01_real": np.real(p01).tolist(), "S01_imag": np.imag(p01).tolist(), "points": {str(point): {"coherence": float(coherence[int(point / 5)]), "cross_phase_rad": float(phase[int(point / 5)])} for point in POINTS}}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--classification", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    classified = json.loads(args.classification.read_text(encoding="utf-8"))
    paths0 = {record_key(path): path for path in noise_paths(args.target_root, "CH0")}
    paths1 = {record_key(path): path for path in noise_paths(args.target_root, "CH1")}
    rows = classified["records"]
    common = [key for key in rows if key in paths0 and key in paths1 and rows[key]["CH0"]["classification"] != "not_accepted" and rows[key]["CH1"]["classification"] != "not_accepted"]
    keys = {"all": common, "pulse_free": [key for key in common if rows[key]["CH0"]["classification"] == rows[key]["CH1"]["classification"] == "pulse_free"], "pulse_contaminated": [key for key in common if rows[key]["CH0"]["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"} or rows[key]["CH1"]["classification"] in {"definitely_pulse_contaminated", "likely_pulse_contaminated"}], "ambiguous": [key for key in common if rows[key]["CH0"]["classification"] == "ambiguous" or rows[key]["CH1"]["classification"] == "ambiguous"]}
    result = {"stage": "partitioned_coherence", "estimator": "per-record mean removal, Hann, rFFT, cross-power average", "frequency_range_Hz": [5, 500], "pairing": "exact event key", "subsets": {name: _calculate([(paths0[key], paths1[key]) for key in values]) for name, values in keys.items()}}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "partitioned_coherence.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

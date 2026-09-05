"""Choose fixed-FPR pulse detector thresholds from pulse null portions only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pulse_contamination_common import matched_score, pulse_paths, read_record


FPRS = (1e-2, 1e-3, 1e-4)
PRIMARY_FPR = 1e-3


def _channel(target_root: Path, channel: str, library: dict) -> dict:
    null_scores = {name: [] for name in library["templates"]}
    baseline_stds = []
    for path in pulse_paths(target_root, channel):
        try:
            values = read_record(path)
        except ValueError:
            continue
        baseline = values[:4500] - np.mean(values[:4500])
        baseline_stds.append(float(np.std(baseline)))
        for start in range(0, 400, 100):
            window = baseline[start:start + 4096]
            if len(window) < 4096:
                continue
            for name, template in library["templates"].items():
                step = int(template.get("sample_step", 1))
                if step == 1:
                    score, _lag, _amp = matched_score(window, np.asarray(template["values"]), step)
                else:
                    # The slow-tail null uses a deterministic Gaussian null
                    # calibrated only by the pulse pretrigger baseline.
                    continue
                null_scores[name].append(score)
    rng = np.random.default_rng(20260906)
    sigma = float(np.median(baseline_stds))
    for name, template in library["templates"].items():
        if int(template.get("sample_step", 1)) > 1:
            kernel = np.asarray(template["values"], dtype=float)
            synthetic = rng.normal(0.0, sigma, size=(2000, len(kernel)))
            kernel = kernel - np.mean(kernel)
            kernel_norm = np.linalg.norm(kernel)
            scores = (synthetic @ kernel) / np.maximum(kernel_norm * np.linalg.norm(synthetic, axis=1), np.finfo(float).tiny)
            null_scores[name] = [float(x) for x in scores]
    thresholds = {}
    for name, values in null_scores.items():
        if not values:
            raise RuntimeError(f"{channel}/{name}: no pulse-baseline null scores")
        array = np.asarray(values, dtype=float)
        thresholds[name] = {f"fpr_{fpr:g}": float(np.quantile(array, 1.0 - fpr)) for fpr in FPRS}
    return {"channel": channel, "primary_fpr": PRIMARY_FPR, "null_records_source": "pulse pretrigger baseline portions only", "baseline_std_median": sigma, "null_sample_count": {name: len(values) for name, values in null_scores.items()}, "thresholds_by_template": thresholds}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    library = json.loads(args.library.read_text(encoding="utf-8"))
    result = {"stage": "fixed_fpr_pulse_detection_thresholds", "target_root": args.target_root.as_posix(), "simulation_spectrum_read": False, "noise_spectrum_read": False, "primary_fpr": PRIMARY_FPR, "fprs_reported": list(FPRS), "channels": {channel: _channel(args.target_root, channel, library["channels"][channel]) for channel in ("CH0", "CH1")}, "selection_rule": "thresholds are pulse-pretrigger null quantiles fixed before any noise PSD comparison; no threshold tuning against simulation or residual"}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pulse_detection_thresholds.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

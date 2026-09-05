"""Build pulse templates from target pulse records only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pulse_contamination_common import normalized_pulse, pulse_paths, sha256


TARGET_ROOT = Path(r"G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2")


def _cross(values: np.ndarray, peak: int, fraction: float, rising: bool) -> int | None:
    if rising:
        found = np.flatnonzero(values[:peak] >= fraction * values[peak])
        return int(found[0]) if len(found) else None
    found = np.flatnonzero(values[peak:] <= fraction * values[peak])
    return int(peak + found[0]) if len(found) else None


def _channel(target_root: Path, channel: str) -> dict:
    rows = []
    files = pulse_paths(target_root, channel)
    for path in files:
        try:
            raw = np.fromfile(path, dtype=np.float64, offset=4)
            if raw.size != 100000 or not np.all(np.isfinite(raw)):
                continue
            row = normalized_pulse(raw)
            oriented = row["full_template"]
            peak = 256
            rise20 = _cross(oriented, peak, 0.2, True)
            rise90 = _cross(oriented, peak, 0.9, True)
            decay90 = _cross(oriented, peak, 0.9, False)
            decay10 = _cross(oriented, peak, 0.1, False)
            row["rise_time_s"] = (rise90 - rise20) / 500000.0 if rise20 is not None and rise90 is not None else None
            row["decay_time_s"] = (decay10 - decay90) / 500000.0 if decay90 is not None and decay10 is not None else None
            row["tail_level"] = float(np.mean(row["tail_template"][-500:]))
            row["source_file"] = path.as_posix()
            rows.append(row)
        except ValueError:
            continue
    if len(rows) < 10:
        raise RuntimeError(f"{channel}: too few pulse records for a template library")
    def median_array(key):
        return np.median(np.asarray([row[key] for row in rows]), axis=0)
    return {
        "channel": channel,
        "records_found": len(files),
        "records_used": len(rows),
        "normalization": "per-record baseline from first 4500 samples, polarity oriented, peak-normalized, aligned at peak",
        "template_length": 4096,
        "sample_rate_Hz": 500000.0,
        "template_source_files": [{"path": row["source_file"], "sha256": sha256(Path(row["source_file"]))} for row in rows],
        "templates": {
            "full_pulse": {"sample_step": 1, "relative_start_samples": -256, "values": median_array("full_template").tolist()},
            "post_peak_tail": {"sample_step": 1, "relative_start_samples": 0, "values": median_array("tail_template").tolist()},
            "slow_tail": {"sample_step": 10, "relative_start_samples": 0, "values": median_array("slow_tail_template").tolist()},
        },
        "rise_time_distribution_s": {key: float(np.quantile([row["rise_time_s"] for row in rows if row["rise_time_s"] is not None], q)) for key, q in (("q05", 0.05), ("q50", 0.5), ("q95", 0.95))},
        "decay_distribution_s": {key: float(np.quantile([row["decay_time_s"] for row in rows if row["decay_time_s"] is not None], q)) for key, q in (("q05", 0.05), ("q50", 0.5), ("q95", 0.95))},
        "tail_distribution": {key: float(np.quantile([row["tail_level"] for row in rows], q)) for key, q in (("q05", 0.05), ("q50", 0.5), ("q95", 0.95))},
        "provenance": "pulse-only; no target noise path, spectrum, residual, or simulation output was read",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target-root", type=Path, default=TARGET_ROOT)
    args = parser.parse_args()
    result = {"stage": "pulse_template_library", "target_root": args.target_root.as_posix(), "noise_records_read": False, "channels": {channel: _channel(args.target_root, channel) for channel in ("CH0", "CH1")}}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pulse_template_library.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

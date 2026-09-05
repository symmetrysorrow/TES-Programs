"""Report fixed-FPR clean-spectrum sensitivity without retuning thresholds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Analyze_Experimental_Data.tes_analysis.noise_utils import estimate_one_sided_asd
from pulse_contamination_common import noise_paths, read_record


RATE_HZ = 500000.0
SAMPLES = 100000
FPRS = (1e-2, 1e-3, 1e-4)
ANCHORS = (10, 20, 50, 100, 200, 500, 1000, 3000, 5000, 7000, 10000)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    classified = json.loads(args.classification.read_text(encoding="utf-8"))
    thresholds = json.loads(args.thresholds.read_text(encoding="utf-8"))
    result = {"stage": "pulse_contamination_threshold_sensitivity", "thresholds_fixed_before_psd": True, "fprs": list(FPRS), "channels": {}}
    for channel in ("CH0", "CH1"):
        paths = noise_paths(args.target_root, channel)
        rows = classified["records"]
        accepted = {key: row[channel] for key, row in rows.items() if row[channel]["classification"] != "not_accepted"}
        output = {}
        for fpr in FPRS:
            keys = {key for key, row in accepted.items() if row["best_template"] is not None and row["matched_score"] < float(thresholds["channels"][channel]["thresholds_by_template"][row["best_template"]][f"fpr_{fpr:g}"])}
            if not keys:
                output[f"fpr_{fpr:g}"] = {"record_count": 0, "anchors_asd": None}
                continue
            selected = [path for path in paths if path.stem.split("_", 1)[1] in keys]
            asd, accepted_count = estimate_one_sided_asd((read_record(path) for path in selected), SAMPLES, RATE_HZ, cutoff=0.0, remove_mean=True)
            frequencies = np.arange(asd.size) * RATE_HZ / SAMPLES
            output[f"fpr_{fpr:g}"] = {"record_count": len(keys), "estimator_accepted_count": accepted_count, "anchors_asd": {str(anchor): float(asd[int(np.argmin(np.abs(frequencies - anchor)))]) for anchor in ANCHORS}}
        result["channels"][channel] = output
    for channel, rows in result["channels"].items():
        values = [row["anchors_asd"]["10"] for row in rows.values() if row["anchors_asd"]]
        result["channels"][channel]["classification"] = "threshold_sensitive" if values and max(values) / min(values) > 1.5 else ("robust" if values else "inconclusive")
    result["interpretation"] = "FPR thresholds are descriptive robustness checks; none was chosen by agreement with a simulation spectrum."
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pulse_contamination_threshold_sensitivity.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

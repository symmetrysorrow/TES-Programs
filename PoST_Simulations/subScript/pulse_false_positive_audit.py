"""False-positive controls for the fixed-FPR pulse detector."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pulse_contamination_common import matched_score, pulse_paths, read_record


def _scores(target_root: Path, channel: str, template: np.ndarray, transform: str) -> list[float]:
    values = []
    rng = np.random.default_rng(20260907)
    for path in pulse_paths(target_root, channel):
        try:
            record = read_record(path)
        except ValueError:
            continue
        baseline = record[:4500] - np.mean(record[:4500])
        for start in rng.integers(0, 405, size=8):
            window = baseline[int(start):int(start) + len(template)]
            if len(window) != len(template):
                continue
            kernel = template[::-1] if transform == "time_reversed" else (-template if transform == "sign_inverted" else template)
            score, _lag, _amp = matched_score(window, kernel)
            values.append(score)
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    library = json.loads(args.library.read_text(encoding="utf-8"))
    thresholds = json.loads(args.thresholds.read_text(encoding="utf-8"))
    controls = {}
    for channel in ("CH0", "CH1"):
        template = np.asarray(library["channels"][channel]["templates"]["full_pulse"]["values"], dtype=float)
        threshold = float(thresholds["channels"][channel]["thresholds_by_template"]["full_pulse"]["fpr_0.001"])
        channel_controls = {}
        for name in ("forward", "time_reversed", "sign_inverted"):
            transform = "time_reversed" if name == "time_reversed" else ("sign_inverted" if name == "sign_inverted" else "forward")
            values = _scores(args.target_root, channel, template, transform)
            channel_controls[name] = {"sample_count": len(values), "threshold_primary": threshold, "false_positive_count": int(np.count_nonzero(np.asarray(values) >= threshold)), "false_positive_rate": float(np.mean(np.asarray(values) >= threshold)) if values else 0.0, "score_q99": float(np.quantile(values, 0.99)) if values else None}
        controls[channel] = channel_controls
    result = {"stage": "pulse_false_positive_audit", "null_source": "pulse pretrigger baseline portions", "simulation_spectrum_read": False, "controls": controls, "controls_included": ["pulse pretrigger baseline", "time-reversed template", "sign-inverted template", "random baseline lag windows"], "interpretation": "The primary threshold remains the fixed 1e-3 pulse-null quantile; these controls do not tune it against any spectrum."}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pulse_false_positive_audit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

"""Compare pulse-free experimental ASD to the frozen pulse-gated physics ensemble."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from proxy_physics import noise_components, operating_point


FREQUENCIES = np.array([10, 20, 50, 100, 200, 500, 1000, 3000, 5000, 7000, 10000], dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-a-dir", type=Path, required=True)
    parser.add_argument("--partitioned-spectra", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    scenarios = json.loads((args.stage_a_dir / "proxy_scenarios.json").read_text(encoding="utf-8"))["pulse_consistent_scenarios"]
    partitioned = json.loads(args.partitioned_spectra.read_text(encoding="utf-8"))
    measured = partitioned["subsets"]["CH0"]["pulse_free"]["asd"]
    measured_frequency = np.asarray(partitioned["subsets"]["CH0"]["pulse_free"]["frequencies_Hz"])
    measured = np.asarray(measured)
    measured_anchor = np.asarray([measured[int(np.argmin(np.abs(measured_frequency - f)))] for f in FREQUENCIES])
    measured_normalized = measured_anchor / measured_anchor[6]
    values = []
    for scenario in scenarios:
        if not operating_point(scenario["parameters"])["stable"]:
            continue
        _components, meta = noise_components(scenario["parameters"], FREQUENCIES)
        values.append(meta["total_asd"] / meta["total_asd"][6])
    values = np.asarray(values)
    median = np.median(values, axis=0)
    low, high = np.min(values, axis=0), np.max(values, axis=0)
    rows = [{"frequency_Hz": float(f), "pulse_free_experimental_normalized": float(v), "proxy_min": float(lo), "proxy_median": float(mid), "proxy_max": float(hi), "inside_envelope": bool(lo <= v <= hi)} for f, v, lo, mid, hi in zip(FREQUENCIES, measured_normalized, low, median, high)]
    per_scenario_error = np.max(np.abs(np.log(values / measured_normalized[None, :])), axis=1)
    single_scenario_within_10pct = bool(
        np.any(np.all(np.abs(values / measured_normalized[None, :] - 1.0) <= 0.10, axis=1))
    )
    result = {"stage": "pulse_free_simulation_comparison", "experimental_subset": "CH0 pulse_free strict clean", "scenario_selection": "frozen pulse_consistent_scenarios only", "parameter_generation_called": False, "frequencies_Hz": FREQUENCIES.tolist(), "rows": rows, "envelope_coverage_count": int(sum(row["inside_envelope"] for row in rows)), "single_scenario_within_10pct_all_anchor": single_scenario_within_10pct, "best_member_not_promoted": True, "best_member_max_abs_log_ratio_descriptive": float(np.min(per_scenario_error)), "strict_target_conclusion": "C — exact target physical case remains unidentified"}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pulse_free_simulation_comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

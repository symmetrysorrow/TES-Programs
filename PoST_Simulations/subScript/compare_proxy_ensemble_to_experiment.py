"""Stage B comparison: compare a completed frozen ensemble to stored spectra.

No parameter construction happens here.  The only model input is the completed
``proxy_noise_envelope.json``; range changes require a new Stage-A run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import general


FREQUENCIES = np.array([10.0, 100.0, 1000.0, 3000.0, 5000.0, 7000.0, 10000.0])
BANDS = {"10-100_Hz": (10.0, 100.0), "100-1000_Hz": (100.0, 1000.0), "1-3_kHz": (1000.0, 3000.0), "3-10_kHz": (3000.0, 10000.0)}


def _dump(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _classification(value: float, low: float, q05: float, q95: float, high: float) -> str:
    if q05 <= value <= q95:
        return "inside_sampled_q05_q95"
    if low <= value <= high:
        return "inside_sampled_min_max"
    return "outside_sampled_envelope"


def _compare(exp: dict, envelope: dict) -> tuple[dict, dict]:
    points = {float(row["frequency_Hz"]): row for row in exp["points"]}
    model = envelope["normalized_total_asd"]
    rows = []
    for index, frequency in enumerate(FREQUENCIES):
        row = points[float(frequency)]
        value = float(row["experimental_pre_analysis_normalized"])
        low, q05, q50, q95, high = (model[key][index] for key in ("min", "q05", "q50", "q95", "max"))
        rows.append({"frequency_Hz": float(frequency), "experimental_pre_analysis_normalized": value, "proxy_min": low, "proxy_q05": q05, "proxy_median": q50, "proxy_q95": q95, "proxy_max": high, "classification": _classification(value, low, q05, q95, high)})
    bands = {}
    for name, (left, right) in BANDS.items():
        band_rows = [row for row in rows if left <= row["frequency_Hz"] <= right]
        bands[name] = {"frequencies_Hz": [row["frequency_Hz"] for row in band_rows], "outside_count": sum(row["classification"] == "outside_sampled_envelope" for row in band_rows), "status": "partially_or_fully_outside" if any(row["classification"] == "outside_sampled_envelope" for row in band_rows) else "covered_by_sampled_min_max"}
    summary = {"stage": "B_comparison", "comparison_kind": "normalized_shape_only", "experimental_source": "existing comparison_summary.json; pre-analysis values", "proxy_source": "proxy_noise_envelope.json", "parameter_generation_called": False, "range_adjustment_performed": False, "sampled_quantiles_are_not_probabilities": True, "rows": rows, "bands": bands, "interpretation": "Envelope membership is a diagnostic of the frozen ensemble, not a parameter estimate."}
    return summary, points


def _post_factor() -> np.ndarray:
    return general.BesselMagnitudeResponse(FREQUENCIES, 500000.0, 10000.0, passes=2)


def _plot(output_dir: Path, exp: dict, envelope: dict, post: bool = False) -> None:
    import matplotlib.pyplot as plt
    model = envelope["normalized_total_asd"]
    factor = _post_factor() if post else np.ones(len(FREQUENCIES))
    if post:
        factor = factor / factor[2]
    key = "experimental_post_analysis_normalized" if post else "experimental_pre_analysis_normalized"
    label = "post-analysis" if post else "pre-analysis"
    rows = exp["points"]
    x = FREQUENCIES
    plt.figure(figsize=(7, 4))
    plt.fill_between(x, np.asarray(model["min"]) * factor, np.asarray(model["max"]) * factor, alpha=0.25, label="frozen proxy min/max")
    plt.plot(x, np.asarray(model["q50"]) * factor, label="sampled median")
    plt.plot(x, [row[key] for row in rows], "o-", label=f"experiment {label}")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized ASD")
    plt.grid(True, which="both")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_dir / "plots" / ("proxy_noise_envelope_post_analysis.png" if post else "proxy_noise_envelope_pre_analysis.png"), dpi=160)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ensemble", type=Path, required=True)
    parser.add_argument("--experimental-summary", type=Path, required=True)
    parser.add_argument("--stage-a-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "plots").mkdir(parents=True, exist_ok=True)
    envelope = json.loads(args.ensemble.read_text(encoding="utf-8"))
    experimental = json.loads(args.experimental_summary.read_text(encoding="utf-8"))
    stage_a = json.loads(args.stage_a_manifest.read_text(encoding="utf-8"))
    if stage_a.get("freeze_status") != "frozen" or envelope.get("input_stage") != "A_frozen_proxy_scenarios_only":
        raise RuntimeError("comparison requires a completed frozen Stage-A input")
    summary, _ = _compare(experimental, envelope)
    post_factor = _post_factor()
    post_factor = post_factor / post_factor[2]
    post_rows = []
    for index, row in enumerate(experimental["points"]):
        value = float(row["experimental_post_analysis_normalized"])
        low = envelope["normalized_total_asd"]["min"][index] * post_factor[index]
        high = envelope["normalized_total_asd"]["max"][index] * post_factor[index]
        q05 = envelope["normalized_total_asd"]["q05"][index] * post_factor[index]
        q95 = envelope["normalized_total_asd"]["q95"][index] * post_factor[index]
        post_rows.append({"frequency_Hz": float(row["frequency_Hz"]), "experimental_post_analysis_normalized": value, "proxy_post_min": low, "proxy_post_q05": q05, "proxy_post_q95": q95, "proxy_post_max": high, "classification": _classification(value, low, q05, q95, high)})
    summary["ensemble_completed_before_experiment_read"] = True
    summary["post_analysis"] = {"filter": "2nd-order 10 kHz digital Bessel; zero-phase ASD factor |H|^2", "rows": post_rows}
    summary["exploratory_conclusion"] = "P4 — proxy parameter space is still too underconstrained to make a useful target reproduction statement"
    summary["strict_target_conclusion"] = "C — exact target physical case remains unidentified"
    summary["why_not_P1_to_P3"] = "All runnable scenarios require simulation_reference_only values for unresolved thermal/electrical parameters; the frozen ensemble is a sensitivity/reference ensemble, not an independently supported target envelope."
    _dump(args.output_dir / "conditional_comparison_summary.json", summary)
    markdown = ["# Conditional proxy ensemble comparison", "", "Strict target conclusion: **C — exact target physical case remains unidentified**.", "", "Exploratory conclusion: **P4 — proxy parameter space is still too underconstrained to make a useful target reproduction statement**.", "", "The Stage-A range was frozen before this comparison. Sampled q05/q95 are descriptive quantiles, not probabilities.", "", "## Pre-analysis shape", "", "| Hz | experiment | proxy min | proxy q05 | proxy q50 | proxy q95 | proxy max | status |", "|---:|---:|---:|---:|---:|---:|---:|---|"]
    for row in summary["rows"]:
        markdown.append(f"| {row['frequency_Hz']:.0f} | {row['experimental_pre_analysis_normalized']:.6g} | {row['proxy_min']:.6g} | {row['proxy_q05']:.6g} | {row['proxy_median']:.6g} | {row['proxy_q95']:.6g} | {row['proxy_max']:.6g} | {row['classification']} |")
    markdown.extend(["", "## Band diagnostics", ""])
    for name, band in summary["bands"].items():
        markdown.append(f"- **{name}**: {band['status']} ({band['outside_count']} sampled-anchor points outside min/max).")
    markdown.extend(["", "Post-analysis comparison applies the same deterministic 2nd-order 10 kHz Bessel zero-phase ASD factor. It is secondary; detector-side shape interpretation uses pre-analysis.", "", "No parameter range was adjusted and no best-looking member was promoted to a parameter estimate."])
    (args.output_dir / "conditional_comparison_summary.md").write_text("\n".join(markdown) + "\n", encoding="utf-8")
    _plot(args.output_dir, experimental, envelope, post=False)
    _plot(args.output_dir, experimental, envelope, post=True)


if __name__ == "__main__":
    main()

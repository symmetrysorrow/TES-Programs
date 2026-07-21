"""Compare the single-pixel current pulse from COMSOL and Elmer.

The default inputs are the COMSOL probe export in ``docs/Single-Pixel.txt``
and the frozen Elmer 3x-refined pulse series.  The response is defined as the
current decrease from its pre-pulse baseline.  Its 10--90 % rise time uses
linear interpolation, avoiding a dependence on either solver's output grid.

Usage:
    python scripts/analysis/compare_single_pixel_current.py
"""
from __future__ import annotations

import csv
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
COMSOL_PATH = ROOT / "docs" / "Single-Pixel.txt"
ELMER_PATH = ROOT / "artifacts" / "series" / "tes_pulse_20ms_3x_series.csv"
OUTDIR = ROOT / "artifacts" / "comparison" / "single_pixel_current"
PULSE_START_MS = 20.02
BASELINE_WINDOW_MS = (19.5, PULSE_START_MS)


def crossing_time(time_ms: np.ndarray, response_uA: np.ndarray, level_uA: float) -> float:
    """First rising crossing of ``level_uA``, linearly interpolated."""
    indices = np.flatnonzero((response_uA[:-1] < level_uA) & (response_uA[1:] >= level_uA))
    if len(indices) == 0:
        raise ValueError(f"response does not cross {level_uA} uA")
    i = int(indices[0])
    fraction = (level_uA - response_uA[i]) / (response_uA[i + 1] - response_uA[i])
    return float(time_ms[i] + fraction * (time_ms[i + 1] - time_ms[i]))


def metrics(name: str, time_ms: np.ndarray, current_uA: np.ndarray) -> dict[str, float | str]:
    # Include the endpoint: the Elmer export has its last millisecond-scale
    # pre-pulse sample at exactly 20.000 ms.
    mask = (time_ms >= BASELINE_WINDOW_MS[0]) & (time_ms <= BASELINE_WINDOW_MS[1])
    if not np.any(mask):
        raise ValueError(f"{name}: no samples in baseline window {BASELINE_WINDOW_MS} ms")
    baseline = float(np.mean(current_uA[mask]))
    response = baseline - current_uA
    # The pulse response is the largest decrease after the pulse begins.
    post = np.flatnonzero(time_ms >= PULSE_START_MS)
    peak_index = int(post[np.argmax(response[post])])
    amplitude = float(response[peak_index])
    t10 = crossing_time(time_ms[:peak_index + 1], response[:peak_index + 1], 0.1 * amplitude)
    t90 = crossing_time(time_ms[:peak_index + 1], response[:peak_index + 1], 0.9 * amplitude)
    return {
        "model": name,
        "baseline_current_uA": baseline,
        "minimum_current_uA": float(current_uA[peak_index]),
        "peak_current_drop_uA": amplitude,
        "peak_time_ms": float(time_ms[peak_index]),
        "peak_delay_from_pulse_ms": float(time_ms[peak_index] - PULSE_START_MS),
        "t10_ms": t10,
        "t90_ms": t90,
        "rise_time_10_90_ms": t90 - t10,
    }


def relative_error(elmer: float, comsol: float) -> float:
    return 100.0 * (elmer - comsol) / comsol


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--elmer", type=Path, default=ELMER_PATH,
                        help="Elmer series CSV (default: frozen 3x-refined pulse series)")
    parser.add_argument("--label", default="Elmer single-pixel (3x refined)",
                        help="label used for the Elmer trace and metric table")
    args = parser.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    comsol = np.loadtxt(COMSOL_PATH, comments="%", encoding="utf-8")
    elmer_path = args.elmer if args.elmer.is_absolute() else ROOT / args.elmer
    elmer = np.genfromtxt(elmer_path, delimiter=",", names=True)
    t_comsol, i_comsol = comsol[:, 0], comsol[:, 4]
    t_elmer, i_elmer = elmer["time_s"] * 1e3, elmer["tes_current_A"] * 1e6

    comsol_metrics = metrics("COMSOL Single-Pixel", t_comsol, i_comsol)
    elmer_metrics = metrics(args.label, t_elmer, i_elmer)
    error_fields = ("baseline_current_uA", "minimum_current_uA", "peak_current_drop_uA",
                    "peak_delay_from_pulse_ms", "rise_time_10_90_ms")
    errors = {f"{field}_error_percent": relative_error(float(elmer_metrics[field]), float(comsol_metrics[field]))
              for field in error_fields}
    report = {"definition": {
        "pulse_start_ms": PULSE_START_MS,
        "baseline_window_ms": BASELINE_WINDOW_MS,
        "response": "pre-pulse baseline current minus current",
        "rise_time": "10% to 90% of the peak current drop; crossings are linearly interpolated",
    }, "comsol": comsol_metrics, "elmer": elmer_metrics, "elmer_minus_comsol": errors}
    (OUTDIR / "metrics.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    fields = ["model", "baseline_current_uA", "minimum_current_uA", "peak_current_drop_uA",
              "peak_time_ms", "peak_delay_from_pulse_ms", "t10_ms", "t90_ms", "rise_time_10_90_ms"]
    with (OUTDIR / "metrics.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows((comsol_metrics, elmer_metrics))

    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.8), sharex=True, layout="constrained")
    for ax in axes:
        ax.axvline(PULSE_START_MS, color="0.45", ls=":", lw=1, label="pulse start (20.02 ms)")
        ax.grid(color="0.88", lw=0.8)
    axes[0].plot(t_comsol, i_comsol, label="COMSOL Single-Pixel", color="#2468a2", lw=1.8)
    axes[0].plot(t_elmer, i_elmer, label="Elmer single-pixel (3x refined)", color="#d86438", lw=1.6, marker="o", ms=3)
    axes[0].set_ylabel("TES current [µA]")
    axes[0].set_title("Single-pixel TES current: raw time series")
    axes[0].legend(frameon=False, ncol=2, fontsize=8)

    axes[1].plot(t_comsol, comsol_metrics["baseline_current_uA"] - i_comsol,
                 label="COMSOL current drop", color="#2468a2", lw=1.8)
    axes[1].plot(t_elmer, elmer_metrics["baseline_current_uA"] - i_elmer,
                 label="Elmer current drop", color="#d86438", lw=1.6, marker="o", ms=3)
    axes[1].set_xlim(19.8, 22.0)
    axes[1].set_xlabel("time [ms]")
    axes[1].set_ylabel("current drop from baseline [µA]")
    axes[1].set_title("Pulse response normalized to each model's pre-pulse baseline")
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(OUTDIR / "current_timeseries_comparison.png", dpi=200)
    plt.close(fig)

    summary = (
        "# Single-pixel TES current comparison\n\n"
        "Reference: COMSOL `docs/Single-Pixel.txt`; Elmer: "
        f"`{elmer_path.relative_to(ROOT) if elmer_path.is_relative_to(ROOT) else elmer_path}`.\n\n"
        "| Metric | COMSOL | Elmer | Elmer error vs COMSOL |\n|---|---:|---:|---:|\n"
        f"| Baseline current [µA] | {comsol_metrics['baseline_current_uA']:.6f} | {elmer_metrics['baseline_current_uA']:.6f} | {errors['baseline_current_uA_error_percent']:+.2f}% |\n"
        f"| Minimum current [µA] | {comsol_metrics['minimum_current_uA']:.6f} | {elmer_metrics['minimum_current_uA']:.6f} | {errors['minimum_current_uA_error_percent']:+.2f}% |\n"
        f"| Peak current drop [µA] | {comsol_metrics['peak_current_drop_uA']:.6f} | {elmer_metrics['peak_current_drop_uA']:.6f} | {errors['peak_current_drop_uA_error_percent']:+.2f}% |\n"
        f"| Peak delay from pulse [ms] | {comsol_metrics['peak_delay_from_pulse_ms']:.6f} | {elmer_metrics['peak_delay_from_pulse_ms']:.6f} | {errors['peak_delay_from_pulse_ms_error_percent']:+.2f}% |\n"
        f"| 10–90% rise time [ms] | {comsol_metrics['rise_time_10_90_ms']:.6f} | {elmer_metrics['rise_time_10_90_ms']:.6f} | {errors['rise_time_10_90_ms_error_percent']:+.2f}% |\n\n"
        f"Baseline window: {BASELINE_WINDOW_MS[0]:.2f}–{BASELINE_WINDOW_MS[1]:.2f} ms. "
        "The peak is the maximum post-pulse current decrease; 10%/90% crossings are linearly interpolated.\n"
    )
    (OUTDIR / "summary.md").write_text(summary, encoding="utf-8")
    print(f"Wrote {OUTDIR}")


if __name__ == "__main__":
    main()

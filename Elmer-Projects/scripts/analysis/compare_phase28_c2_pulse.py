"""Compare the Phase23 rectangular pulse, Phase28 C2 pulse, and COMSOL."""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator


ROOT = Path(__file__).resolve().parents[2]
PULSE_S = 20.020e-3
OUT = ROOT / "artifacts" / "hybrid_prism_diagnostics" / "phase28_c2_pulse" / "comparison"


@dataclass(frozen=True)
class Series:
    label: str
    time_us: np.ndarray
    drop_uA: np.ndarray
    baseline_uA: float


def read_elmer(path: Path, label: str) -> Series:
    table = np.loadtxt(path, delimiter=",", skiprows=1)
    time_s, current_uA = table[:, 0], table[:, 2] * 1.0e6
    baseline = float(np.mean(current_uA[(time_s >= 19.5e-3) & (time_s < PULSE_S)]))
    return Series(label, (time_s - PULSE_S) * 1.0e6, baseline - current_uA, baseline)


def read_comsol() -> Series:
    table = np.loadtxt(ROOT / "docs" / "Single-Pixel.txt", comments="%", usecols=(0, 4), encoding="utf-8")
    time_s, current_uA = table[:, 0] * 1.0e-3, table[:, 1]
    baseline = float(np.mean(current_uA[(time_s >= 19.5e-3) & (time_s < PULSE_S)]))
    return Series("COMSOL", (time_s - PULSE_S) * 1.0e6, baseline - current_uA, baseline)


def peak(series: Series) -> tuple[float, float]:
    mask = (series.time_us >= 0.0) & (series.time_us <= 500.0)
    i = np.flatnonzero(mask)[np.argmax(series.drop_uA[mask])]
    return float(series.drop_uA[i]), float(series.time_us[i])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    series = [
        read_comsol(),
        read_elmer(
            ROOT / "results" / "case_p19_pulse_phase23_tight" / "case_p19_pulse_phase23_tight_series.csv",
            "Elmer Phase23 rectangle",
        ),
        read_elmer(
            ROOT / "results" / "case_p19_pulse_phase28_c2_01ps" / "case_p19_pulse_phase28_c2_01ps_series.csv",
            "Elmer Phase28 C2, 1 ps",
        ),
    ]
    interp = {s.label: PchipInterpolator(s.time_us, s.drop_uA, extrapolate=False) for s in series}
    grid = np.unique(np.concatenate((np.linspace(-2.0, 0.0, 201), np.geomspace(5.0e-4, 500.0, 2600))))
    curves = {s.label: interp[s.label](grid) for s in series}
    p23, p28, comsol = series[1], series[2], series[0]
    p23_curve, p28_curve, comsol_curve = curves[p23.label], curves[p28.label], curves[comsol.label]
    valid = ~np.isnan(p23_curve) & ~np.isnan(p28_curve) & ~np.isnan(comsol_curve)
    max_p23_comsol = float(np.max(np.abs(p23_curve[valid] - comsol_curve[valid])))
    max_p28_comsol = float(np.max(np.abs(p28_curve[valid] - comsol_curve[valid])))
    max_p28_p23 = float(np.max(np.abs(p28_curve[valid] - p23_curve[valid])))

    metrics: dict[str, object] = {"window_us": [-2.0, 500.0], "series": {}, "max_abs_difference_uA": {}}
    for s in series:
        amp, time = peak(s)
        metrics["series"][s.label] = {"baseline_uA": s.baseline_uA, "peak_drop_uA": amp, "peak_time_us": time}
    metrics["max_abs_difference_uA"] = {
        "phase23_minus_comsol": max_p23_comsol,
        "phase28_minus_comsol": max_p28_comsol,
        "phase28_minus_phase23": max_p28_p23,
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    with (OUT / "baseline_corrected_comparison.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["time_from_pulse_us", "comsol_drop_uA", "phase23_rectangle_drop_uA", "phase28_c2_drop_uA", "phase28_minus_phase23_uA"])
        for i, time in enumerate(grid):
            writer.writerow([time, comsol_curve[i], p23_curve[i], p28_curve[i], p28_curve[i] - p23_curve[i]])

    fig, (ax, diff_ax) = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True, height_ratios=(3, 1), constrained_layout=True)
    styles = [(comsol, "#2166ac", ":"), (p23, "#d95f02", "-"), (p28, "#1b9e77", "-")]
    for s, color, ls in styles:
        ax.plot(grid, curves[s.label], label=s.label, color=color, linestyle=ls, linewidth=2.0)
    diff_ax.plot(grid, p28_curve - p23_curve, color="#6a3d9a", linewidth=1.8, label="Phase28 − Phase23")
    for axes in (ax, diff_ax):
        axes.axvline(0.0, color="#555555", linestyle="--", linewidth=0.9)
        axes.grid(True, alpha=0.25)
        axes.set_xscale("symlog", linthresh=0.1, linscale=1.0)
        axes.set_xlim(-2.0, 500.0)
    ax.set_ylabel("Current drop from baseline [µA]")
    ax.set_title("COMSOL comparison: rectangular vs 1 ps C²-smoothed heat pulse")
    ax.legend(frameon=False, loc="upper left")
    diff_ax.set_xlabel("Time from pulse [µs] (symlog)")
    diff_ax.set_ylabel("ΔI [µA]")
    diff_ax.legend(frameon=False, loc="best")
    fig.savefig(OUT / "comsol_phase23_phase28_c2_symlog.png", dpi=240)
    fig.savefig(OUT / "comsol_phase23_phase28_c2_symlog.svg")
    plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

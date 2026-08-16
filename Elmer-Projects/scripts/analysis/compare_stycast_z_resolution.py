"""Compare COMSOL with one- and eight-layer Stycast pulse responses."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PULSE_S = 20.020e-3
BASELINE_START_S = 19.5e-3
COMSOL = ROOT / "docs" / "Single-Pixel.txt"
ONE_LAYER = ROOT / "results" / "case_p19_pulse_phase23_tight" / "case_p19_pulse_phase23_tight_series.csv"
EIGHT_LAYER_RESULT = ROOT / "results" / "case_stycast_z8_pulse_225us_tight" / "case_stycast_z8_pulse_225us_tight_series.csv"
EIGHT_LAYER_RUNNING = ROOT / "case_stycast_z8_pulse_225us_tight_series.csv"
OUT = ROOT / "artifacts" / "comparison" / "stycast_z_resolution"
SAMPLE_US = np.asarray([10.0, 20.0, 40.0, 50.0, 100.0, 150.0, 200.0, 225.0])


def load_comsol() -> tuple[np.ndarray, np.ndarray, float]:
    rows = []
    for line in COMSOL.read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("%"):
            rows.append([float(value) for value in line.split()])
    table = np.asarray(rows)
    time_s, current_uA = table[:, 0] * 1e-3, table[:, 4]
    baseline = float(np.mean(current_uA[(time_s >= BASELINE_START_S) & (time_s < PULSE_S)]))
    return (time_s - PULSE_S) * 1e6, baseline - current_uA, baseline


def load_elmer(path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    table = np.genfromtxt(path, delimiter=",", names=True)
    time_s = np.atleast_1d(table["time_s"])
    current_uA = np.atleast_1d(table["tes_current_A"]) * 1e6
    pre = (time_s >= BASELINE_START_S) & (time_s < PULSE_S)
    baseline = float(np.mean(current_uA[pre]))
    return (time_s - PULSE_S) * 1e6, baseline - current_uA, baseline


def first_crossing(time_us: np.ndarray, response: np.ndarray, level: float) -> float | None:
    post = time_us >= 0.0
    x, y = time_us[post], response[post]
    indexes = np.flatnonzero(y >= level)
    if len(indexes) == 0:
        return None
    i = int(indexes[0])
    if i == 0:
        return float(x[0])
    return float(x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1]))


def main() -> None:
    eight_path = EIGHT_LAYER_RESULT if EIGHT_LAYER_RESULT.exists() else EIGHT_LAYER_RUNNING
    series = {
        "COMSOL": load_comsol(),
        "Elmer Stycast z1": load_elmer(ONE_LAYER),
        "Elmer Stycast z8": load_elmer(eight_path),
    }
    comsol_peak = float(np.max(series["COMSOL"][1]))
    metrics = {"comsol_peak_uA": comsol_peak, "series": {}}
    comsol_time, comsol_drop, _ = series["COMSOL"]
    rows = []
    for label, (time_us, drop_uA, baseline) in series.items():
        item = {
            "baseline_uA": baseline,
            "last_time_us": float(time_us[-1]),
            "t_at_10pct_comsol_peak_us": first_crossing(time_us, drop_uA, 0.1 * comsol_peak),
            "t_at_50pct_comsol_peak_us": first_crossing(time_us, drop_uA, 0.5 * comsol_peak),
            "t_at_90pct_comsol_peak_us": first_crossing(time_us, drop_uA, 0.9 * comsol_peak),
            "samples_uA": {},
        }
        common = (time_us >= 0.0) & (time_us <= min(220.0, float(comsol_time[-1])))
        difference = drop_uA[common] - np.interp(time_us[common], comsol_time, comsol_drop)
        max_index = int(np.argmax(np.abs(difference)))
        item["max_abs_difference_0_220us_uA"] = float(np.abs(difference[max_index]))
        item["max_abs_difference_time_us"] = float(time_us[common][max_index])
        item["max_abs_difference_0_220us_pct_comsol_peak"] = (
            100.0 * item["max_abs_difference_0_220us_uA"] / comsol_peak
        )
        for sample in SAMPLE_US:
            if sample <= time_us[-1]:
                value = float(np.interp(sample, time_us, drop_uA))
                item["samples_uA"][f"{sample:g}"] = value
                rows.append([label, sample, value])
        metrics["series"][label] = item

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    with (OUT / "fixed_time_current_drop.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["series", "time_from_pulse_us", "current_drop_uA"])
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    colors = {"COMSOL": "#2166ac", "Elmer Stycast z1": "#d95f02", "Elmer Stycast z8": "#1b9e77"}
    for label, (time_us, drop_uA, _) in series.items():
        mask = (time_us >= 0.0) & (time_us <= 225.0)
        ax.plot(time_us[mask], drop_uA[mask], label=label, color=colors[label], linewidth=2.0)
    ax.set(xlabel="Time from pulse [us]", ylabel="TES current drop [uA]", xlim=(0.0, 225.0))
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(OUT / "current_rise_comparison.png", dpi=220)
    fig.savefig(OUT / "current_rise_comparison.svg")
    plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

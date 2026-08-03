"""Plot the +Z 20 um source-offset result against Phase28 and COMSOL."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[2]
PULSE_S = 20.020e-3
OUT = ROOT / "artifacts" / "hybrid_prism_diagnostics" / "phase29_zplus20um" / "comparison"


def elmer(path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    t, i = a[:, 0], a[:, 2] * 1e6
    b = float(np.mean(i[(t >= 19.5e-3) & (t < PULSE_S)]))
    return (t - PULSE_S) * 1e6, b - i, b


def comsol() -> tuple[np.ndarray, np.ndarray, float]:
    a = np.loadtxt(ROOT / "docs" / "Single-Pixel.txt", comments="%", usecols=(0, 4), encoding="utf-8")
    t, i = a[:, 0] * 1e-3, a[:, 1]
    b = float(np.mean(i[(t >= 19.5e-3) & (t < PULSE_S)]))
    return (t - PULSE_S) * 1e6, b - i, b


def metric(t: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    m = (t >= 0) & (t <= 500)
    ids = np.flatnonzero(m)
    k = ids[np.argmax(y[m])]
    return float(y[k]), float(t[k])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = {
        "COMSOL": comsol(),
        "Phase23 rectangle": elmer(ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv"),
        "Phase28 C2": elmer(ROOT / "results/case_p19_pulse_phase28_c2_01ps/case_p19_pulse_phase28_c2_01ps_series.csv"),
        "Phase29 C2 +Z20um": elmer(ROOT / "results/case_p19_pulse_phase29_c2_zplus20um/case_p19_pulse_phase29_c2_zplus20um_series.csv"),
    }
    grid = np.unique(np.concatenate((np.linspace(-2, 0, 201), np.geomspace(5e-4, 500, 2600))))
    curves = {label: PchipInterpolator(t, y, extrapolate=False)(grid) for label, (t, y, _) in raw.items()}
    valid = np.all([~np.isnan(v) for v in curves.values()], axis=0)
    metrics = {label: {"baseline_uA": b, "peak_drop_uA": metric(t, y)[0], "peak_time_us": metric(t, y)[1]} for label, (t, y, b) in raw.items()}
    metrics["max_abs_difference_uA"] = {
        "phase28_minus_phase29": float(np.max(np.abs(curves["Phase28 C2"][valid] - curves["Phase29 C2 +Z20um"][valid]))),
        "phase29_minus_comsol": float(np.max(np.abs(curves["Phase29 C2 +Z20um"][valid] - curves["COMSOL"][valid]))),
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    with (OUT / "baseline_corrected_comparison.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["time_from_pulse_us", *[f"{k}_drop_uA" for k in curves], "phase29_minus_phase28_uA"])
        for i, t in enumerate(grid):
            w.writerow([t, *[curves[k][i] for k in curves], curves["Phase29 C2 +Z20um"][i] - curves["Phase28 C2"][i]])
    colors = {"COMSOL": "#2166ac", "Phase23 rectangle": "#d95f02", "Phase28 C2": "#1b9e77", "Phase29 C2 +Z20um": "#6a3d9a"}
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    for label, y in curves.items():
        ax.plot(grid, y, label=label, color=colors[label], linestyle=":" if label == "COMSOL" else "-", linewidth=2)
    ax.set_xscale("symlog", linthresh=0.1, linscale=1)
    ax.set_xlim(-2, 500)
    ax.axvline(0, color="#555", linestyle="--", linewidth=.9)
    ax.grid(True, alpha=.25)
    ax.set_xlabel("Time from pulse [µs] (symlog)")
    ax.set_ylabel("Current drop from baseline [µA]")
    ax.set_title("COMSOL / Phase23 / Phase28 / Phase29 (+Z 20 µm)")
    ax.legend(frameon=False, loc="upper left")
    fig.savefig(OUT / "comsol_phase23_phase28_phase29_zplus20um.png", dpi=240)
    fig.savefig(OUT / "comsol_phase23_phase28_phase29_zplus20um.svg")
    plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

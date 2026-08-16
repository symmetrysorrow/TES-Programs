"""Plot COMSOL/Elmer absorber-temperature responses for the pulse sweep."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PULSE_MS = 20.020
BASELINE_START_MS = 19.5
CASES = [
    ("COMSOL 1 ns", ROOT / "reference" / "SignglePixel.txt", 1),
    ("COMSOL 1 ps", ROOT / "reference" / "SinglePixel_1ps.txt", 1),
    ("COMSOL 100 ns", ROOT / "reference" / "SinglePixel_dt=1e-7.txt", 1),
    ("COMSOL 1 us", ROOT / "reference" / "SinglePixel_dt=1e-6.txt", 1),
    ("COMSOL 1 ms", ROOT / "reference" / "SinglePixel_dt=1e-3.txt", 1),
]
ELMER_DIR = ROOT / "results" / "case_tes_mpi_comsol_grid_full_uniform_continuous"


def read_case(path: Path, temperature_column: int) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("%"):
            continue
        rows.append([float(value) for value in line.split()])
    table = np.asarray(rows, dtype=float)
    time_us = (table[:, 0] - PULSE_MS) * 1.0e3
    temperature_k = table[:, temperature_column]
    baseline = float(np.mean(temperature_k[(table[:, 0] >= BASELINE_START_MS) & (table[:, 0] < PULSE_MS)]))
    delta_mk = (temperature_k - baseline) * 1.0e3
    post = time_us >= 0.0
    peak_index = int(np.argmax(delta_mk[post]))
    post_time = time_us[post]
    post_delta = delta_mk[post]
    peak = float(post_delta[peak_index])

    def crossing(level: float) -> float:
        index = int(np.flatnonzero(post_delta >= level)[0])
        if index == 0:
            return float(post_time[0])
        return float(post_time[index - 1] + (level - post_delta[index - 1]) *
                     (post_time[index] - post_time[index - 1]) /
                     (post_delta[index] - post_delta[index - 1]))

    t10, t90 = crossing(0.1 * peak), crossing(0.9 * peak)
    metric = {
        "baseline_K": baseline,
        "peak_delta_mK": peak,
        "peak_time_us": float(post_time[peak_index]),
        "final_delta_mK": float(delta_mk[-1]),
        "t10_us": t10,
        "t90_us": t90,
        "rise_10_90_us": t90 - t10,
    }
    return time_us, delta_mk, metric


def read_elmer_abs() -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Read the converged AbsorberT line from each Elmer solver time section."""
    import re

    log = (ELMER_DIR / "solver.log").read_text(encoding="utf-8", errors="ignore")
    sections = re.split(r"MAIN: Time:\s*\d+/\d+:\s*", log)[1:]
    absorber = np.asarray([
        float(re.findall(r"AbsorberT=\s*([0-9.E+-]+)", section)[-1])
        for section in sections
    ])
    series = np.genfromtxt(
        ELMER_DIR / "tes_mpi_comsol_grid_full_uniform_continuous_series.csv",
        delimiter=",",
        names=True,
    )
    n = min(len(absorber), len(series))
    time_s = series["time_s"][:n]
    temperature_k = absorber[:n]
    time_us = (time_s - PULSE_MS * 1.0e-3) * 1.0e6
    baseline = float(np.mean(temperature_k[(time_s >= BASELINE_START_MS * 1.0e-3) & (time_s < PULSE_MS * 1.0e-3)]))
    delta_mk = (temperature_k - baseline) * 1.0e3
    post = time_us >= 0.0
    post_time, post_delta = time_us[post], delta_mk[post]
    peak_index = int(np.argmax(post_delta))
    peak = float(post_delta[peak_index])

    def crossing(level: float) -> float:
        index = int(np.flatnonzero(post_delta >= level)[0])
        if index == 0:
            return float(post_time[0])
        return float(post_time[index - 1] + (level - post_delta[index - 1]) *
                     (post_time[index] - post_time[index - 1]) /
                     (post_delta[index] - post_delta[index - 1]))

    t10, t90 = crossing(0.1 * peak), crossing(0.9 * peak)
    return time_us, delta_mk, {
        "baseline_K": baseline,
        "peak_delta_mK": peak,
        "peak_time_us": float(post_time[peak_index]),
        "final_delta_mK": float(delta_mk[-1]),
        "t10_us": t10,
        "t90_us": t90,
        "rise_10_90_us": t90 - t10,
    }


def main() -> None:
    out_dir = ROOT / "artifacts" / "comparison" / "abs_temperature_dt"
    out_dir.mkdir(parents=True, exist_ok=True)
    series = []
    metrics = []
    for label, path, column in CASES:
        time_us, delta_mk, metric = read_case(path, column)
        series.append((label, time_us, delta_mk))
        metrics.append((label, path.name, metric))
    elmer_time_us, elmer_delta_mk, elmer_metric = read_elmer_abs()
    series.append(("Elmer abs (solver.log)", elmer_time_us, elmer_delta_mk))
    metrics.append(("Elmer abs (solver.log)", "solver.log", elmer_metric))

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.3), constrained_layout=True)
    for ax, xlim, title in zip(
        axes,
        ((0.0, 20.0), (0.0, 1300.0)),
        ("Abs temperature: first 20 µs", "Abs temperature: full 1.3 ms window"),
    ):
        for index, (label, time_us, delta_mk) in enumerate(series):
            mask = (time_us >= xlim[0]) & (time_us <= xlim[1])
            style = {"linestyle": "--", "linewidth": 2.2} if label.startswith("Elmer") else {"linewidth": 1.8}
            ax.plot(time_us[mask], delta_mk[mask], label=label, color=colors[index], **style)
        ax.axvline(0.0, color="0.35", linewidth=0.9, linestyle="--")
        ax.set_xlim(*xlim)
        ax.set_xlabel("Time from pulse start [µs]")
        ax.set_ylabel("Abs temperature rise [mK]")
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
    axes[0].legend(loc="upper left", fontsize=8, frameon=False)
    axes[1].text(0.98, 0.04, "Elmer: AbsorberT from converged solver.log section", transform=axes[1].transAxes,
                 ha="right", va="bottom", fontsize=9, color="0.25",
                 bbox={"facecolor": "white", "edgecolor": "0.65", "alpha": 0.85, "pad": 4})
    fig.suptitle("SinglePixel absorber temperature vs heat-application duration", fontsize=14)
    fig.savefig(out_dir / "abs_temperature_dt_comparison.png", dpi=220)
    plt.close(fig)

    with (out_dir / "abs_temperature_dt_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["series", "source", "baseline_K", "peak_delta_mK", "peak_time_us", "t10_us", "t90_us", "rise_10_90_us", "final_delta_mK"])
        for label, source, metric in metrics:
            writer.writerow([label, source, metric["baseline_K"], metric["peak_delta_mK"], metric["peak_time_us"], metric["t10_us"], metric["t90_us"], metric["rise_10_90_us"], metric["final_delta_mK"]])
    print(out_dir / "abs_temperature_dt_comparison.png")


if __name__ == "__main__":
    main()

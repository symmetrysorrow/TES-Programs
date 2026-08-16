"""Compare the normalized rise of absorber T, TES T, and TES current.

COMSOL uses the standard 1 ns reference table. Elmer absorber temperature is
the converged ``AbsorberT`` value from each solver time section; TES temperature
and current come from the corresponding series CSV.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PULSE_MS = 20.020
BASELINE_START_MS = 19.5
TIME_ZERO_MS = 19.05
PULSE_OFFSET_US = (PULSE_MS - TIME_ZERO_MS) * 1.0e3
COMSOL = ROOT / "reference" / "SignglePixel.txt"
COMSOL_1US = ROOT / "reference" / "SinglePixel_dt=1e-6.txt"
ELMER_DIR = ROOT / "results" / "case_tes_mpi_comsol_grid_full_uniform_continuous"


def crossings(time_us: np.ndarray, response: np.ndarray) -> dict[str, float]:
    peak_index = int(np.argmax(response))
    peak = float(response[peak_index])
    rising_time = time_us[: peak_index + 1]
    rising_response = response[: peak_index + 1]

    def at(level: float) -> float:
        indexes = np.flatnonzero(rising_response >= level)
        if len(indexes) == 0:
            return float("nan")
        index = int(indexes[0])
        if index == 0:
            return float(rising_time[0])
        return float(rising_time[index - 1] + (level - rising_response[index - 1]) *
                     (rising_time[index] - rising_time[index - 1]) /
                     (rising_response[index] - rising_response[index - 1]))

    t10, t90 = at(0.1 * peak), at(0.9 * peak)
    return {"peak": peak, "peak_time_us": float(time_us[peak_index]),
            "t10_us": t10, "t90_us": t90, "rise_10_90_us": t90 - t10}


def normalize(time_us: np.ndarray, signal: np.ndarray, baseline_mask: np.ndarray,
              post_mask: np.ndarray, descending: bool = False,
              metric_origin_us: float = 0.0) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    baseline = float(np.mean(signal[baseline_mask]))
    response = baseline - signal if descending else signal - baseline
    time_post, response_post = time_us[post_mask], response[post_mask]
    metric = crossings(time_post - metric_origin_us, response_post)
    amplitude = metric["peak"]
    plot_mask = time_us >= 0.0
    plot_metric = {
        "t10_plot_us": metric["t10_us"] + metric_origin_us,
        "t90_plot_us": metric["t90_us"] + metric_origin_us,
    }
    return time_us[plot_mask], response[plot_mask] / amplitude, {
        "baseline": baseline, "amplitude": amplitude, **metric, **plot_metric,
    }


def read_comsol(path: Path = COMSOL, abs_column: int = 1, tes_column: int = 3,
                current_column: int = 4) -> dict[str, tuple[np.ndarray, np.ndarray, dict[str, float]]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("%"):
            rows.append([float(value) for value in line.split()])
    table = np.asarray(rows, dtype=float)
    time_us = (table[:, 0] - TIME_ZERO_MS) * 1.0e3
    baseline_mask = (table[:, 0] >= BASELINE_START_MS) & (table[:, 0] < PULSE_MS)
    post_mask = table[:, 0] >= PULSE_MS
    return {
        "abs": normalize(time_us, table[:, abs_column], baseline_mask, post_mask,
                          metric_origin_us=PULSE_OFFSET_US),
        "tes": normalize(time_us, table[:, tes_column], baseline_mask, post_mask,
                          metric_origin_us=PULSE_OFFSET_US),
        "current": normalize(time_us, table[:, current_column], baseline_mask, post_mask,
                              descending=True, metric_origin_us=PULSE_OFFSET_US),
    }


def read_elmer() -> dict[str, tuple[np.ndarray, np.ndarray, dict[str, float]]]:
    log = (ELMER_DIR / "solver.log").read_text(encoding="utf-8", errors="ignore")
    sections = re.split(r"MAIN: Time:\s*\d+/\d+:\s*", log)[1:]
    absorber = np.asarray([
        float(re.findall(r"AbsorberT=\s*([0-9.E+-]+)", section)[-1])
        for section in sections
    ])
    series = np.genfromtxt(
        ELMER_DIR / "tes_mpi_comsol_grid_full_uniform_continuous_series.csv",
        delimiter=",", names=True,
    )
    n = min(len(absorber), len(series))
    time_us = (series["time_s"][:n] - TIME_ZERO_MS * 1.0e-3) * 1.0e6
    pre = (time_us >= (BASELINE_START_MS - TIME_ZERO_MS) * 1.0e3) & (time_us < PULSE_OFFSET_US)
    post = time_us >= PULSE_OFFSET_US
    return {
        "abs": normalize(time_us, absorber[:n], pre, post, metric_origin_us=PULSE_OFFSET_US),
        "tes": normalize(time_us, series["tes_temperature_K"][:n], pre, post,
                          metric_origin_us=PULSE_OFFSET_US),
        "current": normalize(time_us, series["tes_current_A"][:n], pre, post,
                              descending=True, metric_origin_us=PULSE_OFFSET_US),
    }


def main() -> None:
    comsol = read_comsol()
    comsol_1us = read_comsol(COMSOL_1US, abs_column=1, tes_column=2, current_column=3)
    elmer = read_elmer()
    out_dir = ROOT / "artifacts" / "comparison" / "rise_abs_tes_current"
    out_dir.mkdir(parents=True, exist_ok=True)
    panels = [("abs", "Absorber temperature rise", 1_000.0),
              ("tes", "TES temperature rise", 1_300.0),
              ("current", "TES current drop", 1_300.0)]
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8), constrained_layout=True)
    sources = [("COMSOL 1 ns", comsol), ("COMSOL 1 us", comsol_1us), ("Elmer", elmer)]
    colors = {"COMSOL 1 ns": "#2166ac", "COMSOL 1 us": "#2ca25f", "Elmer": "#d95f02"}
    linestyles = {"COMSOL 1 ns": "-", "COMSOL 1 us": "--", "Elmer": "-"}
    metrics_rows = []
    for ax, (key, title, xmax) in zip(axes, panels):
        for label, collection in sources:
            source = collection[key]
            time_us, normalized_response, metric = source
            mask = (time_us >= 0.0) & (time_us <= xmax)
            ax.plot(time_us[mask], normalized_response[mask], color=colors[label],
                    linestyle=linestyles[label], linewidth=2.0, label=label)
            ax.axvline(metric["t10_plot_us"], color=colors[label], alpha=.35, linewidth=.8)
            ax.axvline(metric["t90_plot_us"], color=colors[label], alpha=.35, linewidth=.8)
            metrics_rows.append([title, label, metric["baseline"], metric["amplitude"],
                                 metric["peak_time_us"], metric["t10_us"], metric["t90_us"], metric["rise_10_90_us"]])
        ax.set_xlim(0.0, xmax)
        ax.set_ylim(-0.02, 1.05)
        ax.set_title(title)
        ax.set_xlabel("Time from 19.05 ms [µs]")
        ax.set_ylabel("Normalized response")
        ax.grid(True, alpha=.25)
        ax.axhline(.1, color="0.45", linewidth=.7, linestyle=":")
        ax.axhline(.9, color="0.45", linewidth=.7, linestyle=":")
        ax.axvline(PULSE_OFFSET_US, color="0.2", linewidth=1.0, linestyle="--")
    axes[-1].legend(frameon=False, loc="lower right")
    fig.suptitle("Rise comparison (time zero = 19.05 ms; pulse = 20.020 ms)", fontsize=14)
    fig.savefig(out_dir / "rise_abs_tes_current.png", dpi=220)
    plt.close(fig)

    with (out_dir / "rise_abs_tes_current_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "series", "baseline", "response_amplitude", "peak_time_us", "t10_us", "t90_us", "rise_10_90_us"])
        writer.writerows(metrics_rows)
    print(out_dir / "rise_abs_tes_current.png")


if __name__ == "__main__":
    main()

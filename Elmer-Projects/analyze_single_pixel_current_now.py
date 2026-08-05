from __future__ import annotations

from pathlib import Path
import csv
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
COMSOL_PATH = ROOT / "reference" / "SignglePixel.txt"
ELMER_PATH = (
    ROOT
    / "results"
    / "raw"
    / "legacy-root-output"
    / "tes_mpi_comsol_grid_series.csv"
)
OUT = ROOT / "artifacts" / "single_pixel_comparison_now"

PULSE_START_MS = 20.02
BASELINE_WINDOW_MS = (19.5, PULSE_START_MS)


def baseline(time_ms: np.ndarray, current_a: np.ndarray) -> float:
    mask = (time_ms >= BASELINE_WINDOW_MS[0]) & (time_ms <= BASELINE_WINDOW_MS[1])
    if not np.any(mask):
        raise ValueError(f"no samples in baseline window {BASELINE_WINDOW_MS} ms")
    return float(np.mean(current_a[mask]))


def crossing_time(
    time_ms: np.ndarray,
    response_a: np.ndarray,
    level_a: float,
    peak_index: int,
) -> float:
    before_peak = response_a[: peak_index + 1]
    indices = np.flatnonzero(
        (before_peak[:-1] < level_a) & (before_peak[1:] >= level_a)
    )
    if len(indices) == 0:
        return float("nan")
    index = int(indices[-1])
    fraction = (level_a - response_a[index]) / (
        response_a[index + 1] - response_a[index]
    )
    return float(
        time_ms[index] + fraction * (time_ms[index + 1] - time_ms[index])
    )


def pulse_metrics(
    name: str,
    time_ms: np.ndarray,
    current_a: np.ndarray,
    common_end_ms: float,
) -> dict[str, float | str]:
    base_a = baseline(time_ms, current_a)
    response_a = base_a - current_a
    post = np.flatnonzero(
        (time_ms >= PULSE_START_MS) & (time_ms <= common_end_ms)
    )
    if len(post) == 0:
        raise ValueError(f"{name}: no samples after pulse start")
    peak_index = int(post[np.argmax(response_a[post])])
    peak_a = float(response_a[peak_index])
    t10_ms = crossing_time(time_ms, response_a, 0.1 * peak_a, peak_index)
    t90_ms = crossing_time(time_ms, response_a, 0.9 * peak_a, peak_index)
    return {
        "series": name,
        "baseline_uA": base_a * 1e6,
        "peak_drop_uA": peak_a * 1e6,
        "peak_time_ms": float(time_ms[peak_index]),
        "peak_delay_ms": float(time_ms[peak_index] - PULSE_START_MS),
        "rise_time_10_90_ms": float(t90_ms - t10_ms),
        "comparison_end_ms": common_end_ms,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    comsol = np.loadtxt(COMSOL_PATH, comments="%", encoding="utf-8")
    elmer = np.genfromtxt(ELMER_PATH, delimiter=",", names=True)

    # COMSOL columns are: time, absorber T, StyCast T, TES T, current, TES R.
    # The previous version incorrectly used column 3 (TES temperature).
    t_comsol_ms = comsol[:, 0]
    i_comsol_a = comsol[:, 4] * 1e-6
    t_elmer_ms = elmer["time_s"] * 1e3
    i_elmer_a = elmer["tes_current_A"]

    common_end_ms = float(min(t_comsol_ms.max(), t_elmer_ms.max()))
    rows = [
        pulse_metrics(
            "COMSOL single-pixel", t_comsol_ms, i_comsol_a, common_end_ms
        ),
        pulse_metrics(
            "Elmer 4-rank all-tet", t_elmer_ms, i_elmer_a, common_end_ms
        ),
    ]

    with (OUT / "metrics_corrected.csv").open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "metrics_corrected.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    base_comsol_a = baseline(t_comsol_ms, i_comsol_a)
    base_elmer_a = baseline(t_elmer_ms, i_elmer_a)
    response_comsol_a = base_comsol_a - i_comsol_a
    response_elmer_a = base_elmer_a - i_elmer_a
    peak_comsol_a = rows[0]["peak_drop_uA"] * 1e-6
    peak_elmer_a = rows[1]["peak_drop_uA"] * 1e-6

    mask_comsol = (t_comsol_ms >= 19.5) & (t_comsol_ms <= common_end_ms)
    mask_elmer = (t_elmer_ms >= 19.5) & (t_elmer_ms <= common_end_ms)

    grid_ms = np.unique(
        np.concatenate(
            [t_comsol_ms[mask_comsol], t_elmer_ms[mask_elmer]]
        )
    )
    with (OUT / "timeseries_corrected.csv").open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "time_ms",
                "COMSOL_current_uA",
                "Elmer_current_uA",
                "COMSOL_drop_uA",
                "Elmer_drop_uA",
                "COMSOL_normalized",
                "Elmer_normalized",
            ]
        )
        for time_ms in grid_ms:
            comsol_current_a = np.interp(time_ms, t_comsol_ms, i_comsol_a)
            elmer_current_a = np.interp(time_ms, t_elmer_ms, i_elmer_a)
            comsol_drop_a = base_comsol_a - comsol_current_a
            elmer_drop_a = base_elmer_a - elmer_current_a
            writer.writerow(
                [
                    time_ms,
                    comsol_current_a * 1e6,
                    elmer_current_a * 1e6,
                    comsol_drop_a * 1e6,
                    elmer_drop_a * 1e6,
                    comsol_drop_a / peak_comsol_a,
                    elmer_drop_a / peak_elmer_a,
                ]
            )

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.2), sharex=True)
    for axis in axes:
        axis.axvline(PULSE_START_MS, color="0.4", linestyle="--", linewidth=1)
        axis.grid(alpha=0.25)

    axes[0].plot(
        t_comsol_ms[mask_comsol],
        i_comsol_a[mask_comsol] * 1e6,
        label="COMSOL",
        linewidth=2,
    )
    axes[0].plot(
        t_elmer_ms[mask_elmer],
        i_elmer_a[mask_elmer] * 1e6,
        label="Elmer",
        linewidth=1.7,
    )
    axes[0].set_ylabel("TES current [µA]")
    axes[0].set_title("Single-pixel TES current comparison (correct COMSOL column)")
    axes[0].legend()

    axes[1].plot(
        t_comsol_ms[mask_comsol],
        response_comsol_a[mask_comsol] / peak_comsol_a,
        label="COMSOL",
        linewidth=2,
    )
    axes[1].plot(
        t_elmer_ms[mask_elmer],
        response_elmer_a[mask_elmer] / peak_elmer_a,
        label="Elmer",
        linewidth=1.7,
    )
    axes[1].set_xlabel("Time [ms]")
    axes[1].set_ylabel("Baseline-corrected / peak")
    axes[1].set_ylim(-0.05, 1.08)
    axes[1].legend()

    summary = (
        f"COMSOL: peak={rows[0]['peak_drop_uA']:.3f} µA, "
        f"delay={rows[0]['peak_delay_ms']:.3f} ms, "
        f"rise={rows[0]['rise_time_10_90_ms']:.3f} ms\n"
        f"Elmer: peak={rows[1]['peak_drop_uA']:.3f} µA, "
        f"delay={rows[1]['peak_delay_ms']:.3f} ms, "
        f"rise={rows[1]['rise_time_10_90_ms']:.3f} ms"
    )
    fig.text(0.5, 0.01, summary, ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=(0, 0.065, 1, 1))
    fig.savefig(OUT / "current_comparison_corrected.png", dpi=220)
    plt.close(fig)

    print(OUT / "current_comparison_corrected.png")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()

"""Separate Stycast vertical-mesh and time-step effects in the early TES rise."""

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
OUT = ROOT / "artifacts" / "comparison" / "stycast_z16_resolution"
SAMPLE_US = np.asarray([10.0, 20.0, 40.0, 50.0, 75.0, 100.0])
ONE25_SHORT = ROOT / (
    "results/case_stycast_z16_pulse_105us_fine1p25us_tight/"
    "case_stycast_z16_pulse_105us_fine1p25us_tight_series.csv"
)
ONE25_LONG = ROOT / (
    "results/case_stycast_z16_pulse_225us_fine1p25us_tight/"
    "case_stycast_z16_pulse_225us_fine1p25us_tight_series.csv"
)
ELMER_PATHS = {
    "Elmer z1, 10 us": ROOT
    / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv",
    "Elmer z8, 10 us": ROOT
    / "results/case_stycast_z8_pulse_225us_tight/case_stycast_z8_pulse_225us_tight_series.csv",
    "Elmer z16, 10 us": ROOT
    / "results/case_stycast_z16_pulse_105us_tight/case_stycast_z16_pulse_105us_tight_series.csv",
    "Elmer z16, 5 us": ROOT
    / "results/case_stycast_z16_pulse_105us_fine5us_tight/case_stycast_z16_pulse_105us_fine5us_tight_series.csv",
    "Elmer z16, 2.5 us": ROOT
    / "results/case_stycast_z16_pulse_105us_fine2p5us_tight/case_stycast_z16_pulse_105us_fine2p5us_tight_series.csv",
    "Elmer z16, 1.25 us": ONE25_LONG if ONE25_LONG.exists() else ONE25_SHORT,
    "Elmer z32, 1.25 us": ROOT
    / "results/case_stycast_z32_pulse_105us_fine1p25us_tight/case_stycast_z32_pulse_105us_fine1p25us_tight_series.csv",
    "Elmer z32, 0.625 us": ROOT
    / "results/case_stycast_z32_pulse_105us_fine0p625us_tight/case_stycast_z32_pulse_105us_fine0p625us_tight_series.csv",
}


def load_comsol() -> tuple[np.ndarray, np.ndarray, float]:
    table = np.loadtxt(COMSOL, comments="%", encoding="utf-8")
    time_s, current_uA = table[:, 0] * 1e-3, table[:, 4]
    pre = (time_s >= BASELINE_START_S) & (time_s < PULSE_S)
    baseline = float(np.mean(current_uA[pre]))
    return (time_s - PULSE_S) * 1e6, baseline - current_uA, baseline


def load_elmer(path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    table = np.genfromtxt(path, delimiter=",", names=True)
    time_s = np.atleast_1d(table["time_s"])
    current_uA = np.atleast_1d(table["tes_current_A"]) * 1e6
    pre = (time_s >= BASELINE_START_S) & (time_s < PULSE_S)
    baseline = float(np.mean(current_uA[pre]))
    return (time_s - PULSE_S) * 1e6, baseline - current_uA, baseline


def crossing(time_us: np.ndarray, response_uA: np.ndarray, level_uA: float) -> float | None:
    mask = time_us >= 0.0
    x, y = time_us[mask], response_uA[mask]
    indexes = np.flatnonzero(y >= level_uA)
    if not len(indexes):
        return None
    i = int(indexes[0])
    if i == 0:
        return float(x[0])
    return float(x[i - 1] + (level_uA - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1]))


def max_difference(
    left: tuple[np.ndarray, np.ndarray, float],
    right: tuple[np.ndarray, np.ndarray, float],
    end_us: float = 100.0,
) -> dict[str, float]:
    grid = np.linspace(0.0, end_us, 1001)
    delta = np.interp(grid, left[0], left[1]) - np.interp(grid, right[0], right[1])
    i = int(np.argmax(np.abs(delta)))
    return {
        "max_abs_uA": float(abs(delta[i])),
        "time_us": float(grid[i]),
        "signed_left_minus_right_uA": float(delta[i]),
    }


def main() -> None:
    missing = [str(path) for path in ELMER_PATHS.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing completed series:\n" + "\n".join(missing))

    series = {"COMSOL": load_comsol()}
    series.update({label: load_elmer(path) for label, path in ELMER_PATHS.items()})
    peak_uA = float(np.max(series["COMSOL"][1]))
    metrics: dict[str, object] = {"comsol_peak_uA": peak_uA, "series": {}, "pairwise": {}}
    rows: list[list[object]] = []

    for label, (time_us, response_uA, baseline_uA) in series.items():
        item = {
            "baseline_uA": baseline_uA,
            "last_logged_time_us": float(time_us[-1]),
            "t10_us": crossing(time_us, response_uA, 0.1 * peak_uA),
            "t50_us": crossing(time_us, response_uA, 0.5 * peak_uA),
            "t90_us": crossing(time_us, response_uA, 0.9 * peak_uA),
            "samples_uA": {},
        }
        item["t10_to_t90_us"] = (
            None if item["t90_us"] is None else item["t90_us"] - item["t10_us"]
        )
        for sample_us in SAMPLE_US:
            value = float(np.interp(sample_us, time_us, response_uA))
            item["samples_uA"][f"{sample_us:g}"] = value
            rows.append([label, sample_us, value])
        metrics["series"][label] = item

    for name, left, right in (
        ("z8_10us_minus_z16_10us", "Elmer z8, 10 us", "Elmer z16, 10 us"),
        ("z16_10us_minus_z16_5us", "Elmer z16, 10 us", "Elmer z16, 5 us"),
        ("z16_5us_minus_z16_2p5us", "Elmer z16, 5 us", "Elmer z16, 2.5 us"),
        ("z16_2p5us_minus_z16_1p25us", "Elmer z16, 2.5 us", "Elmer z16, 1.25 us"),
        ("z16_1p25us_minus_z32_1p25us", "Elmer z16, 1.25 us", "Elmer z32, 1.25 us"),
        ("z32_1p25us_minus_z32_0p625us", "Elmer z32, 1.25 us", "Elmer z32, 0.625 us"),
        ("z16_10us_minus_comsol", "Elmer z16, 10 us", "COMSOL"),
        ("z16_5us_minus_comsol", "Elmer z16, 5 us", "COMSOL"),
        ("z16_2p5us_minus_comsol", "Elmer z16, 2.5 us", "COMSOL"),
        ("z16_1p25us_minus_comsol", "Elmer z16, 1.25 us", "COMSOL"),
        ("z32_1p25us_minus_comsol", "Elmer z32, 1.25 us", "COMSOL"),
        ("z32_0p625us_minus_comsol", "Elmer z32, 0.625 us", "COMSOL"),
    ):
        item = max_difference(series[left], series[right])
        item["pct_comsol_peak"] = 100.0 * item["max_abs_uA"] / peak_uA
        metrics["pairwise"][name] = item

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    with (OUT / "fixed_time_current_drop.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["series", "time_from_pulse_us", "current_drop_uA"])
        writer.writerows(rows)

    colors = {
        "COMSOL": "#2166ac",
        "Elmer z1, 10 us": "#d95f02",
        "Elmer z8, 10 us": "#1b9e77",
        "Elmer z16, 10 us": "#7570b3",
        "Elmer z16, 5 us": "#e7298a",
        "Elmer z16, 2.5 us": "#66a61e",
        "Elmer z16, 1.25 us": "#e6ab02",
        "Elmer z32, 1.25 us": "#a6761d",
        "Elmer z32, 0.625 us": "#666666",
    }
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 7.0), sharex=True, constrained_layout=True)
    for label, (time_us, response_uA, _) in series.items():
        mask = (time_us >= 0.0) & (time_us <= 105.0)
        style = ":" if label == "COMSOL" else "-"
        axes[0].plot(time_us[mask], response_uA[mask], label=label, color=colors[label],
                     linestyle=style, linewidth=2.0)
    grid = np.linspace(0.0, 100.0, 1001)
    comsol = np.interp(grid, series["COMSOL"][0], series["COMSOL"][1])
    for label in (
        "Elmer z8, 10 us", "Elmer z16, 10 us", "Elmer z16, 5 us", "Elmer z16, 2.5 us",
        "Elmer z16, 1.25 us", "Elmer z32, 1.25 us", "Elmer z32, 0.625 us"
    ):
        curve = np.interp(grid, series[label][0], series[label][1])
        axes[1].plot(grid, curve - comsol, label=label, color=colors[label], linewidth=1.8)
    axes[0].set_ylabel("TES current drop [uA]")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False, fontsize=9)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set(xlabel="Time from pulse [us]", ylabel="Elmer - COMSOL [uA]", xlim=(0.0, 105.0))
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False, fontsize=9)
    fig.savefig(OUT / "current_rise_and_residual.png", dpi=220)
    fig.savefig(OUT / "current_rise_and_residual.svg")
    plt.close(fig)

    s = metrics["series"]
    p = metrics["pairwise"]
    summary = f"""# Stycast vertical-resolution and time-step validation

All crossing levels use the COMSOL peak current drop ({peak_uA:.6f} uA).

| Series | Baseline [uA] | t10 [us] | t50 [us] | Drop at 40 us [uA] | Drop at 100 us [uA] |
|---|---:|---:|---:|---:|---:|
"""
    for label in series:
        item = s[label]
        summary += (
            f"| {label} | {item['baseline_uA']:.6f} | {item['t10_us']:.4f} | "
            f"{item['t50_us']:.4f} | {item['samples_uA']['40']:.6f} | "
            f"{item['samples_uA']['100']:.6f} |\n"
        )
    summary += f"""
## Separation of effects over 0--100 us

- Spatial refinement, z8 to z16 at 10 us: maximum change {p['z8_10us_minus_z16_10us']['max_abs_uA']:.6f} uA ({p['z8_10us_minus_z16_10us']['pct_comsol_peak']:.3f}% of COMSOL peak).
- Time refinement, 10 us to 5 us at z16: maximum change {p['z16_10us_minus_z16_5us']['max_abs_uA']:.6f} uA ({p['z16_10us_minus_z16_5us']['pct_comsol_peak']:.3f}% of COMSOL peak).
- Time refinement, 5 us to 2.5 us at z16: maximum change {p['z16_5us_minus_z16_2p5us']['max_abs_uA']:.6f} uA ({p['z16_5us_minus_z16_2p5us']['pct_comsol_peak']:.3f}% of COMSOL peak).
- Time refinement, 2.5 us to 1.25 us at z16: maximum change {p['z16_2p5us_minus_z16_1p25us']['max_abs_uA']:.6f} uA ({p['z16_2p5us_minus_z16_1p25us']['pct_comsol_peak']:.3f}% of COMSOL peak).
- Spatial refinement, z16 to z32 at 1.25 us: maximum change {p['z16_1p25us_minus_z32_1p25us']['max_abs_uA']:.6f} uA ({p['z16_1p25us_minus_z32_1p25us']['pct_comsol_peak']:.3f}% of COMSOL peak).
- Time refinement, 1.25 us to 0.625 us at z32: maximum change {p['z32_1p25us_minus_z32_0p625us']['max_abs_uA']:.6f} uA ({p['z32_1p25us_minus_z32_0p625us']['pct_comsol_peak']:.3f}% of COMSOL peak).
- z16/10 us versus COMSOL: maximum residual {p['z16_10us_minus_comsol']['max_abs_uA']:.6f} uA ({p['z16_10us_minus_comsol']['pct_comsol_peak']:.3f}% of COMSOL peak).
- z16/5 us versus COMSOL: maximum residual {p['z16_5us_minus_comsol']['max_abs_uA']:.6f} uA ({p['z16_5us_minus_comsol']['pct_comsol_peak']:.3f}% of COMSOL peak).
- z16/2.5 us versus COMSOL: maximum residual {p['z16_2p5us_minus_comsol']['max_abs_uA']:.6f} uA ({p['z16_2p5us_minus_comsol']['pct_comsol_peak']:.3f}% of COMSOL peak).
- z16/1.25 us versus COMSOL: maximum residual {p['z16_1p25us_minus_comsol']['max_abs_uA']:.6f} uA ({p['z16_1p25us_minus_comsol']['pct_comsol_peak']:.3f}% of COMSOL peak).
- z32/1.25 us versus COMSOL: maximum residual {p['z32_1p25us_minus_comsol']['max_abs_uA']:.6f} uA ({p['z32_1p25us_minus_comsol']['pct_comsol_peak']:.3f}% of COMSOL peak).
- z32/0.625 us versus COMSOL: maximum residual {p['z32_0p625us_minus_comsol']['max_abs_uA']:.6f} uA ({p['z32_0p625us_minus_comsol']['pct_comsol_peak']:.3f}% of COMSOL peak).

## Interpretation

The z8-to-z16 spatial change and successive time-step changes show which
discretization is limiting the early response.  The remaining error is not by itself evidence that the
0.16 um TES needs through-thickness refinement.  With the model properties,
the TES diffusivity is about 3.15 m2/s and its thickness diffusion scale is
only 8.1e-15 s.  In contrast, the 20 um Stycast diffusivity is about
5.09e-7 m2/s, giving a slab diffusion scale L2/(pi2 alpha) of 79.6 us.  A
single z16 Stycast element has h2/alpha = 3.07 us, making the 2.5 us and
1.25 us grids physically relevant convergence checks.
"""
    (OUT / "summary.md").write_text(summary, encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

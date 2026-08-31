"""Compare optimized SinglePixel AMGX current with the COMSOL reference."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PULSE_S = 20.020e-3
BASELINE_START_S = 19.5e-3
DEFAULT_ELMER = (
    ROOT
    / "results"
    / "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid_amgx_smoke_177step"
    / "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid_amgx_smoke_177step_series.csv"
)
DEFAULT_OUT = (
    ROOT
    / "artifacts"
    / "comparison"
    / "comsol_amgx_singlepixel_prod_v2_hybrid_100us"
)
SAMPLE_US = (10.0, 20.0, 40.0, 50.0, 100.0)


@dataclass(frozen=True)
class Series:
    time_us: np.ndarray
    current_uA: np.ndarray
    baseline_uA: float

    @property
    def drop_uA(self) -> np.ndarray:
        return self.baseline_uA - self.current_uA


def normalized_series(time_s: np.ndarray, current_uA: np.ndarray) -> Series:
    pre = (time_s >= BASELINE_START_S) & (time_s < PULSE_S)
    if not np.any(pre):
        raise ValueError("series has no samples in the pre-pulse baseline window")
    baseline = float(np.mean(current_uA[pre]))
    return Series((time_s - PULSE_S) * 1.0e6, current_uA, baseline)


def read_comsol(path: Path) -> Series:
    table = np.loadtxt(path, comments="%", encoding="utf-8")
    return normalized_series(table[:, 0] * 1.0e-3, table[:, 4])


def read_elmer(path: Path) -> Series:
    table = np.genfromtxt(path, delimiter=",", names=True)
    time_s = np.atleast_1d(table["time_s"])
    current_uA = np.atleast_1d(table["tes_current_A"]) * 1.0e6
    order = np.argsort(time_s, kind="stable")
    time_s, current_uA = time_s[order], current_uA[order]
    unique_time, starts, counts = np.unique(
        time_s, return_index=True, return_counts=True
    )
    unique_current = np.empty_like(unique_time)
    for index, (start, count) in enumerate(zip(starts, counts)):
        values = current_uA[start : start + count]
        if float(np.ptp(values)) > 1.0e-9:
            raise ValueError(
                f"conflicting Elmer current values at time {unique_time[index]:.12g} s"
            )
        unique_current[index] = float(np.mean(values))
    return normalized_series(
        unique_time,
        unique_current,
    )


def crossing(time_us: np.ndarray, response_uA: np.ndarray, level_uA: float) -> float | None:
    post = time_us >= 0.0
    x, y = time_us[post], response_uA[post]
    indexes = np.flatnonzero(y >= level_uA)
    if not len(indexes):
        return None
    index = int(indexes[0])
    if index == 0:
        return float(x[0])
    return float(
        x[index - 1]
        + (level_uA - y[index - 1])
        * (x[index] - x[index - 1])
        / (y[index] - y[index - 1])
    )


def formatted_time(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.4f}"


def compare(comsol: Series, elmer: Series, end_us: float) -> tuple[dict[str, object], np.ndarray]:
    common_end = min(end_us, float(comsol.time_us[-1]), float(elmer.time_us[-1]))
    if common_end <= 0.0:
        raise ValueError("series do not overlap after the pulse")
    if common_end < end_us - 1.0e-2:
        raise ValueError(
            f"Elmer series ends at {common_end:.3f} us; "
            f"requested comparison requires {end_us:.3f} us"
        )
    grid = np.linspace(0.0, common_end, int(np.ceil(common_end / 0.05)) + 1)
    comsol_drop = np.interp(grid, comsol.time_us, comsol.drop_uA)
    elmer_drop = np.interp(grid, elmer.time_us, elmer.drop_uA)
    difference = elmer_drop - comsol_drop
    max_index = int(np.argmax(np.abs(difference)))
    comsol_peak = float(np.max(comsol.drop_uA[comsol.time_us >= 0.0]))

    samples: dict[str, object] = {}
    for sample in SAMPLE_US:
        if sample <= common_end + 1.0e-9:
            c = float(np.interp(sample, comsol.time_us, comsol.drop_uA))
            e = float(np.interp(sample, elmer.time_us, elmer.drop_uA))
            samples[f"{sample:g}"] = {
                "comsol_drop_uA": c,
                "amgx_drop_uA": e,
                "amgx_minus_comsol_uA": e - c,
                "difference_pct_comsol_peak": 100.0 * (e - c) / comsol_peak,
            }

    metrics: dict[str, object] = {
        "comparison_window_us": [0.0, common_end],
        "comsol_peak_full_trace_uA": comsol_peak,
        "baseline_uA": {
            "COMSOL": comsol.baseline_uA,
            "AMGX": elmer.baseline_uA,
            "AMGX_error_pct": 100.0
            * (elmer.baseline_uA - comsol.baseline_uA)
            / comsol.baseline_uA,
        },
        "crossing_us_at_comsol_peak_fraction": {},
        "max_abs_difference_uA": float(abs(difference[max_index])),
        "max_abs_difference_time_us": float(grid[max_index]),
        "max_abs_difference_pct_comsol_peak": float(
            100.0 * abs(difference[max_index]) / comsol_peak
        ),
        "rmse_uA": float(np.sqrt(np.mean(difference * difference))),
        "rmse_pct_comsol_peak": float(
            100.0 * np.sqrt(np.mean(difference * difference)) / comsol_peak
        ),
        "samples": samples,
    }
    crossings = metrics["crossing_us_at_comsol_peak_fraction"]
    assert isinstance(crossings, dict)
    for fraction in (0.1, 0.5, 0.9):
        level = fraction * comsol_peak
        c = crossing(comsol.time_us, comsol.drop_uA, level)
        e = crossing(elmer.time_us, elmer.drop_uA, level)
        crossings[f"{fraction:g}"] = {
            "COMSOL": c,
            "AMGX": e,
            "AMGX_minus_COMSOL": None if c is None or e is None else e - c,
        }
    aligned = np.column_stack((grid, comsol_drop, elmer_drop, difference))
    return metrics, aligned


def write_outputs(
    out: Path,
    comsol: Series,
    elmer: Series,
    metrics: dict[str, object],
    aligned: np.ndarray,
) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
    )
    np.savetxt(
        out / "aligned_current.csv",
        aligned,
        delimiter=",",
        header=(
            "time_from_pulse_us,comsol_current_drop_uA,"
            "amgx_current_drop_uA,amgx_minus_comsol_uA"
        ),
        comments="",
    )

    end_us = float(metrics["comparison_window_us"][1])  # type: ignore[index]
    fig, axes = plt.subplots(
        2, 1, figsize=(8.6, 7.0), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )
    mask_c = (comsol.time_us >= 0.0) & (comsol.time_us <= end_us)
    mask_e = (elmer.time_us >= 0.0) & (elmer.time_us <= end_us)
    axes[0].plot(
        comsol.time_us[mask_c], comsol.drop_uA[mask_c],
        label="COMSOL", color="#2166ac", linewidth=2.0,
    )
    axes[0].plot(
        elmer.time_us[mask_e], elmer.drop_uA[mask_e],
        label="Elmer AMGX / RTX 3060 Ti", color="#d95f02", linewidth=1.7,
    )
    axes[0].set_ylabel("TES current drop [µA]")
    axes[0].legend(frameon=False)
    axes[0].grid(True, alpha=0.25)
    axes[1].plot(aligned[:, 0], aligned[:, 3], color="#6a3d9a", linewidth=1.5)
    axes[1].axhline(0.0, color="#555555", linewidth=0.8)
    axes[1].set(xlabel="Time from pulse [µs]", ylabel="AMGX − COMSOL [µA]")
    axes[1].grid(True, alpha=0.25)
    fig.savefig(out / "current_comparison.png", dpi=220)
    fig.savefig(out / "current_comparison.svg")
    plt.close(fig)

    baseline = metrics["baseline_uA"]
    crossings = metrics["crossing_us_at_comsol_peak_fraction"]
    assert isinstance(baseline, dict) and isinstance(crossings, dict)
    summary = f"""# Optimized SinglePixel AMGX versus COMSOL

- Comparison window: 0–{end_us:.3f} µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: {baseline['COMSOL']:.6f} µA
- AMGX baseline: {baseline['AMGX']:.6f} µA ({baseline['AMGX_error_pct']:+.3f}%)
- Maximum absolute waveform difference: {metrics['max_abs_difference_uA']:.6f} µA at {metrics['max_abs_difference_time_us']:.3f} µs ({metrics['max_abs_difference_pct_comsol_peak']:.3f}% of COMSOL full-trace peak)
- RMSE: {metrics['rmse_uA']:.6f} µA ({metrics['rmse_pct_comsol_peak']:.3f}% of COMSOL full-trace peak)
- t10 (COMSOL / AMGX): {formatted_time(crossings['0.1']['COMSOL'])} / {formatted_time(crossings['0.1']['AMGX'])} µs
- t50 (COMSOL / AMGX): {formatted_time(crossings['0.5']['COMSOL'])} / {formatted_time(crossings['0.5']['AMGX'])} µs

The traces are compared after subtracting each model's own pre-pulse baseline.
"""
    (out / "summary.md").write_text(summary, encoding="utf-8", newline="\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--elmer", type=Path, default=DEFAULT_ELMER)
    parser.add_argument("--comsol", type=Path, default=ROOT / "docs" / "Single-Pixel.txt")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--end-us", type=float, default=100.0)
    args = parser.parse_args()
    metrics, aligned = compare(
        read_comsol(args.comsol), read_elmer(args.elmer), args.end_us
    )
    write_outputs(args.out, read_comsol(args.comsol), read_elmer(args.elmer), metrics, aligned)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

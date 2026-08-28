"""Plot the COMSOL, serial-direct, and four-rank MPI TES current pulse.

The curves are aligned at the physical 20.020 ms pulse and expressed as a
drop from each solver's own pre-pulse baseline. Shape-preserving PCHIP curves
make the nonuniform time grids readable.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator


ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Series:
    label: str
    time_us: np.ndarray
    current_uA: np.ndarray
    baseline_uA: float

    @property
    def drop_uA(self) -> np.ndarray:
        return self.baseline_uA - self.current_uA


def read_comsol(path: Path, pulse_ms: float, baseline_start_ms: float) -> Series:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("%"):
            continue
        values = line.split()
        rows.append([float(value) for value in values[:6]])
    table = np.asarray(rows)
    baseline = np.mean(
        table[(table[:, 0] >= baseline_start_ms) & (table[:, 0] < pulse_ms), 4]
    )
    return Series(
        "COMSOL",
        (table[:, 0] - pulse_ms) * 1.0e3,
        table[:, 4],
        float(baseline),
    )


def read_elmer(
    path: Path,
    label: str,
    pulse_ms: float,
    baseline_start_ms: float,
) -> Series:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    time_ms = np.asarray([float(row["time_s"]) * 1.0e3 for row in rows])
    current_uA = np.asarray([float(row["tes_current_A"]) * 1.0e6 for row in rows])
    baseline = np.mean(
        current_uA[(time_ms >= baseline_start_ms) & (time_ms < pulse_ms)]
    )
    return Series(
        label,
        (time_ms - pulse_ms) * 1.0e3,
        current_uA,
        float(baseline),
    )


def clipped(series: Series, start_us: float, end_us: float) -> tuple[np.ndarray, np.ndarray]:
    tolerance_us = 1.0e-2
    mask = (
        (series.time_us >= start_us - tolerance_us)
        & (series.time_us <= end_us + tolerance_us)
    )
    return series.time_us[mask], series.drop_uA[mask]


def smooth_grid(start_us: float, end_us: float) -> np.ndarray:
    pieces = [
        np.arange(start_us, min(0.0, end_us), 0.25),
        np.arange(max(start_us, 0.0), min(10.0, end_us), 0.05),
        np.arange(max(start_us, 10.0), min(120.0, end_us), 0.25),
        np.arange(max(start_us, 120.0), min(2000.0, end_us), 1.0),
        # Beyond a couple ms the pulse response is essentially settled; a 1 us
        # grid out to a 100+ ms end_us would be excessive, so widen the step.
        np.arange(max(start_us, 2000.0), end_us + 0.5, 20.0),
    ]
    nonempty = [piece for piece in pieces if piece.size]
    return np.unique(np.concatenate(nonempty + [np.asarray([start_us, end_us])]))


def log_grid(start_us: float, end_us: float) -> np.ndarray:
    if start_us <= 0.0 or end_us <= start_us:
        raise ValueError("log time axis requires 0 < start_us < end_us")
    key_times = np.asarray(
        [value for value in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0)
         if start_us <= value <= end_us],
        dtype=float,
    )
    return np.unique(
        np.concatenate([np.geomspace(start_us, end_us, 4000), key_times])
    )


def interpolator(series: Series, start_us: float, end_us: float) -> PchipInterpolator:
    time_us, drop_uA = clipped(series, start_us, end_us)
    return PchipInterpolator(time_us, drop_uA, extrapolate=False)


def crossing(t: np.ndarray, y: np.ndarray, level: float, rising: bool) -> float:
    """First time *y* crosses *level* (linear interpolation between samples)."""
    idx = np.where(y >= level if rising else y <= level)[0]
    if len(idx) == 0:
        return float("nan")
    i = int(idx[0])
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))


def pulse_metrics(grid: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Peak value plus 10-90% rise time and 90-10% fall time, in us."""
    valid = ~np.isnan(y)
    grid, y = grid[valid], y[valid]
    peak_i = int(np.argmax(np.abs(y)))
    peak = float(y[peak_i])
    rising = peak >= 0
    lo, hi = 0.1 * peak, 0.9 * peak
    t_rise, y_rise = grid[: peak_i + 1], y[: peak_i + 1]
    t_fall, y_fall = grid[peak_i:], y[peak_i:]
    t10 = crossing(t_rise, y_rise, lo, rising)
    t90 = crossing(t_rise, y_rise, hi, rising)
    t90f = crossing(t_fall, y_fall, hi, not rising)
    t10f = crossing(t_fall, y_fall, lo, not rising)
    return {
        "peak_uA": peak,
        "rise_time_10_90_us": t90 - t10,
        "fall_time_90_10_us": t10f - t90f,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comsol", type=Path, default=ROOT / "docs" / "Single-Pixel.txt")
    parser.add_argument(
        "--direct",
        type=Path,
        default=ROOT / "artifacts" / "series" / "tes_pulse_20ms_3x_series.csv",
    )
    parser.add_argument(
        "--mpi",
        type=Path,
        default=ROOT / "results" / "raw" / "legacy-root-output" / "tes_mpi_legacy_regression_series.csv",
    )
    parser.add_argument(
        "--elmer-label",
        default="Elmer MPI (4 ranks)",
        help="legend label for the Elmer series supplied with --mpi",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "artifacts" / "comparison" / "comsol_direct_mpi",
    )
    parser.add_argument(
        "--tag",
        default="",
        help="suffix inserted before the extension of every output file "
        "(e.g. --tag logx -> current_timeseries_comparison_logx.png), so "
        "multiple variants can share one --out directory without overwriting",
    )
    parser.add_argument("--pulse-ms", type=float, default=20.020)
    parser.add_argument("--baseline-start-ms", type=float, default=19.5)
    parser.add_argument("--start-us", type=float, default=-20.0)
    parser.add_argument("--end-us", type=float, default=600.0)
    parser.add_argument(
        "--skip-direct",
        action="store_true",
        help="omit the serial 'Elmer direct' reference curve from the plot/CSV",
    )
    parser.add_argument(
        "--linear-x",
        action="store_true",
        help="force a linear time axis even when end_us > 600 (default: symlog)",
    )
    parser.add_argument(
        "--log-x",
        action="store_true",
        help="use a true logarithmic post-pulse time axis; start_us must be positive",
    )
    args = parser.parse_args()
    if args.linear_x and args.log_x:
        parser.error("--linear-x and --log-x are mutually exclusive")
    if args.log_x and args.start_us <= 0.0:
        parser.error("--log-x requires --start-us > 0")

    series = [read_comsol(args.comsol, args.pulse_ms, args.baseline_start_ms)]
    if not args.skip_direct:
        series.append(
            read_elmer(args.direct, "Elmer direct", args.pulse_ms, args.baseline_start_ms)
        )
    series.append(
        read_elmer(args.mpi, args.elmer_label, args.pulse_ms, args.baseline_start_ms)
    )
    grid = log_grid(args.start_us, args.end_us) if args.log_x else smooth_grid(args.start_us, args.end_us)
    curves = {
        item.label: interpolator(item, args.start_us, args.end_us)(grid)
        for item in series
    }

    def tagged(name: str) -> Path:
        if not args.tag:
            return args.out / name
        stem, _, ext = name.rpartition(".")
        return args.out / f"{stem}_{args.tag}.{ext}"

    args.out.mkdir(parents=True, exist_ok=True)
    aligned_csv = tagged("baseline_corrected_current_smooth.csv")
    with aligned_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        header = ["time_from_pulse_us", "comsol_current_drop_uA"]
        if not args.skip_direct:
            header += ["direct_current_drop_uA", "direct_minus_comsol_uA"]
        header += ["mpi_current_drop_uA", "mpi_minus_comsol_uA"]
        writer.writerow(header)
        for index, time_us in enumerate(grid):
            comsol = curves["COMSOL"][index]
            mpi = curves[args.elmer_label][index]
            row = [time_us, comsol]
            if not args.skip_direct:
                direct = curves["Elmer direct"][index]
                row += [direct, direct - comsol]
            row += [mpi, mpi - comsol]
            writer.writerow(row)

    metrics = {label: pulse_metrics(grid, curve) for label, curve in curves.items()}
    comsol_metrics = metrics["COMSOL"]
    metrics_csv = tagged("pulse_metrics.csv")
    with metrics_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["series", "peak_uA", "rise_time_10_90_us", "fall_time_90_10_us",
             "peak_error_pct", "rise_time_error_pct", "fall_time_error_pct"]
        )
        for label, values in metrics.items():
            if label == "COMSOL":
                writer.writerow([label, values["peak_uA"], values["rise_time_10_90_us"],
                                  values["fall_time_90_10_us"], "", "", ""])
                continue
            errs = {
                key: (values[key] - comsol_metrics[key]) / comsol_metrics[key] * 100.0
                for key in ("peak_uA", "rise_time_10_90_us", "fall_time_90_10_us")
            }
            writer.writerow(
                [label, values["peak_uA"], values["rise_time_10_90_us"],
                 values["fall_time_90_10_us"], errs["peak_uA"],
                 errs["rise_time_10_90_us"], errs["fall_time_90_10_us"]]
            )
    print(f"\n{'series':<20}{'peak_uA':>12}{'rise_us':>12}{'fall_us':>12}"
          f"{'peak_err%':>12}{'rise_err%':>12}{'fall_err%':>12}")
    for label, values in metrics.items():
        if label == "COMSOL":
            print(f"{label:<20}{values['peak_uA']:>12.4f}{values['rise_time_10_90_us']:>12.4f}"
                  f"{values['fall_time_90_10_us']:>12.4f}{'':>12}{'':>12}{'':>12}")
        else:
            errs = {
                key: (values[key] - comsol_metrics[key]) / comsol_metrics[key] * 100.0
                for key in ("peak_uA", "rise_time_10_90_us", "fall_time_90_10_us")
            }
            print(f"{label:<20}{values['peak_uA']:>12.4f}{values['rise_time_10_90_us']:>12.4f}"
                  f"{values['fall_time_90_10_us']:>12.4f}{errs['peak_uA']:>12.4f}"
                  f"{errs['rise_time_10_90_us']:>12.4f}{errs['fall_time_90_10_us']:>12.4f}")

    colors = {
        "COMSOL": "#2166ac",
        "Elmer direct": "#1b9e77",
        args.elmer_label: "#d95f02",
    }
    fig, response_ax = plt.subplots(1, 1, figsize=(10.2, 5.6), constrained_layout=True)
    for item in series:
        response_ax.plot(
            grid,
            curves[item.label],
            label=item.label,
            color=colors[item.label],
            linewidth=2.0,
            linestyle=":" if item.label == "COMSOL" else "-",
        )

    if not args.log_x:
        response_ax.axvline(0.0, color="#555555", linestyle="--", linewidth=1.0)
    response_ax.grid(True, alpha=0.25)
    response_ax.set_xlim(args.start_us, args.end_us)
    # A long tail (e.g. out to COMSOL's full 180 ms) squashes the sub-ms
    # rise/decay into an unreadable sliver under a linear axis; symlog
    # keeps that detail visible while still showing the long relaxation.
    if args.log_x:
        response_ax.set_xscale("log")
    elif args.end_us > 600.0 and not args.linear_x:
        response_ax.set_xscale("symlog", linthresh=100.0, linscale=1.0)
    response_ax.set_ylabel("Current drop from baseline [uA]")
    response_ax.set_xlabel("Time from pulse [us]")
    if args.end_us <= 600.0:
        response_ax.set_xticks([-20.0, 0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0])
    # else: symlog's default locator already places sensible ticks.
    response_ax.legend(loc="upper left", frameon=False, ncols=3)
    response_ax.set_title(
        "TES current pulse: COMSOL vs Elmer MPI"
        if args.skip_direct
        else "TES current pulse: COMSOL vs Elmer direct vs MPI"
    )

    png_path = tagged("current_timeseries_comparison.png")
    svg_path = tagged("current_timeseries_comparison.svg")
    fig.savefig(png_path, dpi=240)
    fig.savefig(svg_path)
    plt.close(fig)

    print(png_path)
    print(svg_path)
    print(aligned_csv)
    print(metrics_csv)


if __name__ == "__main__":
    main()

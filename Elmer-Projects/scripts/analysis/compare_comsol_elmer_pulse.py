"""Baseline-corrected COMSOL/Elmer TES pulse comparison.

The COMSOL table is referenced to its physical pulse time.  The compact Elmer
probe starts at its pulse, so this script aligns the two at ``event_ms`` and
compares the TES-current decrease from each solver's own baseline.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def comsol_rows(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("%"):
            continue
        values = line.split()
        try:
            rows.append([float(value) for value in values[:5]])
        except ValueError:
            continue
    return np.asarray(rows, dtype=float)


def elmer_rows(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return (
        np.asarray([float(row["time_s"]) * 1.0e6 for row in rows]),
        np.asarray([float(row["tes_current_A"]) * 1.0e6 for row in rows]),
    )


def crossing(time_us: np.ndarray, response: np.ndarray, level: float, start: int, stop: int, rising: bool) -> float | None:
    for index in range(start + 1, stop + 1):
        if (response[index - 1] - level) * (response[index] - level) <= 0:
            if (response[index] >= response[index - 1]) == rising:
                fraction = (level - response[index - 1]) / (response[index] - response[index - 1])
                return float(time_us[index - 1] + fraction * (time_us[index] - time_us[index - 1]))
    return None


def metrics(time_us: np.ndarray, response: np.ndarray) -> dict[str, float | None]:
    peak = int(np.argmax(response))
    amplitude = float(response[peak])
    rise_10 = crossing(time_us, response, 0.1 * amplitude, 0, peak, True)
    rise_90 = crossing(time_us, response, 0.9 * amplitude, 0, peak, True)
    decay_1e = crossing(time_us, response, amplitude / np.e, peak, len(response) - 1, False)
    return {
        "peak_uA": amplitude,
        "peak_time_us": float(time_us[peak]),
        "rise_10_90_us": None if rise_10 is None or rise_90 is None else rise_90 - rise_10,
        "decay_1e_us": decay_1e,
        "window_ends_before_decay_1e": decay_1e is None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comsol", type=Path, default=Path("docs/Single-Pixel.txt"))
    parser.add_argument("--elmer", type=Path, default=Path("generated/validation_repart_x_mpi4_86step_series.csv"))
    parser.add_argument("--out", type=Path, default=Path("artifacts/comsol_elmer_pulse"))
    parser.add_argument("--event-ms", type=float, default=20.020)
    parser.add_argument("--baseline-start-ms", type=float, default=19.5)
    parser.add_argument("--elmer-event-ms", type=float, default=0.0,
                        help="Physical Elmer pulse time; use 0 for a pulse-at-start probe.")
    parser.add_argument("--elmer-baseline-uA", type=float, default=None,
                        help="Override Elmer baseline. By default, mean the pre-pulse window.")
    parser.add_argument("--elmer-baseline-start-ms", type=float, default=None,
                        help="Start of automatic Elmer baseline window; defaults to COMSOL baseline start.")
    args = parser.parse_args()

    comsol = comsol_rows(args.comsol)
    elmer_time_us, elmer_current_uA = elmer_rows(args.elmer)
    elmer_time_us = elmer_time_us - args.elmer_event_ms * 1.0e3
    baseline_mask = (comsol[:, 0] >= args.baseline_start_ms) & (comsol[:, 0] < args.event_ms)
    comsol_baseline_uA = float(np.mean(comsol[baseline_mask, 4]))
    comsol_time_us = (comsol[:, 0] - args.event_ms) * 1.0e3
    comsol_response_uA = comsol_baseline_uA - comsol[:, 4]
    if args.elmer_baseline_uA is None:
        baseline_start_ms = args.baseline_start_ms if args.elmer_baseline_start_ms is None else args.elmer_baseline_start_ms
        elmer_baseline_mask = (elmer_time_us >= (baseline_start_ms - args.elmer_event_ms) * 1.0e3) & (elmer_time_us < 0.0)
        if not np.any(elmer_baseline_mask):
            raise ValueError("No Elmer samples in the requested pre-pulse baseline window; supply --elmer-baseline-uA.")
        elmer_baseline_uA = float(np.mean(elmer_current_uA[elmer_baseline_mask]))
    else:
        elmer_baseline_uA = args.elmer_baseline_uA
    valid = (elmer_time_us >= comsol_time_us[0]) & (elmer_time_us <= comsol_time_us[-1])
    elmer_time_us, elmer_current_uA = elmer_time_us[valid], elmer_current_uA[valid]
    comsol_on_elmer_uA = np.interp(elmer_time_us, comsol_time_us, comsol_response_uA)
    elmer_response_uA = elmer_baseline_uA - elmer_current_uA

    comsol_post = comsol_time_us >= 0.0
    comsol_metric = metrics(comsol_time_us[comsol_post], comsol_response_uA[comsol_post])
    elmer_post = elmer_time_us >= 0.0
    elmer_metric = metrics(elmer_time_us[elmer_post], elmer_response_uA[elmer_post])
    comparison = elmer_response_uA - comsol_on_elmer_uA

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "baseline_corrected_current.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_from_pulse_us", "elmer_current_drop_uA", "comsol_current_drop_uA", "elmer_minus_comsol_uA"])
        writer.writerows(zip(elmer_time_us, elmer_response_uA, comsol_on_elmer_uA, comparison))

    result = {
        "baseline_uA": {"comsol": comsol_baseline_uA, "elmer": elmer_baseline_uA},
        "metrics": {"comsol": comsol_metric, "elmer": elmer_metric},
        "matched_window": {
            "start_us": float(elmer_time_us[0]),
            "end_us": float(elmer_time_us[-1]),
            "max_abs_difference_uA": float(np.max(np.abs(comparison))),
            "difference_at_end_uA": float(comparison[-1]),
        },
    }
    (args.out / "metrics.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, constrained_layout=True)
    axes[0].plot(comsol_time_us[comsol_post], comsol_response_uA[comsol_post], label="COMSOL", color="#2166ac", linewidth=2)
    axes[0].plot(elmer_time_us[elmer_post], elmer_response_uA[elmer_post], "o-", label="Elmer MPI (4 ranks)", color="#b2182b", markersize=2)
    axes[0].set_xlim(0.0, float(elmer_time_us[-1]))
    axes[0].set_ylabel("TES current drop from baseline [µA]")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[1].plot(elmer_time_us[elmer_post], comparison[elmer_post], "o-", color="#4d4d4d", markersize=2)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_xlabel("Time from pulse [µs]")
    axes[1].set_ylabel("Elmer − COMSOL [µA]")
    axes[1].grid(alpha=0.3)
    fig.savefig(args.out / "baseline_corrected_current.png", dpi=180)


if __name__ == "__main__":
    main()

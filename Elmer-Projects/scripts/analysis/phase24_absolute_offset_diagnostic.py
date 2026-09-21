"""Decompose Phase24 absolute-current offset from available series outputs.

This is read-only: it does not start Elmer.  It compares the completed native
HYPRE, same-path Direct MUMPS, and historical CPU/MUMPS Phase23 series using
the same pre-pulse baseline window and common post-event time grid.  Current,
temperature, resistance, and power are reported separately so a current-only
extraction offset can be distinguished from a state/solution offset.

Run from the repository root::

    python scripts/analysis/phase24_absolute_offset_diagnostic.py
"""
from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parents[2]
EVENT_S = 0.020020
BASELINE_START_US = -2.0
BASELINE_END_US = 0.0
CHECKPOINTS_US = (0.1, 0.5, 0.9)
DEFAULT_OUT = ROOT / "artifacts/phase24_absolute_offset_diagnostic_1us"
SERIES = {
    "hypre": ROOT / "results/case_phase24_short_hypre_tol5e7_1us/case_phase24_short_hypre_tol5e7_1us_series.csv",
    "mumps": ROOT / "results/case_phase24_short_same_path_mumps_1us/case_phase24_short_same_path_mumps_1us_series.csv",
    "cpu_reference": ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv",
}
LABELS = {"hypre": "HYPRE 5e-7", "mumps": "same Phase24 Direct MUMPS", "cpu_reference": "CPU/MUMPS Phase23"}
FIELDS = ("tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W")


@dataclass
class Series:
    label: str
    time_us: list[float]
    values: dict[str, list[float]]
    baseline: dict[str, float]


def interpolate(x: float, xs: list[float], ys: list[float]) -> float:
    index = bisect.bisect_left(xs, x)
    if index == 0:
        return ys[0]
    if index == len(xs):
        return ys[-1]
    if math.isclose(xs[index], x, rel_tol=0.0, abs_tol=1.0e-12):
        return ys[index]
    fraction = (x - xs[index - 1]) / (xs[index] - xs[index - 1])
    return ys[index - 1] + fraction * (ys[index] - ys[index - 1])


def load(path: Path, label: str) -> Series:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"time_s", *FIELDS}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path}: missing required columns {sorted(required)}")
    parsed = sorted((float(row["time_s"]) * 1.0e6 - EVENT_S * 1.0e6,
                     {field: float(row[field]) for field in FIELDS}) for row in rows)
    time_us: list[float] = []
    values = {field: [] for field in FIELDS}
    for time, row in parsed:
        if time_us and math.isclose(time, time_us[-1], rel_tol=0.0, abs_tol=1.0e-10):
            for field in FIELDS:
                values[field][-1] = row[field]
        else:
            time_us.append(time)
            for field in FIELDS:
                values[field].append(row[field])
    baseline = {
        field: fmean(values[field][i] for i, time in enumerate(time_us)
                     if BASELINE_START_US <= time < BASELINE_END_US)
        for field in FIELDS
    }
    return Series(label, time_us, values, baseline)


def pair_report(left: Series, right: Series) -> dict:
    common_start = max(0.0, left.time_us[0], right.time_us[0])
    common_end = min(left.time_us[-1], right.time_us[-1])
    grid = [time for time in left.time_us if common_start <= time <= common_end]
    metrics = {}
    for field in FIELDS:
        differences = []
        corrected = []
        for time in grid:
            lv = interpolate(time, left.time_us, left.values[field])
            rv = interpolate(time, right.time_us, right.values[field])
            differences.append(lv - rv)
            corrected.append((lv - left.baseline[field]) - (rv - right.baseline[field]))
        index_raw = max(range(len(differences)), key=lambda i: abs(differences[i]))
        index_corrected = max(range(len(corrected)), key=lambda i: abs(corrected[i]))
        metrics[field] = {
            "left_baseline": left.baseline[field],
            "right_baseline": right.baseline[field],
            "baseline_offset": left.baseline[field] - right.baseline[field],
            "max_raw_abs": abs(differences[index_raw]),
            "max_raw_time_us": grid[index_raw],
            "max_corrected_abs": abs(corrected[index_corrected]),
            "max_corrected_signed": corrected[index_corrected],
            "max_corrected_time_us": grid[index_corrected],
            "checkpoints": {
                str(checkpoint): {
                    "left": interpolate(checkpoint, left.time_us, left.values[field]),
                    "right": interpolate(checkpoint, right.time_us, right.values[field]),
                    "raw_difference": interpolate(checkpoint, left.time_us, left.values[field]) - interpolate(checkpoint, right.time_us, right.values[field]),
                    "baseline_corrected_difference": (
                        interpolate(checkpoint, left.time_us, left.values[field]) - left.baseline[field]
                    ) - (
                        interpolate(checkpoint, right.time_us, right.values[field]) - right.baseline[field]
                    ),
                }
                for checkpoint in CHECKPOINTS_US
                if common_start <= checkpoint <= common_end
            },
        }
    return {"left": left.label, "right": right.label, "common_window_us": [common_start, common_end], "metrics": metrics}


def write_report(payload: dict, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = [
        "# Phase24 absolute operating-point offset diagnostic",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        "Baseline window: `-2..0 us` relative to the event; common post-event window is reported per comparison.",
        "",
        "## Baselines",
        "",
        "| series | T [K] | I [uA] | R [ohm] | P [W] |",
        "|---|---:|---:|---:|---:|",
    ]
    for key, series in payload["series"].items():
        b = series["baseline"]
        lines.append(f"| {series['label']} | {b['tes_temperature_K']:.12g} | {b['tes_current_A']*1e6:.12g} | {b['tes_resistance_ohm']:.12g} | {b['tes_power_W']:.12g} |")
    lines += ["", "## Pairwise decomposition", "", "| comparison | quantity | baseline offset | max raw difference | max baseline-corrected difference |", "|---|---|---:|---:|---:|"]
    for name, pair in payload["comparisons"].items():
        for field, metric in pair["metrics"].items():
            unit = {"tes_temperature_K": "T [K]", "tes_current_A": "I [uA]", "tes_resistance_ohm": "R [ohm]", "tes_power_W": "P [W]"}[field]
            scale = 1.0e6 if field == "tes_current_A" else 1.0
            lines.append(f"| {name} | {unit} | {metric['baseline_offset']*scale:.9g} | {metric['max_raw_abs']*scale:.9g} | {metric['max_corrected_abs']*scale:.9g} |")
    lines += ["", "## Current checkpoints", "", "| comparison | 0.1 us corrected [uA] | 0.5 us corrected [uA] | 0.9 us corrected [uA] |", "|---|---:|---:|---:|"]
    for name, pair in payload["comparisons"].items():
        metric = pair["metrics"]["tes_current_A"]["checkpoints"]
        values = [metric[str(checkpoint)]["baseline_corrected_difference"] * 1.0e6 for checkpoint in CHECKPOINTS_US]
        lines.append(f"| {name} | {values[0]:.9g} | {values[1]:.9g} | {values[2]:.9g} |")
    lines += ["", "Interpretation: a large current baseline offset with small corrected T/R/P differences indicates readout or circuit operating-point extraction. A large corrected temperature or resistance difference indicates a state/solution mismatch. A large corrected current difference only, with matching state variables, indicates current extraction or circuit bookkeeping.", ""]
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    series = {}
    for key, path in SERIES.items():
        if not path.is_file():
            raise SystemExit(f"series not found: {path}")
        series[key] = load(path, LABELS[key])
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "event_s": EVENT_S,
        "baseline_window_us": [BASELINE_START_US, BASELINE_END_US],
        "series": {key: {"label": item.label, "path": str(SERIES[key]), "baseline": item.baseline, "samples": len(item.time_us)} for key, item in series.items()},
        "comparisons": {
            "hypre_vs_mumps": pair_report(series["hypre"], series["mumps"]),
            "hypre_vs_cpu_reference": pair_report(series["hypre"], series["cpu_reference"]),
            "mumps_vs_cpu_reference": pair_report(series["mumps"], series["cpu_reference"]),
        },
    }
    out = args.out if args.out.is_absolute() else ROOT / args.out
    write_report(payload, out)
    print(f"summary={out / 'summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

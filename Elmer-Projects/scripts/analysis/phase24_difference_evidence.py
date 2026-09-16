"""Create a reproducible, short-window Phase24 difference-evidence bundle.

The default invocation compares the already completed 1 ms native-HYPRE run
with the previous Phase23 CPU/MUMPS run.  It does not start Elmer.  The common
post-pulse interval is aligned by interpolation, so no manual tail/grep work
is needed.  The output contains raw and baseline-corrected differences,
checkpoint values, first-threshold crossings, input hashes, a Markdown report,
and a plot when matplotlib is available.

Run from the repository root:

    python scripts/analysis/phase24_difference_evidence.py

For another pair, pass ``--candidate`` and ``--reference``.  Paths may be
relative to the repository root.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CANDIDATE = ROOT / "results/case_phase24_hypre_cpu_pulse_1ms_5us/case_phase24_hypre_cpu_pulse_1ms_5us_series.csv"
DEFAULT_REFERENCE = ROOT / "results/case_p19_pulse_phase23_tight/case_p19_pulse_phase23_tight_series.csv"
DEFAULT_OUT = ROOT / "artifacts/comparison/phase24_difference_evidence"
DEFAULT_CHECKPOINTS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 5.0, 10.0)


def resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_series(path: Path) -> tuple[list[float], list[float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"time_s", "tes_current_A"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path}: expected CSV columns {sorted(required)}")
    pairs = sorted(
        (float(row["time_s"]) * 1.0e6, float(row["tes_current_A"]) * 1.0e6)
        for row in rows
    )
    # Duplicate output times do not carry additional information for this
    # comparison.  Keep the last value, matching the solver's final output.
    times: list[float] = []
    currents: list[float] = []
    for time_us, current_uA in pairs:
        if times and math.isclose(time_us, times[-1], rel_tol=0.0, abs_tol=1.0e-10):
            currents[-1] = current_uA
        else:
            times.append(time_us)
            currents.append(current_uA)
    return times, currents


def interpolate(x: float, xs: list[float], ys: list[float]) -> float:
    if x < xs[0] or x > xs[-1]:
        raise ValueError(f"requested time {x:g} is outside {xs[0]:g}..{xs[-1]:g}")
    index = bisect.bisect_left(xs, x)
    if index == 0:
        return ys[0]
    if index == len(xs):
        return ys[-1]
    if math.isclose(xs[index], x, rel_tol=0.0, abs_tol=1.0e-12):
        return ys[index]
    fraction = (x - xs[index - 1]) / (xs[index] - xs[index - 1])
    return ys[index - 1] + fraction * (ys[index] - ys[index - 1])


def mean_baseline(times: list[float], currents: list[float], start_us: float, end_us: float) -> float:
    values = [current for time, current in zip(times, currents) if start_us <= time < end_us]
    if not values:
        raise ValueError(f"no samples in baseline window [{start_us}, {end_us}) us")
    return fmean(values)


def first_crossing(times: list[float], values: list[float], threshold: float) -> float | None:
    for time, value in zip(times, values):
        if abs(value) >= threshold:
            return time
    return None


def manifest_snapshot(series_path: Path) -> dict:
    manifest_path = series_path.parent / "manifest.json"
    snapshot: dict[str, object] = {
        "series": str(series_path.resolve()),
        "series_sha256": sha256(series_path),
    }
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            snapshot["manifest"] = str(manifest_path.resolve())
            snapshot["manifest_sha256"] = sha256(manifest_path)
            for key in (
                "case", "exit_code", "solver", "elmer_prefix", "mesh",
                "started", "finished", "solver_completed", "runtime_artifacts_sha256",
            ):
                if key in manifest:
                    snapshot[key] = manifest[key]
        except (OSError, json.JSONDecodeError) as exc:
            snapshot["manifest_error"] = str(exc)
    return snapshot


def git_revision() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def fmt(value: float | None, digits: int = 9) -> str:
    return "n/a" if value is None else f"{value:.{digits}g}"


def write_plot(out: Path, rows: list[dict[str, float]], first_raw: float | None, first_corrected: float | None) -> str | None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    times = [row["time_from_event_us"] for row in rows]
    candidate = [row["candidate_current_uA"] for row in rows]
    reference = [row["reference_current_uA"] for row in rows]
    raw = [row["raw_difference_uA"] for row in rows]
    corrected = [row["baseline_corrected_difference_uA"] for row in rows]
    figure, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
    axes[0].plot(times, candidate, label="native HYPRE / candidate", linewidth=1.5)
    axes[0].plot(times, reference, label="CPU/MUMPS / reference", linewidth=1.5)
    axes[0].set_ylabel("TES current [µA]")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, raw, label="raw candidate − reference", linewidth=1.4)
    axes[1].plot(times, corrected, label="baseline-corrected difference", linewidth=1.4)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    for crossing, label, color in ((first_raw, "raw 0.1 µA", "#d62728"), (first_corrected, "corrected 0.1 µA", "#9467bd")):
        if crossing is not None:
            axes[1].axvline(crossing, color=color, linestyle="--", alpha=0.7, label=f"{label} @ {crossing:g} µs")
    axes[1].set_xlabel("time from pulse [µs]")
    axes[1].set_ylabel("difference [µA]")
    axes[1].legend()
    axes[1].grid(alpha=0.25)
    figure.savefig(out / "comparison.png", dpi=180)
    plt.close(figure)
    return "comparison.png"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--candidate-label", default="native HYPRE")
    parser.add_argument("--reference-label", default="CPU/MUMPS")
    parser.add_argument("--event-s", type=float, default=0.020020)
    parser.add_argument("--baseline-start-us", type=float, default=-2.0)
    parser.add_argument("--baseline-end-us", type=float, default=0.0)
    parser.add_argument("--threshold-uA", type=float, default=0.1)
    parser.add_argument("--relative-threshold-percent", type=float, default=0.1)
    parser.add_argument("--checkpoint-us", type=float, nargs="*", default=list(DEFAULT_CHECKPOINTS))
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    candidate_path = resolve(args.candidate)
    reference_path = resolve(args.reference)
    out = resolve(args.out)
    if not candidate_path.is_file():
        raise SystemExit(f"candidate CSV not found: {candidate_path}")
    if not reference_path.is_file():
        raise SystemExit(f"reference CSV not found: {reference_path}")
    if args.threshold_uA <= 0 or args.relative_threshold_percent <= 0:
        raise SystemExit("thresholds must be positive")

    candidate_time_abs, candidate_current = load_series(candidate_path)
    reference_time_abs, reference_current = load_series(reference_path)
    candidate_time = [time - args.event_s * 1.0e6 for time in candidate_time_abs]
    reference_time = [time - args.event_s * 1.0e6 for time in reference_time_abs]
    candidate_baseline = mean_baseline(candidate_time, candidate_current, args.baseline_start_us, args.baseline_end_us)
    reference_baseline = mean_baseline(reference_time, reference_current, args.baseline_start_us, args.baseline_end_us)

    common_start = max(0.0, candidate_time[0], reference_time[0])
    common_end = min(candidate_time[-1], reference_time[-1])
    if common_end <= common_start:
        raise SystemExit(f"no common post-event interval: {common_start}..{common_end} us")
    rows: list[dict[str, float]] = []
    for time in candidate_time:
        if common_start <= time <= common_end:
            candidate_value = interpolate(time, candidate_time, candidate_current)
            reference_value = interpolate(time, reference_time, reference_current)
            raw_difference = candidate_value - reference_value
            corrected_difference = (candidate_value - candidate_baseline) - (reference_value - reference_baseline)
            rows.append({
                "time_from_event_us": time,
                "candidate_current_uA": candidate_value,
                "reference_current_uA": reference_value,
                "raw_difference_uA": raw_difference,
                "baseline_corrected_difference_uA": corrected_difference,
                "raw_relative_percent": 100.0 * abs(raw_difference) / max(abs(reference_value), 1.0e-30),
            })

    raw_times = [row["time_from_event_us"] for row in rows]
    raw_values = [row["raw_difference_uA"] for row in rows]
    corrected_values = [row["baseline_corrected_difference_uA"] for row in rows]
    raw_crossing = first_crossing(raw_times, raw_values, args.threshold_uA)
    corrected_crossing = first_crossing(raw_times, corrected_values, args.threshold_uA)
    relative_crossing = first_crossing(
        raw_times,
        [row["raw_relative_percent"] for row in rows],
        args.relative_threshold_percent,
    )

    out.mkdir(parents=True, exist_ok=True)
    with (out / "aligned_difference.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    checkpoint_rows: list[dict[str, object]] = []
    for checkpoint in sorted(set(args.checkpoint_us)):
        if checkpoint < common_start or checkpoint > common_end:
            checkpoint_rows.append({"time_from_event_us": checkpoint, "available": False})
            continue
        candidate_value = interpolate(checkpoint, candidate_time, candidate_current)
        reference_value = interpolate(checkpoint, reference_time, reference_current)
        raw_difference = candidate_value - reference_value
        corrected_difference = (candidate_value - candidate_baseline) - (reference_value - reference_baseline)
        checkpoint_rows.append({
            "time_from_event_us": checkpoint,
            "available": True,
            "candidate_current_uA": candidate_value,
            "reference_current_uA": reference_value,
            "raw_difference_uA": raw_difference,
            "baseline_corrected_difference_uA": corrected_difference,
            "raw_relative_percent": 100.0 * abs(raw_difference) / max(abs(reference_value), 1.0e-30),
        })
    with (out / "checkpoints.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["time_from_event_us", "available", "candidate_current_uA", "reference_current_uA", "raw_difference_uA", "baseline_corrected_difference_uA", "raw_relative_percent"]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(checkpoint_rows)

    available_checkpoints = [row for row in checkpoint_rows if row.get("available")]
    max_raw = max(rows, key=lambda row: abs(row["raw_difference_uA"]))
    max_corrected = max(rows, key=lambda row: abs(row["baseline_corrected_difference_uA"]))
    plot_name = None if args.no_plot else write_plot(out, rows, raw_crossing, corrected_crossing)
    generated = datetime.now(timezone.utc).isoformat()
    summary = {
        "generated_utc": generated,
        "git_revision": git_revision(),
        "comparison": {
            "candidate_label": args.candidate_label,
            "reference_label": args.reference_label,
            "event_s": args.event_s,
            "baseline_window_us": [args.baseline_start_us, args.baseline_end_us],
            "common_post_event_window_us": [common_start, common_end],
            "common_rows": len(rows),
            "candidate_baseline_uA": candidate_baseline,
            "reference_baseline_uA": reference_baseline,
            "baseline_offset_candidate_minus_reference_uA": candidate_baseline - reference_baseline,
            "first_raw_abs_difference_ge_threshold_us": raw_crossing,
            "first_baseline_corrected_abs_difference_ge_threshold_us": corrected_crossing,
            "first_raw_relative_difference_ge_threshold_us": relative_crossing,
            "threshold_uA": args.threshold_uA,
            "relative_threshold_percent": args.relative_threshold_percent,
            "max_abs_raw_difference": {
                "time_from_event_us": max_raw["time_from_event_us"],
                "difference_uA": max_raw["raw_difference_uA"],
            },
            "max_abs_baseline_corrected_difference": {
                "time_from_event_us": max_corrected["time_from_event_us"],
                "difference_uA": max_corrected["baseline_corrected_difference_uA"],
            },
        },
        "inputs": {
            "candidate": manifest_snapshot(candidate_path),
            "reference": manifest_snapshot(reference_path),
        },
        "outputs": {
            "aligned_difference_csv": "aligned_difference.csv",
            "checkpoints_csv": "checkpoints.csv",
            "plot": plot_name,
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Phase24 difference evidence",
        "",
        f"Generated (UTC): `{generated}`",
        f"Git revision: `{summary['git_revision'] or 'unavailable'}`",
        "",
        f"- Candidate: `{args.candidate_label}` — `{candidate_path}`",
        f"- Reference: `{args.reference_label}` — `{reference_path}`",
        f"- Event: `{args.event_s:.15g} s`",
        f"- Common post-event window: `{common_start:.9g}..{common_end:.9g} µs` ({len(rows)} rows)",
        f"- Baselines: candidate `{candidate_baseline:.9f} µA`, reference `{reference_baseline:.9f} µA`, offset `{candidate_baseline - reference_baseline:+.9f} µA`",
        "",
        "## Difference onset",
        "",
        f"- First raw absolute difference ≥ `{args.threshold_uA:g} µA`: **{fmt(raw_crossing)} µs**",
        f"- First baseline-corrected absolute difference ≥ `{args.threshold_uA:g} µA`: **{fmt(corrected_crossing)} µs**",
        f"- First raw relative difference ≥ `{args.relative_threshold_percent:g}%`: **{fmt(relative_crossing)} µs**",
        f"- Maximum raw difference: `{max_raw['raw_difference_uA']:+.9f} µA` at `{max_raw['time_from_event_us']:.9g} µs`",
        f"- Maximum baseline-corrected difference: `{max_corrected['baseline_corrected_difference_uA']:+.9f} µA` at `{max_corrected['time_from_event_us']:.9g} µs`",
        "",
        "## Checkpoints",
        "",
        "| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in available_checkpoints:
        lines.append(
            f"| {row['time_from_event_us']:g} | {row['candidate_current_uA']:.9f} | {row['reference_current_uA']:.9f} | {row['raw_difference_uA']:+.9f} | {row['baseline_corrected_difference_uA']:+.9f} |"
        )
    lines.extend([
        "",
        "The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.",
        "The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.",
        "",
        "Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`" + (", `comparison.png`" if plot_name else ""),
        "",
    ])
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"out={out}")
    print(f"common_window_us={common_start:.9g}..{common_end:.9g}")
    print(f"first_raw_threshold_us={fmt(raw_crossing)}")
    print(f"first_corrected_threshold_us={fmt(corrected_crossing)}")
    print(f"first_relative_threshold_us={fmt(relative_crossing)}")
    print(f"summary={out / 'summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

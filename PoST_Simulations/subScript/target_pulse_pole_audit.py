"""Noise-blind, record-wise audit of target CH0/CH1 pulse time poles."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


TARGET_ROOT = Path(r"G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2")
RATE_HZ = 500000.0
SAMPLES = 100000
PRETRIGGER = 5000
MAX_RECORDS = 300
FIT_SECONDS = 0.08
FIT_STEP = 10


def _read(path: Path) -> np.ndarray:
    values = np.fromfile(path, dtype=np.float64, offset=4)
    if values.size != SAMPLES or not np.all(np.isfinite(values)):
        raise ValueError("invalid pulse record")
    return values


def _crossing(values: np.ndarray, peak: int, fraction: float, rising: bool) -> int | None:
    if rising:
        indices = np.flatnonzero(values[:peak] >= fraction * values[peak])
        return int(indices[0]) if len(indices) else None
    indices = np.flatnonzero(values[peak:] <= fraction * values[peak])
    return int(peak + indices[0]) if len(indices) else None


def _record_features(raw: np.ndarray) -> tuple[dict, np.ndarray] | None:
    baseline_window = raw[:PRETRIGGER - 500]
    baseline = float(np.mean(baseline_window))
    baseline_std = float(np.std(baseline_window))
    oriented = raw - baseline
    search = oriented[PRETRIGGER:]
    if np.max(np.abs(search)) <= 0.0:
        return None
    if abs(float(np.min(search))) > abs(float(np.max(search))):
        oriented = -oriented
    peak = int(PRETRIGGER + np.argmax(oriented[PRETRIGGER:]))
    amplitude = float(oriented[peak])
    if peak <= PRETRIGGER or amplitude <= 10.0 * max(baseline_std, np.finfo(float).tiny):
        return None
    rise20 = _crossing(oriented, peak, 0.2, True)
    rise90 = _crossing(oriented, peak, 0.9, True)
    decay90 = _crossing(oriented, peak, 0.9, False)
    decay10 = _crossing(oriented, peak, 0.1, False)
    if None in (rise20, rise90, decay90, decay10):
        return None
    tail = float(np.mean(oriented[-1000:]))
    return ({
        "peak_index": peak,
        "peak_time_s": peak / RATE_HZ,
        "amplitude_raw": amplitude,
        "baseline_raw": baseline,
        "baseline_std_raw": baseline_std,
        "tail_offset_raw": tail,
        "rise_20_to_90_s": (rise90 - rise20) / RATE_HZ,
        "decay_90_to_10_s": (decay10 - decay90) / RATE_HZ,
    }, oriented)


def _exp_model(order: int):
    def model(t, *parameters):
        amplitudes = np.asarray(parameters[:order])
        taus = np.asarray(parameters[order:2 * order])
        offset = parameters[-1]
        return offset + np.sum(amplitudes[:, None] * np.exp(-t[None, :] / taus[:, None]), axis=0)
    return model


def _fit_models(t: np.ndarray, values: np.ndarray) -> dict:
    fits = {}
    scale = max(float(np.max(values)), np.finfo(float).tiny)
    for order in (1, 2, 3):
        model = _exp_model(order)
        tau_guess = np.geomspace(2.0e-4, 2.0e-2, order)
        p0 = [scale / order] * order + list(tau_guess) + [float(values[-1])]
        lower = [0.0] * order + [1.0e-5] * order + [-scale]
        upper = [scale * 2.0] * order + [1.0] * order + [scale]
        try:
            fitted, _ = curve_fit(model, t, values, p0=p0, bounds=(lower, upper), maxfev=30000)
            residual = values - model(t, *fitted)
            rss = float(np.sum(residual**2))
            n = len(values)
            k = len(fitted)
            fits[str(order)] = {
                "order": order,
                "aic": float(n * np.log(max(rss / n, np.finfo(float).tiny)) + 2 * k),
                "bic": float(n * np.log(max(rss / n, np.finfo(float).tiny)) + k * np.log(n)),
                "rss": rss,
                "residual_lag1_correlation": float(np.corrcoef(residual[:-1], residual[1:])[0, 1]),
                "amplitudes": [float(x) for x in fitted[:order]],
                "time_constants_s": sorted(float(x) for x in fitted[order:2 * order]),
                "offset": float(fitted[-1]),
                "parameters": [float(x) for x in fitted],
            }
        except (RuntimeError, ValueError, FloatingPointError):
            fits[str(order)] = {"order": order, "fit_status": "failed"}
    return fits


def _evaluate_holdout(fits: dict, t: np.ndarray, values: np.ndarray) -> None:
    for row in fits.values():
        if row.get("fit_status") == "failed":
            continue
        model = _exp_model(row["order"])
        residual = values - model(t, *row["parameters"])
        rss = float(np.sum(residual**2))
        row["holdout_rss"] = rss
        row["holdout_rmse"] = float(np.sqrt(np.mean(residual**2)))
        row["holdout_residual_lag1_correlation"] = float(np.corrcoef(residual[:-1], residual[1:])[0, 1])


def _quantiles(rows: list[dict], key: str) -> dict:
    values = np.asarray([row[key] for row in rows], dtype=float)
    return {name: float(np.quantile(values, q)) for name, q in (("q05", 0.05), ("q50", 0.5), ("q95", 0.95))}


def _audit_channel(channel: str, output_dir: Path) -> dict:
    paths = sorted((TARGET_ROOT / f"{channel}_pulse/rawdata").glob(f"{channel}_*.dat"), key=lambda p: int(p.stem.split("_")[1]))[:MAX_RECORDS]
    rows = []
    aligned = []
    failures = 0
    for path in paths:
        try:
            result = _record_features(_read(path))
            if result is None:
                failures += 1
                continue
            features, waveform = result
            rows.append(features)
            relative = (np.arange(SAMPLES) - features["peak_index"]) / RATE_HZ
            aligned.append(np.interp(np.arange(0.0, FIT_SECONDS, FIT_STEP / RATE_HZ), relative, waveform, left=np.nan, right=np.nan))
        except (OSError, ValueError):
            failures += 1
    if len(rows) < 10:
        raise RuntimeError(f"{channel}: too few valid pulse records ({len(rows)})")
    split = max(5, int(0.7 * len(aligned)))
    aligned = np.asarray(aligned, dtype=float)
    train = np.nanmedian(aligned[:split], axis=0)
    holdout = np.nanmedian(aligned[split:], axis=0)
    time = np.arange(train.size, dtype=float) * FIT_STEP / RATE_HZ
    valid = np.isfinite(train) & np.isfinite(holdout)
    fits = _fit_models(time[valid], train[valid])
    _evaluate_holdout(fits, time[valid], holdout[valid])
    successful = [row for row in fits.values() if row.get("fit_status") != "failed" and np.isfinite(row.get("holdout_rss", np.nan))]
    selected = min(successful, key=lambda row: row["holdout_rss"]) if successful else None
    if selected is None:
        raise RuntimeError(f"{channel}: no decay model converged")
    selected_tau = selected["time_constants_s"][-1]
    return {
        "channel": channel,
        "source_directory": (TARGET_ROOT / f"{channel}_pulse/rawdata").as_posix(),
        "source_kind": "raw .dat pulse records only; no experimental noise records opened",
        "records_found": len(paths),
        "records_valid": len(rows),
        "records_rejected": failures,
        "feature_quantiles": {key: _quantiles(rows, key) for key in ("peak_time_s", "amplitude_raw", "baseline_std_raw", "rise_20_to_90_s", "decay_90_to_10_s")},
        "decay_model_fit": {"train_records": split, "heldout_records": len(aligned) - split, "fit_window_s": FIT_SECONDS, "downsample_step_samples": FIT_STEP, "models": fits, "selection_rule": "lowest held-out RSS among one/two/three exponential decay models"},
        "selected_slow_pole": {"time_constant_s": float(selected_tau), "frequency_Hz": float(1.0 / (2.0 * np.pi * selected_tau)), "model_order": int(selected["order"])},
        "record_wise_pulse_constraints": "descriptive time-domain constraints only; no physical parameter fit",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    channels = {channel: _audit_channel(channel, args.output_dir) for channel in ("CH0", "CH1")}
    result = {
        "stage": "A_pulse_audit",
        "noise_spectrum_read": False,
        "noise_records_read": False,
        "target_root": TARGET_ROOT.as_posix(),
        "sample_rate_Hz": RATE_HZ,
        "channels": channels,
        "strict_target_conclusion": "C — exact target physical case remains unidentified",
        "interpretation": "Pulse poles constrain observable time-domain combinations; they do not identify alpha, beta, C, G, L, or the target internal TES model individually.",
    }
    (args.output_dir / "target_pulse_pole_audit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = ["# Target pulse pole audit", "", "Noise records were not opened. The audit uses raw CH0/CH1 pulse `.dat` records only.", "", "| channel | valid records | selected order | slow tau [ms] | slow pole [Hz] |", "|---|---:|---:|---:|---:|"]
    for channel, row in channels.items():
        pole = row["selected_slow_pole"]
        lines.append(f"| {channel} | {row['records_valid']} | {pole['model_order']} | {pole['time_constant_s'] * 1e3:.4g} | {pole['frequency_Hz']:.4g} |")
    lines += ["", "Strict conclusion: **C — exact target physical case remains unidentified**.", "", "The fit is a pulse-shape pole audit, not a noise residual fit or physical-parameter optimizer."]
    (args.output_dir / "target_pulse_pole_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

"""Audit low-frequency stationarity and CH0/CH1 coherence after Stage B."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import signal

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Analyze_Experimental_Data.tes_analysis.noise_utils import accepted_noise_indices


TARGET_ROOT = Path(r"G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2")
RATE_HZ = 500000.0
SAMPLES = 100000
CUTOFF_HZ = 10000.0
POINTS_HZ = (10.0, 20.0, 100.0, 1000.0, 3000.0, 10000.0)


def _paths(channel: str) -> list[Path]:
    root = TARGET_ROOT / f"{channel}_noise/rawdata"
    return sorted(root.glob(f"{channel}_*.dat"), key=lambda p: int(p.stem.split("_")[1]))


def _read(path: Path) -> np.ndarray:
    values = np.fromfile(path, dtype=np.float64, offset=4)
    if values.size != SAMPLES or not np.all(np.isfinite(values)):
        raise ValueError(f"invalid record: {path}")
    return values


def _range_ok(values: np.ndarray) -> bool:
    return float(np.max(values) - np.min(values)) <= 0.04


def _accepted(paths: list[Path]) -> set[int]:
    def records():
        for path in paths:
            try:
                yield _read(path)
            except ValueError:
                yield np.empty(0, dtype=float)
    indices = accepted_noise_indices(
        records(), SAMPLES, RATE_HZ,
        cutoff=CUTOFF_HZ, remove_mean=True,
        accept_raw=_range_ok, accept_processed=_range_ok,
    )
    return set(indices)


def _record_metrics(ch0: np.ndarray, ch1: np.ndarray) -> dict:
    x0 = ch0 - np.mean(ch0)
    x1 = ch1 - np.mean(ch1)
    detrended0 = signal.detrend(x0, type="linear")
    detrended1 = signal.detrend(x1, type="linear")
    fft0 = np.fft.rfft(detrended0)
    fft1 = np.fft.rfft(detrended1)
    scale = RATE_HZ * SAMPLES
    p0 = np.abs(fft0) ** 2 / scale
    p1 = np.abs(fft1) ** 2 / scale
    pxy = fft0 * np.conjugate(fft1) / scale
    p0[1:-1] *= 2.0
    p1[1:-1] *= 2.0
    pxy[1:-1] *= 2.0
    frequencies = np.arange(p0.size) * RATE_HZ / SAMPLES
    baseline_n = 5000
    time = np.arange(baseline_n, dtype=float) / RATE_HZ
    baseline0 = ch0[:baseline_n]
    baseline1 = ch1[:baseline_n]
    slope0 = float(np.polyfit(time, baseline0, 1)[0])
    slope1 = float(np.polyfit(time, baseline1, 1)[0])
    return {"p0": p0, "p1": p1, "pxy": pxy, "frequencies": frequencies, "baseline_mean": (float(np.mean(baseline0)), float(np.mean(baseline1))), "baseline_slope": (slope0, slope1)}


def _point(values: np.ndarray, frequencies: np.ndarray, frequency: float) -> float:
    return float(values[int(np.argmin(np.abs(frequencies - frequency)))])


def _summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    return {"q05": float(np.quantile(values, 0.05)), "q50": float(np.quantile(values, 0.5)), "q95": float(np.quantile(values, 0.95)), "mean": float(np.mean(values))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ensemble", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not args.ensemble.exists():
        raise RuntimeError("stationarity audit requires the completed Stage-B ensemble first")
    ensemble = json.loads(args.ensemble.read_text(encoding="utf-8"))
    if ensemble.get("stage") != "B_model_ensemble":
        raise RuntimeError("stationarity audit requires a completed Stage-B ensemble")
    paths0, paths1 = _paths("CH0"), _paths("CH1")
    accepted0 = _accepted(paths0)
    accepted1 = _accepted(paths1)
    common = sorted(accepted0 & accepted1)
    records = []
    sum_p0 = sum_p1 = sum_pxy = None
    usable_common = []
    for index in common:
        try:
            metrics = _record_metrics(_read(paths0[index]), _read(paths1[index]))
        except ValueError:
            continue
        usable_common.append(index)
        if sum_p0 is None:
            sum_p0 = np.zeros_like(metrics["p0"])
            sum_p1 = np.zeros_like(metrics["p1"])
            sum_pxy = np.zeros_like(metrics["pxy"])
            frequencies = metrics["frequencies"]
        sum_p0 += metrics["p0"]
        sum_p1 += metrics["p1"]
        sum_pxy += metrics["pxy"]
        records.append({
            "index": index,
            "p10_ch0": _point(metrics["p0"], frequencies, 10.0),
            "p100_ch0": _point(metrics["p0"], frequencies, 100.0),
            "p10_ch1": _point(metrics["p1"], frequencies, 10.0),
            "p100_ch1": _point(metrics["p1"], frequencies, 100.0),
            "baseline_mean_ch0": metrics["baseline_mean"][0],
            "baseline_mean_ch1": metrics["baseline_mean"][1],
            "baseline_slope_ch0": metrics["baseline_slope"][0],
            "baseline_slope_ch1": metrics["baseline_slope"][1],
        })
    common = usable_common
    if not records:
        raise RuntimeError("no common accepted CH0/CH1 records")
    if len(records) < 4:
        result = {
            "stage": "C_stationarity_audit",
            "ensemble_completed_before_experiment_read": True,
            "source_kind": "target CH0/CH1 noise raw records using independent production acceptance",
            "accepted_counts": {"CH0_available": len(accepted0), "CH1_available": len(accepted1), "common": len(records)},
            "accepted_indices_common": common,
            "classification": "inconclusive",
            "reason": "independent CH0/CH1 production acceptance leaves fewer than four paired records for a stationarity or coherence estimate",
            "strict_target_conclusion": "C — exact target physical case remains unidentified",
        }
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "low_frequency_stationarity_audit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        (args.output_dir / "low_frequency_stationarity_audit.md").write_text("# Low-frequency stationarity and coherence audit\n\nClassification: **inconclusive**. Independent CH0/CH1 acceptance leaves fewer than four paired records.\n\nStrict conclusion: **C — exact target physical case remains unidentified**.\n", encoding="utf-8")
        return
    mean_p0, mean_p1, mean_pxy = sum_p0 / len(records), sum_p1 / len(records), sum_pxy / len(records)
    coherence = np.abs(mean_pxy) ** 2 / np.maximum(mean_p0 * mean_p1, np.finfo(float).tiny)
    logp0 = np.log10(np.maximum([row["p10_ch0"] for row in records], np.finfo(float).tiny))
    logp1 = np.log10(np.maximum([row["p10_ch1"] for row in records], np.finfo(float).tiny))
    def corr(key, values):
        x = np.asarray([row[key] for row in records], dtype=float)
        return float(np.corrcoef(x, values)[0, 1]) if np.std(x) > 0.0 and np.std(values) > 0.0 else 0.0
    half = len(records) // 2
    full = {"CH0": [_point(mean_p0, frequencies, f) for f in POINTS_HZ], "CH1": [_point(mean_p1, frequencies, f) for f in POINTS_HZ]}
    half_rows = {}
    for label, subset in (("first_half", records[:half]), ("second_half", records[half:])):
        half_rows[label] = {"CH0_p10_median": float(np.median([row["p10_ch0"] for row in subset])), "CH0_p100_median": float(np.median([row["p100_ch0"] for row in subset])), "CH1_p10_median": float(np.median([row["p10_ch1"] for row in subset])), "CH1_p100_median": float(np.median([row["p100_ch1"] for row in subset]))}
    half_ratios = {key: half_rows["second_half"][key] / half_rows["first_half"][key] for key in half_rows["first_half"]}
    coherence_points = {str(int(f)): _point(coherence, frequencies, f) for f in POINTS_HZ}
    low_band = coherence[(frequencies >= 10.0) & (frequencies <= 100.0)]
    low_coherence = float(np.median(low_band))
    low_coherence_range = [float(np.min(low_band)), float(np.max(low_band))]
    drift_ratio = max(abs(np.log(max(value, np.finfo(float).tiny))) for value in half_ratios.values())
    correlations = {"CH0_log10_P10_vs_baseline_mean": corr("baseline_mean_ch0", logp0), "CH0_log10_P10_vs_baseline_slope": corr("baseline_slope_ch0", logp0), "CH1_log10_P10_vs_baseline_mean": corr("baseline_mean_ch1", logp1), "CH1_log10_P10_vs_baseline_slope": corr("baseline_slope_ch1", logp1)}
    stationary = drift_ratio <= np.log(1.5) and max(abs(value) for value in correlations.values()) < 0.3
    if stationary and low_coherence >= 0.5:
        classification = "stationary_detector_like"
    elif stationary and low_coherence < 0.2:
        classification = "stationary_channel_specific"
    else:
        classification = "inconclusive"
    result = {
        "stage": "C_stationarity_audit",
        "ensemble_completed_before_experiment_read": True,
        "source_kind": "target CH0/CH1 noise raw records using the production acceptance predicates",
        "accepted_counts": {"CH0_available": len(accepted0), "CH1_available": len(accepted1), "common": len(records)},
        "accepted_indices_common": common,
        "frequency_points_Hz": list(POINTS_HZ),
        "mean_psd": full,
        "recordwise_psd_summary": {key: _summary([row[key] for row in records]) for key in ("p10_ch0", "p100_ch0", "p10_ch1", "p100_ch1")},
        "half_split": {"rows": half_rows, "second_over_first": half_ratios},
        "baseline_correlation": correlations,
        "coherence_points": coherence_points,
        "coherence_10_to_100_Hz_median": low_coherence,
        "coherence_10_to_100_Hz_range": low_coherence_range,
        "classification": classification,
        "stationarity_metrics": {"half_split_max_abs_log_ratio": drift_ratio, "max_abs_baseline_correlation": max(abs(value) for value in correlations.values()), "stationary_thresholds": {"half_split_ratio": 1.5, "baseline_correlation": 0.3}},
        "classification_rule": "stationary_detector_like if stationarity passes and 10-100 Hz coherence median >=0.5; stationary_channel_specific if stationarity passes and median <0.2; otherwise inconclusive",
        "strict_target_conclusion": "C — exact target physical case remains unidentified",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "low_frequency_stationarity_audit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = ["# Low-frequency stationarity and coherence audit", "", f"Common accepted CH0/CH1 records: **{len(records)}**.", "", f"Classification: **{classification}**.", "", f"Median CH0/CH1 coherence over 10–100 Hz: `{low_coherence:.6g}`.", "", "The audit uses the production acceptance predicates and does not fit noise residuals or adjust any model parameter.", "", "Strict conclusion: **C — exact target physical case remains unidentified**."]
    (args.output_dir / "low_frequency_stationarity_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

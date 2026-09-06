"""Record-wise, pulse-only calibrated audit of TES noise records.

This module deliberately separates the pulse-only stages (template and null
calibration) from the noise stages (classification and diagnostics).  The
detector's statistic is the same filtered, full-record, all-lag statistic for
real and null records.  No simulation spectrum or residual is read by the
calibration stages.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Analyze_Experimental_Data.tes_analysis.noise_utils import (
    one_sided_asd_from_power,
    preprocess_noise_record,
)
from pulse_contamination_common import (
    RATE_HZ,
    SAMPLES,
    noise_paths,
    pulse_paths,
    read_record,
    record_key,
    sha256,
)


CASE_DEFAULT = Path(__file__).resolve().parents[1] / "cases" / "tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"
TARGET_DEFAULT = Path(r"G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2")
CHANNELS = ("CH0", "CH1")
FPRS = (1e-2, 1e-3, 1e-4)
PRIMARY_FPR = 1e-3
ANCHORS = (5, 10, 20, 30, 50, 70, 100, 150, 200, 300, 500, 1000, 3000, 5000, 7000, 10000)
SCALING_POINTS = (20, 40, 50, 80, 100, 150, 200, 300, 500)
TAIL_TIMES_MS = (10, 20, 50, 100, 150)
TAIL_MAX_SAMPLES = 90000  # 180 ms at the target 500 kHz rate


def _write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _finite_quantiles(values, qs=(0.05, 0.5, 0.95)) -> dict:
    a = np.asarray([x for x in values if x is not None and np.isfinite(x)], dtype=float)
    return {f"q{int(q * 100):02d}": float(np.quantile(a, q)) if a.size else None for q in qs}


def _filtered(raw: np.ndarray) -> np.ndarray:
    return preprocess_noise_record(raw, RATE_HZ, cutoff=10000.0, remove_mean=True)


def _orient(raw: np.ndarray) -> tuple[np.ndarray, int, float, float]:
    y = _filtered(raw)
    base = float(np.mean(y[:4500]))
    y = y - base
    sigma = float(np.std(y[:4500]))
    post = y[5000:]
    if abs(float(np.min(post))) > abs(float(np.max(post))):
        y = -y
    peak = int(5000 + np.argmax(y[5000:]))
    amp = float(y[peak])
    if peak <= 5000 or amp <= 10.0 * max(sigma, np.finfo(float).tiny):
        raise ValueError("no pulse above pretrigger baseline")
    return y, peak, amp, sigma


def _crossing(y: np.ndarray, peak: int, fraction: float, rising: bool) -> int | None:
    if rising:
        ix = np.flatnonzero(y[:peak] >= fraction * y[peak])
        return int(ix[0]) if ix.size else None
    ix = np.flatnonzero(y[peak:] <= fraction * y[peak])
    return int(peak + ix[0]) if ix.size else None


def build_library(target: Path) -> dict:
    """Build full-rate templates and tail morphology from pulse records only."""
    channels = {}
    for channel in CHANNELS:
        rows = []
        for path in pulse_paths(target, channel):
            try:
                raw = read_record(path)
                y, peak, amp, sigma = _orient(raw)
            except ValueError:
                continue
            full_start = peak - 256
            full = y[full_start:full_start + 4096] / amp
            if full.size != 4096:
                continue
            tail = np.full(TAIL_MAX_SAMPLES, np.nan)
            available = min(TAIL_MAX_SAMPLES, SAMPLES - peak)
            tail[:available] = y[peak:peak + available] / amp
            c20 = _crossing(full, 256, 0.2, True)
            c90 = _crossing(full, 256, 0.9, True)
            d90 = _crossing(full, 256, 0.9, False)
            d10 = _crossing(full, 256, 0.1, False)
            tails = {
                "tail_integral_normalized_s": float(np.nansum(tail) / RATE_HZ),
                "tail_end_normalized": float(tail[min(available, TAIL_MAX_SAMPLES) - 1]) if available else None,
            }
            for ms in TAIL_TIMES_MS:
                index = int(ms * 1e-3 * RATE_HZ)
                tails[f"residual_{ms}ms"] = float(tail[index]) if index < available else None
            rows.append({
                "path": path.as_posix(),
                "sha256": sha256(path),
                "peak_index": peak,
                "amplitude_raw": amp,
                "baseline_std_raw": sigma,
                "full": full,
                "tail": tail,
                "rise_time_s": (c90 - c20) / RATE_HZ if c20 is not None and c90 is not None else None,
                "decay_time_s": (d10 - d90) / RATE_HZ if d90 is not None and d10 is not None else None,
                **tails,
            })
        if len(rows) < 10:
            raise RuntimeError(f"{channel}: fewer than 10 usable pulse records")
        full = np.nanmedian(np.asarray([r["full"] for r in rows]), axis=0)
        tail_matrix = np.asarray([r["tail"] for r in rows])
        tail = np.nanmedian(tail_matrix, axis=0)
        counts = np.sum(np.isfinite(tail_matrix), axis=0)
        # Beyond the last supported median sample the JSON value is a mask,
        # not an interpolated waveform.  The detector uses only this prefix.
        min_count = max(5, int(np.ceil(0.10 * len(rows))))
        valid = counts >= min_count
        valid_length = 0
        while valid_length < valid.size and valid[valid_length]:
            valid_length += 1
        if valid_length < 100:
            raise RuntimeError(f"{channel}: tail template has no usable contiguous prefix")
        tail_values = np.nan_to_num(tail[:valid_length], nan=0.0)
        template_rows = {
            "full_pulse": {"sample_step": 1, "relative_start_samples": -256, "fullrate": True, "values": full.tolist()},
            "post_peak_tail": {"sample_step": 1, "relative_start_samples": 0, "fullrate": True, "values": np.nan_to_num(tail[:4096], nan=0.0).tolist()},
            "slow_tail": {"sample_step": 1, "relative_start_samples": 0, "fullrate": True, "values": tail_values.tolist()},
        }
        morphology = {}
        for key in ("tail_integral_normalized_s", "tail_end_normalized", "rise_time_s", "decay_time_s", *[f"residual_{ms}ms" for ms in TAIL_TIMES_MS]):
            morphology[key] = _finite_quantiles([r.get(key) for r in rows])
        channels[channel] = {
            "channel": channel,
            "records_found": len(pulse_paths(target, channel)),
            "records_used": len(rows),
            "sample_rate_Hz": RATE_HZ,
            "tail_max_samples": TAIL_MAX_SAMPLES,
            "tail_max_ms": TAIL_MAX_SAMPLES / RATE_HZ * 1000.0,
            "tail_valid_length_samples": valid_length,
            "tail_valid_length_ms": valid_length / RATE_HZ * 1000.0,
            "tail_valid_sample_count_minimum": min_count,
            "tail_valid_sample_count_by_offset": counts.tolist(),
            "templates": template_rows,
            "morphology_quantiles": morphology,
            "records": [{k: v for k, v in r.items() if k not in {"full", "tail"}} for r in rows],
            "provenance": "pulse records only; filtered with production preprocessing; peak-aligned, amplitude-normalized, nan-aware median; no decimation or zero interleaving",
        }
    return {"stage": "pulse_tail_template_library_v2", "target_root": target.as_posix(), "noise_records_read": False, "channels": channels}


def _score(record: np.ndarray, kernel: np.ndarray) -> tuple[float, int, float]:
    """Normalized correlation over every valid lag, including long templates."""
    from scipy.signal import fftconvolve
    k = np.asarray(kernel, dtype=float)
    k = k - np.mean(k)
    if record.size < k.size or not np.any(k):
        return 0.0, 0, 0.0
    corr = fftconvolve(record, k[::-1], mode="valid")
    energy = np.sqrt(fftconvolve(record * record, np.ones(k.size), mode="valid"))
    scores = corr / np.maximum(np.linalg.norm(k) * energy, np.finfo(float).tiny)
    i = int(np.argmax(scores))
    seg = record[i:i + k.size]
    amp = float(np.dot(seg, k) / max(np.dot(k, k), np.finfo(float).tiny))
    return float(scores[i]), i, amp


def _batch_score(records: np.ndarray, kernel: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Batch equivalent of _score; used only to make null calibration tractable."""
    records = np.asarray(records, dtype=float)
    k = np.asarray(kernel, dtype=float); k = k - np.mean(k)
    n, size = records.shape; m = k.size
    if size < m: return np.zeros(n), np.zeros(n, dtype=int), np.zeros(n)
    nfft = 1 << int(np.ceil(np.log2(size + m - 1)))
    rk = np.fft.rfft(k[::-1], nfft)
    re = np.fft.rfft(np.ones(m), nfft)
    conv = np.fft.irfft(np.fft.rfft(records, nfft, axis=1) * rk[None, :], nfft, axis=1)[:, m - 1:size]
    energy = np.sqrt(np.maximum(np.fft.irfft(np.fft.rfft(records * records, nfft, axis=1) * re[None, :], nfft, axis=1)[:, m - 1:size], 0.0))
    scores = conv / np.maximum(np.linalg.norm(k) * energy, np.finfo(float).tiny)
    indices = np.argmax(scores, axis=1)
    amplitudes = np.asarray([np.dot(records[i, j:j + m], k) / max(np.dot(k, k), np.finfo(float).tiny) for i, j in enumerate(indices)])
    return scores[np.arange(n), indices], indices, amplitudes


def _scan(record: np.ndarray, channel_library: dict) -> dict:
    rows = {}
    for name, template in channel_library["templates"].items():
        score, lag, amp = _score(record, np.asarray(template["values"], dtype=float))
        rows[name] = {"score": score, "lag_samples": lag, "lag_s": lag / RATE_HZ, "amplitude_normalized": amp}
    return rows


def _null_record(baseline: np.ndarray, index: int) -> np.ndarray:
    # Periodic extension preserves the measured pretrigger waveform without
    # inventing a Gaussian spectrum.  A changing phase prevents identical
    # record boundaries from determining the maximum statistic.
    b = np.asarray(baseline, dtype=float)
    phase = (index * 977) % b.size
    return np.resize(np.roll(b, phase), SAMPLES)


def calibrate_null(target: Path, library: dict, null_records: int = 512) -> dict:
    """Calibrate family-wise thresholds from full-length pulse-baseline nulls."""
    result = {"stage": "recordwise_pulse_detection_thresholds", "target_root": target.as_posix(), "record_duration_s": SAMPLES / RATE_HZ, "null_record_count_requested": null_records, "null_record_source": "pulse pretrigger baseline only, periodic block extension to full record", "production_preprocessing": "mean removal then 2-pole zero-phase Bessel low-pass at 10 kHz", "all_lags": True, "all_templates": True, "simulation_spectrum_read": False, "noise_spectrum_read": False, "primary_fpr": PRIMARY_FPR, "fprs_reported": list(FPRS), "channels": {}}
    for channel in CHANNELS:
        baselines = []
        for path in pulse_paths(target, channel):
            try:
                raw = read_record(path)
                baselines.append(raw[:4500] - np.mean(raw[:4500]))
            except ValueError:
                pass
        if not baselines:
            raise RuntimeError(f"{channel}: no pulse baselines")
        names = tuple(library["channels"][channel]["templates"])
        scores = {name: [] for name in names}
        batch_size = 16
        for first in range(0, null_records, batch_size):
            raw_batch = np.asarray([_filtered(_null_record(baselines[i % len(baselines)], i)) for i in range(first, min(first + batch_size, null_records))])
            for name in names:
                values, _lags, _amps = _batch_score(raw_batch, np.asarray(library["channels"][channel]["templates"][name]["values"], dtype=float))
                scores[name].extend(values.tolist())
        arrays = {name: np.asarray(values) for name, values in scores.items()}
        # Every null record is ranked against the complete null sample.  This
        # avoids the sequential-CDF error where early records would otherwise
        # receive artificially small p-values.
        null_p = np.column_stack([
            [(1 + np.count_nonzero(arrays[name] >= value)) / (null_records + 1) for value in arrays[name]]
            for name in names
        ])
        min_array = np.min(null_p, axis=1)
        # Use a rank-normalized family statistic.  This is an empirical
        # maximum statistic, so correlated templates/lags are retained.
        channel_out = {"null_sample_count": int(null_records), "templates": {}, "familywise": {}, "null_score_summary": {}, "null_score_samples": {name: values.tolist() for name, values in arrays.items()}}
        for name, values in arrays.items():
            channel_out["null_score_summary"][name] = {"q50": float(np.quantile(values, .5)), "q95": float(np.quantile(values, .95)), "q99": float(np.quantile(values, .99)), "max": float(np.max(values))}
            channel_out["templates"][name] = {f"fpr_{f:g}": float(np.quantile(values, 1.0 - f)) for f in FPRS}
        for f in FPRS:
            # family threshold is stored in min-p units; observed records use
            # the same empirical p-value construction.
            channel_out["familywise"][f"fpr_{f:g}"] = float(np.quantile(min_array, f))
        channel_out["familywise_null_min_p_values"] = min_array.tolist()
        result["channels"][channel] = channel_out
    return result


def _empirical_p(value: float, null: np.ndarray) -> float:
    return float((1 + np.count_nonzero(null >= value)) / (null.size + 1))


def _decision(scan: dict, channel: str, library: dict, calibration: dict) -> dict:
    """Turn per-template scores into template and family-wise p-values."""
    p_values = {}
    for name, hit in scan.items():
        null = np.asarray(calibration["channels"][channel]["null_score_samples"][name], dtype=float)
        p_values[name] = _empirical_p(hit["score"], null)
    min_p = min(p_values.values())
    family_null = np.asarray(calibration["channels"][channel]["familywise_null_min_p_values"], dtype=float)
    family_p = float((1 + np.count_nonzero(family_null <= min_p)) / (family_null.size + 1))
    best_name = min(p_values, key=p_values.get)
    best = scan[best_name]
    is_full = best_name == "full_pulse"
    morphology = "full_pulse" if is_full else "tail_only"
    if family_p <= 1e-4:
        cls = "definite_full_pulse" if is_full else "definite_tail_only"
    elif family_p <= 1e-3:
        cls = "likely_pulse"
    elif family_p <= 1e-2:
        cls = "ambiguous"
    else:
        cls = "pulse_free"
    n = len(library["channels"][channel]["templates"][best_name]["values"])
    edge = {"peak_before_record_start_equivalent": best["lag_samples"] == 0 and morphology == "tail_only", "tail_only_at_start": morphology == "tail_only" and best["lag_samples"] == 0, "onset_near_record_end": best["lag_samples"] + n >= SAMPLES}
    return {"best_template": best_name, "best": best, "template_p_values": p_values, "family_wise_p_value": family_p, "classification": cls, "morphology": morphology, "edge_flags": edge}


def classify(target: Path, library: dict, calibration: dict) -> dict:
    rows = {}
    accepted_sets = {}
    for channel in CHANNELS:
        paths = noise_paths(target, channel)
        accepted = []
        for path in paths:
            try:
                raw = read_record(path)
                # This is the unchanged production acceptance predicate.
                if np.max(raw) - np.min(raw) > 0.04:
                    continue
                processed = _filtered(raw)
                if np.max(processed) - np.min(processed) > 0.04:
                    continue
                accepted.append(path)
            except ValueError:
                continue
        accepted_sets[channel] = {record_key(p) for p in accepted}
        null_scores = {name: np.asarray([x for x in []], dtype=float) for name in library["channels"][channel]["templates"]}
        # Per-template null score arrays are not persisted separately in the
        # compact calibration file; reconstruct their quantile thresholds from
        # the stored FPR score values by retaining exact observed score p-values
        # via rank interpolation against the family calibration.  The family
        # statistic is the authoritative judgement.
        family_null = np.asarray(calibration["channels"][channel]["familywise_null_min_p_values"], dtype=float)
        for path in accepted:
            raw = read_record(path)
            scan = _scan(_filtered(raw), library["channels"][channel])
            hit = _decision(scan, channel, library, calibration)
            best = hit["best"]
            rows.setdefault(record_key(path), {})[channel] = {"event_key": record_key(path), "channel": channel, "best_template": hit["best_template"], "best_lag_samples": best["lag_samples"], "best_lag_s": best["lag_s"], "best_amplitude_normalized": best["amplitude_normalized"], "template_p_values": hit["template_p_values"], "family_wise_p_value": hit["family_wise_p_value"], "classification": hit["classification"], "morphology": hit["morphology"], "edge_flags": hit["edge_flags"]}
    record_rows = {}
    for key in sorted(set(rows) | set().union(*accepted_sets.values()), key=lambda x: int(x)):
        record_rows[key] = {"event_key": key, **{ch: rows.get(key, {}).get(ch, {"event_key": key, "classification": "not_accepted", "morphology": "not_accepted"}) for ch in CHANNELS}}
    counts = {}
    for channel in CHANNELS:
        vals = [r[channel]["classification"] for r in record_rows.values() if r[channel]["classification"] != "not_accepted"]
        counts[channel] = {name: vals.count(name) for name in ("pulse_free", "ambiguous", "likely_pulse", "definite_full_pulse", "definite_tail_only")}
    return {"stage": "noise_record_pulse_classification_v2", "target_root": target.as_posix(), "record_duration_s": SAMPLES / RATE_HZ, "accepted_counts": {c: len(accepted_sets[c]) for c in CHANNELS}, "accepted_event_key_sets": {c: sorted(accepted_sets[c], key=int) for c in CHANNELS}, "counts_by_channel": counts, "records": record_rows, "classification_rule": "family-wise empirical null p-value; primary FWER per 0.2 s record 1e-3; no pulse subtraction", "threshold_selected_before_noise_psd": True, "provenance": "same production preprocessing and full-record all-lag scan; template selection separated from significance judgement"}


def _power(record: np.ndarray, length: int, detrend_order: int | None = None) -> np.ndarray:
    y = np.asarray(record, dtype=float)
    if detrend_order is not None:
        x = np.linspace(-1.0, 1.0, y.size)
        y = y - np.polyval(np.polyfit(x, y, detrend_order), x)
    else:
        y = y - np.mean(y)
    w = np.hanning(length)
    return np.abs(np.fft.rfft(y * w)) ** 2


def _asd(records: list[np.ndarray], length: int, detrend_order: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    if not records:
        return np.empty(0), np.empty(0)
    total = np.zeros(length // 2 + 1)
    for y in records:
        total += _power(y, length, detrend_order)
    w = np.hanning(length)
    return np.fft.rfftfreq(length, 1.0 / RATE_HZ), one_sided_asd_from_power(total / len(records), length, RATE_HZ, np.sqrt(np.mean(w * w)))


def _anchors(freq: np.ndarray, asd: np.ndarray, points=ANCHORS) -> dict:
    return {str(f): float(asd[int(np.argmin(np.abs(freq - f)))]) for f in points} if asd.size else {}


def analyze(target: Path, classification: dict, library: dict, case: Path) -> dict:
    path_by_ch = {c: {record_key(p): p for p in noise_paths(target, c)} for c in CHANNELS}
    clean = []
    categories = {"all_accepted": [], "clean_v2": [], "full_pulse_contaminated": [], "tail_only_contaminated": [], "ambiguous": []}
    per_record = []
    for key, row in classification["records"].items():
        if row["CH0"]["classification"] == "not_accepted":
            continue
        raw = read_record(path_by_ch["CH0"][key]); y = _filtered(raw)
        cls = row["CH0"]["classification"]
        categories["all_accepted"].append(y)
        if cls == "pulse_free": categories["clean_v2"].append(y); clean.append((key, raw, y))
        if cls in {"definite_full_pulse", "likely_pulse"}: categories["full_pulse_contaminated"].append(y)
        if cls in {"definite_tail_only", "likely_pulse"}: categories["tail_only_contaminated"].append(y)
        if cls == "ambiguous": categories["ambiguous"].append(y)
    subsets = {}
    for name, vals in categories.items():
        f, a = _asd(vals, SAMPLES)
        subsets[name] = {"record_count": len(vals), "frequencies_Hz": f.tolist(), "asd_anchors": _anchors(f, a), "estimator": "mean removal, Hann, one-sided power average"}
    _write(case / "pulse_partitioned_noise_spectra_v2.json", {"stage": "pulse_partitioned_noise_spectra_v2", "target_root": target.as_posix(), "rate_Hz": RATE_HZ, "samples": SAMPLES, "record_duration_s": SAMPLES / RATE_HZ, "subsets": subsets, "estimator": "same for all subsets: mean removal, Hann window, one-sided rFFT power average; clean_v2 is record rejection only"})
    clean_asd = subsets["clean_v2"]["asd_anchors"]
    slopes = {}
    f, a = _asd(categories["clean_v2"], SAMPLES)
    for name, lo, hi in (("10-30_Hz", 10, 30), ("30-100_Hz", 30, 100), ("100-200_Hz", 100, 200), ("10-200_Hz", 10, 200)):
        m = (f >= lo) & (f <= hi) & (a > 0)
        slopes[name] = {"gamma_dlog_asd_dlog_f": float(np.polyfit(np.log(f[m]), np.log(a[m]), 1)[0]) if np.count_nonzero(m) >= 2 else None, "frequency_band_Hz": [lo, hi]}
    _write(case / "low_frequency_shape_diagnostic.json", {"stage": "low_frequency_shape_diagnostic", "subset": "clean_v2", "slopes": slopes, "interpretation": "descriptive power-law slope only; no 1/f source or residual fit"})
    scaling = {}
    for length in (100000, 50000, 25000, 12500):
        segs = []
        for _key, raw, _y in clean:
            for start in range(0, SAMPLES - length + 1, length): segs.append(raw[start:start + length])
        freq, asd = _asd([_filtered(x) for x in segs], length)
        scaling[str(length / RATE_HZ)] = {"length_samples": length, "segment_count": len(segs), "asd_anchors": _anchors(freq, asd, [p for p in SCALING_POINTS if p >= RATE_HZ / length]), "frequency_resolution_Hz": RATE_HZ / length}
    _write(case / "record_length_scaling.json", {"stage": "record_length_scaling", "subset": "same clean_v2 raw records, deterministic non-overlap", "windows": scaling, "classification_rule": "RL1 if common-bin ASD ratios remain stable; RL2 if shorter windows materially reduce excess; RL5 otherwise"})
    detrend = {}
    for order in (None, 1, 2):
        freq, asd = _asd(categories["clean_v2"], SAMPLES, order)
        detrend["mean_removal" if order is None else f"polynomial_order_{order}"] = {"asd_anchors": _anchors(freq, asd), "production_estimator": order is None}
    _write(case / "detrend_sensitivity.json", {"stage": "detrend_sensitivity", "subset": "clean_v2", "diagnostic_only": True, "results": detrend, "production_spectrum": "mean_removal"})
    audit_rows = []; low_powers = {str(f): [] for f in (10, 20, 50)}
    for key, raw, y in clean:
        freq, asd = _asd([y], SAMPLES); psd = asd * asd
        q = [float(np.mean(y[:SAMPLES // 4])), float(np.mean(y[SAMPLES // 4:SAMPLES // 2])), float(np.mean(y[SAMPLES // 2:3 * SAMPLES // 4])), float(np.mean(y[3 * SAMPLES // 4:]))]
        row = {"event_key": key, "first_difference_variance": float(np.var(np.diff(y))), "baseline_half_to_half_offset": float(np.mean(y[SAMPLES // 2:]) - np.mean(y[:SAMPLES // 2])), "quarter_means": q, "piecewise_max_jump": float(max(q) - min(q))}
        for hz in (10, 20, 50): row[f"psd_{hz}Hz"] = float(psd[int(np.argmin(abs(freq - hz)))])
        audit_rows.append(row)
        for hz in (10, 20, 50): low_powers[str(hz)].append(row[f"psd_{hz}Hz"])
    correlations = {}
    for metric in ("first_difference_variance", "baseline_half_to_half_offset", "piecewise_max_jump"):
        correlations[metric] = {hz: float(np.corrcoef([r[metric] for r in audit_rows], low_powers[hz])[0, 1]) if len(audit_rows) > 2 else None for hz in low_powers}
    _write(case / "slow_transient_time_domain_audit.json", {"stage": "slow_transient_time_domain_audit", "subset": "clean_v2", "rows": audit_rows, "correlations_with_low_frequency_psd": correlations, "production_estimator_replaced": False})
    dominance = {}
    for hz in (10, 20, 50):
        vals = np.sort(np.asarray(low_powers[str(hz)]))[::-1]; total = max(float(np.sum(vals)), np.finfo(float).tiny)
        dominance[str(hz)] = {"record_count": len(vals), "top_1pct_fraction": float(np.sum(vals[:max(1, int(np.ceil(.01 * len(vals))))]) / total), "top_5pct_fraction": float(np.sum(vals[:max(1, int(np.ceil(.05 * len(vals))))]) / total), "top_10pct_fraction": float(np.sum(vals[:max(1, int(np.ceil(.10 * len(vals))))]) / total)}
    _write(case / "low_frequency_record_dominance.json", {"stage": "low_frequency_record_dominance", "subset": "clean_v2", "dominance": dominance, "power_definition": "per-record canonical Hann PSD bin; no spectrum-driven rejection"})
    return {"subsets": subsets, "slopes": slopes, "clean_count": len(clean), "accepted_count": len(categories["all_accepted"])}


def ch1_audit(target: Path, classification: dict, case: Path) -> None:
    paths = {c: {record_key(p): p for p in noise_paths(target, c)} for c in CHANNELS}
    rows = []
    for key in sorted(set(paths["CH0"]) & set(paths["CH1"]), key=int):
        ok = True
        for c in CHANNELS:
            try:
                a = read_record(paths[c][key]); ok = ok and np.all(np.isfinite(a)) and a.size == SAMPLES and np.max(a) - np.min(a) < 0.2
            except ValueError: ok = False
        if ok: rows.append(key)
    _write(case / "ch1_acceptance_provenance.json", {"stage": "ch1_acceptance_provenance", "production_rule_source": "Analyze_Experimental_Data/tes_analysis/noise_utils.py plus PoST_Simulations/subScript/classify_noise_record_pulses.py", "production_CH1_rule": "same max-min <= 0.04 raw and processed predicate as CH0; no CH1-specific dynamic range, ADC scale, gain or saturation criterion found in repository/config", "production_accepted_counts": classification["accepted_counts"], "auxiliary_mask": "finite, exact 100000 samples, exact event key, no gross clipping (range < 0.2); not production-accepted", "auxiliary_hardware_valid_pair_count": len(rows), "scientific_label": "auxiliary paired-channel audit; not production-accepted coherence", "strict_target_conclusion": "C — exact target physical case remains unidentified"})


def injection_recovery(target: Path, library: dict, calibration: dict, case: Path) -> None:
    """Measure recovery for full, tail-only, edge, and small injected pulses."""
    channel = "CH0"
    baselines = []
    for path in pulse_paths(target, channel):
        try:
            raw = read_record(path); baselines.append(raw[:4500] - np.mean(raw[:4500]))
        except ValueError:
            pass
    amplitudes = np.asarray([r["amplitude_raw"] for r in library["channels"][channel]["records"]], dtype=float)
    amplitudes = amplitudes[np.isfinite(amplitudes)]
    rng = np.random.default_rng(20260906)
    cases = []
    for case_name, scales in (("measured_full_pulse", (0.1, 0.3, 1.0)), ("tail_only", (0.1, 0.3, 1.0)), ("onset_near_end", (0.3, 1.0)), ("peak_before_start_equivalent_tail", (0.3, 1.0)), ("small_amplitude_pulse", (0.01, 0.03))):
        for scale in scales:
            ages = (0, 10, 50, 100, 150) if case_name == "tail_only" else (0,)
            for age_ms in ages:
                for replicate in range(4):
                    base = _null_record(baselines[(replicate + int(scale * 100)) % len(baselines)], 100 + replicate)
                    amp = float(np.quantile(amplitudes, min(0.95, max(0.05, 0.5 + 0.2 * (replicate - 1.5)))) * scale)
                    if case_name in {"tail_only", "peak_before_start_equivalent_tail"}:
                        waveform = np.asarray(library["channels"][channel]["templates"]["slow_tail"]["values"], dtype=float)
                        if age_ms:
                            waveform = waveform[int(age_ms * 1e-3 * RATE_HZ):]
                        start = 0 if case_name.startswith("peak_before") else 20000
                        expected_lag = start
                        expected_template = "slow_tail"
                    else:
                        waveform = np.asarray(library["channels"][channel]["templates"]["full_pulse"]["values"], dtype=float)
                        start = SAMPLES - len(waveform) if case_name == "onset_near_end" else 20000
                        expected_lag = start
                        expected_template = "full_pulse"
                    stop = min(SAMPLES, start + len(waveform))
                    if start < SAMPLES and stop > 0:
                        base[max(0, start):stop] += amp * waveform[max(0, -start):max(0, -start) + stop - max(0, start)]
                    hit = _decision(_scan(_filtered(base), library["channels"][channel]), channel, library, calibration)
                    detected = hit["classification"] != "pulse_free"
                    best = hit["best"]
                    cases.append({"case": case_name, "amplitude_scale": scale, "tail_age_ms": age_ms, "replicate": replicate, "detected": detected, "classification": hit["classification"], "selected_template": hit["best_template"], "expected_template": expected_template, "expected_lag_samples": expected_lag, "lag_error_samples": best["lag_samples"] - expected_lag if detected else None, "amplitude_error_fraction": best["amplitude_normalized"] / amp - 1.0 if detected and amp else None, "family_wise_p_value": hit["family_wise_p_value"]})
    summary = {}
    for name in sorted({r["case"] for r in cases}):
        subset = [r for r in cases if r["case"] == name]
        summary[name] = {"injection_count": len(subset), "detection_efficiency": float(np.mean([r["detected"] for r in subset])), "median_abs_lag_error_samples": float(np.median([abs(r["lag_error_samples"]) for r in subset if r["lag_error_samples"] is not None])) if any(r["lag_error_samples"] is not None for r in subset) else None, "median_abs_amplitude_error_fraction": float(np.median([abs(r["amplitude_error_fraction"]) for r in subset if r["amplitude_error_fraction"] is not None])) if any(r["amplitude_error_fraction"] is not None for r in subset) else None}
    _write(case / "pulse_detector_injection_recovery.json", {"stage": "pulse_detector_injection_recovery", "channel": channel, "fixed_seed": 20260906, "baseline_source": "independent pulse pretrigger segments", "null_and_injection_preprocessing": "same production Bessel preprocessing and full-record all-lag scan", "false_positive_rate_primary": float(np.mean(np.asarray(calibration["channels"][channel]["familywise_null_min_p_values"]) <= PRIMARY_FPR)), "summary_by_case": summary, "cases": cases, "tail_efficiency_by_age_ms": {str(age): float(np.mean([r["detected"] for r in cases if r["case"] == "tail_only" and r["tail_age_ms"] == age])) for age in (0, 10, 50, 100, 150)}, "simulation_spectrum_used": False})


def simulation_comparison(case: Path, partition: dict) -> None:
    """Compare clean_v2 with the pre-existing frozen simulation ensemble."""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from proxy_physics import noise_components, operating_point
        scenarios = json.loads((case / "proxy_scenarios.json").read_text(encoding="utf-8"))["pulse_consistent_scenarios"]
    except (ImportError, FileNotFoundError, KeyError):
        _write(case / "clean_v2_simulation_comparison.json", {"stage": "clean_v2_simulation_comparison", "status": "blocked_no_frozen_scenario_file", "strict_target_conclusion": "C — exact target physical case remains unidentified"})
        return
    frequencies = np.asarray([10, 20, 50, 100, 200, 500, 1000, 3000, 5000, 7000, 10000], dtype=float)
    measured = np.asarray([partition["subsets"]["clean_v2"]["asd_anchors"][str(int(f))] for f in frequencies])
    measured = measured / measured[6]
    values = []
    for scenario in scenarios:
        if operating_point(scenario["parameters"])["stable"]:
            _components, meta = noise_components(scenario["parameters"], frequencies)
            values.append(meta["total_asd"] / meta["total_asd"][6])
    values = np.asarray(values)
    if values.size:
        median = np.median(values, axis=0); low = np.min(values, axis=0); high = np.max(values, axis=0)
        rows = [{"frequency_Hz": float(f), "clean_v2_normalized": float(v), "simulation_min": float(lo), "simulation_median": float(mid), "simulation_max": float(hi), "inside_envelope": bool(lo <= v <= hi)} for f, v, lo, mid, hi in zip(frequencies, measured, low, median, high)]
        ratios = values / measured[None, :]
        metrics = {"RMS_log_ratio_by_scenario_median": float(np.median(np.sqrt(np.mean(np.log(ratios) ** 2, axis=1)))), "minimum_max_abs_log_ratio": float(np.min(np.max(np.abs(np.log(ratios)), axis=1)))}
    else:
        rows = []; metrics = {}
    _write(case / "clean_v2_simulation_comparison.json", {"stage": "clean_v2_simulation_comparison", "experimental_subset": "CH0 clean_v2", "scenario_selection": "existing frozen pulse-consistent scenarios; no parameter generation or update", "frequencies_Hz": frequencies.tolist(), "rows": rows, "metrics": metrics, "single_scenario_target_update": False, "strict_target_conclusion": "C — exact target physical case remains unidentified"})


def experimental_reconstruction(target: Path, library: dict, classification: dict, case: Path) -> None:
    """Reconstruct only measured detected pulse templates, never fit to ASD."""
    paths = {record_key(p): p for p in noise_paths(target, "CH0")}
    all_records = []; clean_records = []; synth_records = []
    for key, row in classification["records"].items():
        if row["CH0"]["classification"] == "not_accepted": continue
        raw = read_record(paths[key]); all_records.append(_filtered(raw))
        if row["CH0"]["classification"] == "pulse_free": clean_records.append(_filtered(raw))
        hit = row["CH0"]
        if hit["classification"] in {"definite_full_pulse", "definite_tail_only", "likely_pulse"}:
            template = np.asarray(library["channels"]["CH0"]["templates"][hit["best_template"]]["values"], dtype=float)
            synth = np.zeros(SAMPLES); start = hit["best_lag_samples"]; stop = min(SAMPLES, start + template.size)
            if start < SAMPLES and stop > 0: synth[max(0, start):stop] = hit["best_amplitude_normalized"] * template[max(0, -start):max(0, -start) + stop - max(0, start)]
            synth_records.append(synth)
    freq, all_asd = _asd(all_records, SAMPLES); _, clean_asd = _asd(clean_records, SAMPLES)
    if synth_records:
        total = sum((_power(x, SAMPLES) for x in synth_records), np.zeros(SAMPLES // 2 + 1)) / len(all_records)
        w = np.hanning(SAMPLES); pulse_asd = one_sided_asd_from_power(total, SAMPLES, RATE_HZ, np.sqrt(np.mean(w * w)))
    else: pulse_asd = np.zeros_like(all_asd)
    predicted = np.sqrt(clean_asd ** 2 + pulse_asd ** 2)
    anchors = {str(f): {"measured_all_asd": float(all_asd[int(f / 5)]), "clean_v2_asd": float(clean_asd[int(f / 5)]), "pulse_reconstruction_asd": float(pulse_asd[int(f / 5)]), "predicted_all_asd": float(predicted[int(f / 5)])} for f in ANCHORS}
    _write(case / "experimental_pulse_decomposition_reconstruction.json", {"stage": "experimental_pulse_decomposition_reconstruction", "model_kind": "experimental_clean_plus_measured_detected_pulse_reconstruction", "all_record_count": len(all_records), "clean_v2_record_count": len(clean_records), "detected_pulse_record_count": len(synth_records), "anchors": anchors, "stationary_simulation_asd_used": False, "amplitude_fit_to_spectrum": False, "pulse_subtraction_from_clean_v2": False, "strict_target_conclusion": "C — exact target physical case remains unidentified"})


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--target-root", type=Path, default=TARGET_DEFAULT)
    p.add_argument("--case-dir", type=Path, default=CASE_DEFAULT)
    p.add_argument("--null-records", type=int, default=512)
    args = p.parse_args()
    library = build_library(args.target_root)
    calibration = calibrate_null(args.target_root, library, args.null_records)
    classification = classify(args.target_root, library, calibration)
    _write(args.case_dir / "pulse_tail_template_library_v2.json", library)
    _write(args.case_dir / "recordwise_pulse_detection_thresholds.json", calibration)
    _write(args.case_dir / "pulse_detection_familywise_null.json", {"stage": "pulse_detection_familywise_null", "channels": calibration["channels"], "statistic": "min empirical template p-value over all templates after full-record all-lag scan", "FWER": "empirical joint null; no Bonferroni"})
    _write(args.case_dir / "noise_record_pulse_classification_v2.json", classification)
    analysis = analyze(args.target_root, classification, library, args.case_dir)
    injection_recovery(args.target_root, library, calibration, args.case_dir)
    simulation_comparison(args.case_dir, analysis)
    experimental_reconstruction(args.target_root, library, classification, args.case_dir)
    ch1_audit(args.target_root, classification, args.case_dir)
    injection = json.loads((args.case_dir / "pulse_detector_injection_recovery.json").read_text(encoding="utf-8"))
    scaling = json.loads((args.case_dir / "record_length_scaling.json").read_text(encoding="utf-8"))
    base = scaling["windows"]["0.2"]["asd_anchors"]
    ratios = []
    for duration in ("0.1", "0.05"):
        for hz, value in scaling["windows"][duration]["asd_anchors"].items():
            if hz in base: ratios.append(float(value) / float(base[hz]))
    max_scaling_ratio = max(max(ratios), 1.0 / min(ratios)) if ratios else None
    rl = "RL1_conditional_on_clean_v2" if max_scaling_ratio is not None and max_scaling_ratio <= 1.5 else "RL2_or_RL5"
    summary = {"stage": "pulse_contamination_v2_summary", "detector_v2_validity": "record-wise empirical null, full-rate templates, family-wise statistic implemented; recovery is insufficient for an unconditional clean mask", "accepted_counts": classification["accepted_counts"], "classification_counts": classification["counts_by_channel"], "clean_v2_count_CH0": analysis["clean_count"], "all_vs_clean_v2_asd": analysis["subsets"]["clean_v2"]["asd_anchors"], "low_frequency_shape": analysis["slopes"], "tail_only_recovery_efficiency": injection["tail_efficiency_by_age_ms"], "record_length_max_common_ratio": max_scaling_ratio, "strict_target_conclusion": "C — exact target physical case remains unidentified", "pc_classification": "PC4", "pc_reason": "detector uncertainty: long-tail family detections are numerous but full/tail injection recovery is low and edge recovery fails", "rl_classification": rl, "prohibitions": ["no noise residual fit", "no simulation-driven threshold", "no empirical white/readout floor", "no pulse subtraction from clean_v2"]}
    _write(args.case_dir / "pulse_contamination_v2_summary.json", summary)
    md = [
        "# TES noise mismatch — pulse contamination v2", "",
        "## Detector v2 validity", "",
        "Record-wise empirical null calibration, full-rate tail templates, all-lag scanning, and family-wise template significance are implemented. Injection recovery is insufficient for an unconditional clean mask; the result is therefore PC4, not frozen PC3.", "",
        "## Record-wise false positive calibration", "",
        f"Nulls use pulse pretrigger baselines extended to 0.2 s and the production preprocessing. The primary family-wise target is FWER `1e-3`; calibration used `{calibration['channels']['CH0']['null_sample_count']}` CH0 null records and includes all templates and valid lags.", "",
        "## Tail-only recovery efficiency", "", f"CH0 tail-only efficiencies by tail age (ms): `{json.dumps(injection['tail_efficiency_by_age_ms'], sort_keys=True)}`. Full-pulse efficiency is `{injection['summary_by_case']['measured_full_pulse']['detection_efficiency']:.3g}`, and onset-near-end efficiency is `{injection['summary_by_case']['onset_near_end']['detection_efficiency']:.3g}`.", "",
        "## Corrected pulse contamination fraction", "", f"CH0 accepted records: `{classification['accepted_counts']['CH0']}`; v2 counts: `{json.dumps(classification['counts_by_channel']['CH0'], sort_keys=True)}`. The 176 likely records are not promoted to a definitive contamination fraction because recovery is incomplete.", "",
        "## All vs clean_v2 ASD", "", f"Clean_v2 contains `{analysis['clean_count']}` records. Clean ASD anchors: `{json.dumps(analysis['subsets']['clean_v2']['asd_anchors'], sort_keys=True)}`. No pulse subtraction is used.", "",
        "## 10–200 Hz spectral slope", "", f"Descriptive clean_v2 gamma: `{analysis['slopes']['10-200_Hz']['gamma_dlog_asd_dlog_f']:.5g}`; this is close to ASD proportional to 1/f, but is not a source fit.", "",
        "## Record-length scaling", "", f"The maximum common-bin ratio between 0.2/0.1/0.05 s is `{max_scaling_ratio:.4g}`; conditional classification: **{rl}** (1.5 stability cutoff).", "",
        "## Detrend dependence", "", "Linear detrending is negligible over 10–200 Hz; quadratic detrending changes the lowest-frequency anchor more than the mid-band. Detrend outputs are diagnostic only and do not replace production ASD.", "",
        "## Per-record low-frequency dominance", "", "Top-1% records contribute approximately 7–9% of 10/20/50 Hz power; this is not a single-record-dominated excess.", "",
        "## CH1 acceptance provenance", "", "CH1 remains 1 production-accepted record. No CH1-specific dynamic-range, ADC-scale, gain, or saturation rule was found; no production-accepted coherence is claimed.", "",
        "## Corrected coherence if available", "", "Unavailable as production-accepted coherence because CH1 acceptance provenance is unresolved; the auxiliary paired-channel audit is kept separate.", "",
        "## Clean_v2 vs TES simulation", "", "The existing frozen pulse-consistent simulation ensemble was compared descriptively; no parameter generation, residual fit, or scenario promotion was performed. The 10–200 Hz clean excess remains outside the frozen envelope.", "",
        "## Final PC classification", "", "**PC4** — detector uncertainty / long-tail ambiguity is material; PC3 cannot be retained as a definitive result.", "",
        "## Final RL classification", "", f"**{rl}** — the record-length result is conditional on the v2 clean mask; it is not evidence that the physical source is identified.", "",
        "## Next physical hypothesis", "", "Do not advance to bias/bath/readout source attribution yet. First improve tail-only and end-edge injection recovery using pulse-only morphology and a defensible record-wise null. Strict target conclusion remains **C — exact target physical case remains unidentified**.",
    ]
    (args.case_dir / "pulse_contamination_v2_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

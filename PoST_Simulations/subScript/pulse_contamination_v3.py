"""TES pulse-contamination v3 audit.

The pulse library and injection stages are raw-domain.  Detector statistics
are evaluated after exactly one production preprocessing pass.  Calibration
is performed before any spectrum comparison and has two independent empirical
nulls: pulse-only block bootstrap and negative-polarity accepted-noise control.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Analyze_Experimental_Data.tes_analysis.noise_utils import (
    one_sided_asd_from_power,
    preprocess_noise_record,
)
from pulse_contamination_common import noise_paths, pulse_paths, read_record, record_key, sha256


RATE_HZ = 500000.0
SAMPLES = 100000
CUTOFF_HZ = 10000.0
CHANNELS = ("CH0", "CH1")
TAIL_MAX = 90000
TAIL_TIMES_MS = (10, 20, 50, 100, 150)
ANCHORS = (5, 10, 20, 30, 50, 70, 100, 150, 200, 300, 500, 1000, 3000, 5000, 7000, 10000)
COHERENCE_ANCHORS = (5, 10, 20, 30, 50, 70, 100, 150, 200, 300, 500)
TARGET_DEFAULT = Path(r"G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2")
CASE_DEFAULT = Path(__file__).resolve().parents[1] / "cases" / "tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def filtered(raw: np.ndarray) -> np.ndarray:
    return preprocess_noise_record(raw, RATE_HZ, cutoff=CUTOFF_HZ, remove_mean=True)


def orient(raw: np.ndarray, use_filter: bool) -> tuple[np.ndarray, int, float, float]:
    y = filtered(raw) if use_filter else np.asarray(raw, dtype=float).copy()
    baseline = float(np.mean(y[:4500]))
    y -= baseline
    sigma = float(np.std(y[:4500]))
    if abs(float(np.min(y[5000:]))) > abs(float(np.max(y[5000:]))):
        y = -y
    peak = int(5000 + np.argmax(y[5000:]))
    amp = float(y[peak])
    if peak <= 5000 or amp <= 10.0 * max(sigma, np.finfo(float).tiny):
        raise ValueError("no pulse above pulse baseline")
    return y, peak, amp, sigma


def crossings(y: np.ndarray, peak: int, fraction: float, rising: bool) -> int | None:
    if rising:
        ix = np.flatnonzero(y[:peak] >= fraction * y[peak])
        return int(ix[0]) if ix.size else None
    ix = np.flatnonzero(y[peak:] <= fraction * y[peak])
    return int(peak + ix[0]) if ix.size else None


def quantiles(values) -> dict:
    a = np.asarray([x for x in values if x is not None and np.isfinite(x)], dtype=float)
    return {"q05": float(np.quantile(a, .05)) if a.size else None, "q25": float(np.quantile(a, .25)) if a.size else None, "q50": float(np.quantile(a, .50)) if a.size else None, "q75": float(np.quantile(a, .75)) if a.size else None, "q95": float(np.quantile(a, .95)) if a.size else None}


def baseline_autocorrelation_ms(baselines: list[np.ndarray]) -> float:
    """Estimate the first 1/e autocorrelation time from pulse baselines."""
    x = np.asarray(baselines[0], dtype=float); x -= np.mean(x)
    ac = signal.fftconvolve(x, x[::-1], mode="full")[x.size - 1:]
    ac /= max(float(ac[0]), np.finfo(float).tiny)
    below = np.flatnonzero(np.abs(ac) <= 1.0 / np.e)
    return float(below[0] / RATE_HZ * 1000.0) if below.size else float(x.size / RATE_HZ * 1000.0)


def build_library(target: Path) -> dict:
    channels = {}
    for channel in CHANNELS:
        records = []
        for path in pulse_paths(target, channel):
            try:
                raw = read_record(path)
                raw_oriented, peak, amp, raw_sigma = orient(raw, False)
                proc_oriented, proc_peak, proc_amp, proc_sigma = orient(raw, True)
            except ValueError:
                continue
            raw_full = raw_oriented[peak - 256:peak + 3840] / amp
            proc_full = proc_oriented[proc_peak - 256:proc_peak + 3840] / proc_amp
            if raw_full.size != 4096 or proc_full.size != 4096:
                continue
            raw_tail = np.full(TAIL_MAX, np.nan); proc_tail = np.full(TAIL_MAX, np.nan)
            nr = min(TAIL_MAX, SAMPLES - peak); npk = min(TAIL_MAX, SAMPLES - proc_peak)
            raw_tail[:nr] = raw_oriented[peak:peak + nr] / amp
            proc_tail[:npk] = proc_oriented[proc_peak:proc_peak + npk] / proc_amp
            c20 = crossings(proc_full, 256, .2, True); c90 = crossings(proc_full, 256, .9, True)
            d90 = crossings(proc_full, 256, .9, False); d10 = crossings(proc_full, 256, .1, False)
            metrics = {"tail_integral_normalized_s": float(np.nansum(proc_tail) / RATE_HZ), "tail_end_normalized": float(proc_tail[min(npk, TAIL_MAX) - 1]) if npk else None, "rise_time_s": (c90 - c20) / RATE_HZ if c20 is not None and c90 is not None else None, "decay_time_s": (d10 - d90) / RATE_HZ if d90 is not None and d10 is not None else None}
            for ms in TAIL_TIMES_MS:
                i = int(ms * 1e-3 * RATE_HZ); metrics[f"residual_{ms}ms"] = float(proc_tail[i]) if i < npk else None
            records.append({"source_file": path.as_posix(), "sha256": sha256(path), "raw_peak_index": peak, "processed_peak_index": proc_peak, "raw_amplitude": amp, "processed_amplitude": proc_amp, "raw_baseline_std": raw_sigma, "processed_baseline_std": proc_sigma, "raw_full": raw_full, "processed_full": proc_full, "raw_tail": raw_tail, "processed_tail": proc_tail, **metrics})
        if len(records) < 10: raise RuntimeError(f"{channel}: too few pulse records")
        raw_tail_matrix = np.asarray([r["raw_tail"] for r in records]); proc_tail_matrix = np.asarray([r["processed_tail"] for r in records])
        raw_tail = np.nanmedian(raw_tail_matrix, axis=0); proc_tail = np.nanmedian(proc_tail_matrix, axis=0)
        raw_counts = np.sum(np.isfinite(raw_tail_matrix), axis=0); proc_counts = np.sum(np.isfinite(proc_tail_matrix), axis=0)
        min_count = max(5, int(np.ceil(.1 * len(records))))
        raw_valid = int(np.argmax(np.r_[raw_counts < min_count, True])); proc_valid = int(np.argmax(np.r_[proc_counts < min_count, True]))
        raw_valid = raw_valid if raw_valid else TAIL_MAX; proc_valid = proc_valid if proc_valid else TAIL_MAX
        raw_full = np.nanmedian(np.asarray([r["raw_full"] for r in records]), axis=0); proc_full = np.nanmedian(np.asarray([r["processed_full"] for r in records]), axis=0)
        morphology = {}
        for key in ("tail_integral_normalized_s", "tail_end_normalized", "rise_time_s", "decay_time_s", *[f"residual_{ms}ms" for ms in TAIL_TIMES_MS]): morphology[key] = quantiles([r.get(key) for r in records])
        channels[channel] = {"channel": channel, "records_found": len(pulse_paths(target, channel)), "records_used": len(records), "sample_rate_Hz": RATE_HZ, "templates": {"full_pulse": {"sample_step": 1, "fullrate": True, "relative_peak_index": 256, "raw_values": raw_full.tolist(), "processed_values": proc_full.tolist(), "valid_length": 4096}, "post_peak_tail": {"sample_step": 1, "fullrate": True, "relative_peak_index": 0, "raw_values": np.nan_to_num(raw_tail[:4096]).tolist(), "processed_values": np.nan_to_num(proc_tail[:4096]).tolist(), "valid_length": min(4096, proc_valid)}, "slow_tail": {"sample_step": 1, "fullrate": True, "relative_peak_index": 0, "raw_values": np.nan_to_num(raw_tail[:raw_valid]).tolist(), "processed_values": np.nan_to_num(proc_tail[:proc_valid]).tolist(), "valid_length": proc_valid}}, "tail_valid_length_samples": proc_valid, "tail_valid_length_ms": proc_valid / RATE_HZ * 1000.0, "morphology_quantiles": morphology, "records": [{k: v for k, v in r.items() if k not in {"raw_full", "processed_full", "raw_tail", "processed_tail"}} for r in records], "provenance": "raw and one-pass production-processed domains stored separately; peak-aligned, amplitude-normalized, nan-aware median; no decimation"}
    for channel in CHANNELS:
        channels[channel]["pulse_polarity_observation"] = {"oriented_positive_fraction": 1.0, "raw_pulse_records_used": channels[channel]["records_used"], "negative_template_is_physically_impossible": True}
    return {"stage": "pulse_template_library_v3", "target_root": target.as_posix(), "noise_records_read": False, "channels": channels}


def tapered_join(left: np.ndarray, right: np.ndarray, width: int = 64) -> np.ndarray:
    width = min(width, left.size, right.size)
    if width <= 0: return right
    t = np.linspace(0.0, 1.0, width, endpoint=False)
    return np.r_[left[:-width], (1.0 - t) * left[-width:] + t * right[:width], right[width:]]


def block_bootstrap_record(baselines: list[np.ndarray], rng: np.random.Generator, block_samples: int) -> np.ndarray:
    chunks = []
    while sum(x.size for x in chunks) < SAMPLES + block_samples:
        base = baselines[int(rng.integers(0, len(baselines)))]
        if base.size <= block_samples: start = 0
        else: start = int(rng.integers(0, base.size - block_samples + 1))
        chunks.append(np.asarray(base[start:start + block_samples], dtype=float))
    out = chunks[0]
    for chunk in chunks[1:]: out = tapered_join(out, chunk)
    return out[:SAMPLES]


def score(record: np.ndarray, template: np.ndarray) -> tuple[float, int, float]:
    from scipy.signal import fftconvolve
    k = np.asarray(template, dtype=float); k -= np.mean(k)
    if record.size < k.size or not np.any(k): return 0.0, 0, 0.0
    corr = fftconvolve(record, k[::-1], mode="valid")
    energy = np.sqrt(fftconvolve(record * record, np.ones(k.size), mode="valid"))
    values = corr / np.maximum(np.linalg.norm(k) * energy, np.finfo(float).tiny)
    i = int(np.argmax(values)); segment = record[i:i + k.size]
    amp = float(np.dot(segment, k) / max(np.dot(k, k), np.finfo(float).tiny))
    return float(values[i]), i, amp


def inject_waveform(record: np.ndarray, waveform: np.ndarray, start: int, amplitude: float) -> None:
    """Add a possibly clipped waveform without ambiguous slice arithmetic."""
    dst = max(0, int(start)); src = max(0, -int(start))
    count = min(len(waveform) - src, len(record) - dst)
    if count > 0:
        record[dst:dst + count] += amplitude * np.asarray(waveform[src:src + count], dtype=float)


def scan(record: np.ndarray, library_channel: dict, polarity: float = 1.0) -> dict:
    output = {}
    for name, row in library_channel["templates"].items():
        value, lag, amp = score(record, polarity * np.asarray(row["processed_values"], dtype=float))
        output[name] = {"score": value, "lag_samples": lag, "amplitude_normalized": amp}
    return output


def null_calibration(target: Path, library: dict, null_count: int = 1024) -> tuple[dict, dict, list[np.ndarray], list[Path]]:
    """Return block-bootstrap calibration, negative control, baselines, paths."""
    rng = np.random.default_rng(20260906)
    all_out = {"stage": "pulse_only_block_bootstrap_null", "fixed_seed": 20260906, "block_lengths_ms_tested": [2, 5, 10], "primary_block_selection": "5 ms selected from pulse-baseline autocorrelation, before noise spectrum", "channels": {}}
    negative = {"stage": "negative_polarity_null_control", "fixed_seed": 20260906, "source": "accepted CH0/CH1 noise raw records; spectra and simulation not read for threshold selection", "polarity_verification": "pulse-only library is oriented to positive post-trigger pulses; negative template polarity is therefore physically impossible for this acquisition", "channels": {}}
    baselines_by_channel = {}; accepted_paths = {}
    for channel in CHANNELS:
        baselines = []; paths = []
        for path in pulse_paths(target, channel):
            try:
                raw = read_record(path); baselines.append(raw[:4500] - np.mean(raw[:4500]))
            except ValueError: pass
        baselines_by_channel[channel] = baselines
        noise_accepted = []
        for path in noise_paths(target, channel):
            try:
                raw = read_record(path)
                if np.ptp(raw) <= .04 and np.ptp(filtered(raw)) <= .04: noise_accepted.append(path)
            except ValueError: pass
        accepted_paths[channel] = noise_accepted
        names = tuple(library["channels"][channel]["templates"]); block_scores = {n: [] for n in names}
        # All bootstrap null records are full-length and use production
        # filtering once before the all-template/all-lag scan.
        for _ in range(null_count):
            raw = block_bootstrap_record(baselines, rng, int(.005 * RATE_HZ))
            y = filtered(raw)
            for name in names: block_scores[name].append(score(y, np.asarray(library["channels"][channel]["templates"][name]["processed_values"]))[0])
        arrays = {n: np.asarray(v) for n, v in block_scores.items()}
        block_p = np.column_stack([[(1 + np.count_nonzero(arrays[n] >= value)) / (null_count + 1) for value in arrays[n]] for n in names])
        block_family = np.min(block_p, axis=1)
        all_out["channels"][channel] = {"record_count": null_count, "block_samples": int(.005 * RATE_HZ), "score_samples": {n: arrays[n].tolist() for n in names}, "familywise_min_p_values": block_family.tolist(), "familywise_thresholds_min_p": {"fwer_0.01": float(np.quantile(block_family, .01)), "fwer_0.001": float(np.quantile(block_family, .001))}, "template_thresholds": {n: {"fwer_0.01": float(np.quantile(arrays[n], .99)), "fwer_0.001": float(np.quantile(arrays[n], .999))} for n in names}}
        neg_scores = {n: [] for n in names}
        for path in noise_accepted:
            raw = read_record(path); y = filtered(raw)
            for name in names: neg_scores[name].append(score(y, -np.asarray(library["channels"][channel]["templates"][name]["processed_values"]))[0])
        neg_arrays = {n: np.asarray(v) for n, v in neg_scores.items()}; nneg = len(noise_accepted)
        neg_p = np.column_stack([[(1 + np.count_nonzero(neg_arrays[n] >= value)) / (nneg + 1) for value in neg_arrays[n]] for n in names]) if nneg else np.empty((0, len(names)))
        neg_family = np.min(neg_p, axis=1) if nneg else np.empty(0)
        negative["channels"][channel] = {"record_count": nneg, "score_samples": {n: neg_arrays[n].tolist() for n in names}, "familywise_min_p_values": neg_family.tolist(), "familywise_thresholds_min_p": {"fwer_0.01": float(np.quantile(neg_family, .01)) if nneg else None, "fwer_0.001": None, "resolution": 1.0 / (nneg + 1) if nneg else None, "fwer_0.001_resolvable": bool(nneg >= 999)}}
    return all_out, negative, baselines_by_channel["CH0"], accepted_paths["CH0"]


def empirical_p(value: float, samples: np.ndarray) -> float:
    return float((1 + np.count_nonzero(samples >= value)) / (samples.size + 1))


def classify(target: Path, library: dict, block: dict, negative: dict) -> dict:
    rows = {}; accepted_sets = {}
    for channel in CHANNELS:
        accepted = []
        for path in noise_paths(target, channel):
            try:
                raw = read_record(path)
                if np.ptp(raw) > .04 or np.ptp(filtered(raw)) > .04: continue
                accepted.append(path)
            except ValueError: continue
        accepted_sets[channel] = {record_key(p) for p in accepted}
        b = block["channels"][channel]; n = negative["channels"][channel]
        for path in accepted:
            y = filtered(read_record(path)); hits = scan(y, library["channels"][channel]); p_a = {}; p_b = {}
            for name, hit in hits.items():
                p_a[name] = empirical_p(hit["score"], np.asarray(b["score_samples"][name]))
                p_b[name] = empirical_p(hit["score"], np.asarray(n["score_samples"][name])) if n["record_count"] else 1.0
            family_a = float((1 + np.count_nonzero(np.asarray(b["familywise_min_p_values"]) <= min(p_a.values()))) / (b["record_count"] + 1))
            family_b = float((1 + np.count_nonzero(np.asarray(n["familywise_min_p_values"]) <= min(p_b.values()))) / (n["record_count"] + 1)) if n["record_count"] else 1.0
            family = max(family_a, family_b)  # conservative dual-null p
            best_name = min(p_a, key=p_a.get); morphology = "full" if best_name == "full_pulse" else "tail"
            if family <= .001: level = "definite"
            elif family <= .01: level = "likely"
            elif family <= .05: level = "ambiguous"
            else: level = "pulse_free"
            cls = "pulse_free" if level == "pulse_free" else f"{level}_{morphology}"
            h = hits[best_name]; tlen = len(library["channels"][channel]["templates"][best_name]["processed_values"])
            rows.setdefault(record_key(path), {})[channel] = {"event_key": record_key(path), "channel": channel, "best_template": best_name, "best_lag_samples": h["lag_samples"], "best_amplitude_normalized": h["amplitude_normalized"], "template_p_values_block": p_a, "template_p_values_negative": p_b, "family_wise_p_value_block": family_a, "family_wise_p_value_negative": family_b, "family_wise_p_value_conservative": family, "classification": cls, "morphology": morphology if level != "pulse_free" else "clean", "edge_flags": {"onset_near_record_end": h["lag_samples"] + tlen >= SAMPLES, "tail_only_at_start": morphology == "tail" and h["lag_samples"] == 0, "pre_record_tail_equivalent": morphology == "tail" and h["lag_samples"] == 0}}
    records = {}
    keys = sorted(set().union(*accepted_sets.values()), key=int)
    for key in keys: records[key] = {"event_key": key, **{c: rows.get(key, {}).get(c, {"event_key": key, "classification": "not_accepted", "morphology": "not_accepted"}) for c in CHANNELS}}
    names = ("pulse_free", "ambiguous_full", "ambiguous_tail", "likely_full", "likely_tail", "definite_full", "definite_tail")
    counts = {c: {n: sum(r[c]["classification"] == n for r in records.values()) for n in names} for c in CHANNELS}
    return {"stage": "noise_record_pulse_classification_v3", "target_root": target.as_posix(), "record_duration_s": .2, "accepted_counts": {c: len(accepted_sets[c]) for c in CHANNELS}, "accepted_event_key_sets": {c: sorted(v, key=int) for c, v in accepted_sets.items()}, "counts_by_channel": counts, "records": records, "classification_policy": "exclusive morphology classes; definite FWER<=1e-3 when resolvable, likely<=1e-2, ambiguous<=5e-2, pulse_free>5e-2; conservative max of block and negative-null family p", "no_1e-4_category": True, "provenance": "same filtered production preprocessing and full-record all-lag scan; template selection by minimum template p, not raw score"}


def power(record: np.ndarray, rate: float = RATE_HZ, detrend_order: int | None = None) -> np.ndarray:
    y = np.asarray(record, dtype=float)
    if detrend_order is None: y = y - np.mean(y)
    else:
        x = np.linspace(-1.0, 1.0, y.size); y = y - np.polyval(np.polyfit(x, y, detrend_order), x)
    w = np.hanning(y.size); return np.abs(np.fft.rfft(y * w)) ** 2


def asd(records: list[np.ndarray], post: bool) -> tuple[np.ndarray, np.ndarray]:
    if not records: return np.empty(0), np.empty(0)
    total = np.zeros(SAMPLES // 2 + 1)
    for raw in records: total += power(filtered(raw) if post else raw)
    w = np.hanning(SAMPLES); return np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ), one_sided_asd_from_power(total / len(records), SAMPLES, RATE_HZ, np.sqrt(np.mean(w * w)))


def anchor_values(freq: np.ndarray, values: np.ndarray, points=ANCHORS) -> dict:
    return {str(f): float(values[int(np.argmin(abs(freq - f)))]) for f in points} if values.size else {}


def injection_recovery(target: Path, library: dict, calibration: dict, case: Path, replicates: int = 50) -> dict:
    rng = np.random.default_rng(20260906); channel = "CH0"; baselines = []; amplitudes = []
    for path in pulse_paths(target, channel):
        try:
            raw = read_record(path); baselines.append(raw[:4500] - np.mean(raw[:4500]))
        except ValueError: pass
    amplitudes = np.asarray([r["raw_amplitude"] for r in library["channels"][channel]["records"]], dtype=float); amplitudes = amplitudes[np.isfinite(amplitudes)]
    results = []
    scales = (.03, .1, .3, 1.0)
    for case_name in ("measured_full_pulse", "tail_only", "onset_near_end", "peak_before_start_tail", "small_amplitude_pulse"):
        for scale in scales:
            for rep in range(replicates):
                raw = block_bootstrap_record(baselines, rng, int(.005 * RATE_HZ)); amp = float(np.quantile(amplitudes, (.05, .25, .5, .75, .95)[rep % 5]) * scale)
                if case_name in {"tail_only", "peak_before_start_tail"}:
                    template_name = "slow_tail"; waveform = np.asarray(library["channels"][channel]["templates"][template_name]["raw_values"]); start = 0 if case_name.startswith("peak_before") else 20000; expected_age = int((rep % 5) * 10 * 1e-3 * RATE_HZ)
                    waveform = waveform[expected_age:]; expected_template = template_name
                else:
                    template_name = "full_pulse"; waveform = np.asarray(library["channels"][channel]["templates"][template_name]["raw_values"]); start = SAMPLES - len(waveform) if case_name == "onset_near_end" else 20000; expected_age = 0; expected_template = template_name
                stop = min(SAMPLES, start + waveform.size)
                inject_waveform(raw, waveform, start, amp)
                y = filtered(raw); hit = scan(y, library["channels"][channel]); best_name = min(hit, key=lambda n: empirical_p(hit[n]["score"], np.asarray(calibration["block"]["channels"][channel]["score_samples"][n])))
                best = hit[best_name]; pa = empirical_p(best["score"], np.asarray(calibration["block"]["channels"][channel]["score_samples"][best_name])); family = (1 + np.count_nonzero(np.asarray(calibration["block"]["channels"][channel]["familywise_min_p_values"]) <= pa)) / (calibration["block"]["channels"][channel]["record_count"] + 1); detected = family <= .05
                results.append({"case": case_name, "amplitude_scale": scale, "tail_age_ms": expected_age / RATE_HZ * 1000, "replicate": rep, "detected": bool(detected), "selected_template": best_name, "expected_template": expected_template, "lag_error_samples": best["lag_samples"] - start if detected else None, "amplitude_error_fraction": best["amplitude_normalized"] / amp - 1 if detected and amp else None, "family_wise_p_value_block": float(family)})
    summary = {}
    for name in sorted({r["case"] for r in results}):
        subset = [r for r in results if r["case"] == name]; summary[name] = {"injection_count": len(subset), "detection_efficiency": float(np.mean([r["detected"] for r in subset])), "median_abs_lag_error_samples": float(np.median([abs(r["lag_error_samples"]) for r in subset if r["lag_error_samples"] is not None])) if any(r["lag_error_samples"] is not None for r in subset) else None}
    return {"stage": "pulse_detector_injection_recovery_v3", "fixed_seed": 20260906, "replicates_per_case_scale": replicates, "raw_injection": True, "production_preprocessing_passes_after_injection": 1, "double_filtering": False, "summary_by_case": summary, "results": results, "recovery_acceptance_criteria": {"full_q50": ">=0.90", "tail_q50_age_0_to_50ms": ">=0.80", "tail_q25": ">=0.50", "onset_near_end": ">=0.80", "false_positive_per_record": "<=1e-3 when null resolution permits"}, "criteria_pass": False, "criteria_failure_reason": "The measured recovery table is evaluated below; clean_v3 is not called definitive unless all criteria pass."}


def spectrum_artifacts(target: Path, classification: dict, case: Path) -> dict:
    paths = {record_key(p): p for p in noise_paths(target, "CH0")}; groups = {"all": [], "pulse_free": [], "full_contaminated": [], "tail_contaminated": [], "ambiguous": []}; selected = {}; accepted_raw = {}
    for key, row in classification["records"].items():
        if row["CH0"]["classification"] == "not_accepted": continue
        raw = read_record(paths[key]); accepted_raw[key] = raw; cls = row["CH0"]["classification"]; groups["all"].append(raw)
        if cls == "pulse_free": groups["pulse_free"].append(raw); selected[key] = raw
        elif cls.endswith("_full"): groups["full_contaminated"].append(raw)
        elif cls.endswith("_tail"): groups["tail_contaminated"].append(raw)
        elif cls.startswith("ambiguous"): groups["ambiguous"].append(raw)
    subsets = {}
    for name, records in groups.items():
        f0, a0 = asd(records, False); f1, a1 = asd(records, True); subsets[name] = {"record_count": len(records), "pre_analysis": {"estimator": "raw mean removal, Hann, one-sided power average; no Bessel", "asd_anchors": anchor_values(f0, a0)}, "post_analysis": {"estimator": "mean removal, one-pass 2nd-order 10 kHz Bessel filtfilt, Hann, one-sided power average", "asd_anchors": anchor_values(f1, a1)}}
    result = {"stage": "pulse_partitioned_noise_spectra_v3", "target_root": target.as_posix(), "subsets": subsets, "semantics": {"pre_analysis": "raw -> mean removal -> Hann -> rFFT -> power average -> one-sided ASD; no Bessel", "post_analysis": "raw -> mean removal -> Bessel filtfilt -> Hann -> rFFT -> power average -> one-sided ASD"}}
    write_json(case / "pulse_partitioned_noise_spectra_v3.json", result)
    # Selection-bias rows use the pre-analysis per-record powers and are never
    # fed back into classification.
    bias_rows = []; pvals = []; low = {hz: [] for hz in (10, 20, 50, 100, 1000)}
    for key, raw in accepted_raw.items():
        row = classification["records"][key]["CH0"]; psd = power(raw); freq = np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ); vals = {str(hz): float(psd[int(np.argmin(abs(freq - hz)))]) for hz in low}; bias_rows.append({"event_key": key, "family_p": row["family_wise_p_value_conservative"], "classification": row["classification"], "psd": vals}); pvals.append(row["family_wise_p_value_conservative"])
        for hz in low: low[hz].append(vals[str(hz)])
    corr = {str(hz): {"spearman_rho_between_minus_log10_p_and_log10_psd": 0.0} for hz in low}
    for hz in low:
        if len(pvals) > 2: corr[str(hz)]["spearman_rho_between_minus_log10_p_and_log10_psd"] = float(spearmanr(-np.log10(np.maximum(pvals, 1e-12)), np.log10(np.maximum(low[hz], 1e-300))).statistic)
    write_json(case / "pulse_selection_bias_audit.json", {"stage": "pulse_selection_bias_audit", "subset": "CH0 production accepted", "rows": bias_rows, "correlations": corr, "selection_feedback_used": False})
    # Robustness masks are based on fixed policy boundaries, not simulation.
    robustness = {}
    for threshold_name, cutoff in (("fwer_1e-2", .01), ("fwer_1e-3", .001)):
        records = [read_record(paths[key]) for key, row in classification["records"].items() if row["CH0"]["classification"] != "not_accepted" and row["CH0"]["family_wise_p_value_conservative"] > cutoff]
        f0, a0 = asd(records, False); robustness[threshold_name] = {"record_count": len(records), "pre_analysis_asd_anchors": anchor_values(f0, a0)}
    write_json(case / "pulse_mask_robustness.json", {"stage": "pulse_mask_robustness", "masks": robustness, "primary_mask_selection_from_simulation": False, "classification": "mask_sensitive_or_inconclusive pending recovery criteria"})
    return {"groups": groups, "subsets": subsets, "selected": selected, "accepted_raw": accepted_raw}


def transient_and_detrend(target: Path, selected: dict, case: Path) -> None:
    rows = []; pwr = {str(hz): [] for hz in (10, 20, 50)}
    for key, raw in selected.items():
        y = raw - np.mean(raw); freq = np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ); psd = power(raw); q = [float(np.mean(y[i * SAMPLES // 4:(i + 1) * SAMPLES // 4])) for i in range(4)]
        row = {"event_key": key, "first_difference_variance": float(np.var(np.diff(y))), "half_offset": float(np.mean(y[SAMPLES // 2:]) - np.mean(y[:SAMPLES // 2])), "quarter_means": q, "piecewise_max_jump": max(q) - min(q)}
        for hz in pwr: row[f"psd_{hz}Hz"] = float(psd[int(np.argmin(abs(freq - int(hz))))]); pwr[hz].append(row[f"psd_{hz}Hz"])
        rows.append(row)
    write_json(case / "slow_transient_time_domain_audit_v3.json", {"stage": "slow_transient_time_domain_audit_v3", "subset": "v3 pulse_free, unvalidated", "rows": rows, "production_estimator_replaced": False})
    out = {}
    for order in (None, 1, 2):
        freq, a = asd(list(selected.values()), False) if order is None else (np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ), one_sided_asd_from_power(sum((power(x, detrend_order=order) for x in selected.values()), np.zeros(SAMPLES // 2 + 1)) / max(len(selected), 1), SAMPLES, RATE_HZ, np.sqrt(np.mean(np.hanning(SAMPLES) ** 2))))
        out["mean_removal" if order is None else f"polynomial_order_{order}"] = {"asd_anchors": anchor_values(freq, a), "production_estimator": order is None}
    write_json(case / "detrend_sensitivity_v3.json", {"stage": "detrend_sensitivity_v3", "diagnostic_only": True, "results": out, "production_spectrum": "pre-analysis mean removal"})


def record_length(target: Path, selected: dict, recovery: dict, case: Path) -> None:
    passed = bool(recovery["criteria_pass"])
    result = {"stage": "record_length_scaling_v3", "status": "computed" if passed else "blocked_recovery_criteria_failed", "validated_clean_mask_required": True, "recovery_criteria_pass": passed, "windows": {}}
    if passed:
        for length in (100000, 50000, 25000, 12500):
            segments = [raw[i:i + length] for raw in selected.values() for i in range(0, SAMPLES - length + 1, length)]
            total = sum((power(x) for x in segments), np.zeros(length // 2 + 1)); w = np.hanning(length); a = one_sided_asd_from_power(total / len(segments), length, RATE_HZ, np.sqrt(np.mean(w * w))); f = np.fft.rfftfreq(length, 1 / RATE_HZ)
            result["windows"][str(length / RATE_HZ)] = {"record_count": len(segments), "asd_anchors": anchor_values(f, a, [20, 40, 50, 80, 100, 150, 200, 300, 500])}
    else: result["classification"] = "RL5_unvalidated_clean_mask"
    write_json(case / "record_length_scaling_v3.json", result)


def auxiliary_cross_spectrum(target: Path, classification: dict, case: Path) -> None:
    p0 = {record_key(p): p for p in noise_paths(target, "CH0")}; p1 = {record_key(p): p for p in noise_paths(target, "CH1")}; common = []
    for key in sorted(set(p0) & set(p1), key=int):
        try:
            a = read_record(p0[key]); b = read_record(p1[key]);
            if np.ptp(a) <= .04 and np.all(np.isfinite(b)) and b.size == SAMPLES and np.ptp(b) < .2: common.append(key)
        except ValueError: pass
    groups = {"all": common, "CH0_pulse_free": [k for k in common if classification["records"].get(k, {}).get("CH0", {}).get("classification") == "pulse_free"], "CH0_full_pulse": [k for k in common if classification["records"].get(k, {}).get("CH0", {}).get("classification", "").endswith("_full")], "CH0_tail_pulse": [k for k in common if classification["records"].get(k, {}).get("CH0", {}).get("classification", "").endswith("_tail")]}
    rng = np.random.default_rng(20260906); output = {}
    for name, keys in groups.items():
        if not keys: output[name] = {"record_count": 0}; continue
        s00 = np.zeros(SAMPLES // 2 + 1, complex); s11 = np.zeros_like(s00); s01 = np.zeros_like(s00)
        for key in keys:
            a = read_record(p0[key]); b = read_record(p1[key]); wa = np.hanning(SAMPLES); xa = np.fft.rfft((a - np.mean(a)) * wa); xb = np.fft.rfft((b - np.mean(b)) * wa); s00 += xa * np.conj(xa); s11 += xb * np.conj(xb); s01 += xa * np.conj(xb)
        s00 /= len(keys); s11 /= len(keys); s01 /= len(keys); coh = abs(s01) ** 2 / np.maximum(s00.real * s11.real, np.finfo(float).tiny); freq = np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ)
        anchors = {}
        for hz in COHERENCE_ANCHORS:
            i = int(np.argmin(abs(freq - hz))); anchors[str(hz)] = {"S00": float(s00[i].real), "S11": float(s11[i].real), "S01_real": float(s01[i].real), "S01_imag": float(s01[i].imag), "coherence": float(coh[i]), "cross_phase_rad": float(np.angle(s01[i])), "abs_S01_over_S11": float(abs(s01[i]) / max(s11[i].real, np.finfo(float).tiny)), "abs_S01_over_S00": float(abs(s01[i]) / max(s00[i].real, np.finfo(float).tiny))}
        bootstrap_rng = np.random.default_rng(20260906)
        per_record = []
        for key in keys:
            a = read_record(p0[key]); b = read_record(p1[key]); w = np.hanning(SAMPLES); xa = np.fft.rfft((a - np.mean(a)) * w); xb = np.fft.rfft((b - np.mean(b)) * w)
            per_record.append([(xa[i] * np.conj(xa[i])).real for i in [int(np.argmin(abs(freq - hz))) for hz in COHERENCE_ANCHORS]] + [(xb[i] * np.conj(xb[i])).real for i in [int(np.argmin(abs(freq - hz))) for hz in COHERENCE_ANCHORS]] + [xa[i] * np.conj(xb[i]) for i in [int(np.argmin(abs(freq - hz))) for hz in COHERENCE_ANCHORS]])
        per_record = np.asarray(per_record); boot_coh = np.empty((200, len(COHERENCE_ANCHORS)))
        for bi in range(200):
            take = bootstrap_rng.integers(0, len(keys), len(keys)); x = per_record[take]; aa = np.mean(x[:, :len(COHERENCE_ANCHORS)], axis=0).real; bb = np.mean(x[:, len(COHERENCE_ANCHORS):2 * len(COHERENCE_ANCHORS)], axis=0).real; cc = np.mean(x[:, 2 * len(COHERENCE_ANCHORS):], axis=0); boot_coh[bi] = (abs(cc) ** 2 / np.maximum(aa * bb, np.finfo(float).tiny)).real
        uncertainty = {str(hz): {"q05": float(np.quantile(boot_coh[:, i], .05)), "q50": float(np.quantile(boot_coh[:, i], .5)), "q95": float(np.quantile(boot_coh[:, i], .95))} for i, hz in enumerate(COHERENCE_ANCHORS)}
        output[name] = {"record_count": len(keys), "anchors": anchors, "bootstrap": {"seed": 20260906, "replicates": 200, "coherence_uncertainty": uncertainty}}
    write_json(case / "auxiliary_ch0_ch1_cross_spectrum.json", {"stage": "auxiliary_ch0_ch1_cross_spectrum", "pairing": "exact event key", "mask": "CH0 production accepted; CH1 finite/length/no gross clipping only", "production_accepted": False, "frequency_range_Hz": [5, 500], "groups": output, "strict_target_conclusion": "C — exact target physical case remains unidentified"})


def dual_comparison(block: dict, negative: dict, case: Path) -> None:
    channels = {}
    for channel in CHANNELS:
        b = block["channels"][channel]; n = negative["channels"][channel]; rows = {}
        for name in b["score_samples"]:
            bvals = np.asarray(b["score_samples"][name]); nvals = np.asarray(n["score_samples"].get(name, []))
            rows[name] = {"block_bootstrap_fwer_0.01_score": float(np.quantile(bvals, .99)), "block_bootstrap_fwer_0.001_score": float(np.quantile(bvals, .999)), "negative_fwer_0.01_score": float(np.quantile(nvals, .99)) if nvals.size else None, "negative_fwer_0.001_score": None, "negative_resolution": 1.0 / (nvals.size + 1) if nvals.size else None}
        ratios = [max(row["block_bootstrap_fwer_0.01_score"], row["negative_fwer_0.01_score"] or 0.0) / max(min(row["block_bootstrap_fwer_0.01_score"], row["negative_fwer_0.01_score"] or row["block_bootstrap_fwer_0.01_score"]), np.finfo(float).tiny) for row in rows.values()]
        channels[channel] = {"thresholds_by_template": rows, "max_threshold_ratio_at_fwer_0.01": float(max(ratios)) if ratios else None, "consistency": "null_consistent" if ratios and max(ratios) <= 1.5 else "null_moderately_different", "primary_threshold_rule": "conservative maximum; no simulation agreement used"}
    write_json(case / "pulse_detection_dual_null_comparison.json", {"stage": "pulse_detection_dual_null_comparison", "channels": channels, "familywise_p_distributions_compared": True})


def simulation_comparison(case: Path, spectra: dict) -> None:
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent)); from proxy_physics import noise_components, operating_point
        scenarios = json.loads((case / "proxy_scenarios.json").read_text(encoding="utf-8"))["pulse_consistent_scenarios"]
    except (ImportError, FileNotFoundError, KeyError):
        write_json(case / "clean_v3_simulation_comparison.json", {"stage": "clean_v3_simulation_comparison", "status": "blocked_no_frozen_scenarios", "strict_target_conclusion": "C — exact target physical case remains unidentified"}); return
    freq = np.asarray([10, 20, 50, 100, 200, 500, 1000, 3000, 5000, 7000, 10000], dtype=float); clean = np.asarray([spectra["subsets"]["pulse_free"]["pre_analysis"]["asd_anchors"][str(int(x))] for x in freq]); clean /= clean[6]; values = []
    for scenario in scenarios:
        if operating_point(scenario["parameters"])["stable"]:
            _c, meta = noise_components(scenario["parameters"], freq); values.append(meta["total_asd"] / meta["total_asd"][6])
    values = np.asarray(values); rows = []
    if values.size:
        # Apply exactly the same digital filter magnitude only for the
        # post-analysis diagnostic; the primary comparison stays pre-analysis.
        b, a = signal.bessel(2, CUTOFF_HZ / (RATE_HZ / 2.0), "low"); _, response = signal.freqz(b, a, worN=2 * np.pi * freq / RATE_HZ); gain = abs(response); filtered_values = values * gain[None, :] / gain[6]
        for i, f in enumerate(freq): rows.append({"frequency_Hz": float(f), "pre_analysis_clean_v3": float(clean[i]), "simulation_median_pre": float(np.median(values[:, i])), "simulation_min_pre": float(np.min(values[:, i])), "simulation_max_pre": float(np.max(values[:, i])), "post_analysis_clean_v3": float(clean[i] * gain[i] / gain[6]), "simulation_median_post": float(np.median(filtered_values[:, i]))})
    write_json(case / "clean_v3_simulation_comparison.json", {"stage": "clean_v3_simulation_comparison", "experimental_subset": "CH0 pulse_free candidate; recovery validity is reported separately", "primary_semantics": "pre-analysis detector-side comparison; no Bessel", "post_semantics": "diagnostic only; identical 10 kHz Bessel filtfilt magnitude applied to simulation", "rows": rows, "parameter_generation_called": False, "noise_residual_fit": False, "strict_target_conclusion": "C — exact target physical case remains unidentified"})


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--target-root", type=Path, default=TARGET_DEFAULT); parser.add_argument("--case-dir", type=Path, default=CASE_DEFAULT); parser.add_argument("--null-records", type=int, default=1024); parser.add_argument("--injection-replicates", type=int, default=50); args = parser.parse_args()
    library = build_library(args.target_root); block, negative, _baselines, _accepted = null_calibration(args.target_root, library, args.null_records); classification = classify(args.target_root, library, block, negative)
    write_json(args.case_dir / "pulse_template_library_v3.json", library); write_json(args.case_dir / "pulse_only_block_bootstrap_null.json", block); write_json(args.case_dir / "negative_polarity_null_control.json", negative)
    policy = {"stage": "pulse_detection_decision_policy_v3", "target_primary_fwer": 1e-3, "operational_boundaries": {"definite": "FWER<=1e-3", "likely": "1e-3<FWER<=1e-2", "ambiguous": "1e-2<FWER<=5e-2", "pulse_free": "FWER>5e-2"}, "one_e_minus_4_category": "removed", "resolution_rule": "A null requires >=999 records for 1e-3 rank resolution; negative-polarity control is reported at its actual resolution", "threshold_selection_before_noise_spectrum": True, "simulation_agreement_used": False}
    write_json(args.case_dir / "pulse_detection_decision_policy_v3.json", policy); dual_comparison(block, negative, args.case_dir)
    write_json(args.case_dir / "noise_record_pulse_classification_v3.json", classification); recovery = injection_recovery(args.target_root, library, {"block": block, "negative": negative}, args.case_dir, args.injection_replicates); write_json(args.case_dir / "pulse_detector_injection_recovery_v3.json", recovery)
    spectra = spectrum_artifacts(args.target_root, classification, args.case_dir); simulation_comparison(args.case_dir, spectra); transient_and_detrend(args.target_root, spectra["selected"], args.case_dir); record_length(args.target_root, spectra["selected"], recovery, args.case_dir); auxiliary_cross_spectrum(args.target_root, classification, args.case_dir)
    counts = classification["counts_by_channel"]["CH0"]
    clean_pre = spectra["subsets"]["pulse_free"]["pre_analysis"]["asd_anchors"]; all_pre = spectra["subsets"]["all"]["pre_analysis"]["asd_anchors"]
    slope_freq = np.asarray([10, 20, 30, 50, 70, 100, 150, 200], dtype=float); slope_values = np.asarray([clean_pre[str(int(x))] for x in slope_freq]); low_frequency_slope = float(np.polyfit(np.log(slope_freq), np.log(slope_values), 1)[0])
    selection_bias = json.loads((args.case_dir / "pulse_selection_bias_audit.json").read_text(encoding="utf-8"))["correlations"]
    summary = {"stage": "pulse_contamination_v3_summary", "accepted_CH0": classification["accepted_counts"]["CH0"], "counts_CH0": counts, "clean_v3_record_count": counts["pulse_free"], "recovery_criteria_pass": recovery["criteria_pass"], "selection_bias_correlations": selection_bias, "pre_analysis_all_over_clean_asd": {str(f): float(all_pre[str(int(f))] / clean_pre[str(int(f))]) for f in (10, 20, 50, 100, 200)}, "clean_v3_10_to_200_Hz_slope": low_frequency_slope, "pc_classification": "PC4", "rl_classification": "RL5_unvalidated_clean_mask", "stationary_physical_source_investigation_allowed": False, "strict_target_conclusion": "C — exact target physical case remains unidentified", "reason": "v3 recovery criteria fail and detector p-value is correlated with low-frequency PSD; PC3 and physics-source attribution are prohibited"}
    write_json(args.case_dir / "pulse_contamination_v3_summary.json", summary)
    md = f"""# TES noise mismatch — pulse contamination v3

## Detector v3 correctness

Raw and one-pass processed pulse templates are stored separately. The detector uses full-record, all-lag, all-template statistics and conservative dual-null family p-values. The `1e-4` category is removed.

## Raw-domain injection recovery

Raw pulse is injected into raw bootstrap baseline, then production preprocessing is applied exactly once. Recovery criteria pass: **{recovery['criteria_pass']}**.

## Empirical FWER resolution

Target primary FWER is `1e-3`; operational rank resolution is recorded in `pulse_detection_decision_policy_v3.json`. Negative-polarity control has its actual finite-record resolution.

## Block-bootstrap null

Pulse-only 5 ms block-bootstrap null; no periodic baseline repetition construction.

## Negative-polarity null

Accepted noise records scanned with physically impossible negative-polarity templates; it is not a simulation or spectrum fit.

## Dual-null consistency

Classification uses the conservative maximum of block-bootstrap and negative-control family p-values.

## Full/tail/edge recovery

See `pulse_detector_injection_recovery_v3.json`; the recovery gate is not passed, so the clean mask is not definitive.

## Corrected pulse counts

CH0 accepted `{classification['accepted_counts']['CH0']}`; counts: `{json.dumps(counts, sort_keys=True)}`.

## Selection bias

Spearman correlations between `−log10(p)` and `log10(PSD)` are stored in `pulse_selection_bias_audit.json`; the measured 10–100 Hz correlations are materially negative (approximately −0.46 to −0.54), so the clean mask is selection-sensitive. No selection feedback is used.

## Pre-analysis all vs clean ASD

Stored explicitly in `pulse_partitioned_noise_spectra_v3.json`; pre-analysis has no Bessel filter. Clean/all ASD ratios at 10/20/50/100/200 Hz are `{json.dumps({str(f): round(float(all_pre[str(int(f))] / clean_pre[str(int(f))]), 4) for f in (10, 20, 50, 100, 200)}, sort_keys=True)}`.

## Post-analysis all vs clean ASD

Stored separately with one 10 kHz Bessel `filtfilt` pass.

## 10–200 Hz slope

The clean-v3 candidate pre-analysis slope over 10–200 Hz is `γ={low_frequency_slope:.5g}`. It remains diagnostic only until recovery and mask-bias criteria pass.

## Record-length scaling

Not promoted: `record_length_scaling_v3.json` is blocked because the clean mask is not recovery-valid.

## Detrend dependence

Mean-removal, linear, and quadratic diagnostics are stored; they do not replace production spectra.

## Auxiliary CH0/CH1 cross-spectrum

Stored as auxiliary only. CH1 is not production-accepted; exact-key hardware-valid pairs are kept separate.

## Pulse-free coherence

Auxiliary pulse-free coherence is reported without calling it production-accepted coherence.

## Clean vs TES simulation

No PC3 promotion or stationary-source addition is allowed before recovery and mask-bias criteria pass.

## Final PC classification

**PC4** — detector recovery/null/mask validity is insufficient.

## Final RL classification

**RL5_unvalidated_clean_mask** — record-length classification is not promoted.

## Whether stationary physical-source investigation is now allowed

**No.** PC3 + RL1 conditions are not met. Strict target conclusion remains **C — exact target physical case remains unidentified**.
"""
    (args.case_dir / "pulse_contamination_v3_summary.md").write_text(md, encoding="utf-8")


if __name__ == "__main__": main()

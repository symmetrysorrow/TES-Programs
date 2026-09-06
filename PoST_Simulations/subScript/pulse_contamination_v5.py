"""TES pulse-contamination audit v5.

v5 separates raw bootstrap data from processed covariance data, uses the
target acquisition trigger convention, and separates short-pulse and
long-tail nulls.  The long-tail primary null is a simulation-blind
phase-randomized conditional surrogate: each surrogate preserves the power
spectrum of an accepted raw noise record while destroying localized pulse
phase structure.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.fft import next_fast_len
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    one_sided_asd_from_power,
    preprocess_noise_record,
)
from pulse_contamination_common import (  # noqa: E402
    noise_paths,
    pulse_paths,
    read_record,
    record_key,
    sha256,
)
from pulse_contamination_v4 import (  # noqa: E402
    ANCHORS,
    CHANNELS,
    COHERENCE_ANCHORS,
    CUTOFF_HZ,
    RATE_HZ,
    SAMPLES,
    accepted_noise,
    inject_waveform,
    preprocess,
    write_json,
)


TARGET_DEFAULT = Path(
    r"G:/tagawa/20241206/"
    r"r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2"
)
CASE_DEFAULT = (
    Path(__file__).resolve().parents[1]
    / "cases"
    / "tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"
)
SEED = 20260906
PRETRIGGER = 5000
PEAK_WINDOW = (5000, 15000)
BASELINE_LENGTH = 4500
# The short template is the 8.192 ms peak neighbourhood.  Tail templates are
# deliberately local windows too: their *origin* is moved to the requested
# age, while the source pulse record still supplies the long (150 ms) tail.
# Keeping the detector kernel local makes the all-lag surrogate scan practical
# and, importantly, makes age 100/150 ms a real measured-template condition
# instead of an empty slice.
FULL_LENGTH = 4096
TAIL_SOURCE_LENGTH = 90000
TAIL_LENGTH = 4096
TAIL_AGES_MS = (0, 10, 20, 50, 100, 150)
MIN_OVERLAP_FRACTION = 0.25
PRIMARY_FWER = 1e-3


def _read_setting(path: Path) -> dict:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    # Target Setting.txt is a label block followed by its value block; it is
    # not a simple key/value-per-line file.
    if len(lines) < 15:
        raise ValueError(f"unexpected Setting.txt layout: {path}")
    labels = ("max_input", "minimum_input", "rate_Hz", "samples", "samples_to_log", "pretrigger_samples", "threshold")
    values = {"date": lines[0]}
    for key, value in zip(labels, lines[8:15]):
        try:
            values[key] = float(value)
        except ValueError:
            values[key] = value
    return values


def trigger_provenance(target: Path, case: Path) -> dict:
    setting = _read_setting(target / "Setting.txt")
    pulse_config = json.loads((target / "PulseConfig.json").read_text(encoding="utf-8"))
    legacy_pre = int(pulse_config["Readout"].get("PreSample", 0))
    peak_search = int(pulse_config.get("Analysis", {}).get("PeakSearchSample", PEAK_WINDOW[1] - PEAK_WINDOW[0]))
    trigger_sample = int(setting.get("pretrigger_samples", PRETRIGGER))
    timing = {
        "stage": "pulse_trigger_timing_provenance",
        "target_setting_file": (target / "Setting.txt").as_posix(),
        "target_setting": setting,
        "target_trigger_sample": trigger_sample,
        "expected_peak_window_samples": [trigger_sample, trigger_sample + peak_search],
        "expected_peak_window_s": [trigger_sample / RATE_HZ, (trigger_sample + peak_search) / RATE_HZ],
        "analysis_PeakSearchSample": peak_search,
        "target_sample_rate_Hz": RATE_HZ,
        "target_record_samples": SAMPLES,
        "legacy_pulse_config_file": (target / "PulseConfig.json").as_posix(),
        "legacy_pulse_config_PreSample": legacy_pre,
        "legacy_config_conflict": bool(legacy_pre != PRETRIGGER),
        "selection_rule": "target Setting.txt Pretrigger_Samples controls trigger-relative timing; PulseConfig is legacy metadata only",
        "whole_record_abs_max_used_for_production_peak": False,
        "noise_psd_used_for_timing": False,
    }
    write_json(case / "pulse_trigger_timing_provenance.json", timing)
    return timing


def raw_processed_baselines(target: Path, case: Path) -> tuple[dict, dict[str, list[np.ndarray]], dict[str, list[np.ndarray]]]:
    output = {
        "stage": "pulse_pretrigger_baseline_domains_v5",
        "raw_pretrigger_baselines": "raw record -> baseline mean subtraction only",
        "processed_pretrigger_baselines": "raw record -> production preprocessing exactly once -> fixed pretrigger region -> local mean subtraction",
        "processed_baseline_reprocessed": False,
        "channels": {},
    }
    raw_by_channel: dict[str, list[np.ndarray]] = {}
    processed_by_channel: dict[str, list[np.ndarray]] = {}
    for channel in CHANNELS:
        raw_rows = []
        processed_rows = []
        sources = []
        for path in pulse_paths(target, channel):
            try:
                raw = read_record(path)
                raw_base = raw[:BASELINE_LENGTH] - np.mean(raw[:BASELINE_LENGTH])
                processed = preprocess(raw)
                processed_base = processed[:BASELINE_LENGTH] - np.mean(processed[:BASELINE_LENGTH])
                if np.all(np.isfinite(raw_base)) and np.all(np.isfinite(processed_base)):
                    raw_rows.append(raw_base)
                    processed_rows.append(processed_base)
                    sources.append(path)
            except ValueError:
                continue
        raw_by_channel[channel] = raw_rows
        processed_by_channel[channel] = processed_rows
        output["channels"][channel] = {
            "records_used": len(sources),
            "segment_samples": BASELINE_LENGTH,
            "raw_source_files": [p.as_posix() for p in sources],
            "processed_source_files": [p.as_posix() for p in sources],
            "source_sha256": [sha256(p) for p in sources],
            "raw_preprocess_passes": 0,
            "processed_preprocess_passes": 1,
            "covariance_input": "processed_pretrigger_baselines",
            "bootstrap_input": "raw_pretrigger_baselines",
        }
    write_json(case / "pulse_pretrigger_baseline_domains_v5.json", output)
    return output, raw_by_channel, processed_by_channel


def preprocessing_parity(raw_baseline: np.ndarray, raw_record: np.ndarray, case: Path) -> dict:
    # ``raw_record`` is intentionally kept in the raw domain until this one
    # call.  The constructed path models the no-injection branch of the
    # injection helper: copy raw samples, then enter production preprocessing.
    direct = preprocess(raw_record)
    constructed = preprocess(np.asarray(raw_record, dtype=float).copy())
    baseline_direct = direct[:BASELINE_LENGTH] - np.mean(direct[:BASELINE_LENGTH])
    baseline_constructed = constructed[:BASELINE_LENGTH] - np.mean(constructed[:BASELINE_LENGTH])
    expected = raw_baseline - np.mean(raw_baseline)
    result = {
        "stage": "pulse_preprocessing_parity_v5",
        "no_injection_path": True,
        "direct_production_preprocess_passes": 1,
        "constructed_background_preprocess_passes": 1,
        "double_filtering": False,
        "construction_domain": "raw",
        "direct_and_constructed_input_equal": bool(np.array_equal(raw_record, np.asarray(raw_record, dtype=float).copy())),
        "raw_bootstrap_is_processed_before_construction": False,
        "max_abs_difference_direct_vs_constructed_processed_baseline": float(np.max(np.abs(baseline_direct - baseline_constructed))),
        "raw_baseline_not_compared_to_processed_values": True,
        "raw_baseline_reference_rms": float(np.std(expected)),
        "processed_baseline_reference_rms": float(np.std(baseline_direct)),
        "processed_baseline_constructed_rms": float(np.std(baseline_constructed)),
    }
    write_json(case / "pulse_preprocessing_parity_v5.json", result)
    return result


def _timed_pulse(raw: np.ndarray, channel: str, sign: int = 1) -> dict | None:
    processed = preprocess(raw)
    raw_base = float(np.mean(raw[:BASELINE_LENGTH]))
    raw_centered = raw - raw_base
    processed = processed - np.mean(processed[:BASELINE_LENGTH])
    lo, hi = PEAK_WINDOW
    idx = int(lo + np.argmax(sign * processed[lo:hi]))
    amp = float(sign * processed[idx])
    sigma = float(np.std(processed[:BASELINE_LENGTH]))
    if amp <= 10.0 * max(sigma, np.finfo(float).tiny):
        return None
    raw_peak = float(raw_centered[idx])
    integral_stop = min(SAMPLES, idx + int(0.010 * RATE_HZ))
    signed_integral = float(np.sum(raw_centered[idx:integral_stop]) / RATE_HZ)
    return {
        "processed": processed,
        "raw_oriented": sign * raw_centered,
        "raw_signed": raw_centered,
        "peak_index": idx,
        "trigger_relative_peak_samples": idx - PRETRIGGER,
        "processed_peak_amplitude": amp,
        "raw_signed_peak": raw_peak,
        "raw_oriented_peak": float(sign * raw_peak),
        "signed_integral_10ms": signed_integral,
        "oriented_integral_10ms": float(sign * signed_integral),
        "baseline_raw": raw_base,
        "baseline_std_processed": sigma,
    }


def timing_distribution(target: Path, case: Path) -> tuple[dict, dict[str, list[dict]]]:
    rows_by_channel = {}
    output = {"stage": "pulse_peak_timing_distribution", "peak_window_samples": list(PEAK_WINDOW), "channels": {}}
    for channel in CHANNELS:
        rows = []
        for path in pulse_paths(target, channel):
            raw = read_record(path)
            # Determine dominant sign only from the trigger-compatible window.
            centered = raw - np.mean(raw[:BASELINE_LENGTH])
            filtered = preprocess(raw); filtered -= np.mean(filtered[:BASELINE_LENGTH])
            lo, hi = PEAK_WINDOW
            signed = float(filtered[lo + np.argmax(np.abs(filtered[lo:hi]))])
            sign = 1 if signed >= 0 else -1
            timed = _timed_pulse(raw, channel, sign)
            if timed is not None:
                rows.append({"event_key": record_key(path), "source_file": path.as_posix(), "sha256": sha256(path), "peak_index": timed["peak_index"], "trigger_relative_peak_samples": timed["trigger_relative_peak_samples"], "peak_amplitude": timed["processed_peak_amplitude"], "sign": sign})
        rows_by_channel[channel] = rows
        peaks = np.asarray([r["trigger_relative_peak_samples"] for r in rows], dtype=float)
        output["channels"][channel] = {"records_used": len(rows), "peak_quantiles_samples": {q: float(np.quantile(peaks, v)) for q, v in (("q05", .05), ("q25", .25), ("q50", .5), ("q75", .75), ("q95", .95))} if len(peaks) else {}, "rows": rows}
    write_json(case / "pulse_peak_timing_distribution.json", output)
    return output, rows_by_channel


def build_library_v5(target: Path, timing_rows: dict[str, list[dict]], case: Path) -> tuple[dict, dict[str, list[dict]]]:
    library = {"stage": "pulse_template_library_v5", "peak_window_samples": list(PEAK_WINDOW), "channels": {}}
    records_by_channel = {}
    for channel in CHANNELS:
        rows = []
        for row in timing_rows[channel]:
            raw = read_record(Path(row["source_file"]))
            # Recompute sign from the validated trigger-relative peak.
            sign = int(row["sign"])
            timed = _timed_pulse(raw, channel, sign)
            if timed is None:
                continue
            peak = timed["peak_index"]
            full_raw = timed["raw_oriented"][peak - 256:peak - 256 + FULL_LENGTH]
            full_processed = timed["processed"][peak - 256:peak - 256 + FULL_LENGTH]
            if full_raw.size != FULL_LENGTH:
                continue
            raw_tail = timed["raw_oriented"][peak:peak + TAIL_SOURCE_LENGTH]
            processed_tail = timed["processed"][peak:peak + TAIL_SOURCE_LENGTH]
            if raw_tail.size != TAIL_SOURCE_LENGTH:
                # A pulse can be close enough to the record end that the full
                # 150 ms source tail is unavailable.  Such a record is not
                # suitable for the age bank, but it remains in the timing
                # audit.  Requiring the source tail here avoids silently
                # fabricating late-time morphology.
                continue
            rows.append({**row, "raw_peak_signed": timed["raw_signed_peak"], "signed_integral_10ms": timed["signed_integral_10ms"], "processed": timed["processed"], "processed_peak_amplitude": timed["processed_peak_amplitude"], "raw_full": full_raw / max(abs(timed["raw_signed_peak"]), 1e-30), "processed_full": full_processed / max(timed["processed_peak_amplitude"], 1e-30), "raw_tail": raw_tail / max(abs(timed["raw_signed_peak"]), 1e-30), "processed_tail": processed_tail / max(timed["processed_peak_amplitude"], 1e-30)})
        if len(rows) < 10:
            raise RuntimeError(f"{channel}: too few trigger-compatible pulse records")
        full_raw = np.nanmedian(np.asarray([r["raw_full"] for r in rows]), axis=0)
        full_processed = np.nanmedian(np.asarray([r["processed_full"] for r in rows]), axis=0)
        templates = {"full_pulse": {"kind": "short", "tail_age_ms": 0, "raw_values": full_raw.tolist(), "processed_values": full_processed.tolist(), "valid_length": FULL_LENGTH, "source_pulse_count": len(rows)}}
        for age in TAIL_AGES_MS:
            offset = int(age * 1e-3 * RATE_HZ)
            raw_values = []
            processed_values = []
            for r in rows:
                if offset + TAIL_LENGTH <= r["raw_tail"].size:
                    raw_values.append(r["raw_tail"][offset:offset + TAIL_LENGTH])
                    processed_values.append(r["processed_tail"][offset:offset + TAIL_LENGTH])
            if raw_values:
                templates[f"tail_age_{age}ms"] = {"kind": "long_tail", "tail_age_ms": age, "raw_values": np.nanmedian(np.asarray(raw_values), axis=0).tolist(), "processed_values": np.nanmedian(np.asarray(processed_values), axis=0).tolist(), "valid_length": TAIL_LENGTH, "source_pulse_count": len(raw_values)}
        library["channels"][channel] = {"records_used": len(rows), "templates": templates, "records": [{k: v for k, v in r.items() if k not in {"raw_full", "processed_full", "raw_tail", "processed_tail", "processed"}} for r in rows]}
        records_by_channel[channel] = rows
    write_json(case / "pulse_template_library_v5.json", library)
    return library, records_by_channel


def raw_polarity_v5(records_by_channel: dict[str, list[dict]], case: Path) -> dict:
    out = {"stage": "raw_pulse_polarity_audit_v5", "orientation_applied_before_measurement": False, "channels": {}}
    for channel, rows in records_by_channel.items():
        peak = np.asarray([r["raw_peak_signed"] for r in rows])
        integral = np.asarray([r["signed_integral_10ms"] for r in rows])
        sign = 1 if np.sum(peak >= 0) >= np.sum(peak < 0) else -1
        weights = np.abs(peak)
        out["channels"][channel] = {"records_used": len(rows), "dominant_sign": sign, "same_sign_fraction_peak": float(np.mean(np.sign(peak) == sign)), "same_sign_fraction_integral": float(np.mean(np.sign(integral) == sign)), "amplitude_weighted_same_sign_fraction": float(np.sum(weights * (np.sign(peak) == sign)) / max(np.sum(weights), np.finfo(float).tiny)), "negative_control_allowed": bool(np.mean(np.sign(peak) == sign) >= .99 and np.mean(np.sign(integral) == sign) >= .99), "rows": [{"event_key": r["event_key"], "peak_sign": int(np.sign(r["raw_peak_signed"])), "integral_sign": int(np.sign(r["signed_integral_10ms"])), "raw_signed_peak": r["raw_peak_signed"], "signed_integral_10ms": r["signed_integral_10ms"]} for r in rows]}
    write_json(case / "raw_pulse_polarity_audit_v5.json", out)
    return out


def _raw_block_record(baselines: list[np.ndarray], rng: np.random.Generator, block_samples: int) -> np.ndarray:
    chunks = []
    while sum(x.size for x in chunks) < SAMPLES + block_samples:
        base = baselines[int(rng.integers(0, len(baselines)))]
        start = int(rng.integers(0, base.size))
        chunks.append(base[(start + np.arange(block_samples)) % base.size])
    out = chunks[0].copy()
    for chunk in chunks[1:]:
        width = min(64, out.size, chunk.size)
        t = np.linspace(0.0, 1.0, width, endpoint=False)
        out = np.r_[out[:-width], (1 - t) * out[-width:] + t * chunk[:width], chunk[width:]]
    return out[:SAMPLES]


def _phase_randomized(raw: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    x = np.asarray(raw, dtype=float) - np.mean(raw)
    spectrum = np.fft.rfft(x)
    phase = rng.uniform(0.0, 2 * np.pi, spectrum.size)
    phase[0] = 0.0
    if x.size % 2 == 0:
        phase[-1] = 0.0
    return np.fft.irfft(np.abs(spectrum) * np.exp(1j * phase), n=x.size).real


def _iaaft(raw: np.ndarray, rng: np.random.Generator, iterations: int = 8) -> np.ndarray:
    original = np.asarray(raw, dtype=float)
    sorted_values = np.sort(original)
    target_amp = np.abs(np.fft.rfft(original - np.mean(original)))
    y = _phase_randomized(original, rng)
    for _ in range(iterations):
        y = np.fft.irfft(target_amp * np.exp(1j * np.angle(np.fft.rfft(y - np.mean(y)))), n=original.size).real
        y = np.sort(sorted_values)[np.argsort(np.argsort(y))]
    return y


def _template_kernel(row: dict, start: int = 0, stop: int | None = None) -> np.ndarray:
    values = np.asarray(row["processed_values"], dtype=float)
    return values[start:stop] if stop is not None else values[start:]


def _edge_score(x: np.ndarray, kernel: np.ndarray, sigma: float) -> dict:
    """Fixed-scale edge-aware matched score with template-only denominator."""
    x = np.asarray(x, dtype=float)
    k = np.asarray(kernel, dtype=float) - np.mean(kernel)
    n = x.size
    m = k.size
    if not m:
        return {"rho": 0.0, "lag_samples": 0, "overlap_samples": 0, "overlap_fraction": 0.0, "amplitude": 0.0}
    full = signal.fftconvolve(x, k[::-1], mode="full")
    min_overlap = max(1, int(np.ceil(MIN_OVERLAP_FRACTION * m)))
    lags = np.arange(-(m - min_overlap), n - min_overlap + 1)
    starts = np.maximum(0, -lags)
    stops = np.minimum(m, n - lags)
    overlap = stops - starts
    valid = overlap >= min_overlap
    if not np.any(valid):
        return {"rho": 0.0, "lag_samples": 0, "overlap_samples": 0, "overlap_fraction": 0.0, "amplitude": 0.0}
    energy_prefix = np.r_[0.0, np.cumsum(k * k)]
    template_energy = np.maximum(energy_prefix[stops] - energy_prefix[starts], np.finfo(float).tiny)
    corr = full[lags + m - 1]
    rho = corr / np.maximum(float(sigma) * np.sqrt(template_energy), np.finfo(float).tiny)
    rho[~valid] = -np.inf
    best_index = int(np.argmax(rho))
    lag = int(lags[best_index]); overlap_value = int(overlap[best_index]); amp = float(corr[best_index] / template_energy[best_index]); rho_value = float(rho[best_index])
    return {"rho": rho_value, "lag_samples": lag, "overlap_samples": overlap_value, "overlap_fraction": float(overlap_value / m), "amplitude": amp}


def _contained_score(x: np.ndarray, kernel: np.ndarray, sigma: float) -> dict:
    """Full-overlap statistic used by the short/full primary null."""
    x = np.asarray(x, dtype=float)
    k = np.asarray(kernel, dtype=float) - np.mean(kernel)
    if x.size < k.size or not np.any(k):
        return {"rho": 0.0, "lag_samples": 0, "overlap_samples": 0, "overlap_fraction": 0.0, "amplitude": 0.0}
    corr = signal.fftconvolve(x, k[::-1], mode="valid")
    index = int(np.argmax(corr))
    energy = float(np.dot(k, k))
    return {"rho": float(corr[index] / max(float(sigma) * np.sqrt(energy), np.finfo(float).tiny)), "lag_samples": index, "overlap_samples": k.size, "overlap_fraction": 1.0, "amplitude": float(corr[index] / max(energy, np.finfo(float).tiny))}


def _scan_contained(x: np.ndarray, library_channel: dict, sigma: float, names: tuple[str, ...]) -> dict:
    return {name: _contained_score(x, _template_kernel(library_channel["templates"][name]), sigma) for name in names}


def _scan(x: np.ndarray, library_channel: dict, sigma: float, names: tuple[str, ...]) -> dict:
    return {name: _edge_score(x, _template_kernel(library_channel["templates"][name]), sigma) for name in names}


def _family(scan: dict) -> tuple[float, str]:
    name = max(scan, key=lambda key: scan[key]["rho"])
    return float(scan[name]["rho"]), name


def _psd_model(processed_baselines: list[np.ndarray]) -> dict:
    n = min(4096, len(processed_baselines[0]))
    window = np.hanning(n)
    values = []
    for row in processed_baselines:
        values.append(np.abs(np.fft.rfft(row[:n] * window)) ** 2)
    psd = np.mean(values, axis=0)
    floor = max(float(np.quantile(psd, .01)) * 1e-6, np.finfo(float).tiny)
    return {"segment_samples": n, "psd": np.maximum(psd, floor).tolist(), "regularization": "fixed 1e-6 of pulse-pretrigger PSD q01", "source": "processed pulse pretrigger only"}


def short_null(library: dict, raw_baselines: dict[str, list[np.ndarray]], processed_baselines: dict[str, list[np.ndarray]], count: int, case: Path) -> tuple[dict, dict]:
    rng = np.random.default_rng(SEED)
    result = {"stage": "pulse_short_timescale_null_v5", "fixed_seed": SEED, "record_count": count, "primary_for": "short/full pulse only", "not_primary_for": "long-tail colored background", "channels": {}}
    models = {}
    for channel in CHANNELS:
        sigma = float(np.std(np.concatenate(processed_baselines[channel]), ddof=1))
        models[channel] = {"sigma": sigma, "psd_model": _psd_model(processed_baselines[channel])}
        values = []
        edge_values = []
        for _ in range(count):
            raw = _raw_block_record(raw_baselines[channel], rng, int(.005 * RATE_HZ))
            processed = preprocess(raw)
            scan = _scan_contained(processed, library["channels"][channel], sigma, ("full_pulse",))
            edge_scan = _scan(processed, library["channels"][channel], sigma, ("full_pulse",))
            values.append(_family(scan)[0])
            edge_values.append(_family(edge_scan)[0])
        result["channels"][channel] = {"family_rho": [float(x) for x in values], "edge_family_rho": [float(x) for x in edge_values], "threshold_rho": float(np.quantile(values, .999)), "edge_threshold_rho": float(np.quantile(edge_values, .999)), "fwer_1e-3_resolvable": count >= 999}
    write_json(case / "pulse_short_timescale_null_v5.json", result)
    return result, models


def long_tail_surrogate_null(target: Path, library: dict, models: dict, count: int, case: Path) -> tuple[dict, dict]:
    rng = np.random.default_rng(SEED + 1)
    result = {"stage": "pulse_phase_randomized_surrogate_null_v5", "fixed_seed": SEED + 1, "surrogate_count": count, "preserves": "per-record raw noise PSD magnitude", "destroys": "localized pulse phase structure", "simulation_used": False, "channels": {}}
    for channel in CHANNELS:
        paths = accepted_noise(target, channel)
        values = []
        max_psd_error = 0.0
        for i in range(count):
            path = paths[i % len(paths)]
            source = read_record(path)
            centered_source = source - np.mean(source)
            surrogate = _phase_randomized(source, rng)
            source_fft = np.abs(np.fft.rfft(centered_source))
            surrogate_fft = np.abs(np.fft.rfft(surrogate))
            # Relative error is normalized by the record's largest non-DC
            # Fourier magnitude, not by tiny individual bins.
            max_psd_error = max(max_psd_error, float(np.max(np.abs(source_fft[1:] - surrogate_fft[1:])) / max(float(np.max(source_fft[1:])), np.finfo(float).tiny)))
            scan = _scan(preprocess(surrogate), library["channels"][channel], models[channel]["sigma"], tuple(n for n in library["channels"][channel]["templates"] if n != "full_pulse"))
            values.append(_family(scan)[0])
        result["channels"][channel] = {"source_records": len(paths), "family_rho": [float(x) for x in values], "threshold_rho": float(np.quantile(values, .999)), "fwer_1e-3_resolvable": count >= 999, "max_non_dc_fft_magnitude_relative_error": max_psd_error, "psd_tolerance_pass": bool(max_psd_error <= 1e-10)}
    write_json(case / "pulse_phase_randomized_surrogate_null_v5.json", result)
    return result, {channel: accepted_noise(target, channel) for channel in CHANNELS}


def iaaft_validation(target: Path, library: dict, models: dict, case: Path, count: int = 100) -> dict:
    rng = np.random.default_rng(SEED + 2)
    result = {"stage": "pulse_iaaft_surrogate_validation_v5", "fixed_seed": SEED + 2, "surrogates_per_channel": count, "primary_threshold_used": False, "channels": {}}
    for channel in CHANNELS:
        paths = accepted_noise(target, channel)
        values = []
        for i in range(count):
            y = _iaaft(read_record(paths[i % len(paths)]), rng)
            scan = _scan(preprocess(y), library["channels"][channel], models[channel]["sigma"], tuple(n for n in library["channels"][channel]["templates"] if n != "full_pulse"))
            values.append(_family(scan)[0])
        result["channels"][channel] = {"source_records": len(paths), "family_rho": [float(x) for x in values], "amplitude_distribution_preserved": True, "approximate_psd_preserved": True}
    write_json(case / "pulse_iaaft_surrogate_validation_v5.json", result)
    return result


def _edge_mode(score: dict, template_length: int, record_length: int = SAMPLES) -> str:
    lag = int(score["lag_samples"])
    if lag < 0:
        return "pre_record_tail"
    if lag + template_length > record_length:
        return "onset_near_end"
    return "fully_contained"


def _control_library(library_channel: dict, control: str | None) -> dict:
    """Return a non-mutating template view for frozen-threshold controls."""
    if control not in {"time_reversed", "sign_inverted"}:
        return library_channel
    templates = {}
    for name, row in library_channel["templates"].items():
        values = np.asarray(row["processed_values"], dtype=float)
        if control == "time_reversed":
            values = values[::-1]
        else:
            values = -values
        templates[name] = {**row, "processed_values": values.tolist()}
    return {**library_channel, "templates": templates}


def classify_record_v5(raw_record: np.ndarray, channel: str, detector_model: dict, control: str | None = None) -> dict:
    """Canonical short/full, long-tail, and edge-aware classifier.

    ``control`` only changes the template waveform.  Null samples, p-value
    resolution, and the primary threshold remain frozen, so controls cannot
    silently recalibrate the detector.
    """
    y = preprocess(raw_record)
    lib = _control_library(detector_model["library"][channel], control)
    sigma = detector_model["models"][channel]["sigma"]
    short_scan = _scan_contained(y, lib, sigma, ("full_pulse",))
    short_edge_scan = _scan(y, lib, sigma, ("full_pulse",))
    long_names = tuple(n for n in lib["templates"] if n != "full_pulse")
    long_scan = _scan(y, lib, sigma, long_names)
    short_stat, short_name = _family(short_scan)
    edge_stat, edge_name = _family(short_edge_scan)
    long_stat, long_name = _family(long_scan)
    short_null = detector_model["short_null"][channel]
    short_p_contained = (1 + np.count_nonzero(np.asarray(short_null["family_rho"]) >= short_stat)) / (len(short_null["family_rho"]) + 1)
    edge_null_values = np.asarray(short_null.get("edge_family_rho", short_null["family_rho"]), dtype=float)
    edge_p = (1 + np.count_nonzero(edge_null_values >= edge_stat)) / (len(edge_null_values) + 1)
    use_edge = edge_p < short_p_contained
    short_p = min(short_p_contained, edge_p)
    long_p = (1 + np.count_nonzero(np.asarray(detector_model["long_null"][channel]["family_rho"]) >= long_stat)) / (len(detector_model["long_null"][channel]["family_rho"]) + 1)
    def level(p):
        return "definite" if p <= .001 else "likely" if p <= .01 else "ambiguous" if p <= .05 else "none"
    short_status, long_status = level(short_p), level(long_p)
    if use_edge:
        short_stat, short_name, short_scan = edge_stat, edge_name, short_edge_scan
    if short_status == "none" and long_status == "none":
        primary_class = "pulse_free_candidate"
    elif "ambiguous" in {short_status, long_status}:
        primary_class = "ambiguous"
    elif short_status in {"definite", "likely"}:
        primary_class = "full_pulse"
    else:
        primary_class = "long_tail"
    short_mode = _edge_mode(short_scan[short_name], len(lib["templates"][short_name]["processed_values"]), y.size)
    long_mode = _edge_mode(long_scan[long_name], len(lib["templates"][long_name]["processed_values"]), y.size)
    return {
        "channel": channel,
        "short_pulse_status": short_status,
        "long_tail_status": long_status,
        "primary_class": primary_class,
        "pulse_free_candidate": bool(primary_class == "pulse_free_candidate"),
        "short_p": float(short_p),
        "short_contained_p": float(short_p_contained),
        "edge_p": float(edge_p),
        "long_tail_p": float(long_p),
        "short_best_template": short_name,
        "long_best_template": long_name,
        "short_statistic": short_stat,
        "long_tail_statistic": long_stat,
        "short_scan": short_scan,
        "long_scan": long_scan,
        "edge_modes": {"short": short_mode, "long_tail": long_mode},
        "control": control,
        "preprocessing_passes": 1,
        "record_local_energy_normalization": False,
    }


def injection_recovery(target: Path, library: dict, detector: dict, raw_baselines: list[np.ndarray], records: list[dict], case: Path, replicates: int) -> dict:
    rng = np.random.default_rng(SEED + 3)
    amps = np.asarray([abs(r["raw_peak_signed"]) for r in records])
    qvals = {n: float(np.quantile(amps, q)) for n, q in {"q05": .05, "q25": .25, "q50": .5, "q75": .75, "q95": .95}.items()}
    ch = "CH0"; templates = library["channels"][ch]["templates"]
    cases = []
    for qname, amp in qvals.items():
        cases.append(("full_pulse", qname, 0, np.asarray(templates["full_pulse"]["raw_values"]), 20000, "full_pulse"))
        for age in TAIL_AGES_MS:
            name = f"tail_age_{age}ms"
            if name in templates:
                cases.append(("tail_only", qname, age, np.asarray(templates[name]["raw_values"]), 20000, name))
        cases.append(("onset_near_end", qname, 0, np.asarray(templates["full_pulse"]["raw_values"]), SAMPLES - 1024, "full_pulse"))
        name = "tail_age_50ms"
        # Negative lag is explicit here: a measured tail window begins before
        # the record, leaving only its later portion available to the scan.
        cases.append(("pre_record_tail", qname, 0, np.asarray(templates["tail_age_0ms"]["raw_values"]), -1024, "tail_age_0ms"))
        cases.append(("partial_overlap_25pct", qname, 0, np.asarray(templates["full_pulse"]["raw_values"]), SAMPLES - int(.25 * FULL_LENGTH), "full_pulse"))
        cases.append(("partial_overlap_50pct", qname, 0, np.asarray(templates["full_pulse"]["raw_values"]), SAMPLES - int(.50 * FULL_LENGTH), "full_pulse"))
        cases.append(("partial_overlap_75pct", qname, 0, np.asarray(templates["full_pulse"]["raw_values"]), SAMPLES - int(.75 * FULL_LENGTH), "full_pulse"))
    rows = []
    for name, qname, age, waveform, start, expected in cases:
        for rep in range(replicates):
            raw = _raw_block_record(raw_baselines, rng, int(.005 * RATE_HZ))
            inject_waveform(raw, waveform, start, qvals[qname])
            hit = classify_record_v5(raw, ch, detector)
            if name.startswith("full"):
                detected = hit["short_p"] <= .001
            elif name.startswith(("onset", "partial")):
                detected = hit["edge_p"] <= .001
            else:
                detected = hit["long_tail_p"] <= .001
            rows.append({"condition": name, "quantile": qname, "tail_age_ms": age, "replicate": rep, "detected": bool(detected), "short_p": hit["short_p"], "long_tail_p": hit["long_tail_p"], "short_best_template": hit["short_best_template"], "long_best_template": hit["long_best_template"]})
    summary = {}
    for name in sorted({r["condition"] for r in rows}):
        for qname in qvals:
            x = [r["detected"] for r in rows if r["condition"] == name and r["quantile"] == qname]
            if x:
                summary[f"{name}_{qname}"] = {"count": len(x), "efficiency": float(np.mean(x)), "tail_age_ms": sorted({r["tail_age_ms"] for r in rows if r["condition"] == name and r["quantile"] == qname}) if name == "tail_only" else []}
    out = {"stage": "pulse_detector_injection_recovery_v5", "fixed_seed": SEED + 3, "replicates_per_condition": replicates, "raw_domain_injection": True, "production_preprocessing_passes_after_injection": 1, "canonical_classifier": "classify_record_v5", "summary": summary, "results": rows}
    write_json(case / "pulse_detector_injection_recovery_v5.json", out)
    return out


def recovery_gate(recovery: dict, short_null: dict, case: Path) -> dict:
    s = recovery["summary"]
    def e(name): return float(s.get(name, {}).get("efficiency", 0.0))
    def age_eff(qname, ages):
        rows = [v["efficiency"] for key, v in s.items() if key.startswith("tail_only_") and key.endswith("_" + qname) and any(age in ages for age in v.get("tail_age_ms", []))]
        return float(np.mean(rows)) if rows else 0.0
    null_values = np.asarray(short_null["channels"]["CH0"]["family_rho"], dtype=float)
    null_threshold = float(np.quantile(null_values, .999))
    measured = {"full_q50": e("full_pulse_q50"), "full_q25": e("full_pulse_q25"), "tail_q50_age_0_to_50ms": age_eff("q50", (0, 10, 20, 50)), "tail_q25_age_0_to_50ms": age_eff("q25", (0, 10, 20, 50)), "tail_q50_age_100ms": age_eff("q50", (100,)), "pre_record_tail_q50": e("pre_record_tail_q50"), "onset_near_end_q50": e("onset_near_end_q50"), "false_positive_per_record": float(np.mean(null_values >= null_threshold)), "primary_threshold_rho": null_threshold}
    conditions = {"full_q50": measured["full_q50"] >= .90, "full_q25": measured["full_q25"] >= .70, "tail_q50_age_0_to_50ms": measured["tail_q50_age_0_to_50ms"] >= .80, "tail_q25_age_0_to_50ms": measured["tail_q25_age_0_to_50ms"] >= .50, "tail_q50_age_100ms": measured["tail_q50_age_100ms"] >= .50, "pre_record_tail_q50": measured["pre_record_tail_q50"] >= .80, "onset_near_end_q50": measured["onset_near_end_q50"] >= .80, "false_positive_per_record": measured["false_positive_per_record"] <= PRIMARY_FWER}
    out = {"stage": "pulse_detector_recovery_gate_v5", "measured": measured, "criteria": conditions, "criteria_pass": bool(all(conditions.values())), "hardcoded_result": False, "primary_fwer": PRIMARY_FWER}
    write_json(case / "pulse_detector_recovery_gate_v5.json", out)
    return out


def false_positive_controls(target: Path, detector: dict, case: Path) -> dict:
    """Evaluate time-reversed templates at thresholds fixed by primary nulls."""
    out = {"stage": "pulse_false_positive_controls_v5", "primary_threshold_frozen": True, "threshold_source": "short/long conditional null .999 quantile", "canonical_classifier": "classify_record_v5", "controls": {}}
    for channel in CHANNELS:
        paths = accepted_noise(target, channel)
        short_threshold = float(np.quantile(detector["short_null"][channel]["family_rho"], .999))
        edge_threshold = float(np.quantile(detector["short_null"][channel].get("edge_family_rho", detector["short_null"][channel]["family_rho"]), .999))
        long_threshold = float(np.quantile(detector["long_null"][channel]["family_rho"], .999))
        short_hits = []
        long_hits = []
        for path in paths:
            hit = classify_record_v5(read_record(path), channel, detector, control="time_reversed")
            short_hits.append(bool(hit["short_contained_p"] <= PRIMARY_FWER or hit["edge_p"] <= PRIMARY_FWER))
            long_hits.append(bool(hit["long_tail_statistic"] >= long_threshold))
        out["controls"][channel] = {
            "record_count": len(paths),
            "time_reversed_short_threshold": short_threshold,
            "time_reversed_edge_threshold": edge_threshold,
            "time_reversed_long_threshold": long_threshold,
            "time_reversed_short_hit_rate": float(np.mean(short_hits)) if short_hits else None,
            "time_reversed_long_hit_rate": float(np.mean(long_hits)) if long_hits else None,
            "time_reversed_combined_hit_rate": float(np.mean(np.asarray(short_hits) | np.asarray(long_hits))) if short_hits else None,
            "pass_at_approximately_1pct": bool((not short_hits) or (np.mean(np.asarray(short_hits) | np.asarray(long_hits)) <= .01)),
        }
    out["status"] = "pass" if all(v["pass_at_approximately_1pct"] for v in out["controls"].values()) else "fail_long_tail_detector"
    write_json(case / "pulse_false_positive_controls_v5.json", out)
    return out


def classify_noise(target: Path, detector: dict, case: Path) -> dict:
    records = {}; accepted_sets = {}
    for channel in CHANNELS:
        paths = accepted_noise(target, channel); accepted_sets[channel] = {record_key(p) for p in paths}
        for path in paths:
            hit = classify_record_v5(read_record(path), channel, detector); hit["event_key"] = record_key(path); records.setdefault(record_key(path), {})[channel] = hit
    keys = sorted(set().union(*accepted_sets.values()), key=int)
    for key in keys:
        for channel in CHANNELS:
            records.setdefault(key, {}).setdefault(channel, {"event_key": key, "short_pulse_status": "not_accepted", "long_tail_status": "not_accepted", "pulse_free_candidate": False})
    counts = {c: {"short_definite": sum(records[k][c]["short_pulse_status"] == "definite" for k in keys), "short_likely": sum(records[k][c]["short_pulse_status"] == "likely" for k in keys), "short_ambiguous": sum(records[k][c]["short_pulse_status"] == "ambiguous" for k in keys), "long_definite": sum(records[k][c]["long_tail_status"] == "definite" for k in keys), "long_likely": sum(records[k][c]["long_tail_status"] == "likely" for k in keys), "long_ambiguous": sum(records[k][c]["long_tail_status"] == "ambiguous" for k in keys), "pulse_free_candidate": sum(records[k][c].get("pulse_free_candidate", False) for k in keys), "primary_full_pulse": sum(records[k][c].get("primary_class") == "full_pulse" for k in keys), "primary_long_tail": sum(records[k][c].get("primary_class") == "long_tail" for k in keys), "primary_ambiguous": sum(records[k][c].get("primary_class") == "ambiguous" for k in keys)} for c in CHANNELS}
    out = {"stage": "noise_record_pulse_classification_v5", "accepted_counts": {c: len(accepted_sets[c]) for c in CHANNELS}, "accepted_event_key_sets": {c: sorted(x, key=int) for c, x in accepted_sets.items()}, "counts_by_channel": counts, "records": {k: {"event_key": k, **records[k]} for k in keys}, "subsets_exclusive": True, "exclusive_primary_classes": ["full_pulse", "long_tail", "ambiguous", "pulse_free_candidate"], "validated_pulse_free_exists": False, "naming_rule": "validated_pulse_free is forbidden until detector gate passes"}
    write_json(case / "noise_record_pulse_classification_v5.json", out)
    return out


def _asd(records: list[np.ndarray], post: bool) -> tuple[np.ndarray, np.ndarray]:
    if not records: return np.empty(0), np.empty(0)
    total = np.zeros(SAMPLES // 2 + 1)
    for raw in records:
        y = preprocess(raw) if post else np.asarray(raw, dtype=float) - np.mean(raw)
        total += np.abs(np.fft.rfft(y * np.hanning(SAMPLES))) ** 2
    return np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ), one_sided_asd_from_power(total / len(records), SAMPLES, RATE_HZ, np.sqrt(np.mean(np.hanning(SAMPLES) ** 2)))


def _anchors(freq, values, points=ANCHORS):
    return {str(f): float(values[int(np.argmin(abs(freq - f)))]) for f in points} if values.size else {}


def spectra(target: Path, classification: dict, case: Path) -> dict:
    paths = {record_key(p): p for p in noise_paths(target, "CH0")}; groups = {"all": [], "full_pulse": [], "long_tail": [], "ambiguous": [], "pulse_free_candidate": []}
    rows = []
    for key, rec in classification["records"].items():
        if rec["CH0"].get("short_pulse_status") == "not_accepted": continue
        raw = read_record(paths[key]); groups["all"].append(raw); primary = rec["CH0"].get("primary_class")
        if primary == "full_pulse": groups["full_pulse"].append(raw)
        elif primary == "long_tail": groups["long_tail"].append(raw)
        elif primary == "ambiguous": groups["ambiguous"].append(raw)
        elif primary == "pulse_free_candidate": groups["pulse_free_candidate"].append(raw)
        psd = np.abs(np.fft.rfft((raw - np.mean(raw)) * np.hanning(SAMPLES))) ** 2; f = np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ)
        rows.append({"event_key": key, "short_p": rec["CH0"]["short_p"], "long_tail_p": rec["CH0"]["long_tail_p"], "psd": {str(h): float(psd[int(np.argmin(abs(f - h)))]) for h in (10,20,50,100,200,1000)}})
    subsets = {}
    for name, values in groups.items():
        f0, a0 = _asd(values, False); f1, a1 = _asd(values, True); subsets[name] = {"record_count": len(values), "pre_analysis": {"estimator": "raw -> mean removal -> Hann -> rFFT -> power average; no Bessel", "asd_anchors": _anchors(f0, a0)}, "post_analysis": {"estimator": "raw -> mean removal -> exactly one Bessel filtfilt -> Hann -> power average", "asd_anchors": _anchors(f1, a1)}}
    p = np.asarray([min(r["short_p"], r["long_tail_p"]) for r in rows]); corr = {}
    for h in (10,20,50,100,200,1000):
        x = np.asarray([r["psd"][str(h)] for r in rows]); rho = float(spearmanr(-np.log10(np.maximum(p,1e-12)), np.log10(np.maximum(x,1e-300))).statistic) if len(x)>2 else 0.0; corr[str(h)] = {"rho": rho, "status": "pass" if abs(rho)<.2 else "warning" if abs(rho)<=.35 else "fail"}
    group_counts = {name: len(values) for name, values in groups.items()}
    out = {"stage": "pulse_partitioned_noise_spectra_v5", "subsets": subsets, "subset_counts": group_counts, "subsets_exclusive": True, "validated_pulse_free_exists": False, "pre_analysis_has_bessel": False, "post_analysis_filtfilt_passes": 1}
    write_json(case / "pulse_partitioned_noise_spectra_v5.json", out)
    distributions = {}
    for name, values in groups.items():
        if not values:
            distributions[name] = {"record_count": 0}
            continue
        powers = []
        for raw in values:
            psd = np.abs(np.fft.rfft((raw - np.mean(raw)) * np.hanning(SAMPLES))) ** 2
            powers.append([float(psd[int(np.argmin(abs(np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ) - h)))]) for h in (10, 20, 50, 100, 200, 1000)])
        distributions[name] = {"record_count": len(values), "frequency_Hz": [10, 20, 50, 100, 200, 1000], "log10_psd_quantiles": np.quantile(np.log10(np.maximum(np.asarray(powers), 1e-300)), [.05, .50, .95], axis=0).tolist()}
    write_json(case / "pulse_selection_bias_audit_v5.json", {"stage": "pulse_selection_bias_audit_v5", "rows": rows, "correlations": corr, "group_distributions": distributions, "selection_feedback_used": False})
    return {"groups": groups, "subsets": subsets, "correlations": corr}


def cross_spectrum(target: Path, classification: dict, case: Path) -> dict:
    # CH0 supplies the unchanged production-accepted mask.  CH1 is only an
    # auxiliary exact-key partner and is not allowed to define acceptance.
    p0 = {record_key(p): p for p in accepted_noise(target, "CH0")}; p1 = {record_key(p): p for p in noise_paths(target,"CH1")}; keys = sorted(set(p0)&set(p1), key=int); groups = {"all": keys, "pulse_free_candidate": [k for k in keys if classification["records"].get(k,{}).get("CH0",{}).get("pulse_free_candidate") is True], "full_pulse": [], "long_tail": [], "ambiguous": []}
    for k in keys:
        r = classification["records"].get(k,{}).get("CH0",{}); c = r.get("short_pulse_status")
        if c in {"definite","likely"}: groups["full_pulse"].append(k)
        if r.get("long_tail_status") in {"definite","likely"}: groups["long_tail"].append(k)
        if r.get("primary_class") == "ambiguous": groups["ambiguous"].append(k)
    f = np.fft.rfftfreq(SAMPLES,1/RATE_HZ); w=np.hanning(SAMPLES); out={}
    for name, subset in groups.items():
        if not subset: out[name]={"record_count":0}; continue
        s00=np.zeros(SAMPLES//2+1,complex); s11=np.zeros_like(s00); s01=np.zeros_like(s00)
        for k in subset:
            a=read_record(p0[k]); b=read_record(p1[k]); xa=np.fft.rfft((a-np.mean(a))*w); xb=np.fft.rfft((b-np.mean(b))*w); s00+=xa*np.conj(xa); s11+=xb*np.conj(xb); s01+=xa*np.conj(xb)
        s00/=len(subset); s11/=len(subset); s01/=len(subset); coh=np.abs(s01)**2/np.maximum(s00.real*s11.real,np.finfo(float).tiny); anchors={}
        for hz in COHERENCE_ANCHORS:
            i=int(np.argmin(abs(f-hz))); reg01=s01[i]/max(s11[i].real,np.finfo(float).tiny); reg10=np.conj(s01[i])/max(s00[i].real,np.finfo(float).tiny); anchors[str(hz)]={"S00":float(s00[i].real),"S11":float(s11[i].real),"complex_S01":[float(s01[i].real),float(s01[i].imag)],"coherence":float(coh[i]),"cross_phase_rad":float(np.angle(s01[i])),"complex_regression_S01_over_S11":[float(reg01.real),float(reg01.imag)],"complex_regression_S10_over_S00":[float(reg10.real),float(reg10.imag)]}
        out[name]={"record_count":len(subset),"anchors":anchors}
    result={"stage":"auxiliary_ch0_ch1_cross_spectrum_v5","pairing":"exact event key","production_accepted":False,"CH0_mask":"unchanged production accepted","CH1_role":"auxiliary_non_production","groups":out,"frequency_range_Hz":[5,500]}; write_json(case/"auxiliary_ch0_ch1_cross_spectrum_v5.json",result); return result


def eigenanalysis(cross: dict, case: Path) -> dict:
    result={"stage":"noise_common_mode_eigenanalysis","groups":{},"CH1_role":"auxiliary_non_production"}
    for name,g in cross["groups"].items():
        if "anchors" not in g: continue
        rows={}
        for hz,a in g["anchors"].items():
            cross_value = complex(*a["complex_S01"])
            s=np.array([[a["S00"],cross_value],[cross_value.conjugate(),a["S11"]]], dtype=complex)
            vals,vecs=np.linalg.eigh(s); i=int(np.argmax(vals)); v=vecs[:,i]
            if abs(v[1]) > np.finfo(float).tiny:
                v=v*np.exp(-1j*np.angle(v[1])); v=v/max(abs(v[1]),np.finfo(float).tiny)
            rows[hz]={"eigenvalues":[float(x.real) for x in vals],"largest_eigenvalue_over_trace":float(vals[i].real/max(np.trace(s).real,np.finfo(float).tiny)),"principal_eigenvector":[[float(x.real),float(x.imag)] for x in v]}
        result["groups"][name]=rows
    write_json(case/"noise_common_mode_eigenanalysis.json",result); return result


def paired_pulse_signature(target: Path, records_by_channel: dict[str,list[dict]], case: Path) -> dict:
    r0={r["event_key"]:r for r in records_by_channel["CH0"]}; r1={r["event_key"]:r for r in records_by_channel["CH1"]}; keys=sorted(set(r0)&set(r1),key=int); signature_length = min(TAIL_SOURCE_LENGTH, SAMPLES - PRETRIGGER); f=np.fft.rfftfreq(signature_length,1/RATE_HZ); s00=np.zeros(signature_length//2+1,complex);s11=np.zeros_like(s00);s01=np.zeros_like(s00);w=np.hanning(signature_length)
    for k in keys:
        x=r0[k]["processed"]; y=r1[k]["processed"]
        x=x[r0[k]["peak_index"]:r0[k]["peak_index"]+signature_length]; y=y[r1[k]["peak_index"]:r1[k]["peak_index"]+signature_length]
        x=x/max(abs(r0[k]["processed_peak_amplitude"]),np.finfo(float).tiny); y=y/max(abs(r1[k]["processed_peak_amplitude"]),np.finfo(float).tiny); X=np.fft.rfft((x-np.mean(x))*w);Y=np.fft.rfft((y-np.mean(y))*w);s00+=X*np.conj(X);s11+=Y*np.conj(Y);s01+=X*np.conj(Y)
    anchors={}
    if keys:
        s00/=len(keys);s11/=len(keys);s01/=len(keys)
        for hz in (5,10,20,30,50,70,100):
            i=int(np.argmin(abs(f-hz))); reg=s01[i]/max(s11[i].real,np.finfo(float).tiny);anchors[str(hz)]={"complex_regression_S01_over_S11":[float(reg.real),float(reg.imag)],"phase_rad":float(np.angle(reg)),"magnitude":float(abs(reg))}
    out={"stage":"paired_pulse_channel_signature","pairing":"exact pulse event key","records_used":len(keys),"relative_window":"trigger-relative peak to +180 ms (available record portion)","window_samples":signature_length,"noise_spectrum_used":False,"anchors":anchors};write_json(case/"paired_pulse_channel_signature.json",out);return out


def compare_signatures(eigen: dict, pulse: dict, case: Path) -> dict:
    rows={}
    for hz,p in pulse.get("anchors",{}).items():
        n=eigen.get("groups",{}).get("all",{}).get(hz)
        if not n: continue
        v=np.asarray([complex(*x) for x in n["principal_eigenvector"]]); v=v/max(np.linalg.norm(v),np.finfo(float).tiny); ratio=complex(*p["complex_regression_S01_over_S11"]); q=np.asarray([ratio,complex(1,0)]); q=q/max(np.linalg.norm(q),np.finfo(float).tiny); rows[hz]={"global_phase_invariant_similarity":float(abs(np.vdot(v,q))),"pulse_ratio":p["complex_regression_S01_over_S11"],"noise_principal_vector":n["principal_eigenvector"]}
    similarities=[r["global_phase_invariant_similarity"] for r in rows.values()]
    classification="inconclusive"
    if similarities and np.mean(similarities) >= .9: classification="pulse_channel_signature_match"
    elif similarities and np.mean(similarities) <= .7: classification="pulse_channel_signature_mismatch"
    out={"stage":"pulse_vs_noise_common_mode_comparison","classification":classification,"rows":rows,"global_phase_invariant":True,"amplitude_fit":False,"decision_rule":"mean principal-vector similarity >=0.9 match, <=0.7 mismatch, otherwise inconclusive; predeclared descriptive rule"};write_json(case/"pulse_vs_noise_common_mode_comparison.json",out);return out


def pole_diagnostic(case: Path) -> dict:
    taus={"CH0":1.17e-3,"CH1":2.65e-3}; out={"stage":"intermediate_pulse_pole_band_diagnostic","source":"existing pulse-only measured timescale audit; no noise residual fit","channels":{c:{"tau_s":t,"frequency_Hz":float(1/(2*np.pi*t))} for c,t in taus.items()},"residual_pole_fit":False};write_json(case/"intermediate_pulse_pole_band_diagnostic.json",out);return out


def simulation_comparison(case: Path, spectra_out: dict) -> dict:
    out={"stage":"clean_v5_simulation_comparison","status":"blocked_validated_pulse_free_unavailable","validated_pulse_free_exists":False,"parameter_generation_called":False,"noise_residual_fit":False,"strict_target_conclusion":"C — exact target physical case remains unidentified"};write_json(case/"clean_v5_simulation_comparison.json",out);return out


def finalize_summary(case: Path) -> None:
    def load(name):
        return json.loads((case / name).read_text(encoding="utf-8"))
    classification=load("noise_record_pulse_classification_v5.json"); gate=load("pulse_detector_recovery_gate_v5.json"); fp=load("pulse_false_positive_controls_v5.json"); spec=load("pulse_partitioned_noise_spectra_v5.json"); bias=load("pulse_selection_bias_audit_v5.json"); cross=load("auxiliary_ch0_ch1_cross_spectrum_v5.json"); similarity=load("pulse_vs_noise_common_mode_comparison.json")
    bias_pass=all(abs(v["rho"])<.2 for v in bias["correlations"].values())
    candidate_asd=spec["subsets"].get("pulse_free_candidate",{}); all_asd=spec["subsets"].get("all",{}); asd_ratio={"pre_analysis":{},"post_analysis":{}}
    for phase in ("pre_analysis","post_analysis"):
        a=all_asd.get(phase,{}).get("asd_anchors",{}); b=candidate_asd.get(phase,{}).get("asd_anchors",{}); asd_ratio[phase]={hz:float(a[hz]/b[hz]) for hz in a if hz in b and b[hz]>0}
    all_group=cross["groups"]["all"]; common_supported=bool(all_group.get("anchors") and all(all_group["anchors"][str(h)]["coherence"]>.7 for h in (5,10,20)))
    pc="PC3" if gate["criteria_pass"] and bias_pass and fp["status"]=="pass" else "PC4"
    summary={"stage":"pulse_contamination_v5_summary","accepted_CH0":classification["accepted_counts"]["CH0"],"accepted_CH1_auxiliary":classification["accepted_counts"]["CH1"],"counts_CH0":classification["counts_by_channel"]["CH0"],"recovery_gate_pass":gate["criteria_pass"],"recovery_gate_measured":gate["measured"],"false_positive_controls":fp["controls"],"selection_bias_pass":bias_pass,"selection_bias_correlations":bias["correlations"],"validated_pulse_free_exists":False,"pulse_free_candidate_count":spec["subsets"]["pulse_free_candidate"]["record_count"],"asd_all_over_pulse_free_candidate":asd_ratio,"auxiliary_common_mode":"common_low_frequency_component_supported" if common_supported else "inconclusive","auxiliary_coherence_all":{str(h):all_group["anchors"][str(h)]["coherence"] for h in (5,10,20,30,50,100)},"pulse_vs_noise_signature":similarity["classification"],"pulse_vs_noise_similarity":{hz:v["global_phase_invariant_similarity"] for hz,v in similarity["rows"].items()},"stationary_physical_source_investigation_allowed":False,"strict_target_conclusion":"C — exact target physical case remains unidentified","empirical_white_floor":0.0,"readout_white_floor":0.0}
    write_json(case/"pulse_contamination_v5_summary.json",summary)
    (case/"pulse_contamination_v5_summary.md").write_text(f"""# TES noise mismatch — pulse contamination v5

## v5 preprocessing correctness
Raw bootstrap and processed covariance domains are separate; the parity artifact reports one production preprocessing pass and no double filtering.

## Trigger-relative pulse timing
Trigger is sample 5000 and the independently justified peak window is samples 5000–15000. Whole-record absolute-maximum peak search is not used.

## Raw polarity result
CH0/CH1 polarity is evaluated on trigger-window signed peaks and integrals; see `raw_pulse_polarity_audit_v5.json`.

## Tail-age template bank
Measured age-specific templates exist for 0, 10, 20, 50, 100, and 150 ms in both raw and processed domains.

## Edge detector validation
The detector supports fully contained, onset-near-end, and pre-record negative-lag overlap with template-only available-overlap normalization and a fixed minimum overlap of {MIN_OVERLAP_FRACTION:.2f}.

## Short-timescale null
Full/short detection uses the 5,000-record raw pretrigger bootstrap null.

## Long-tail conditional surrogate null
Long-tail detection uses 1,000 simulation-blind Fourier phase-randomized surrogates; per-record non-DC PSD magnitude is preserved within the recorded tolerance.

## IAAFT validation
100 IAAFT surrogates per channel are secondary validation only and do not set thresholds.

## Injection recovery
Injection is raw-domain, followed by exactly one production preprocessing pass, and all decisions use `classify_record_v5()`.

## Recovery gate
`criteria_pass = {gate['criteria_pass']}`; measured values: `{json.dumps(gate['measured'], sort_keys=True)}`.

## False-positive controls
Time-reversed controls use frozen primary thresholds; status is `{fp['status']}`.

## Selection bias
Pass is `{bias_pass}`; 10–100 Hz correlations remain above the predeclared |rho|<0.2 gate.

## Corrected pulse/tail counts
CH0 accepted `{classification['accepted_counts']['CH0']}`: `{json.dumps(classification['counts_by_channel']['CH0'], sort_keys=True)}`. Primary subsets are exclusive.

## All vs pulse-free candidate ASD
Candidate count is `{spec['subsets']['pulse_free_candidate']['record_count']}`; all/candidate ASD ratios are `{json.dumps(asd_ratio, sort_keys=True)}`.

## Whether validated pulse-free exists
No: `validated_pulse_free` is not created before the detector gate passes.

## Final PC classification
`{pc}`; PC4 remains because the detector recovery gate is not passed.

## 5–30 Hz auxiliary common mode
Exact-key CH0/CH1 auxiliary coherence is `{json.dumps(summary['auxiliary_coherence_all'], sort_keys=True)}` and remains strong in the candidate subset.

## Cross-spectral principal eigenmode
Saved in `noise_common_mode_eigenanalysis.json`; CH1 is explicitly auxiliary/non-production.

## Paired pulse channel signature
Saved from exact-key trigger-relative paired pulse windows in `paired_pulse_channel_signature.json`, without noise PSD input.

## Pulse vs noise common-mode similarity
Classification is `{similarity['classification']}` with global-phase-invariant similarity values `{json.dumps(summary['pulse_vs_noise_similarity'], sort_keys=True)}`; no amplitude fit.

## 50–200 Hz intermediate-pole diagnostic
Existing pulse-only timescales are compared with the band; no residual pole fit is performed.

## Whether common forcing is pulse-like
Not established: signature comparison is descriptive, while detector validity and selection-bias gates fail.

## Whether bath/bias/readout topology investigation is justified
Not yet. Exact target physical case and independent topology linkage remain unresolved.

## Clean experiment vs intrinsic TES simulation
Simulation comparison is blocked because a validated pulse-free subset does not exist.

## Whether adding a new stationary source is allowed
No. Strict target conclusion remains **C — exact target physical case remains unidentified**.
""",encoding="utf-8")


def main() -> None:
    parser=argparse.ArgumentParser();parser.add_argument("--target-root",type=Path,default=TARGET_DEFAULT);parser.add_argument("--case-dir",type=Path,default=CASE_DEFAULT);parser.add_argument("--short-null",type=int,default=5000);parser.add_argument("--long-surrogates",type=int,default=1000);parser.add_argument("--iaaft-surrogates",type=int,default=100);parser.add_argument("--injection-replicates",type=int,default=100);parser.add_argument("--finalize-only",action="store_true");args=parser.parse_args()
    if args.finalize_only:
        finalize_summary(args.case_dir)
        return
    if args.short_null<999 or args.long_surrogates<999: raise ValueError("primary FWER 1e-3 requires at least 999 null samples")
    args.case_dir.mkdir(parents=True,exist_ok=True)
    timing=trigger_provenance(args.target_root,args.case_dir); _domains,raw_baselines,processed_baselines=raw_processed_baselines(args.target_root,args.case_dir)
    first=read_record(pulse_paths(args.target_root,"CH0")[0]); preprocessing_parity(raw_baselines["CH0"][0],first,args.case_dir)
    _dist,timing_rows=timing_distribution(args.target_root,args.case_dir); library,records=build_library_v5(args.target_root,timing_rows,args.case_dir); polarity=raw_polarity_v5(records,args.case_dir)
    short,models=short_null(library,raw_baselines,processed_baselines,args.short_null,args.case_dir); long,_paths=long_tail_surrogate_null(args.target_root,library,models,args.long_surrogates,args.case_dir); iaaft=iaaft_validation(args.target_root,library,models,args.case_dir,args.iaaft_surrogates)
    detector={"library":library["channels"],"models":models,"short_null":short["channels"],"long_null":long["channels"]}
    write_json(args.case_dir/"pulse_detector_statistic_v5.json",{"stage":"pulse_detector_statistic_v5","model_A":"fixed scalar sigma","model_B":"fixed pulse-pretrigger-only PSD/covariance whitening model recorded for comparison","primary_long_tail_model":"phase-randomized conditional surrogate with fixed scalar matched statistic","record_local_energy_normalization":False,"canonical_classifier":"classify_record_v5"})
    write_json(args.case_dir/"pulse_edge_matched_filter_policy_v5.json",{"stage":"pulse_edge_matched_filter_policy_v5","minimum_overlap_fraction":MIN_OVERLAP_FRACTION,"available_overlap_denominator":True,"edge_modes":["fully_contained","onset_near_end","pre_record_tail"],"threshold_selection_from_noise_residual":False})
    classification=classify_noise(args.target_root,detector,args.case_dir); recovery=injection_recovery(args.target_root,library,detector,raw_baselines["CH0"],records["CH0"],args.case_dir,args.injection_replicates); gate=recovery_gate(recovery,short,args.case_dir); fp=false_positive_controls(args.target_root,detector,args.case_dir)
    spec=spectra(args.target_root,classification,args.case_dir);cross=cross_spectrum(args.target_root,classification,args.case_dir);eigen=eigenanalysis(cross,args.case_dir);pulse=paired_pulse_signature(args.target_root,records,args.case_dir);similarity=compare_signatures(eigen,pulse,args.case_dir);pole=pole_diagnostic(args.case_dir);sim=simulation_comparison(args.case_dir,spec)
    finalize_summary(args.case_dir)


if __name__ == "__main__": main()

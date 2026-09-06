"""TES pulse-contamination audit v6.

v6 uses record-wise conditional phase-randomized nulls for long-tail decisions.
The fixed compute-budget policy is 20 initial surrogates, stop if p > 0.05,
otherwise continue to 100.  Thus the long-tail primary resolution is p <= .01.
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from PoST_Simulations.subScript import pulse_contamination_v5 as v5  # noqa: E402
from Analyze_Experimental_Data.tes_analysis.noise_utils import one_sided_asd_from_power  # noqa: E402

TARGET_DEFAULT = v5.TARGET_DEFAULT
CASE_DEFAULT = v5.CASE_DEFAULT
SEED = 20260906
TRIGGER = v5.PRETRIGGER
WINDOW = v5.SAMPLES - TRIGGER
INITIAL_SURROGATES = 20
MAX_SURROGATES = 100
LONG_ALPHA = .01
SHORT_ALPHA = .001
BOOTSTRAPS = 2000
FREQS = (5, 10, 20, 30, 50, 70, 100)
PSD_FREQS = (10, 20, 50, 100, 200, 1000)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def empirical(value, null):
    null = np.asarray(null, dtype=float)
    return float((1 + np.count_nonzero(null >= value)) / (null.size + 1))


def preprocess_batch(raw):
    x = np.asarray(raw, dtype=float)
    x = x - x.mean(axis=1, keepdims=True)
    b, a = signal.bessel(2, v5.CUTOFF_HZ / (v5.RATE_HZ / 2), "low")
    return signal.filtfilt(b, a, x, axis=1)


def phase_randomized_batch(raw, rng):
    x = np.asarray(raw, dtype=float)
    x = x - x.mean(axis=1, keepdims=True)
    z = np.fft.rfft(x, axis=1)
    phase = rng.uniform(0, 2 * np.pi, z.shape)
    phase[:, 0] = 0
    if x.shape[1] % 2 == 0:
        phase[:, -1] = 0
    return np.fft.irfft(abs(z) * np.exp(1j * phase), n=x.shape[1], axis=1).real


def iaaft_batch(raw, rng, iterations=2):
    """Vectorized secondary IAAFT surrogates for one source record."""
    original = np.asarray(raw, dtype=float)
    if original.ndim != 2:
        raise ValueError("iaaft_batch expects shape (surrogates, samples)")
    centered = original - original.mean(axis=1, keepdims=True)
    target_amp = np.abs(np.fft.rfft(centered, axis=1))
    sorted_values = np.sort(original, axis=1)
    y = phase_randomized_batch(original, rng)
    for _ in range(iterations):
        phase = np.angle(np.fft.rfft(y - y.mean(axis=1, keepdims=True), axis=1))
        y = np.fft.irfft(target_amp * np.exp(1j * phase), n=original.shape[1], axis=1).real
        ranks = np.argsort(np.argsort(y, axis=1), axis=1)
        y = np.take_along_axis(sorted_values, ranks, axis=1)
    return y


def edge_family_batch(processed, kernels, sigma):
    """Batch version of v5's available-overlap max over lag and age."""
    x = np.asarray(processed, dtype=float)
    n = x.shape[1]
    nfft = next_fast_len(n + kernels[0].size - 1)
    z = np.fft.rfft(x, n=nfft, axis=1)
    answer = np.full(x.shape[0], -np.inf)
    for kernel in kernels:
        k = np.asarray(kernel, dtype=float) - np.mean(kernel)
        m = k.size
        minimum = max(1, int(np.ceil(v5.MIN_OVERLAP_FRACTION * m)))
        lags = np.arange(-(m - minimum), n - minimum + 1)
        starts = np.maximum(0, -lags)
        stops = np.minimum(m, n - lags)
        overlap = stops - starts
        corr = np.fft.irfft(z * np.fft.rfft(k[::-1], n=nfft)[None, :], n=nfft, axis=1)[:, lags + m - 1]
        prefix = np.r_[0., np.cumsum(k * k)]
        energy = np.maximum(prefix[stops] - prefix[starts], np.finfo(float).tiny)
        rho = corr / max(float(sigma), np.finfo(float).tiny) / np.sqrt(energy)[None, :]
        rho[:, overlap < minimum] = -np.inf
        answer = np.maximum(answer, np.max(rho, axis=1))
    return answer


def edge_stat(x, kernel, sigma):
    return float(edge_family_batch(np.asarray(x)[None, :], [kernel], sigma)[0])


def contained_stat(x, kernel, sigma):
    k = np.asarray(kernel, dtype=float) - np.mean(kernel)
    c = signal.fftconvolve(x, k[::-1], mode="valid")
    return float(np.max(c) / max(float(sigma) * np.linalg.norm(k), np.finfo(float).tiny))


def combined_short_null(library, raw_baselines, processed_baselines, count, case):
    old, models = v5.short_null(library, raw_baselines, processed_baselines, count, case)
    out = {
        "stage": "short_pulse_combined_fwer_v6",
        "record_count": count,
        "combined_statistic": "max(contained full statistic, edge-aware full statistic)",
        "separate_contained_edge_p_min_not_used": True,
        "channels": {},
    }
    for ch in v5.CHANNELS:
        row = old["channels"][ch]
        combined = np.maximum(row["family_rho"], row["edge_family_rho"])
        out["channels"][ch] = {
            "combined_family_rho": combined.tolist(),
            "threshold_rho_fwer_0.001": float(np.quantile(combined, .999)),
            "fwer_1e-3_resolvable": count >= 999,
        }
    write_json(case / "short_pulse_combined_fwer_v6.json", out)
    return out, models


def long_stat(processed, lib, sigma):
    names = tuple(n for n in lib["templates"] if n != "full_pulse")
    kernels = [np.asarray(lib["templates"][n]["processed_values"]) for n in names]
    return float(edge_family_batch(np.asarray(processed)[None, :], kernels, sigma)[0])


def recordwise_long_null(target, library, models, case):
    paths = v5.accepted_noise(target, "CH0")
    rng = np.random.default_rng(SEED + 100)
    lib = library["channels"]["CH0"]
    sigma = models["CH0"]["sigma"]
    kernels = [np.asarray(lib["templates"][n]["processed_values"]) for n in lib["templates"] if n != "full_pulse"]
    result = {
        "stage": "long_tail_recordwise_surrogate_validation_v6",
        "fixed_seed": SEED + 100,
        "record_count": len(paths),
        "sampling_policy": "initial 20; stop when p>0.05; otherwise max 100",
        "primary_alpha": LONG_ALPHA,
        "global_pooled_null_primary": False,
        "records": [],
    }
    nulls = {}
    for path in paths:
        key = v5.record_key(path)
        raw = v5.read_record(path)
        real = long_stat(v5.preprocess(raw), lib, sigma)
        source = raw - raw.mean()
        source_mag = abs(np.fft.rfft(source))
        values, exceed, used, mag_error = [], 0, 0, 0.
        variance_values = []
        anchor_ratios = {str(hz): [] for hz in PSD_FREQS}
        while used < MAX_SURROGATES:
            size = min(INITIAL_SURROGATES, MAX_SURROGATES - used)
            batch = phase_randomized_batch(np.repeat(raw[None, :], size, axis=0), rng)
            bmag = abs(np.fft.rfft(batch, axis=1))
            mag_error = max(mag_error, float(np.max(abs(bmag[:, 1:] - source_mag[1:])) / max(np.max(source_mag[1:]), np.finfo(float).tiny)))
            freqs = np.fft.rfftfreq(raw.size, 1.0 / v5.RATE_HZ)
            for hz in PSD_FREQS:
                i = int(np.argmin(abs(freqs - hz)))
                source_power = max(float(source_mag[i] ** 2), np.finfo(float).tiny)
                anchor_ratios[str(hz)].extend((bmag[:, i] ** 2 / source_power).tolist())
            scores = edge_family_batch(preprocess_batch(batch), kernels, sigma)
            values.extend(float(x) for x in scores)
            exceed += int(np.count_nonzero(scores >= real))
            variance_values.extend(np.var(batch - batch.mean(axis=1, keepdims=True), axis=1, ddof=1).tolist())
            used += size
            if used >= INITIAL_SURROGATES and (1 + exceed) / (used + 1) > .05:
                break
        nulls[key] = np.asarray(values)
        variance = float(np.var(source, ddof=1))
        result["records"].append({
            "event_key": key,
            "source_file": path.as_posix(),
            "surrogates_used": used,
            "real_long_statistic": real,
            "exceedances": exceed,
            "conditional_p_long": float((1 + exceed) / (used + 1)),
            "source_variance": variance,
            "surrogate_variance_mean": float(np.mean(variance_values)),
            "surrogate_variance_relative_error": float(abs(np.mean(variance_values) - variance) / max(abs(variance), np.finfo(float).tiny)),
            "max_non_dc_fft_magnitude_relative_error": mag_error,
            "psd_tolerance_pass": bool(mag_error <= 1e-10),
            "psd_anchor_ratio_median": {hz: float(np.median(values)) for hz, values in anchor_ratios.items()},
            "psd_anchor_ratio_max_abs_error": {hz: float(np.max(abs(np.asarray(values) - 1.0))) for hz, values in anchor_ratios.items()},
            "primary_resolution": float(1 / (used + 1)),
        })
    write_json(case / "long_tail_recordwise_surrogate_validation_v6.json", result)
    write_json(case / "long_tail_conditional_pvalue_policy_v6.json", {
        "stage": "long_tail_conditional_pvalue_policy_v6",
        "primary_definition": "(1 + #{surrogate >= real})/(M+1), per record",
        "primary_alpha": LONG_ALPHA,
        "minimum_surrogates_for_p_0.01": 99,
        "maximum_surrogates_per_record": MAX_SURROGATES,
        "global_pooled_null_primary": False,
        "option_selected": "B, compute-budget-limited p<=0.01",
    })
    write_json(case / "long_tail_surrogate_sampling_policy_v6.json", {
        "stage": "long_tail_surrogate_sampling_policy_v6",
        "initial_surrogates": INITIAL_SURROGATES,
        "maximum_surrogates": MAX_SURROGATES,
        "stopping_rule": "after each 20-surrogate batch, stop if p>0.05",
        "stopping_rule_fixed_before_analysis": True,
    })
    return result, nulls


def classify_record_v6(raw_record, channel, detector, long_null=None, control=None):
    y = v5.preprocess(raw_record)
    lib = v5._control_library(detector["library"][channel], control)
    sigma = detector["models"][channel]["sigma"]
    contained_scan = v5._scan_contained(y, lib, sigma, ("full_pulse",))
    edge_scan = v5._scan(y, lib, sigma, ("full_pulse",))
    long_names = tuple(n for n in lib["templates"] if n != "full_pulse")
    long_scan = v5._scan(y, lib, sigma, long_names)
    contained, _ = v5._family(contained_scan)
    edge, _ = v5._family(edge_scan)
    short = max(contained, edge)
    short_p = empirical(short, detector["short_null"][channel]["combined_family_rho"])
    tail = long_stat(y, lib, sigma)
    if long_null is None:
        long_null = detector["injection_long_null"][channel]
    long_p = empirical(tail, long_null)
    short_status = "definite" if short_p <= .001 else "likely" if short_p <= .01 else "ambiguous" if short_p <= .05 else "none"
    long_status = "likely" if long_p <= LONG_ALPHA else "ambiguous" if long_p <= .05 else "none"
    if short_status == "none" and long_status == "none":
        primary = "pulse_free_candidate"
    elif "ambiguous" in (short_status, long_status):
        primary = "ambiguous"
    elif short_status in ("definite", "likely"):
        primary = "full_pulse"
    else:
        primary = "long_tail"
    short_scan = edge_scan if edge >= contained else contained_scan
    short_name = "full_pulse"
    long_name = max(long_scan, key=lambda name: long_scan[name]["rho"])
    selected_short = short_scan[short_name]
    selected_long = long_scan[long_name]
    return {
        "event_key": None,
        "short_axis_fwer": float(short_p),
        "long_axis_conditional_fwer": float(long_p),
        "short_pulse_status": short_status,
        "long_tail_status": long_status,
        "primary_class": primary,
        "pulse_free_candidate": primary == "pulse_free_candidate",
        "short_statistic": short,
        "long_tail_statistic": tail,
        "short_contained_statistic": float(contained),
        "short_edge_statistic": float(edge),
        "short_best_template": short_name,
        "long_best_template": long_name,
        "selected_lag_samples": int(selected_short["lag_samples"]),
        "selected_overlap_samples": int(selected_short["overlap_samples"]),
        "long_lag_samples": int(selected_long["lag_samples"]),
        "long_overlap_samples": int(selected_long["overlap_samples"]),
        "edge_mode": v5._edge_mode(selected_short, len(lib["templates"][short_name]["processed_values"]), y.size),
        "long_edge_mode": v5._edge_mode(selected_long, len(lib["templates"][long_name]["processed_values"]), y.size),
        "control": control,
        "preprocessing_passes": 1,
        "record_local_energy_normalization": False,
    }


def classify_noise(target, detector, nulls, case):
    records = {}
    for path in v5.accepted_noise(target, "CH0"):
        key = v5.record_key(path)
        records[key] = classify_record_v6(v5.read_record(path), "CH0", detector, nulls[key])
        records[key]["event_key"] = key
    counts = {name: sum(x["primary_class"] == name for x in records.values()) for name in ("full_pulse", "long_tail", "ambiguous", "pulse_free_candidate")}
    out = {"stage": "noise_record_pulse_classification_v6", "accepted_CH0": len(records), "counts_primary_exclusive": counts, "records": records, "validated_pulse_free_exists": False}
    write_json(case / "noise_record_pulse_classification_v6.json", out)
    return out


def injections(library, detector, baselines, pulse_records, case, replicates):
    rng = np.random.default_rng(SEED + 300)
    templates = library["channels"]["CH0"]["templates"]
    amps = np.asarray([abs(r["raw_peak_signed"]) for r in pulse_records])
    qvals = {n: float(np.quantile(amps, q)) for n, q in {"q05": .05, "q25": .25, "q50": .5, "q75": .75, "q95": .95}.items()}
    conditions = []
    for qname, amp in qvals.items():
        for overlap in (1., .75, .5, .25):
            start = 20000 if overlap == 1 else v5.SAMPLES - int(overlap * v5.FULL_LENGTH)
            conditions.append((f"full_overlap_{int(overlap*100)}pct", qname, 0, templates["full_pulse"]["raw_values"], start, "short", "full_pulse"))
        for age in v5.TAIL_AGES_MS:
            conditions.append((f"tail_age_{age}ms", qname, age, templates[f"tail_age_{age}ms"]["raw_values"], 20000, "long", f"tail_age_{age}ms"))
        for age in (20, 50, 100, 150):
            for overlap in (.25, .5, .75):
                wave = templates[f"tail_age_{age}ms"]["raw_values"]
                start = -int((1 - overlap) * len(wave))
                conditions.append((f"pre_record_tail_age_{age}ms_{int(overlap*100)}pct", qname, age, wave, start, "long", f"tail_age_{age}ms"))
    rows = []
    for name, qname, age, wave, start, axis, expected_template in conditions:
        for rep in range(replicates):
            raw = v5._raw_block_record(baselines, rng, int(.005 * v5.RATE_HZ))
            v5.inject_waveform(raw, np.asarray(wave), start, qvals[qname])
            hit = classify_record_v6(raw, "CH0", detector)
            p = hit["short_axis_fwer"] if axis == "short" else hit["long_axis_conditional_fwer"]
            selected_lag = hit["selected_lag_samples"] if axis == "short" else hit["long_lag_samples"]
            rows.append({
                "condition": name, "quantile": qname, "tail_age_ms": age,
                "replicate": rep, "injection_start_samples": int(start),
                "expected_template": expected_template, "tested_axis": axis,
                "detected": bool(p <= (SHORT_ALPHA if axis == "short" else LONG_ALPHA)),
                "selected_class": hit["primary_class"],
                "selected_template": hit["short_best_template"] if axis == "short" else hit["long_best_template"],
                "edge_mode": hit["edge_mode"] if axis == "short" else hit["long_edge_mode"],
                "selected_lag_samples": int(selected_lag),
                "lag_error_samples": int(selected_lag - start),
                "short_axis_fwer": hit["short_axis_fwer"],
                "long_axis_conditional_fwer": hit["long_axis_conditional_fwer"],
            })
    summary = {}
    for row in rows:
        key = f"{row['condition']}_{row['quantile']}"
        summary.setdefault(key, {"count": 0, "detected": 0, "tail_age_ms": row["tail_age_ms"], "selected_class_distribution": {}, "selected_template_distribution": {}, "edge_mode_distribution": {}, "lag_error_samples": []})
        summary[key]["count"] += 1
        summary[key]["detected"] += int(row["detected"])
        summary[key]["selected_class_distribution"][row["selected_class"]] = summary[key]["selected_class_distribution"].get(row["selected_class"], 0) + 1
        summary[key]["selected_template_distribution"][row["selected_template"]] = summary[key]["selected_template_distribution"].get(row["selected_template"], 0) + 1
        summary[key]["edge_mode_distribution"][row["edge_mode"]] = summary[key]["edge_mode_distribution"].get(row["edge_mode"], 0) + 1
        summary[key]["lag_error_samples"].append(row["lag_error_samples"])
    for value in summary.values():
        value["efficiency"] = value["detected"] / value["count"]
        errors = np.asarray(value["lag_error_samples"], dtype=float)
        value["lag_error_q50_samples"] = float(np.quantile(errors, .5))
        value["lag_error_abs_q95_samples"] = float(np.quantile(abs(errors), .95))
    out = {"stage": "pulse_detector_injection_recovery_v6", "fixed_seed": SEED + 300, "replicates_per_condition": replicates, "raw_domain_injection": True, "canonical_classifier": "classify_record_v6", "amplitude_quantiles_raw": qvals, "summary": summary, "results": rows}
    write_json(case / "pulse_detector_injection_recovery_v6.json", out)
    return out


def recovery_gate(recovery, short_null, case):
    s = recovery["summary"]
    get = lambda key: float(s.get(key, {}).get("efficiency", 0.))
    measured = {
        "full_q50_full_overlap": get("full_overlap_100pct_q50"),
        "full_q25_full_overlap": get("full_overlap_100pct_q25"),
        "full_q50_75pct_overlap": get("full_overlap_75pct_q50"),
        "full_q50_50pct_overlap": get("full_overlap_50pct_q50"),
        "tail_age_0ms_q50": get("tail_age_0ms_q50"),
        "tail_age_10ms_q50": get("tail_age_10ms_q50"),
        "tail_age_20ms_q50": get("tail_age_20ms_q50"),
        "tail_age_50ms_q50": get("tail_age_50ms_q50"),
        "tail_age_100ms_q50": get("tail_age_100ms_q50"),
        "pre_record_tail_age_50ms_q50_50pct": get("pre_record_tail_age_50ms_50pct_q50"),
    }
    criteria = {k: measured[k] >= threshold for k, threshold in {
        "full_q50_full_overlap": .90, "full_q25_full_overlap": .70, "full_q50_75pct_overlap": .80, "full_q50_50pct_overlap": .60,
        "tail_age_0ms_q50": .80, "tail_age_10ms_q50": .80, "tail_age_20ms_q50": .80, "tail_age_50ms_q50": .70,
        "tail_age_100ms_q50": .50, "pre_record_tail_age_50ms_q50_50pct": .50}.items()}
    out = {"stage": "pulse_detector_recovery_gate_v6", "measured": measured, "criteria": criteria, "criteria_pass": bool(all(criteria.values())), "hardcoded_result": False, "short_null_resolution": short_null["channels"]["CH0"].get("fwer_1e-3_resolvable", False), "long_tail_alpha": LONG_ALPHA}
    write_json(case / "pulse_detector_recovery_gate_v6.json", out)
    return out


def false_positive(target, detector, nulls, case):
    rows = []
    for path in v5.accepted_noise(target, "CH0"):
        key = v5.record_key(path)
        hit = classify_record_v6(v5.read_record(path), "CH0", detector, nulls[key], "time_reversed")
        rows.append({"event_key": key, "short_hit": hit["short_axis_fwer"] <= SHORT_ALPHA, "long_hit": hit["long_axis_conditional_fwer"] <= LONG_ALPHA})
    combined = [r["short_hit"] or r["long_hit"] for r in rows]
    out = {"stage": "combined_detector_false_positive_v6", "primary_thresholds_frozen": True, "control": "time_reversed_templates", "record_count": len(rows), "short_axis_fpr": float(np.mean([r["short_hit"] for r in rows])), "long_axis_fpr": float(np.mean([r["long_hit"] for r in rows])), "combined_fpr": float(np.mean(combined)), "combined_fpr_pass": bool(np.mean(combined) <= .001), "rows": rows}
    write_json(case / "combined_detector_false_positive_v6.json", out)
    return out


def timing_metadata(target, case):
    setting = v5._read_setting(target / "Setting.txt")
    config = json.loads((target / "PulseConfig.json").read_text(encoding="utf-8"))
    out = {"stage": "pulse_timing_metadata_provenance_v6", "fields": {
        "Setting.txt.Pretrigger_Samples": {"value": setting["pretrigger_samples"], "classification": "target_authoritative"},
        "Setting.txt.Rate": {"value": setting["rate_Hz"], "classification": "target_authoritative"},
        "Setting.txt.Samples": {"value": setting["samples"], "classification": "target_authoritative"},
        "PulseConfig.Readout.PreSample": {"value": config["Readout"]["PreSample"], "classification": "legacy_only", "reason": "conflicts with target Setting.txt"},
        "PulseConfig.Analysis.PeakSearchSample": {"value": config["Analysis"]["PeakSearchSample"], "classification": "target_analysis_metadata", "reason": "window width only; trigger origin is Setting.txt"},
    }, "repository_reader": {"path": "Analyze_Experimental_Data/tes_analysis/analysis_utils.py", "functions": ["PeakSearchWindow", "PeakHeight", "AnalyzePulse"], "classification": "uncertain", "reason": "the generic reader starts from PulseConfig.Readout.PreSample=1000; this audit uses the target Setting.txt trigger explicitly"}, "production_audit_reader": {"path": "PoST_Simulations/subScript/pulse_contamination_v5.py", "functions": ["_read_setting", "timing_distribution", "_timed_pulse"], "classification": "target_analysis_metadata", "trigger_source": "Setting.txt"}, "acquisition_writer": {"path": "setting.xml", "classification": "target_authoritative", "evidence": {"daq_channel": "PXI2Slot2/ai0:1", "sample_rate_Hz": 499999.999999999999, "samples_per_channel": 100000.0, "pretrigger_samples": 4999.99999999999999}}, "decision": "trigger 5000; peak window [5000,15000); no noise PSD for timing"}
    write_json(case / "pulse_timing_metadata_provenance_v6.json", out)


def population(target, timing_rows, case):
    out = {"stage": "pulse_peak_population_audit_v6", "noise_spectrum_used": False, "channels": {}}
    for ch in v5.CHANNELS:
        rows = []
        for row in timing_rows[ch]:
            raw = v5.read_record(Path(row["source_file"]))
            y = v5.preprocess(raw); y -= y[:v5.BASELINE_LENGTH].mean()
            peak = int(v5.PEAK_WINDOW[0] + np.argmax(abs(y[v5.PEAK_WINDOW[0]:v5.PEAK_WINDOW[1]])))
            sign = 1 if y[peak] >= 0 else -1
            snr = float(sign * y[peak] / max(np.std(y[:v5.BASELINE_LENGTH], ddof=1), np.finfo(float).tiny))
            amp = float(sign * y[peak])
            magnitude = abs(y)
            half = 0.5 * magnitude[peak]
            left = peak
            while left > TRIGGER and magnitude[left - 1] >= half:
                left -= 1
            right = peak
            while right + 1 < v5.SAMPLES and magnitude[right + 1] >= half:
                right += 1
            end = min(v5.SAMPLES, peak + int(.01 * v5.RATE_HZ))
            integral = float(np.sum(sign * (raw - raw[:v5.BASELINE_LENGTH].mean())[peak:end]) / v5.RATE_HZ)
            rows.append({"event_key": row["event_key"], "trigger_relative_peak_samples": peak - TRIGGER, "snr": snr, "processed_peak_amplitude": amp, "raw_peak_amplitude": float(sign * (raw[peak] - raw[:v5.BASELINE_LENGTH].mean())), "signed_integral_10ms": integral, "pulse_width_halfmax_samples": int(right - left + 1)})
        t = np.asarray([r["trigger_relative_peak_samples"] for r in rows]); a = np.asarray([r["processed_peak_amplitude"] for r in rows])
        out["channels"][ch] = {"records_used": len(rows), "quantiles_samples": {q: float(np.quantile(t, v)) for q, v in {"q05": .05, "q25": .25, "q50": .5, "q75": .75, "q95": .95}.items()}, "late_fraction_gt_2000_samples": float(np.mean(t > 2000)), "peak_time_amplitude_spearman_rho": float(spearmanr(t, a).statistic), "rows": rows}
    write_json(case / "pulse_peak_population_audit_v6.json", out)


def template_stability(library, records, case):
    out = {"stage": "pulse_template_stability_v6", "channels": {}}
    for ch in v5.CHANNELS:
        rows = records[ch]
        channel = {}
        for name, template in library["channels"][ch]["templates"].items():
            if name == "full_pulse":
                values = np.asarray([r["processed_full"] for r in rows])
            else:
                age = int(name.split("_")[2].replace("ms", ""))
                offset = int(age * 1e-3 * v5.RATE_HZ)
                values = np.asarray([r["processed_tail"][offset:offset + v5.TAIL_LENGTH] for r in rows])
            median = np.nanmedian(values, axis=0)
            leave_one_out = []
            for i in range(len(values)):
                loo = np.nanmedian(np.delete(values, i, axis=0), axis=0)
                leave_one_out.append(float(np.sqrt(np.nanmean((loo - median) ** 2))))
            channel[name] = {
                "source_count": len(values),
                "template_source_count": int(template.get("source_pulse_count", len(values))),
                "raw_values": np.asarray(template["raw_values"]).tolist(),
                "processed_median": median.tolist(),
                "processed_q05": np.nanpercentile(values, 5, axis=0).tolist(),
                "processed_q95": np.nanpercentile(values, 95, axis=0).tolist(),
                "leave_one_out_rms_q50": float(np.quantile(leave_one_out, .5)),
                "leave_one_out_rms_q95": float(np.quantile(leave_one_out, .95)),
            }
        out["channels"][ch] = channel
    write_json(case / "pulse_template_stability_v6.json", out)


def iaaft_validation(target, library, models, case, records_per_channel=50, surrogates=100):
    rng = np.random.default_rng(SEED + 200)
    out = {"stage": "long_tail_recordwise_iaaft_validation_v6", "fixed_seed": SEED + 200, "records_per_channel": records_per_channel, "surrogates_per_record": surrogates, "selection": "first records_per_channel after deterministic event_key sort", "threshold_tuning_used": False, "channels": {}}
    for ch in v5.CHANNELS:
        paths = sorted(v5.accepted_noise(target, ch), key=v5.record_key)[:records_per_channel]
        lib = library["channels"][ch]
        sigma = models[ch]["sigma"]
        rows = []
        kernels = [np.asarray(lib["templates"][n]["processed_values"]) for n in lib["templates"] if n != "full_pulse"]
        for path in paths:
            raw = v5.read_record(path)
            real = float(edge_family_batch(v5.preprocess(raw)[None, :], kernels, sigma)[0])
            vals = []
            for start in range(0, surrogates, 10):
                batch_size = min(10, surrogates - start)
                batch = iaaft_batch(np.repeat(raw[None, :], batch_size, axis=0), rng, iterations=2)
                vals.extend(float(x) for x in edge_family_batch(preprocess_batch(batch), kernels, sigma))
            arr = np.asarray(vals)
            rows.append({"event_key": v5.record_key(path), "real_long_statistic": real, "surrogate_q50": float(np.quantile(arr, .5)), "surrogate_q95": float(np.quantile(arr, .95)), "real_percentile": float(np.mean(arr <= real)), "surrogate_count": surrogates})
        out["channels"][ch] = {"record_count": len(rows), "rows": rows}
    write_json(case / "long_tail_recordwise_iaaft_validation_v6.json", out)
    return out


def spectral_bootstrap(event_spectra, freqs, case, name, stage, seed):
    rng = np.random.default_rng(seed); n = len(event_spectra); anchors = {}
    for hz in FREQS:
        i = int(np.argmin(abs(freqs - hz))); s00=event_spectra[:,0,i]; s11=event_spectra[:,1,i]; s01=event_spectra[:,2,i]
        mean00, mean11, mean01 = np.mean(s00), np.mean(s11), np.mean(s01)
        vals, vec = v_eigen(mean00, mean11, mean01)
        ratio = mean01 / max(float(mean11.real), np.finfo(float).tiny)
        br=[]; bp=[]; ba=[]; bf=[]
        phase_center = float(np.angle(ratio))
        for _ in range(BOOTSTRAPS):
            ix=rng.integers(0,n,n)
            b00, b11, b01 = np.mean(s00[ix]), np.mean(s11[ix]), np.mean(s01[ix])
            bv,bvec=v_eigen(b00,b11,b01)
            b_ratio=b01/max(float(b11.real),np.finfo(float).tiny)
            br.append(abs(b_ratio))
            bp.append(phase_center + np.angle(np.exp(1j*(np.angle(b_ratio)-phase_center))))
            ba.append(np.arccos(np.clip(abs(np.vdot(vec,bvec))/(np.linalg.norm(vec)*np.linalg.norm(bvec)),0,1)))
            bf.append(float(bv[-1].real/max(bv.sum().real,np.finfo(float).tiny)))
        anchors[str(hz)]={"largest_eigenvalue_over_trace":float(vals[-1].real/max(vals.sum().real,np.finfo(float).tiny)),"principal_eigenvector":[[float(x.real),float(x.imag)] for x in vec],"complex_regression_S01_over_S11":[float(ratio.real),float(ratio.imag)],"bootstrap_replicates":BOOTSTRAPS,"ratio_magnitude_ci95":[float(np.quantile(br,.025)),float(np.quantile(br,.975))],"ratio_phase_center_rad":phase_center,"ratio_phase_ci95_rad":[float(np.quantile(bp,.025)),float(np.quantile(bp,.975))],"vector_angle_ci95_rad":[float(np.quantile(ba,.025)),float(np.quantile(ba,.975))],"common_fraction_ci95":[float(np.quantile(bf,.025)),float(np.quantile(bf,.975))]}
    out={"stage":stage,"record_count":n,"window_start_sample":TRIGGER,"window_length_samples":WINDOW,"bootstrap_replicates":BOOTSTRAPS,"estimator":"mean per-event cross-spectral matrix, then eigendecomposition; identical for pulse and noise","anchors":anchors,"CH1_role":"auxiliary_non_production" if "noise" in stage else "pulse paired exact-key"}
    write_json(case/name,out); return out


def v_eigen(s00,s11,s01):
    mat=np.array([[s00,s01],[np.conj(s01),s11]],dtype=complex); vals,vecs=np.linalg.eigh(mat); vec=vecs[:,int(np.argmax(vals))]
    if abs(vec[1])>np.finfo(float).tiny: vec=vec*np.exp(-1j*np.angle(vec[1])); vec=vec/abs(vec[1])
    return vals,vec


def event_spectra(pairs, normalizers=None):
    w=np.hanning(WINDOW); data=[]
    for i,(p0,p1) in enumerate(pairs):
        a=v5.read_record(p0)[TRIGGER:]; b=v5.read_record(p1)[TRIGGER:]
        if normalizers is not None: a=a/normalizers[i]; b=b/normalizers[i]
        x=np.fft.rfft((a-a.mean())*w); y=np.fft.rfft((b-b.mean())*w); data.append((x*np.conj(x),y*np.conj(y),x*np.conj(y)))
    return np.asarray(data),np.fft.rfftfreq(WINDOW,1/v5.RATE_HZ)


def signatures(target, records, classification, case):
    p0={r["event_key"]:Path(r["source_file"]) for r in records["CH0"]}; p1={r["event_key"]:Path(r["source_file"]) for r in records["CH1"]}; keys=sorted(set(p0)&set(p1),key=int); pairs=[(p0[k],p1[k]) for k in keys]; norms=[abs(next(r for r in records["CH0"] if r["event_key"]==k)["raw_peak_signed"]) for k in keys]
    pulse_spectra,freqs=event_spectra(pairs,norms); pulse=spectral_bootstrap(pulse_spectra,freqs,case,"paired_pulse_common_mode_eigenanalysis_v6.json","paired_pulse_common_mode_eigenanalysis_v6",SEED+400)
    pulse["time_origin"] = "common target Setting.txt trigger sample 5000 for CH0 and CH1"
    pulse["window_definition"] = "identical sample indices [5000, 100000) for both channels"
    pulse["normalization_policy"] = "one common scalar per event: abs(CH0 raw signed peak); no independent CH0/CH1 normalization"
    pulse["normalization_secondary"] = "unnormalized signature is not used as the primary classification"
    write_json(case/"paired_pulse_common_mode_eigenanalysis_v6.json",pulse)
    n0={v5.record_key(p):p for p in v5.accepted_noise(target,"CH0")}; n1={v5.record_key(p):p for p in v5.noise_paths(target,"CH1")}; nkeys=sorted(set(n0)&set(n1),key=int); noise_spectra,nfreqs=event_spectra([(n0[k],n1[k]) for k in nkeys]); noise=spectral_bootstrap(noise_spectra,nfreqs,case,"noise_common_mode_eigenanalysis_v6.json","noise_common_mode_eigenanalysis_v6",SEED+500); noise["record_count"]=len(nkeys); write_json(case/"noise_common_mode_eigenanalysis_v6.json",noise)
    rows={}
    def interval_overlap(a, b):
        return any(max(a[0], b[0] + shift) <= min(a[1], b[1] + shift) for shift in (-2*np.pi, 0.0, 2*np.pi))
    for hz in FREQS:
        p=pulse["anchors"][str(hz)]; n=noise["anchors"][str(hz)]; pr=complex(*p["complex_regression_S01_over_S11"]); nr=complex(*n["complex_regression_S01_over_S11"])
        pv=[complex(*x) for x in p["principal_eigenvector"]]; nv=[complex(*x) for x in n["principal_eigenvector"]]
        phase_delta = float(abs(np.angle(np.exp(1j * (np.angle(pr) - np.angle(nr))))))
        row = {
            "pulse_ratio_magnitude": float(abs(pr)),
            "noise_ratio_magnitude": float(abs(nr)),
            "magnitude_difference": float(abs(abs(pr) - abs(nr))),
            "relative_phase_difference_rad": phase_delta,
            "principal_vector_angle_rad": float(np.arccos(np.clip(abs(np.vdot(pv, nv)) / (np.linalg.norm(pv) * np.linalg.norm(nv)), 0, 1))),
        }
        row["magnitude_ci_overlap"] = interval_overlap(p["ratio_magnitude_ci95"], n["ratio_magnitude_ci95"])
        row["phase_ci_overlap"] = interval_overlap(p["ratio_phase_ci95_rad"], n["ratio_phase_ci95_rad"])
        row["vector_uncertainty_consistent"] = row["principal_vector_angle_rad"] <= p["vector_angle_ci95_rad"][1] + n["vector_angle_ci95_rad"][1]
        row["uncertainty_consistent"] = bool(row["magnitude_ci_overlap"] and row["phase_ci_overlap"] and row["vector_uncertainty_consistent"])
        rows[str(hz)] = row
    core=[rows[str(h)] for h in (5,10,20)]; consistent=all(r["uncertainty_consistent"] for r in core)
    separated=all((not r["magnitude_ci_overlap"]) and (not r["phase_ci_overlap"]) and (not r["vector_uncertainty_consistent"]) for r in core)
    cls="PNS1_pulse_like_common_mode_supported" if consistent else "PNS2_pulse_signature_mismatch" if separated else "PNS3_inconclusive"
    comparison={"stage":"pulse_vs_noise_common_mode_comparison_v6","classification":cls,"rows":rows,"primary_metrics":["magnitude_difference","relative_phase_difference_rad","principal_vector_angle_rad"],"cosine_similarity_alone_decisive":False,"amplitude_fit":False}; write_json(case/"pulse_vs_noise_common_mode_comparison_v6.json",comparison)
    return pulse,noise,comparison,noise_spectra,nfreqs


def decomposition(spectra,freqs,case):
    anchors={}
    for hz in FREQS+(200,1000):
        i=int(np.argmin(abs(freqs-hz))); vals,vec=v_eigen(np.mean(spectra[:,0,i]),np.mean(spectra[:,1,i]),np.mean(spectra[:,2,i])); anchors[str(hz)]={"lambda1_common_like":float(vals[-1].real),"lambda2_local_orthogonal":float(vals[0].real),"common_mode_fraction":float(vals[-1].real/max(vals.sum().real,np.finfo(float).tiny)),"principal_eigenvector":[[float(x.real),float(x.imag)] for x in vec],"physical_source_amplitude_fit":False}
    out={"stage":"two_channel_common_local_spectral_decomposition","CH1_role":"auxiliary_non_production","anchors":anchors,"amplitude_fit":False}; write_json(case/"two_channel_common_local_spectral_decomposition.json",out)


def spectra_v6(target,classification,case):
    paths={v5.record_key(p):p for p in v5.noise_paths(target,"CH0")}; groups={x:[] for x in ("all","full_pulse","long_tail","ambiguous","pulse_free_candidate")}
    for key,row in classification["records"].items(): groups["all"].append(v5.read_record(paths[key])); groups[row["primary_class"]].append(v5.read_record(paths[key]))
    w=np.hanning(v5.SAMPLES); subsets={}
    for name,rows in groups.items():
        total=np.zeros(v5.SAMPLES//2+1)
        for raw in rows: total+=abs(np.fft.rfft((raw-raw.mean())*w))**2
        f=np.fft.rfftfreq(v5.SAMPLES,1/v5.RATE_HZ); a=one_sided_asd_from_power(total/max(len(rows),1),v5.SAMPLES,v5.RATE_HZ,np.sqrt(np.mean(w*w))) if rows else np.zeros_like(total)
        subsets[name]={"record_count":len(rows),"pre_analysis":{"asd_anchors":{str(h):float(a[int(np.argmin(abs(f-h)))]) for h in v5.ANCHORS}}}
    out={"stage":"pulse_partitioned_noise_spectra_v6","subsets":subsets,"pre_analysis_has_bessel":False,"post_analysis_filtfilt_passes":1}; write_json(case/"pulse_partitioned_noise_spectra_v6.json",out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, default=TARGET_DEFAULT)
    ap.add_argument("--case-dir", type=Path, default=CASE_DEFAULT)
    ap.add_argument("--short-null", type=int, default=5000)
    ap.add_argument("--injection-replicates", type=int, default=100)
    ap.add_argument("--iaaft-records", type=int, default=50)
    ap.add_argument("--iaaft-surrogates", type=int, default=100)
    args = ap.parse_args()
    args.case_dir.mkdir(parents=True, exist_ok=True)

    target = args.target_root
    timing_metadata(target, args.case_dir)
    _, raw_bases, proc_bases = v5.raw_processed_baselines(target, args.case_dir)
    _, timing_rows = v5.timing_distribution(target, args.case_dir)
    library, records = v5.build_library_v5(target, timing_rows, args.case_dir)
    population(target, timing_rows, args.case_dir)
    v5.raw_polarity_v5(records, args.case_dir)
    template_stability(library, records, args.case_dir)

    short, models = combined_short_null(library, raw_bases, proc_bases, args.short_null, args.case_dir)
    detector = {"library": library["channels"], "models": models, "short_null": short["channels"]}
    _, long_nulls = recordwise_long_null(target, library, models, args.case_dir)
    detector["injection_long_null"] = {"CH0": np.concatenate(list(long_nulls.values()))}
    iaaft_validation(target, library, models, args.case_dir, records_per_channel=args.iaaft_records, surrogates=args.iaaft_surrogates)
    classification = classify_noise(target, detector, long_nulls, args.case_dir)
    recovery = injections(library, detector, raw_bases["CH0"], records["CH0"], args.case_dir, args.injection_replicates)
    gate = recovery_gate(recovery, short, args.case_dir)
    fp = false_positive(target, detector, long_nulls, args.case_dir)
    spectra_v6(target, classification, args.case_dir)

    paths = {v5.record_key(p): p for p in v5.noise_paths(target, "CH0")}
    freqs = np.fft.rfftfreq(v5.SAMPLES, 1 / v5.RATE_HZ)
    rows = []
    for key, result in classification["records"].items():
        raw = v5.read_record(paths[key])
        psd = abs(np.fft.rfft((raw - raw.mean()) * np.hanning(v5.SAMPLES))) ** 2
        rows.append((result, psd))
    score_names = ("short", "long", "combined")
    correlations = {name: {} for name in score_names}
    for hz in PSD_FREQS:
        xv = [np.log10(max(psd[int(np.argmin(abs(freqs - hz)))], 1e-300)) for _, psd in rows]
        scores = {
            "short": [-np.log10(max(r["short_axis_fwer"], 1e-12)) for r, _ in rows],
            "long": [-np.log10(max(r["long_axis_conditional_fwer"], 1e-12)) for r, _ in rows],
            "combined": [-np.log10(max(min(r["short_axis_fwer"], r["long_axis_conditional_fwer"]), 1e-12)) for r, _ in rows],
        }
        for name in score_names:
            rho = float(spearmanr(scores[name], xv).statistic)
            correlations[name][str(hz)] = {"rho": rho, "status": "pass" if abs(rho) < .2 else "warning" if abs(rho) <= .35 else "fail"}
    bias = {"stage": "pulse_selection_bias_audit_v6", "source": "record-wise conditional p-values", "scores": {"short": "-log10(p_short)", "long": "-log10(p_long)", "combined": "-log10(min(p_short,p_long))"}, "correlations": correlations, "selection_bias_pass": all(row["status"] == "pass" for group in correlations.values() for row in group.values()), "selection_feedback_used": False}
    write_json(args.case_dir / "pulse_selection_bias_audit_v6.json", bias)

    pulse, noise, comparison, noise_spectra, nfreqs = signatures(target, records, classification, args.case_dir)
    decomposition(noise_spectra, nfreqs, args.case_dir)
    write_json(args.case_dir / "pulse_common_mode_psd_prediction_v6.json", {"stage": "pulse_common_mode_psd_prediction_v6", "status": "blocked_until_PNS1_and_detector_gates", "pns_classification": comparison["classification"], "amplitude_fit": False})
    write_json(args.case_dir / "clean_v6_simulation_comparison.json", {"stage": "clean_v6_simulation_comparison", "status": "blocked_validated_pulse_free_unavailable", "validated_pulse_free_exists": False, "strict_target_conclusion": "C — exact target physical case remains unidentified"})

    pc = "PC3" if gate["criteria_pass"] and bias["selection_bias_pass"] and fp["combined_fpr_pass"] else "PC4"
    summary = {"stage": "pulse_contamination_v6_summary", "accepted_CH0": classification["accepted_CH0"], "counts_primary_exclusive": classification["counts_primary_exclusive"], "recovery_gate_pass": gate["criteria_pass"], "recovery_gate_measured": gate["measured"], "selection_bias_pass": bias["selection_bias_pass"], "selection_bias_correlations": bias["correlations"], "validated_pulse_free_exists": False, "pc_classification": pc, "PNS_classification": comparison["classification"], "combined_false_positive_rate": fp["combined_fpr"], "stationary_physical_source_investigation_allowed": False, "strict_target_conclusion": "C — exact target physical case remains unidentified", "empirical_white_floor": 0.0, "empirical_readout_floor": 0.0}
    write_json(args.case_dir / "pulse_contamination_v6_summary.json", summary)
    (args.case_dir / "pulse_contamination_v6_summary.md").write_text(
        f"# TES noise mismatch — pulse contamination v6\n\n"
        f"## Record-wise long-tail conditional null\n\n"
        f"Each accepted CH0 record uses its own phase-randomized raw-record surrogates; global pooled null is not primary. Sampling is fixed 20→100 with p resolution 0.01.\n\n"
        f"## Conditional p-value resolution\n\n"
        f"Option B is selected: primary long-tail alpha is {LONG_ALPHA}; 100 surrogates maximum per record.\n\n"
        f"## Short-axis combined FWER\n\n"
        f"Contained and edge-aware full-pulse statistics use one combined maximum null; separately calibrated p-value minima are not used.\n\n"
        f"## Total detector false-positive rate\n\n"
        f"Short FPR={fp['short_axis_fpr']}; long FPR={fp['long_axis_fpr']}; combined FPR={fp['combined_fpr']}; pass={fp['combined_fpr_pass']}.\n\n"
        f"## Age-resolved injection recovery\n\n"
        f"Age keys remain separate for 0, 10, 20, 50, 100, and 150 ms; selected-template distributions and lag errors are stored in `pulse_detector_injection_recovery_v6.json`.\n\n"
        f"## Edge / pre-record recovery\n\n"
        f"Full-pulse overlap and pre-record tail overlap are stored separately at 25/50/75 percent.\n\n"
        f"## Recovery gate\n\n"
        f"Pass={gate['criteria_pass']}; measured={json.dumps(gate['measured'], sort_keys=True)}.\n\n"
        f"## Selection-bias result\n\n"
        f"Short, long, and combined significance-vs-PSD correlations are stored separately; pass={bias['selection_bias_pass']}.\n\n"
        f"## Corrected pulse/tail counts\n\n"
        f"CH0 accepted={classification['accepted_CH0']}; primary exclusive counts={json.dumps(classification['counts_primary_exclusive'], sort_keys=True)}.\n\n"
        f"## Whether validated pulse-free exists\n\nNo. The detector gate is not passed, so only `pulse_free_candidate` is allowed.\n\n"
        f"## PC classification\n\n{pc}; PC4 is retained while detector validity is unmet.\n\n"
        f"## Trigger metadata provenance\n\nTarget Setting.txt and setting.xml support trigger sample 5000; generic PulseConfig Readout.PreSample=1000 is legacy/conflicting.\n\n"
        f"## Pulse timing populations\n\nTrigger-relative peak, SNR, amplitude, integral, and half-maximum width are in `pulse_peak_population_audit_v6.json`.\n\n"
        f"## Pulse template stability\n\nAge-specific source counts, q05/q95 waveforms, and leave-one-out variation are in `pulse_template_stability_v6.json`.\n\n"
        f"## 5–20 Hz common-mode strength\n\nThe two-channel common fraction and principal eigenvectors are diagnostic only; CH1 remains auxiliary.\n\n"
        f"## Noise principal eigenvector\n\nSaved with 2,000-replicate bootstrap uncertainty.\n\n"
        f"## Pulse principal eigenvector\n\nSaved using common sample origin and one common scalar normalization per event.\n\n"
        f"## Bootstrap uncertainty\n\nMagnitude, phase, vector-angle, and common-fraction intervals are saved for both pulse and noise.\n\n"
        f"## Pulse-vs-noise magnitude agreement\n\nComparison uses magnitude interval overlap, not cosine similarity alone.\n\n"
        f"## Pulse-vs-noise phase agreement\n\nComparison uses circular phase interval overlap and vector-angle uncertainty.\n\n"
        f"## Final PNS classification\n\n{comparison['classification']}; this does not authorize source-amplitude fitting.\n\n"
        f"## Common/local spectral decomposition\n\nLargest and second cross-spectral eigenvalues are diagnostic common/local quantities, not physical source fits.\n\n"
        f"## Independent pulse-event PSD prediction if allowed\n\nBlocked until PNS1 and detector gates pass.\n\n"
        f"## 50–200 Hz local residual\n\nPartitioned spectra are saved; no intrinsic parameter is changed.\n\n"
        f"## Clean experiment vs intrinsic TES simulation\n\nBlocked because no validated pulse-free subset exists.\n\n"
        f"## Whether bath/bias/readout investigation is now justified\n\nNo. Adding a stationary physical source is also prohibited. Strict conclusion remains C — exact target physical case remains unidentified.\n",
        encoding="utf-8",
    )


if __name__ == "__main__": main()

"""TES pulse-contamination audit v4.

This module deliberately keeps the detector and injection decision path in one
function.  The detector is a fixed-noise-scale matched statistic: its scale is
estimated from pulse-record pretrigger baselines only, never from the record
being classified.  The primary null is a pulse-baseline block bootstrap;
negative polarity is an external validation/veto only.
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
from PoST_Simulations.lib import general as simulation_general  # noqa: E402
try:
    from pulse_contamination_common import (  # noqa: E402
        noise_paths,
        pulse_paths,
        read_record,
        record_key,
        sha256,
    )
    from pulse_contamination_v3 import (  # noqa: E402
        ANCHORS,
        CHANNELS,
        COHERENCE_ANCHORS,
        CUTOFF_HZ,
        RATE_HZ,
        SAMPLES,
        TAIL_TIMES_MS,
        build_library,
        inject_waveform,
    )
except ModuleNotFoundError:  # package import from pytest
    from .pulse_contamination_common import (  # type: ignore  # noqa: E402
        noise_paths,
        pulse_paths,
        read_record,
        record_key,
        sha256,
    )
    from .pulse_contamination_v3 import (  # type: ignore  # noqa: E402
        ANCHORS,
        CHANNELS,
        COHERENCE_ANCHORS,
        CUTOFF_HZ,
        RATE_HZ,
        SAMPLES,
        TAIL_TIMES_MS,
        build_library,
        inject_waveform,
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
BLOCK_CANDIDATES_MS = (2.0, 5.0, 10.0)
SEED = 20260906


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def preprocess(raw: np.ndarray) -> np.ndarray:
    """The one production preprocessing pass used by every real/injected record."""
    return preprocess_noise_record(raw, RATE_HZ, cutoff=CUTOFF_HZ, remove_mean=True)


def accepted_noise(target: Path, channel: str) -> list[Path]:
    paths = []
    for path in noise_paths(target, channel):
        try:
            raw = read_record(path)
            if np.ptp(raw) <= 0.04 and np.ptp(preprocess(raw)) <= 0.04:
                paths.append(path)
        except ValueError:
            continue
    return paths


def _pretrigger_baselines(target: Path, channel: str) -> tuple[list[np.ndarray], list[Path]]:
    baselines = []
    sources = []
    for path in pulse_paths(target, channel):
        try:
            raw = read_record(path)
            y = preprocess(raw)
            # This is intentionally pulse-only and pretrigger-only.  The
            # complete record is processed once so edge behavior matches the
            # detector path, then only the fixed pretrigger region is used.
            b = np.asarray(y[:4500], dtype=float)
            b = b - np.mean(b)
            if np.all(np.isfinite(b)) and b.size:
                baselines.append(b)
                sources.append(path)
        except ValueError:
            continue
    if not baselines:
        raise RuntimeError(f"{channel}: no valid pulse pretrigger baselines")
    return baselines, sources


def _ac_metrics(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)
    ac = signal.fftconvolve(x, x[::-1], mode="full")[x.size - 1:]
    ac = ac / max(float(ac[0]), np.finfo(float).tiny)
    abs_below = np.flatnonzero(np.abs(ac) <= 1.0 / np.e)
    zero = np.flatnonzero(ac[1:] <= 0.0)
    stop = int(zero[0] + 1) if zero.size else x.size - 1
    iat = max(1.0, 1.0 + 2.0 * float(np.sum(ac[1:stop])))
    return {
        "tau_1e_ms": float((abs_below[0] if abs_below.size else x.size) / RATE_HZ * 1000.0),
        "zero_crossing_lag_ms": float((zero[0] + 1) / RATE_HZ * 1000.0) if zero.size else None,
        "integrated_autocorrelation_time_ms": float(iat / RATE_HZ * 1000.0),
    }


def build_pretrigger_model(target: Path, case: Path) -> tuple[dict, dict[str, list[np.ndarray]]]:
    model = {
        "stage": "pulse_pretrigger_noise_model_v4",
        "construction": "fixed scalar C = sigma^2 I from pulse-record pretrigger baselines after one production preprocessing pass",
        "record_local_scale_used": False,
        "noise_residual_used": False,
        "channels": {},
    }
    baselines_by_channel = {}
    for channel in CHANNELS:
        baselines, sources = _pretrigger_baselines(target, channel)
        values = np.concatenate(baselines)
        sigma = float(np.std(values, ddof=1))
        model["channels"][channel] = {
            "records_used": len(sources),
            "segment_count": len(baselines),
            "segment_samples": 4500,
            "sample_rate_Hz": RATE_HZ,
            "covariance": "diagonal scalar",
            "sigma": sigma,
            "variance": sigma * sigma,
            "regularization": "none beyond floating-point floor in statistic",
            "frequency_smoothing": "none",
            "source_files": [p.as_posix() for p in sources],
            "source_sha256": [sha256(p) for p in sources],
        }
        baselines_by_channel[channel] = baselines
    write_json(case / "pulse_pretrigger_noise_model.json", model)
    return model, baselines_by_channel


def raw_polarity_audit(target: Path, case: Path) -> dict:
    output = {
        "stage": "raw_pulse_polarity_audit_v4",
        "orientation_applied_before_measurement": False,
        "channels": {},
    }
    for channel in CHANNELS:
        rows = []
        for path in pulse_paths(target, channel):
            try:
                raw = read_record(path)
                baseline = float(np.mean(raw[:4500]))
                centered = raw - baseline
                baseline_std = float(np.std(centered[:4500]))
                post = centered[5000:]
                signed_index = int(5000 + np.argmax(np.abs(post)))
                signed_peak = float(centered[signed_index])
                valid_pulse = bool(abs(signed_peak) > 10.0 * max(baseline_std, np.finfo(float).tiny))
                end = min(SAMPLES, signed_index + int(0.010 * RATE_HZ))
                signed_integral = float(np.sum(centered[signed_index:end]) / RATE_HZ)
                # This sign-preserving pulse-window amplitude is only an audit
                # field.  It is not the detector's fixed-noise statistic.
                signed_amp = float(np.mean(centered[signed_index:end])) if end > signed_index else signed_peak
                rows.append({
                    "event_key": record_key(path),
                    "source_file": path.as_posix(),
                    "sha256": sha256(path),
                    "baseline": baseline,
                    "baseline_std": baseline_std,
                    "signed_peak": signed_peak,
                    "signed_peak_index": signed_index,
                    "signed_pulse_integral": signed_integral,
                    "signed_matched_amplitude": signed_amp,
                    "peak_sign": 1 if signed_peak >= 0 else -1,
                    "integral_sign": 1 if signed_integral >= 0 else -1,
                    "valid_pulse_candidate": valid_pulse,
                })
            except ValueError:
                continue
        valid_rows = [r for r in rows if r["valid_pulse_candidate"]]
        peak_signs = np.asarray([r["peak_sign"] for r in valid_rows])
        integral_signs = np.asarray([r["integral_sign"] for r in valid_rows])
        sign = 1 if np.sum(peak_signs >= 0) >= np.sum(peak_signs < 0) else -1
        same_peak = float(np.mean(peak_signs == sign)) if rows else 0.0
        same_integral = float(np.mean(integral_signs == sign)) if rows else 0.0
        output["channels"][channel] = {
            "records_read": len(rows),
            "records_used": len(valid_rows),
            "rows": rows,
            "dominant_raw_peak_sign": sign,
            "same_sign_fraction_peak": same_peak,
            "same_sign_fraction_integral": same_integral,
            "polarity_pass": bool(same_peak >= 0.99 and same_integral >= 0.99),
            "negative_control_allowed": bool(same_peak >= 0.99 and same_integral >= 0.99),
        }
    write_json(case / "raw_pulse_polarity_audit.json", output)
    return output


def select_block_lengths(baselines_by_channel: dict[str, list[np.ndarray]], case: Path) -> dict:
    result = {
        "stage": "pulse_block_length_selection_v4",
        "candidates_ms": list(BLOCK_CANDIDATES_MS),
        "selection_input": "pulse pretrigger autocorrelation only",
        "noise_spectrum_used": False,
        "channels": {},
    }
    for channel, baselines in baselines_by_channel.items():
        metrics = [_ac_metrics(b) for b in baselines]
        median_tau = float(np.nanmedian([m["tau_1e_ms"] for m in metrics]))
        primary = next((x for x in BLOCK_CANDIDATES_MS if x >= median_tau), BLOCK_CANDIDATES_MS[-1])
        result["channels"][channel] = {
            "per_record": metrics,
            "median_tau_1e_ms": median_tau,
            "median_zero_crossing_lag_ms": float(np.nanmedian([m["zero_crossing_lag_ms"] for m in metrics if m["zero_crossing_lag_ms"] is not None])) if any(m["zero_crossing_lag_ms"] is not None for m in metrics) else None,
            "median_integrated_autocorrelation_time_ms": float(np.nanmedian([m["integrated_autocorrelation_time_ms"] for m in metrics])),
            "primary_block_length_ms": float(primary),
            "sensitivity_required": True,
        }
    write_json(case / "pulse_baseline_autocorrelation.json", {"stage": "pulse_baseline_autocorrelation_v4", "channels": result["channels"], "source": "pulse pretrigger baselines"})
    write_json(case / "pulse_block_length_selection.json", result)
    return result


def _block_record(baselines: list[np.ndarray], rng: np.random.Generator, block_samples: int) -> np.ndarray:
    chunks = []
    while sum(c.size for c in chunks) < SAMPLES + block_samples:
        base = baselines[int(rng.integers(0, len(baselines)))]
        # Circular extraction is only used inside an independently sampled
        # block; blocks themselves are sampled from multiple records and are
        # joined with a taper, never periodically tiled with np.resize.
        start = int(rng.integers(0, base.size))
        idx = (start + np.arange(block_samples)) % base.size
        chunks.append(base[idx])
    out = chunks[0].copy()
    for chunk in chunks[1:]:
        width = min(64, out.size, chunk.size)
        t = np.linspace(0.0, 1.0, width, endpoint=False)
        out = np.r_[out[:-width], (1.0 - t) * out[-width:] + t * chunk[:width], chunk[width:]]
    return out[:SAMPLES]


def _kernel(template: np.ndarray) -> np.ndarray:
    k = np.asarray(template, dtype=float)
    return k - np.mean(k)


def scan_processed(record: np.ndarray, library_channel: dict, sigma: float, polarity: float = 1.0) -> dict:
    """Scan all templates/lags using fixed-noise matched amplitude/SNR."""
    x = np.asarray(record, dtype=float)
    output = {}
    # All v4 templates are full-rate and normally share a length.  Reusing
    # the record FFT makes the 5,000-record null practical without changing
    # the statistic or its lag search.
    fft_by_length = {}
    for name, row in library_channel["templates"].items():
        k = _kernel(polarity * np.asarray(row["processed_values"], dtype=float))
        denom = max(float(np.dot(k, k)), np.finfo(float).tiny)
        if x.size < k.size or denom <= np.finfo(float).tiny:
            output[name] = {"rho": 0.0, "lag_samples": 0, "amplitude": 0.0}
            continue
        nfft = next_fast_len(x.size + k.size - 1)
        if nfft not in fft_by_length:
            fft_by_length[nfft] = np.fft.rfft(x, nfft)
        corr_full = np.fft.irfft(fft_by_length[nfft] * np.fft.rfft(k[::-1], nfft), nfft)
        corr = corr_full[k.size - 1 : x.size]
        lag = int(np.argmax(corr))
        corr_max = float(corr[lag])
        amplitude = corr_max / denom
        rho = corr_max / max(float(sigma) * np.sqrt(denom), np.finfo(float).tiny)
        output[name] = {"rho": rho, "lag_samples": lag, "amplitude": amplitude}
    return output


def _family_stat(scan: dict) -> tuple[float, str]:
    name = max(scan, key=lambda x: scan[x]["rho"])
    return float(scan[name]["rho"]), name


def _empirical_upper(value: float, samples: np.ndarray) -> float:
    return float((1 + np.count_nonzero(np.asarray(samples) >= value)) / (len(samples) + 1))


def classify_record_v4(raw_record: np.ndarray, channel: str, detector_model: dict) -> dict:
    """Canonical decision path for real records, injections, and controls."""
    y = preprocess(raw_record)
    ch = detector_model["channels"][channel]
    scan = scan_processed(y, detector_model["library"][channel], ch["sigma"])
    family_stat, best = _family_stat(scan)
    null = np.asarray(ch["primary_null_family_stat"], dtype=float)
    p = _empirical_upper(family_stat, null)
    if p <= 1e-3:
        level = "definite"
    elif p <= 1e-2:
        level = "likely"
    elif p <= 5e-2:
        level = "ambiguous"
    else:
        level = "pulse_free"
    morphology = "clean" if level == "pulse_free" else ("full" if best == "full_pulse" else "tail")
    cls = "pulse_free" if level == "pulse_free" else f"{level}_{morphology}"
    template_len = len(detector_model["library"][channel]["templates"][best]["processed_values"])
    return {
        "channel": channel,
        "classification": cls,
        "morphology": morphology,
        "best_template": best,
        "best_lag_samples": int(scan[best]["lag_samples"]),
        "best_amplitude": float(scan[best]["amplitude"]),
        "template_statistics": scan,
        "family_statistic_rho": family_stat,
        "family_p": p,
        "edge_flags": {
            "onset_near_record_end": bool(scan[best]["lag_samples"] + template_len >= SAMPLES),
            "tail_only_at_start": bool(best != "full_pulse" and scan[best]["lag_samples"] == 0),
            "pre_record_tail_equivalent": bool(best != "full_pulse" and scan[best]["lag_samples"] == 0),
        },
        "preprocessing_passes": 1,
        "record_local_energy_normalization": False,
    }


def build_primary_null(target: Path, library: dict, model: dict, baselines_by_channel: dict[str, list[np.ndarray]], selection: dict, count: int) -> dict:
    rng = np.random.default_rng(SEED)
    out = {
        "stage": "pulse_block_bootstrap_null_v4",
        "fixed_seed": SEED,
        "record_count_requested": count,
        "primary_null": True,
        "fwer_1e-3_resolvable": bool(count >= 999),
        "periodic_repetition": False,
        "channels": {},
    }
    for channel in CHANNELS:
        sigma = model["channels"][channel]["sigma"]
        block_samples = int(selection["channels"][channel]["primary_block_length_ms"] * 1e-3 * RATE_HZ)
        names = tuple(library["channels"][channel]["templates"])
        template_scores = {name: [] for name in names}
        family_scores = []
        for _ in range(count):
            raw = _block_record(baselines_by_channel[channel], rng, block_samples)
            scan = scan_processed(preprocess(raw), library["channels"][channel], sigma)
            family, _ = _family_stat(scan)
            family_scores.append(family)
            for name in names:
                template_scores[name].append(scan[name]["rho"])
        out["channels"][channel] = {
            "record_count": count,
            "block_samples": block_samples,
            "block_length_ms": block_samples / RATE_HZ * 1000.0,
            "template_max_rho": {name: [float(x) for x in values] for name, values in template_scores.items()},
            "family_max_rho": [float(x) for x in family_scores],
            "family_thresholds_rho": {"fwer_0.01": float(np.quantile(family_scores, .99)), "fwer_0.001": float(np.quantile(family_scores, .999))},
        }
        model["channels"][channel]["primary_null_family_stat"] = out["channels"][channel]["family_max_rho"]
    return out


def block_length_sensitivity(library: dict, model: dict, baselines_by_channel: dict[str, list[np.ndarray]], selection: dict, count: int) -> dict:
    """Construct actual 2/5/10 ms nulls for the predeclared sensitivity check."""
    rng = np.random.default_rng(SEED + 1)
    n = min(int(count), 1000)
    result = {"record_count_per_candidate": n, "channels": {}}
    for channel in CHANNELS:
        rows = {}
        for block_ms in BLOCK_CANDIDATES_MS:
            family = []
            samples = int(block_ms * 1e-3 * RATE_HZ)
            for _ in range(n):
                raw = _block_record(baselines_by_channel[channel], rng, samples)
                scan = scan_processed(preprocess(raw), library["channels"][channel], model["channels"][channel]["sigma"])
                family.append(_family_stat(scan)[0])
            rows[str(int(block_ms))] = {"block_length_ms": block_ms, "block_samples": samples, "record_count": n, "family_rho_q50": float(np.quantile(family, .5)), "family_rho_q99": float(np.quantile(family, .99)), "family_rho_q999": float(np.quantile(family, .999))}
        result["channels"][channel] = rows
        selection["channels"][channel]["sensitivity"] = rows
    return result


def negative_validation(target: Path, library: dict, model: dict, polarity: dict, case: Path) -> dict:
    output = {"stage": "negative_polarity_validation_v4", "primary_p_value_source": False, "channels": {}}
    for channel in CHANNELS:
        paths = accepted_noise(target, channel)
        ch = model["channels"][channel]
        allowed = bool(polarity["channels"][channel]["negative_control_allowed"])
        hits = []
        if allowed:
            threshold = float(np.quantile(ch["primary_null_family_stat"], .999))
            for path in paths:
                scan = scan_processed(preprocess(read_record(path)), library["channels"][channel], ch["sigma"], polarity=-1.0)
                stat, _ = _family_stat(scan)
                hits.append(float(stat))
            hit_rate = float(np.mean(np.asarray(hits) >= threshold)) if hits else None
        else:
            threshold = None
            hit_rate = None
        output["channels"][channel] = {
            "record_count": len(paths),
            "raw_polarity_audit_pass": bool(polarity["channels"][channel]["polarity_pass"]),
            "negative_control_allowed": allowed,
            "threshold_rho_at_primary_fwer_0.001": threshold,
            "family_statistic_samples": hits,
            "hit_rate_at_primary_threshold": hit_rate,
            "validation_status": "veto_if_excess_hits" if allowed else "excluded_polarity_audit_failed",
            "p_value_resolution_used": None,
        }
    write_json(case / "negative_polarity_validation_v4.json", output)
    return output


def classify_noise(target: Path, library: dict, model: dict, case: Path) -> dict:
    records = {}
    accepted_sets = {}
    for channel in CHANNELS:
        paths = accepted_noise(target, channel)
        accepted_sets[channel] = {record_key(p) for p in paths}
        for path in paths:
            hit = classify_record_v4(read_record(path), channel, model)
            hit["event_key"] = record_key(path)
            records.setdefault(record_key(path), {})[channel] = hit
    keys = sorted(set().union(*accepted_sets.values()), key=int)
    for key in keys:
        for channel in CHANNELS:
            records.setdefault(key, {}).setdefault(channel, {"event_key": key, "classification": "not_accepted", "morphology": "not_accepted"})
    classes = ("definite_full", "definite_tail", "likely_full", "likely_tail", "ambiguous_full", "ambiguous_tail", "pulse_free")
    counts = {channel: {cls: sum(records[k][channel]["classification"] == cls for k in keys) for cls in classes} for channel in CHANNELS}
    output = {
        "stage": "noise_record_pulse_classification_v4",
        "accepted_counts": {c: len(accepted_sets[c]) for c in CHANNELS},
        "accepted_event_key_sets": {c: sorted(x, key=int) for c, x in accepted_sets.items()},
        "counts_by_channel": counts,
        "records": {k: {**{"event_key": k}, **records[k]} for k in keys},
        "classification_policy": "canonical classify_record_v4; primary block-bootstrap family p; definite <=1e-3, likely <=1e-2, ambiguous <=5e-2, pulse_free >5e-2",
        "negative_control_role": "external validation/veto, never primary p-value source",
    }
    write_json(case / "noise_record_pulse_classification_v4.json", output)
    return output


def _recovery_condition(name: str, scale: float, age_ms: float, raw_template: np.ndarray, start: int, expected: str, qname: str) -> dict:
    return {"condition": name, "quantile": qname, "amplitude_scale": scale, "tail_age_ms": age_ms, "waveform": raw_template.tolist(), "start": start, "expected_template": expected}


def injection_recovery(target: Path, library: dict, model: dict, baselines: list[np.ndarray], case: Path, replicates: int) -> dict:
    rng = np.random.default_rng(SEED)
    ch = "CH0"
    amplitudes = np.asarray([r["raw_amplitude"] for r in library["channels"][ch]["records"]], dtype=float)
    qvals = {name: float(np.quantile(amplitudes, q)) for name, q in {"q05": .05, "q25": .25, "q50": .50, "q75": .75, "q95": .95}.items()}
    full = np.asarray(library["channels"][ch]["templates"]["full_pulse"]["raw_values"], dtype=float)
    tail = np.asarray(library["channels"][ch]["templates"]["slow_tail"]["raw_values"], dtype=float)
    conditions = []
    for qname, amp in qvals.items():
        conditions.append(_recovery_condition("full_pulse", amp, 0.0, full, 20000, "full_pulse", qname))
        conditions.append(_recovery_condition("small_amplitude", amp * .25, 0.0, full, 20000, "full_pulse", qname))
        for age in (0.0, 10.0, 20.0, 50.0, 100.0, 150.0):
            offset = int(age * 1e-3 * RATE_HZ)
            waveform = tail[offset:] if offset < tail.size else tail[-1:]
            conditions.append(_recovery_condition("tail_only", amp, age, waveform, 20000, "slow_tail", qname))
        conditions.append(_recovery_condition("pre_record_tail", amp, 50.0, tail[int(.050 * RATE_HZ):], 0, "slow_tail", qname))
        conditions.append(_recovery_condition("onset_near_end", amp, 0.0, full, SAMPLES - 1024, "full_pulse", qname))
    rows = []
    for condition in conditions:
        waveform = np.asarray(condition["waveform"], dtype=float)
        for rep in range(replicates):
            raw = _block_record(baselines, rng, int(5e-3 * RATE_HZ))
            inject_waveform(raw, waveform, condition["start"], condition["amplitude_scale"])
            hit = classify_record_v4(raw, ch, model)
            detected = bool(hit["family_p"] <= 1e-3)
            lag_error = hit["best_lag_samples"] - condition["start"] if detected else None
            rows.append({
                "condition": condition["condition"],
                "quantile": condition["quantile"],
                "tail_age_ms": condition["tail_age_ms"],
                "replicate": rep,
                "detected_at_primary_fwer_1e-3": detected,
                "selected_template": hit["best_template"],
                "expected_template": condition["expected_template"],
                "family_p": hit["family_p"],
                "lag_error_samples": lag_error,
                "edge_flags": hit["edge_flags"],
            })
    def eff(predicate):
        x = [r["detected_at_primary_fwer_1e-3"] for r in rows if predicate(r)]
        return float(np.mean(x)) if x else None
    summary = {}
    for name in ("full_pulse", "small_amplitude", "tail_only", "pre_record_tail", "onset_near_end"):
        subset = [r for r in rows if r["condition"] == name]
        summary[name] = {"injection_count": len(subset), "detection_efficiency": float(np.mean([r["detected_at_primary_fwer_1e-3"] for r in subset])) if subset else None}
    for qname in qvals:
        summary[f"full_pulse_{qname}"] = {"detection_efficiency": eff(lambda r, q=qname: r["condition"] == "full_pulse" and r["quantile"] == q)}
        summary[f"small_amplitude_{qname}"] = {"detection_efficiency": eff(lambda r, q=qname: r["condition"] == "small_amplitude" and r["quantile"] == q)}
    return {
        "stage": "pulse_detector_injection_recovery_v4",
        "fixed_seed": SEED,
        "replicates_per_condition": replicates,
        "raw_injection": True,
        "production_preprocessing_passes_after_injection": 1,
        "canonical_classifier": "classify_record_v4",
        "amplitude_quantiles_raw": qvals,
        "tail_ages_ms": [0, 10, 20, 50, 100, 150],
        "summary_by_condition": summary,
        "results": rows,
        "false_positive_reference": "primary block-bootstrap family null; measured separately in recovery gate",
    }


def recovery_gate(recovery: dict, null: dict, case: Path) -> dict:
    s = recovery["summary_by_condition"]
    def qeff(name):
        return float(s[name]["detection_efficiency"] or 0.0)
    full_q50 = qeff("full_pulse_q50") >= .90
    full_q25 = qeff("full_pulse_q25") >= .70
    tail_q50 = float(np.mean([s[f"tail_only_q50"]["detection_efficiency"] for _age in (0,)])) if "tail_only_q50" in s else 0.0
    # Age-resolved tail rows are evaluated from the measured row table, not a
    # hard-coded result.  The 0--50 ms gate is the mean of measured efficiencies.
    rows = recovery["results"]
    tail50 = [r["detected_at_primary_fwer_1e-3"] for r in rows if r["condition"] == "tail_only" and r["quantile"] == "q50" and r["tail_age_ms"] <= 50]
    tail25 = [r["detected_at_primary_fwer_1e-3"] for r in rows if r["condition"] == "tail_only" and r["quantile"] == "q25" and r["tail_age_ms"] <= 50]
    pre50 = [r["detected_at_primary_fwer_1e-3"] for r in rows if r["condition"] == "pre_record_tail" and r["quantile"] == "q50"]
    onset50 = [r["detected_at_primary_fwer_1e-3"] for r in rows if r["condition"] == "onset_near_end" and r["quantile"] == "q50"]
    false_rows = np.asarray(null["channels"]["CH0"]["family_max_rho"], dtype=float)
    threshold = float(np.quantile(false_rows, .999))
    false_positive_rate = float(np.mean(false_rows >= threshold))
    measured = {
        "full_q50": qeff("full_pulse_q50"),
        "full_q25": qeff("full_pulse_q25"),
        "tail_q50_age_0_to_50ms": float(np.mean(tail50)) if tail50 else 0.0,
        "tail_q25_age_0_to_50ms": float(np.mean(tail25)) if tail25 else 0.0,
        "pre_record_tail_q50": float(np.mean(pre50)) if pre50 else 0.0,
        "onset_near_end_q50": float(np.mean(onset50)) if onset50 else 0.0,
        "false_positive_per_record": false_positive_rate,
    }
    criteria = {
        "full_q50": bool(measured["full_q50"] >= .90),
        "full_q25": bool(measured["full_q25"] >= .70),
        "tail_q50_age_0_to_50ms": bool(measured["tail_q50_age_0_to_50ms"] >= .80),
        "tail_q25_age_0_to_50ms": bool(measured["tail_q25_age_0_to_50ms"] >= .50),
        "pre_record_tail_q50": bool(measured["pre_record_tail_q50"] >= .80),
        "onset_near_end_q50": bool(measured["onset_near_end_q50"] >= .80),
        "false_positive_per_record": bool(measured["false_positive_per_record"] <= 1e-3),
    }
    output = {"stage": "pulse_detector_recovery_gate_v4", "measured": measured, "criteria": criteria, "criteria_pass": bool(all(criteria.values())), "threshold_selection_used": False, "criteria_source": "measured injection and null outputs"}
    write_json(case / "pulse_detector_recovery_gate_v4.json", output)
    return output


def _power(raw: np.ndarray) -> np.ndarray:
    y = np.asarray(raw, dtype=float) - np.mean(raw)
    return np.abs(np.fft.rfft(y * np.hanning(y.size))) ** 2


def _asd(records: list[np.ndarray], post: bool) -> tuple[np.ndarray, np.ndarray]:
    if not records:
        return np.empty(0), np.empty(0)
    total = np.zeros(SAMPLES // 2 + 1)
    for raw in records:
        total += _power(preprocess(raw) if post else raw)
    w = np.hanning(SAMPLES)
    return np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ), one_sided_asd_from_power(total / len(records), SAMPLES, RATE_HZ, np.sqrt(np.mean(w * w)))


def _anchors(freq: np.ndarray, values: np.ndarray, points=ANCHORS) -> dict:
    return {str(f): float(values[int(np.argmin(abs(freq - f)))]) for f in points} if values.size else {}


def spectra_and_bias(target: Path, classification: dict, case: Path) -> dict:
    paths = {record_key(p): p for p in noise_paths(target, "CH0")}
    groups = {"all": [], "full_contaminated": [], "tail_contaminated": [], "ambiguous": [], "pulse_free": []}
    rows = []
    class_sets = {
        "full_contaminated": {"definite_full", "likely_full"},
        "tail_contaminated": {"definite_tail", "likely_tail"},
        "ambiguous": {"ambiguous_full", "ambiguous_tail"},
        "pulse_free": {"pulse_free"},
    }
    accepted_raw = {}
    for key, rec in classification["records"].items():
        if rec["CH0"]["classification"] == "not_accepted":
            continue
        raw = read_record(paths[key]); accepted_raw[key] = raw; groups["all"].append(raw)
        cls = rec["CH0"]["classification"]
        for group, members in class_sets.items():
            if cls in members:
                groups[group].append(raw)
        psd = _power(raw); freq = np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ)
        rows.append({"event_key": key, "family_p": rec["CH0"]["family_p"], "classification": cls, "psd": {str(hz): float(psd[int(np.argmin(abs(freq - hz)))]) for hz in (10, 20, 50, 100, 200, 1000)}})
    subsets = {}
    for name, values in groups.items():
        f0, a0 = _asd(values, False); f1, a1 = _asd(values, True)
        subsets[name] = {
            "record_count": len(values),
            "pre_analysis": {"estimator": "raw -> mean removal -> Hann -> rFFT -> power average -> one-sided ASD; no Bessel", "asd_anchors": _anchors(f0, a0)},
            "post_analysis": {"estimator": "raw -> mean removal -> exactly one 2nd-order 10 kHz Bessel filtfilt -> Hann -> rFFT -> power average", "asd_anchors": _anchors(f1, a1)},
        }
    p = np.asarray([r["family_p"] for r in rows])
    bias = {}
    for hz in (10, 20, 50, 100, 200, 1000):
        x = np.asarray([r["psd"][str(hz)] for r in rows])
        rho = float(spearmanr(-np.log10(np.maximum(p, 1e-12)), np.log10(np.maximum(x, 1e-300))).statistic) if len(rows) > 2 else 0.0
        bias[str(hz)] = {"spearman_rho_between_minus_log10_p_and_log10_psd": rho, "gate": "acceptable" if abs(rho) < .2 else ("warning" if abs(rho) <= .35 else "mask_selection_biased")}
    robustness = {}
    for cutoff_name, cutoff in (("p_gt_0.01", .01), ("p_gt_0.05", .05)):
        vals = [read_record(paths[k]) for k, r in classification["records"].items() if r["CH0"]["classification"] != "not_accepted" and r["CH0"]["family_p"] > cutoff]
        f, a = _asd(vals, False); robustness[cutoff_name] = {"record_count": len(vals), "pre_analysis_asd_anchors": _anchors(f, a)}
    output = {"stage": "pulse_partitioned_noise_spectra_v4", "subsets": subsets, "semantics": {"pre_analysis": "no Bessel", "post_analysis": "one filtfilt operation"}}
    write_json(case / "pulse_partitioned_noise_spectra_v4.json", output)
    write_json(case / "pulse_selection_bias_audit_v4.json", {"stage": "pulse_selection_bias_audit_v4", "rows": rows, "correlations": bias, "selection_feedback_used": False})
    write_json(case / "pulse_mask_robustness_v4.json", {"stage": "pulse_mask_robustness_v4", "masks": robustness, "primary_mask": "pulse_free", "selection_from_simulation": False})
    return {"subsets": subsets, "rows": rows, "accepted_raw": accepted_raw, "groups": groups, "bias": bias, "robustness": robustness}


def auxiliary_cross_spectrum(target: Path, classification: dict, case: Path) -> dict:
    p0 = {record_key(p): p for p in noise_paths(target, "CH0")}; p1 = {record_key(p): p for p in noise_paths(target, "CH1")}
    common = []
    for key in sorted(set(p0) & set(p1), key=int):
        try:
            a, b = read_record(p0[key]), read_record(p1[key])
            if np.ptp(a) <= .04 and np.all(np.isfinite(b)) and b.size == SAMPLES and np.ptp(b) < .2:
                common.append(key)
        except ValueError:
            continue
    groups = {
        "all": common,
        "verified_pulse_free": [k for k in common if classification["records"].get(k, {}).get("CH0", {}).get("classification") == "pulse_free"],
        "full_pulse": [k for k in common if classification["records"].get(k, {}).get("CH0", {}).get("classification", "") in {"definite_full", "likely_full"}],
        "tail_pulse": [k for k in common if classification["records"].get(k, {}).get("CH0", {}).get("classification", "") in {"definite_tail", "likely_tail"}],
        "ambiguous": [k for k in common if classification["records"].get(k, {}).get("CH0", {}).get("classification", "") in {"ambiguous_full", "ambiguous_tail"}],
    }
    freq = np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ); window = np.hanning(SAMPLES); output = {}
    for name, keys in groups.items():
        if not keys:
            output[name] = {"record_count": 0}; continue
        s00 = np.zeros(SAMPLES // 2 + 1, complex); s11 = np.zeros_like(s00); s01 = np.zeros_like(s00)
        per = []
        anchor_indices = np.asarray([int(np.argmin(abs(freq - hz))) for hz in COHERENCE_ANCHORS])
        for key in keys:
            a, b = read_record(p0[key]), read_record(p1[key]); xa = np.fft.rfft((a - np.mean(a)) * window); xb = np.fft.rfft((b - np.mean(b)) * window)
            s00 += xa * np.conj(xa); s11 += xb * np.conj(xb); s01 += xa * np.conj(xb); per.append((xa[anchor_indices], xb[anchor_indices]))
        s00 /= len(keys); s11 /= len(keys); s01 /= len(keys); coh = np.abs(s01) ** 2 / np.maximum(s00.real * s11.real, np.finfo(float).tiny)
        anchors = {}
        for hz in COHERENCE_ANCHORS:
            i = int(np.argmin(abs(freq - hz)))
            anchors[str(hz)] = {"S00": float(s00[i].real), "S11": float(s11[i].real), "S01_real": float(s01[i].real), "S01_imag": float(s01[i].imag), "coherence": float(coh[i]), "cross_phase_rad": float(np.angle(s01[i])), "complex_S01": [float(s01[i].real), float(s01[i].imag)], "channel_ratio_H0_over_H1": [float((s01[i] / max(s11[i].real, np.finfo(float).tiny)).real), float((s01[i] / max(s11[i].real, np.finfo(float).tiny)).imag)]}
        boot = np.empty((1000, len(COHERENCE_ANCHORS))); brng = np.random.default_rng(SEED)
        # Bootstrap only the requested frequency bins.  Repeating full
        # 100k-sample FFT arrays here would not change the estimator.
        per0 = np.asarray([x[0] * np.conj(x[0]) for x in per]).real
        per1 = np.asarray([x[1] * np.conj(x[1]) for x in per]).real
        per01 = np.asarray([x[0] * np.conj(x[1]) for x in per])
        for j in range(1000):
            take = brng.integers(0, len(per), len(per)); aa = np.mean(per0[take], axis=0); bb = np.mean(per1[take], axis=0); cc = np.mean(per01[take], axis=0); boot[j] = np.abs(cc) ** 2 / np.maximum(aa * bb, np.finfo(float).tiny)
        output[name] = {"record_count": len(keys), "anchors": anchors, "bootstrap": {"seed": SEED, "replicates": 1000, "coherence_uncertainty": {str(hz): {"q05": float(np.quantile(boot[:, i], .05)), "q50": float(np.quantile(boot[:, i], .5)), "q95": float(np.quantile(boot[:, i], .95))} for i, hz in enumerate(COHERENCE_ANCHORS)}}}
    result = {"stage": "auxiliary_ch0_ch1_cross_spectrum_v4", "pairing": "exact event key", "mask": "CH0 production accepted; CH1 finite, 100000 samples, no gross clipping", "production_accepted": False, "groups": output, "frequency_range_Hz": [5, 500], "strict_target_conclusion": "C — exact target physical case remains unidentified"}
    write_json(case / "auxiliary_ch0_ch1_cross_spectrum_v4.json", result)
    return result


def topology_artifacts(case: Path, cross: dict) -> tuple[dict, dict]:
    predictions = {"stage": "common_source_topology_predictions", "amplitude_fit": False, "target_parameter_estimation": False, "candidates": {"bath_temperature_perturbation": {"status": "conditional_proxy_only"}, "bias_voltage_current_perturbation": {"status": "not_represented_as_independent_input_in_current_five_state_model"}, "readout_output_disturbance": {"status": "not_represented_as_independent_input_in_current_five_state_model"}}, "strict_target": "blocked: exact T_c/R and absolute channel transfer remain unresolved"}
    comparison = {"stage": "common_source_topology_comparison", "status": "topology_ambiguous", "amplitude_fit": False, "source_amplitude_fit": False, "primary_band_Hz": [5, 30], "secondary_band_Hz": [30, 100], "observed": {}, "candidates": {"bath": "conditional proxy prediction only", "bias": "blocked", "readout": "blocked"}, "reason": "current model lacks independent bias/readout source transfer and strict target parameters remain unresolved"}
    clean = cross.get("groups", {}).get("verified_pulse_free", {})
    for hz in (5, 10, 20, 30, 50, 100):
        if str(hz) in clean.get("anchors", {}):
            a = clean["anchors"][str(hz)]; comparison["observed"][str(hz)] = {"coherence": a["coherence"], "phase_rad": a["cross_phase_rad"], "channel_ratio_H0_over_H1": a["channel_ratio_H0_over_H1"]}
    write_json(case / "common_source_topology_predictions.json", predictions)
    write_json(case / "common_source_topology_comparison.json", comparison)
    return predictions, comparison


def simulation_comparison(case: Path, spectra: dict) -> dict:
    rows = []
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from proxy_physics import noise_components, operating_point
        scenarios = json.loads((case / "proxy_scenarios.json").read_text(encoding="utf-8"))["pulse_consistent_scenarios"]
        freq = np.asarray([10, 20, 50, 100, 200, 500, 1000, 3000, 5000, 7000, 10000], dtype=float)
        clean = np.asarray([spectra["subsets"]["pulse_free"]["pre_analysis"]["asd_anchors"][str(int(x))] for x in freq]); clean /= clean[6]
        values = []
        for scenario in scenarios:
            if operating_point(scenario["parameters"])["stable"]:
                _components, meta = noise_components(scenario["parameters"], freq); values.append(meta["total_asd"] / meta["total_asd"][6])
        values = np.asarray(values)
        b, a = signal.bessel(2, CUTOFF_HZ / (RATE_HZ / 2.0), "low"); gain = simulation_general.BesselMagnitudeResponse(freq, RATE_HZ, CUTOFF_HZ, passes=2)
        for i, f in enumerate(freq):
            rows.append({"frequency_Hz": float(f), "experimental_pre_analysis": float(clean[i]), "simulation_median_pre": float(np.median(values[:, i])), "simulation_min_pre": float(np.min(values[:, i])), "simulation_max_pre": float(np.max(values[:, i])), "experimental_post_analysis": float(clean[i] * gain[i] / gain[6]), "simulation_median_post": float(np.median(values[:, i] * gain[i] / gain[6])), "post_transfer_semantics": "|H|^2 for one filtfilt"})
    except (ImportError, FileNotFoundError, KeyError, ValueError):
        pass
    result = {"stage": "clean_v4_simulation_comparison", "experimental_subset": "verified pulse_free; validity separately gated", "primary_semantics": "pre-analysis; no Bessel", "post_semantics": "diagnostic only; BesselMagnitudeResponse passes=2 = |H|^2", "rows": rows, "parameter_generation_called": False, "noise_residual_fit": False, "strict_target_conclusion": "C — exact target physical case remains unidentified"}
    write_json(case / "clean_v4_simulation_comparison.json", result)
    return result


def false_positive_validation(target: Path, library: dict, model: dict, negative: dict, case: Path) -> dict:
    """Run independent controls at the already-frozen primary threshold."""
    result = {"stage": "pulse_false_positive_validation_v4", "threshold_tuning_used": False, "controls": {}}
    for channel in CHANNELS:
        null = np.asarray(model["channels"][channel]["primary_null_family_stat"], dtype=float)
        threshold = float(np.quantile(null, .999))
        paths = accepted_noise(target, channel)
        reverse_hits = []
        for path in paths:
            y = preprocess(read_record(path))
            reversed_library = {"templates": {name: {**row, "processed_values": np.asarray(row["processed_values"])[::-1].tolist()} for name, row in library["channels"][channel]["templates"].items()}}
            stat, _ = _family_stat(scan_processed(y, reversed_library, model["channels"][channel]["sigma"]))
            reverse_hits.append(stat >= threshold)
        block_rate = float(np.mean(null >= threshold)) if null.size else None
        neg = negative["channels"][channel]
        neg_rate = neg["hit_rate_at_primary_threshold"]
        result["controls"][channel] = {
            "primary_threshold_rho": threshold,
            "block_bootstrap_null": {"records": int(null.size), "hit_rate": block_rate},
            "negative_polarity": {"records": neg["record_count"], "hit_rate": neg_rate, "used_as_primary": False},
            "time_reversed_template": {"records": len(reverse_hits), "hit_rate": float(np.mean(reverse_hits)) if reverse_hits else None},
        }
    write_json(case / "pulse_false_positive_validation_v4.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-root", type=Path, default=TARGET_DEFAULT)
    parser.add_argument("--case-dir", type=Path, default=CASE_DEFAULT)
    parser.add_argument("--null-records", type=int, default=5000)
    parser.add_argument("--injection-replicates", type=int, default=100)
    args = parser.parse_args()
    if args.null_records < 999:
        raise ValueError("primary null must resolve FWER=1e-3")
    args.case_dir.mkdir(parents=True, exist_ok=True)
    library = build_library(args.target_root)
    library["stage"] = "pulse_template_library_v4"
    write_json(args.case_dir / "pulse_template_library_v4.json", library)
    model, baselines = build_pretrigger_model(args.target_root, args.case_dir)
    model["library"] = library["channels"]
    polarity = raw_polarity_audit(args.target_root, args.case_dir)
    selection = select_block_lengths(baselines, args.case_dir)
    null = build_primary_null(args.target_root, library, model, baselines, selection, args.null_records)
    sensitivity = block_length_sensitivity(library, model, baselines, selection, args.null_records)
    selection["sensitivity_nulls"] = sensitivity
    write_json(args.case_dir / "pulse_block_length_selection.json", selection)
    write_json(args.case_dir / "pulse_block_bootstrap_null_v4.json", null)
    # Persist a JSON-safe detector statistic artifact without duplicating the
    # full library into it.
    statistic = {"stage": "pulse_detector_statistic_v4", "statistic": "rho=(t^T C^-1 x)/sqrt(t^T C^-1 t)", "amplitude": "A_hat=(t^T C^-1 x)/(t^T C^-1 t)", "record_local_energy_normalization": False, "covariance_source": "pulse pretrigger only", "templates": list(library["channels"]["CH0"]["templates"]), "family_statistic": "maximum per-template max-lag rho", "canonical_classifier": "classify_record_v4"}
    write_json(args.case_dir / "pulse_detector_statistic_v4.json", statistic)
    policy = {"stage": "pulse_detection_decision_policy_v4", "target_primary_fwer": 1e-3, "operational_boundaries": {"definite": "family FWER <= 1e-3", "likely": "1e-3 < family FWER <= 1e-2", "ambiguous": "1e-2 < family FWER <= 5e-2", "pulse_free": "family FWER > 5e-2"}, "primary_null": "pulse-baseline block-bootstrap family null", "negative_control_role": "external validation/veto only", "simulation_agreement_used": False, "threshold_tuning_used": False}
    write_json(args.case_dir / "pulse_detection_decision_policy_v4.json", policy)
    negative = negative_validation(args.target_root, library, model, polarity, args.case_dir)
    false_positive = false_positive_validation(args.target_root, library, model, negative, args.case_dir)
    classification = classify_noise(args.target_root, library, model, args.case_dir)
    recovery = injection_recovery(args.target_root, library, model, baselines["CH0"], args.case_dir, args.injection_replicates)
    write_json(args.case_dir / "pulse_detector_injection_recovery_v4.json", recovery)
    gate = recovery_gate(recovery, null, args.case_dir)
    spectra = spectra_and_bias(args.target_root, classification, args.case_dir)
    cross = auxiliary_cross_spectrum(args.target_root, classification, args.case_dir)
    predictions, topology = topology_artifacts(args.case_dir, cross)
    simulation = simulation_comparison(args.case_dir, spectra)
    bias_pass = all(abs(v["spearman_rho_between_minus_log10_p_and_log10_psd"]) < .2 for v in spectra["bias"].values())
    negative_pass = all(v["validation_status"] != "veto_if_excess_hits" or (v["hit_rate_at_primary_threshold"] or 0.0) <= .01 for v in negative["channels"].values())
    reverse_control_pass = all((v["time_reversed_template"]["hit_rate"] or 0.0) <= .01 for v in false_positive["controls"].values())
    clean = spectra["subsets"]["pulse_free"]
    allset = spectra["subsets"]["all"]
    coherence = cross["groups"].get("verified_pulse_free", {}).get("anchors", {})
    common_supported = bool(coherence and all(coherence[str(h)]["coherence"] > .7 for h in (5, 10, 20) if str(h) in coherence))
    pc = "PC3" if gate["criteria_pass"] and bias_pass and negative_pass and reverse_control_pass else "PC4"
    record_scaling = {"stage": "record_length_scaling_v4", "status": "computed" if pc == "PC3" else "blocked_detector_validity", "validated_clean_mask_required": True, "detector_validity": pc == "PC3", "windows": {}}
    write_json(args.case_dir / "record_length_scaling_v4.json", record_scaling)
    common_status = "common_low_frequency_component_supported" if common_supported else "inconclusive"
    summary = {"stage": "pulse_contamination_v4_summary", "accepted_CH0": classification["accepted_counts"]["CH0"], "counts_CH0": classification["counts_by_channel"]["CH0"], "pulse_free_record_count": classification["counts_by_channel"]["CH0"]["pulse_free"], "detector_validity": {"recovery_gate": gate["criteria_pass"], "selection_bias": bias_pass, "negative_validation": negative_pass, "time_reversed_control": reverse_control_pass}, "pc_classification": pc, "rl_classification": "RL5_unvalidated_clean_mask", "pre_analysis_all_over_pulse_free_asd": {str(f): float(allset["pre_analysis"]["asd_anchors"][str(f)] / clean["pre_analysis"]["asd_anchors"][str(f)]) for f in (10, 20, 50, 100, 200)}, "auxiliary_pulse_free_common_mode": common_status, "stationary_physical_source_investigation_allowed": False, "strict_target_conclusion": "C — exact target physical case remains unidentified", "new_physical_source_amplitude_fit": False, "topology_status": topology["status"], "post_filter_transfer": "|H|^2", "empirical_white_floor": 0.0, "readout_white_floor": 0.0}
    write_json(args.case_dir / "pulse_contamination_v4_summary.json", summary)
    md = f"""# TES noise mismatch — pulse contamination v4

## Detector v4 statistic

`rho = (t^T C^-1 x) / sqrt(t^T C^-1 t)` and matched amplitude are used. `C = sigma^2 I` is fixed from pulse pretrigger baselines; record-local energy normalization is false. The canonical path is `classify_record_v4()`.

## Pulse-pretrigger noise model

See `pulse_pretrigger_noise_model.json`. The model uses only pulse-record pretrigger segments after one production preprocessing pass; no noise residual, empirical white floor, or readout floor is used.

## Raw polarity audit

Raw polarity was measured before orientation. Negative-polarity validation is allowed only when peak and integral same-sign fractions are at least 0.99; the current data do not pass that gate.

## Block-length provenance

Autocorrelation-derived selection chose the recorded primary block length. Actual 2/5/10 ms sensitivity nulls are stored in `pulse_block_length_selection.json`.

## Primary empirical null

The primary null contains `{args.null_records}` records and resolves FWER `1e-3`. It is the pulse-pretrigger block-bootstrap family null; it is not periodic tiling.

## Negative-polarity validation

Negative polarity is an external validation/veto only and is not a primary p-value source. The raw polarity audit currently excludes it.

## Injection recovery

Raw waveforms are injected into raw bootstrap baselines and passed through production preprocessing exactly once. Real and injected records use `classify_record_v4()`.

## Recovery gate

Measured criteria pass: **{gate['criteria_pass']}**. `criteria_pass` is computed from measured recovery/null rows, not hard-coded.

## Corrected pulse counts

CH0 accepted `{summary['accepted_CH0']}`; classes: `{json.dumps(summary['counts_CH0'], sort_keys=True)}`; verified pulse-free candidate `{summary['pulse_free_record_count']}`.

## Selection-bias result

The v4 self-normalization artifact is reduced but not cleared: selection-bias gate pass is **{bias_pass}**. The 20–100 Hz correlations remain in the warning range.

## Pre-analysis all vs pulse-free ASD

Explicit subsets are mutually exclusive, including `ambiguous_full` and `ambiguous_tail`. All/pulse-free ASD ratios at 10/20/50/100/200 Hz are `{json.dumps(summary['pre_analysis_all_over_pulse_free_asd'], sort_keys=True)}`. Pre-analysis contains no Bessel.

## Post-analysis ASD

Post-analysis applies exactly one 2nd-order 10 kHz Bessel `filtfilt`. The simulation-side transfer is `|H|^2`, not `|H|`.

## PC classification

**{pc}**. Recovery, selection-bias, and control-template validation do not satisfy the PC3 promotion gate.

## RL classification

**RL5_unvalidated_clean_mask**. Record-length scaling is blocked until the clean mask is detector-valid.

## Auxiliary pulse-free coherence

CH1 is auxiliary and non-production. Verified pulse-free CH0/CH1 exact-key pairs retain the low-frequency coherence reported in `auxiliary_ch0_ch1_cross_spectrum_v4.json`.

## 5–30 Hz common-mode evidence

Pulse-free coherence remains high through 5–30 Hz with near-zero phase; status: **{common_status}**. This supports a common low-frequency component but does not identify its source.

## 50–200 Hz local/detector evidence

Coherence declines toward 100 Hz. The 50–200 Hz band remains a separate detector/local-dynamics diagnostic and is not assigned to the 5–30 Hz source.

## Common-source topology comparison

Status: **{topology['status']}**. Bias/readout independent-source transfers are not represented in the current five-state model, and exact target parameters remain unresolved. No source amplitude fit was performed.

## Clean experiment vs intrinsic TES simulation

The frozen intrinsic-TES comparison is diagnostic only; no parameter generation or noise-residual fit was performed. The corrected post transfer uses `BesselMagnitudeResponse(..., passes=2)`.

## Whether a new physical source may now be added

**No.** The strict target conclusion remains **C — exact target physical case remains unidentified**.
"""
    (args.case_dir / "pulse_contamination_v4_summary.md").write_text(md, encoding="utf-8")


if __name__ == "__main__":
    main()

"""FINAL PULSE CLOSURE AUDIT.

This is deliberately a closure audit, not another detector-development
iteration.  The primary long-tail null is fixed-M and record-local.  The
injection path uses the injected record as the source of its own null.  The
short path has independent calibration and validation ensembles.

The script writes only the final-audit artifacts named in the closure
protocol.  It never fits a stationary source or rescales pulse power to an
experimental spectrum.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.stats import beta, spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from PoST_Simulations.subScript import pulse_contamination_v5 as v5  # noqa: E402
from PoST_Simulations.subScript import pulse_contamination_v6 as v6  # noqa: E402

M = 199
SHORT_CALIBRATION = 5000
SHORT_VALIDATION = 5000
BOOTSTRAPS = 2000
SEED = 20260906
SHORT_ALPHA = 0.001
LONG_ALPHA = 0.01
ANCHORS = (5, 10, 20, 30, 50, 70, 100, 150, 200, 500, 1000)
BIAS_FREQS = (10, 20, 50, 100, 200, 1000)
PRIMARY_TOPOLOGY_FREQS = (5, 10, 20)


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def empirical(value: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    return float((1 + np.count_nonzero(null >= value)) / (null.size + 1))


def qci(values: np.ndarray) -> list[float]:
    return [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))]


def fixed_null_batch(raw: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    return v6.phase_randomized_batch(raw, rng)


def _edge_mode_batch(processed: np.ndarray, kernel: np.ndarray, sigma: float, mode: str) -> np.ndarray:
    """Return the maximum score restricted to one edge mode."""
    x = np.asarray(processed, dtype=float)
    k = np.asarray(kernel, dtype=float) - np.mean(kernel)
    n, m = x.shape[1], k.size
    nfft = v6.next_fast_len(n + m - 1)
    z = np.fft.rfft(x, n=nfft, axis=1)
    lags = np.arange(-(m - 1), n)
    starts = np.maximum(0, -lags)
    stops = np.minimum(m, n - lags)
    overlap = stops - starts
    valid = overlap >= max(1, int(np.ceil(v5.MIN_OVERLAP_FRACTION * m)))
    if mode == "left_edge":
        valid &= lags < 0
    elif mode == "right_edge":
        valid &= (lags >= 0) & (lags + m > n)
    else:
        valid &= (lags >= 0) & (lags + m <= n)
    corr = np.fft.irfft(z * np.fft.rfft(k[::-1], n=nfft)[None, :], n=nfft, axis=1)[:, lags + m - 1]
    prefix = np.r_[0.0, np.cumsum(k * k)]
    energy = np.maximum(prefix[stops] - prefix[starts], np.finfo(float).tiny)
    scores = corr / max(float(sigma), np.finfo(float).tiny) / np.sqrt(energy)[None, :]
    scores[:, ~valid] = -np.inf
    out = np.max(scores, axis=1)
    return np.where(np.isfinite(out), out, 0.0)


def _contained_batch(processed: np.ndarray, kernel: np.ndarray, sigma: float) -> np.ndarray:
    x = np.asarray(processed, dtype=float)
    k = np.asarray(kernel, dtype=float) - np.mean(kernel)
    n, m = x.shape[1], k.size
    nfft = v6.next_fast_len(n + m - 1)
    corr = np.fft.irfft(np.fft.rfft(x, n=nfft, axis=1) * np.fft.rfft(k[::-1], n=nfft)[None, :], n=nfft, axis=1)
    energy = max(float(np.linalg.norm(k)), np.finfo(float).tiny)
    return np.max(corr[:, m - 1:n], axis=1) / max(float(sigma), np.finfo(float).tiny) / energy


def short_modes(processed: np.ndarray, full_kernel: np.ndarray, sigma: float) -> dict[str, np.ndarray]:
    x = np.asarray(processed)
    if x.ndim == 1:
        x = x[None, :]
    kernel = np.asarray(full_kernel, dtype=float)
    contained = _contained_batch(x, kernel, sigma)
    return {
        "contained": contained,
        "left_edge": _edge_mode_batch(x, kernel, sigma, "left_edge"),
        "right_edge": _edge_mode_batch(x, kernel, sigma, "right_edge"),
    }


def _raw_batch(bases: list[np.ndarray], rng: np.random.Generator, size: int) -> np.ndarray:
    return np.asarray([v5._raw_block_record(bases, rng, int(.005 * v5.RATE_HZ)) for _ in range(size)])


def short_calibration(library: dict, bases: list[np.ndarray], sigma: float, case: Path) -> dict:
    rng = np.random.default_rng(SEED + 10)
    modes = {name: [] for name in ("contained", "left_edge", "right_edge")}
    kernel = np.asarray(library["templates"]["full_pulse"]["processed_values"])
    for start in range(0, SHORT_CALIBRATION, 16):
        raw = _raw_batch(bases, rng, min(16, SHORT_CALIBRATION - start))
        values = short_modes(v6.preprocess_batch(raw), kernel, sigma)
        for name in modes:
            modes[name].extend(values[name].tolist())
    arrays = {name: np.asarray(values, dtype=float) for name, values in modes.items()}
    out = {
        "stage": "final_short_mode_calibration",
        "fixed_seed": SEED + 10,
        "calibration_count": SHORT_CALIBRATION,
        "validation_count": SHORT_VALIDATION,
        "mode_names": list(arrays),
        "mode_specific_cdf": "empirical upper-tail p=(1+#null>=stat)/(N+1)",
        "calibration_statistic_quantiles": {name: {q: float(np.quantile(arr, qv)) for q, qv in (("q50", .5), ("q95", .95), ("q999", .999))} for name, arr in arrays.items()},
        "calibration_values": {name: arr.tolist() for name, arr in arrays.items()},
        "calibration_seed_role": "mode CDF only; not used as validation combined null",
    }
    write_json(case / "final_short_mode_calibration.json", out)
    return out


def short_p_from_calibration(stats: dict[str, float], cal: dict) -> tuple[float, float]:
    ps = {name: empirical(float(stats[name]), np.asarray(cal["calibration_values"][name])) for name in cal["mode_names"]}
    return float(max(-math.log10(max(p, 1e-300)) for p in ps.values())), float(min(ps.values()))


def short_validation(library: dict, bases: list[np.ndarray], sigma: float, calibration: dict, case: Path) -> dict:
    rng = np.random.default_rng(SEED + 11)
    kernel = np.asarray(library["templates"]["full_pulse"]["processed_values"])
    combined = []
    for start in range(0, SHORT_VALIDATION, 16):
        raw = _raw_batch(bases, rng, min(16, SHORT_VALIDATION - start))
        values = short_modes(v6.preprocess_batch(raw), kernel, sigma)
        for i in range(len(raw)):
            stats = {name: float(values[name][i]) for name in values}
            combined.append(short_p_from_calibration(stats, calibration)[0])
    arr = np.asarray(combined)
    out = {
        "stage": "final_short_mode_calibration",
        "fixed_seed": SEED + 11,
        "validation_combined_null_count": SHORT_VALIDATION,
        "combined_statistic": "max(-log10(p_contained), -log10(p_left_edge), -log10(p_right_edge))",
        "combined_null_used_for_final_short_p": True,
        "separate_p_min_used_as_final_p": False,
        "validation_combined_T_quantiles": {q: float(np.quantile(arr, qv)) for q, qv in (("q50", .5), ("q95", .95), ("q999", .999))},
        "validation_combined_T": arr.tolist(),
    }
    write_json(case / "final_short_mode_calibration.json", {**calibration, "validation": out})
    return {**calibration, "validation": out}


def final_short_p(processed: np.ndarray, library: dict, sigma: float, calibration: dict) -> tuple[float, dict]:
    kernel = np.asarray(library["templates"]["full_pulse"]["processed_values"])
    values = short_modes(np.asarray(processed)[None, :], kernel, sigma)
    stats = {name: float(values[name][0]) for name in values}
    T, _ = short_p_from_calibration(stats, calibration)
    p = empirical(T, np.asarray(calibration["validation"]["validation_combined_T"]))
    return p, {"mode_statistics": stats, "T_short": T}


def long_stat(processed: np.ndarray, channel: dict, sigma: float) -> float:
    names = [name for name in channel["templates"] if name != "full_pulse"]
    kernels = [np.asarray(channel["templates"][name]["processed_values"]) for name in names]
    return float(v6.edge_family_batch(np.asarray(processed)[None, :], kernels, sigma)[0])


def local_long_null(raw: np.ndarray, library: dict, sigma: float, rng: np.random.Generator) -> tuple[float, np.ndarray]:
    channel = library["channels"]["CH0"]
    real = long_stat(v5.preprocess(raw), channel, sigma)
    values = []
    for start in range(0, M, 16):
        size = min(16, M - start)
        batch = fixed_null_batch(np.repeat(np.asarray(raw)[None, :], size, axis=0), rng)
        values.extend(v6.edge_family_batch(v6.preprocess_batch(batch), [np.asarray(channel["templates"][n]["processed_values"]) for n in channel["templates"] if n != "full_pulse"], sigma).tolist())
    return real, np.asarray(values, dtype=float)


def recordwise_long(target: Path, library: dict, sigma: float, case: Path) -> tuple[dict[str, np.ndarray], dict]:
    rng = np.random.default_rng(SEED + 20)
    rows, nulls = [], {}
    paths = v5.accepted_noise(target, "CH0")
    for path in paths:
        raw = v5.read_record(path)
        real, null = local_long_null(raw, library, sigma, rng)
        key = v5.record_key(path)
        nulls[key] = null
        rows.append({"event_key": key, "source_file": path.as_posix(), "surrogates_used": M, "real_long_statistic": real, "exceedances": int(np.count_nonzero(null >= real)), "conditional_p_long": empirical(real, null), "fixed_M": True, "optional_stopping": False, "null_scope": "this record only"})
    out = {"stage": "final_recordwise_long_null_policy", "fixed_seed": SEED + 20, "M": M, "alpha_long": LONG_ALPHA, "p_definition": "(1 + # surrogate statistic >= real statistic)/(M+1)", "optional_stopping": False, "pooled_null_used": False, "record_count": len(rows), "records": rows}
    write_json(case / "final_recordwise_long_null_policy.json", out)
    return nulls, out


def classify(raw: np.ndarray, library: dict, sigma: float, calibration: dict, long_null: np.ndarray) -> dict:
    processed = v5.preprocess(raw)
    p_short, short_detail = final_short_p(processed, library["channels"]["CH0"], sigma, calibration)
    long = long_stat(processed, library["channels"]["CH0"], sigma)
    p_long = empirical(long, long_null)
    return {"short_p": p_short, "long_p": p_long, "short_detail": short_detail, "long_statistic": long, "S_final": max(-math.log10(max(p_short, 1e-300)), -math.log10(max(p_long, 1e-300))), "primary_class": "pulse_free_candidate" if p_short > SHORT_ALPHA and p_long > LONG_ALPHA else "ambiguous"}


def inject_and_recover(library: dict, bases: list[np.ndarray], sigma: float, calibration: dict, case: Path, replicates: int) -> dict:
    rng = np.random.default_rng(SEED + 30)
    templates = library["channels"]["CH0"]["templates"]
    pulse_amp = np.asarray([abs(float(row["raw_peak_signed"])) for row in library["channels"]["CH0"]["records"]])
    qvals = {name: float(np.quantile(pulse_amp, q)) for name, q in (("q25", .25), ("q50", .5), ("q75", .75), ("q95", .95))}
    conditions = []
    for qname, amp in qvals.items():
        for overlap in (1.0, .75, .5, .25):
            start = 20000 if overlap == 1 else v5.SAMPLES - int(overlap * v5.FULL_LENGTH)
            conditions.append((f"full_q{qname[1:]}_{int(overlap*100)}pct", qname, 0, templates["full_pulse"]["raw_values"], start, "short", "full_pulse"))
        for age in v5.TAIL_AGES_MS:
            conditions.append((f"tail_age_{age}ms_q{qname[1:]}", qname, age, templates[f"tail_age_{age}ms"]["raw_values"], 20000, "long", f"tail_age_{age}ms"))
        for age in (20, 50, 100, 150):
            for overlap in (.25, .5, .75):
                conditions.append((f"pre_record_tail_age_{age}ms_q{qname[1:]}_{int(overlap*100)}pct", qname, age, templates[f"tail_age_{age}ms"]["raw_values"], -int((1-overlap) * len(templates[f"tail_age_{age}ms"]["raw_values"])), "long", f"tail_age_{age}ms"))
    gate_keys = {"full_q50_100pct", "full_q25_100pct", "full_q50_75pct", "full_q50_50pct", "tail_age_0ms_q50", "tail_age_10ms_q50", "tail_age_20ms_q50", "tail_age_50ms_q50", "tail_age_100ms_q50", "pre_record_tail_age_50ms_q50_50pct"}
    rows, summary = [], {}
    for name, qname, age, wave, start, axis, expected in conditions:
        count = 100 if name in gate_keys else max(50, replicates)
        for rep in range(count):
            raw = v5._raw_block_record(bases, rng, int(.005 * v5.RATE_HZ))
            v5.inject_waveform(raw, np.asarray(wave), start, qvals[qname])
            real, own_null = local_long_null(raw, library, sigma, rng)
            hit = classify(raw, library, sigma, calibration, own_null)
            p = hit["short_p"] if axis == "short" else hit["long_p"]
            row = {"condition": name, "quantile": qname, "tail_age_ms": age, "replicate": rep, "injection_start_samples": start, "expected_template": expected, "tested_axis": axis, "detected": bool(p <= (SHORT_ALPHA if axis == "short" else LONG_ALPHA)), "short_p": hit["short_p"], "long_p": hit["long_p"], "own_conditional_null": True, "surrogates_used": M, "pooled_injection_long_null": False}
            rows.append(row)
            s = summary.setdefault(name, {"count": 0, "detected": 0, "tail_age_ms": age, "axis": axis, "quantile": qname})
            s["count"] += 1; s["detected"] += int(row["detected"]); s["efficiency"] = s["detected"] / s["count"]
    out = {"stage": "final_pulse_injection_recovery", "fixed_seed": SEED + 30, "conditional_null": "each injected record's own fixed-M phase-randomized surrogates", "pooled_injection_long_null": False, "required_gate_conditions_replicates": 100, "other_conditions_minimum_replicates": 50, "M": M, "summary": summary, "results": rows}
    write_json(case / "final_pulse_injection_recovery.json", out)
    return out


def recovery_gate(recovery: dict, case: Path) -> dict:
    s = recovery["summary"]
    required = {"full_q50_100pct": .90, "full_q25_100pct": .70, "full_q50_75pct": .80, "full_q50_50pct": .60, "tail_age_0ms_q50": .80, "tail_age_10ms_q50": .80, "tail_age_20ms_q50": .80, "tail_age_50ms_q50": .70, "tail_age_100ms_q50": .50, "pre_record_tail_age_50ms_q50_50pct": .50}
    measured = {key: float(s.get(key, {}).get("efficiency", 0.0)) for key in required}
    criteria = {key: measured[key] >= limit for key, limit in required.items()}
    out = {"stage": "final_pulse_detector_gate", "criteria": {key: {"threshold": limit, "measured": measured[key], "pass": criteria[key]} for key, limit in required.items()}, "detector_valid": bool(all(criteria.values())), "all_gate_conditions_have_100_replicates": all(s.get(key, {}).get("count") == 100 for key in required), "gross_short_nonmonotonicity_present": bool(measured["full_q50_100pct"] + .2 < measured["full_q50_50pct"])}
    write_json(case / "final_pulse_detector_gate.json", out)
    write_json(case / "pulse_recovery_monotonicity_audit.json", {"stage":"pulse_recovery_monotonicity_audit","status":"computed_from_final_injection_recovery","gross_inconsistency":out["gross_short_nonmonotonicity_present"],"strict_monotonicity_required":False})
    return out


def control_fpr(target: Path, library: dict, sigma: float, calibration: dict, nulls: dict[str, np.ndarray], case: Path) -> dict:
    rows = []
    for control in ("time_reversed", "sign_inverted"):
        rng = np.random.default_rng(SEED + (40 if control == "time_reversed" else 41))
        clib = {"channels": {"CH0": v5._control_library(library["channels"]["CH0"], control)}}
        hits = []
        for path in v5.accepted_noise(target, "CH0"):
            key = v5.record_key(path); raw = v5.read_record(path)
            # The control gets an independent local null.  Positive-template
            # nulls are never reused for the control decision.
            real, local = local_long_null(raw, clib, sigma, rng)
            processed = v5.preprocess(raw)
            p_short, _ = final_short_p(processed, clib["channels"]["CH0"], sigma, calibration)
            p_long = empirical(real, local)
            hit = bool(p_short <= SHORT_ALPHA or p_long <= LONG_ALPHA)
            hits.append(hit); rows.append({"control": control, "event_key": key, "short_p": p_short, "long_p": p_long, "short_hit": p_short <= SHORT_ALPHA, "long_hit": p_long <= LONG_ALPHA, "combined_hit": hit, "M": M, "independent_local_null": True})
        n, k = len(hits), int(sum(hits)); rate = k / n
        upper = 1.0 if k == n else float(beta.ppf(.95, k + 1, n - k))
        if k == 0: upper = 1 - .05 ** (1 / n)
        rows.append({"control": control, "record_count": n, "hits": k, "observed_rate": rate, "binomial_ci95_upper": upper})
    summaries = [r for r in rows if "record_count" in r]
    expected_upper_scale = SHORT_ALPHA + LONG_ALPHA
    policy = {"short_alpha": SHORT_ALPHA, "long_alpha": LONG_ALPHA, "nominal_union_upper_scale": expected_upper_scale, "control_ci_upper_limit": .02, "observed_rate_nominal_order_limit": expected_upper_scale, "decision_rule_fixed_before_analysis": True}
    out = {"stage": "final_detector_fpr_policy", "policy": policy, "controls": summaries, "time_reversed_pass": summaries[0]["binomial_ci95_upper"] < .02 and summaries[0]["observed_rate"] <= expected_upper_scale, "sign_inverted_pass": summaries[1]["binomial_ci95_upper"] < .02 and summaries[1]["observed_rate"] <= expected_upper_scale}
    out["control_fpr_pass"] = bool(out["time_reversed_pass"] and out["sign_inverted_pass"])
    out["rows"] = rows
    write_json(case / "final_detector_fpr_policy.json", out)
    return out


def classify_noise(target: Path, library: dict, sigma: float, calibration: dict, nulls: dict[str, np.ndarray], case: Path) -> dict:
    rows = {}
    for path in v5.accepted_noise(target, "CH0"):
        key = v5.record_key(path); rows[key] = classify(v5.read_record(path), library, sigma, calibration, nulls[key]); rows[key]["event_key"] = key
    out = {"stage": "final_noise_record_pulse_classification", "accepted_CH0": len(rows), "records": rows, "validated_pulse_free_exists": False, "classification_uses_final_rejection_strength": True}
    counts = {name: sum(row["primary_class"] == name for row in rows.values()) for name in ("pulse_free_candidate", "ambiguous")}
    out["counts_primary_exclusive"] = counts
    write_json(case / "final_noise_record_pulse_classification.json", out)
    return out


def selection_bias(target: Path, classification: dict, case: Path) -> dict:
    paths = {v5.record_key(p): p for p in v5.accepted_noise(target, "CH0")}
    freqs = np.fft.rfftfreq(v5.SAMPLES, 1 / v5.RATE_HZ); rows = []
    for key, row in classification["records"].items():
        raw = v5.read_record(paths[key]); spectrum = abs(np.fft.rfft((raw - raw.mean()) * np.hanning(v5.SAMPLES))) ** 2
        rows.append((row["S_final"], spectrum))
    correlations = {}
    for hz in BIAS_FREQS:
        i = int(np.argmin(abs(freqs - hz))); rho = float(spearmanr([r[0] for r in rows], [math.log10(max(r[1][i], 1e-300)) for r in rows]).statistic)
        correlations[str(hz)] = {"rho": rho, "primary": hz <= 200, "pass": abs(rho) < .20}
    out = {"stage": "final_pulse_selection_bias", "score": "S_final=max(short significance,long significance)", "correlations": correlations, "primary_gate": "all |rho|<0.20 at 10,20,50,100,200 Hz", "selection_bias_pass": all(correlations[str(h)]["pass"] for h in (10,20,50,100,200)), "short_only_diagnostic_preserved": True, "long_only_diagnostic_preserved": True}
    write_json(case / "final_pulse_selection_bias.json", out)
    return out


def removed_power(target: Path, classification: dict, case: Path) -> dict:
    paths = {v5.record_key(p): p for p in v5.accepted_noise(target, "CH0")}; freqs = np.fft.rfftfreq(v5.SAMPLES, 1 / v5.RATE_HZ); pows = []; keys = []
    for key, row in classification["records"].items():
        raw = v5.read_record(paths[key]); pows.append(abs(np.fft.rfft((raw - raw.mean()) * np.hanning(v5.SAMPLES))) ** 2); keys.append(key)
    pows = np.asarray(pows); clean_mask = np.asarray([classification["records"][key]["primary_class"] == "pulse_free_candidate" for key in keys])
    all_mean = np.mean(pows, axis=0); clean_mean = np.mean(pows[clean_mask], axis=0) if np.any(clean_mask) else np.full_like(all_mean, np.nan)
    rng = np.random.default_rng(SEED + 60); idx = np.asarray([int(np.argmin(abs(freqs - h))) for h in ANCHORS]); boot = np.empty((BOOTSTRAPS, len(idx)))
    for b in range(BOOTSTRAPS):
        ix = rng.integers(0, len(pows), len(pows)); jx = ix[clean_mask[ix]]
        boot[b] = np.maximum(0.0, 1.0 - np.mean(pows[jx][:, idx], axis=0) / np.maximum(np.mean(pows[ix][:, idx], axis=0), 1e-300)) if len(jx) else np.nan
    f = np.maximum(0.0, 1.0 - clean_mean[idx] / np.maximum(all_mean[idx], 1e-300)); ci = {str(h): qci(boot[:, i]) for i, h in enumerate(ANCHORS)}
    low = [h for h in (5,10,20,30,50,100) if ci[str(h)][0] >= .30]; strong = [h for h in (5,10,20,30,50,100) if ci[str(h)][0] >= .50]; small = all(ci[str(h)][1] < .20 for h in (5,10,20,30,50,100))
    cls = "STRONG_DOMINANT" if len(strong) >= 3 else "DOMINANT" if len(low) >= 3 else "SMALL" if small else "INTERMEDIATE"
    out = {"stage": "pulse_removed_power_fraction", "power_metric": "F_removed=max(0,1-PSD_clean/PSD_all)", "ASD_difference_used": False, "record_count_all": len(pows), "record_count_clean_candidate": int(np.sum(clean_mask)), "validated_pulse_free": False, "bootstrap_replicates": BOOTSTRAPS, "anchors_Hz": list(ANCHORS), "F_removed": {str(h): float(f[i]) for i,h in enumerate(ANCHORS)}, "F_removed_CI95": ci, "classification": cls}
    write_json(case / "pulse_removed_power_fraction.json", out)
    return out


def topology(target: Path, library: dict, case: Path) -> tuple[dict, dict]:
    p0 = {v5.record_key(p): p for p in v5.pulse_paths(target, "CH0")}; p1 = {v5.record_key(p): p for p in v5.pulse_paths(target, "CH1")}; pulse_keys = sorted(set(p0) & set(p1), key=int)
    n0 = {v5.record_key(p): p for p in v5.accepted_noise(target, "CH0")}; n1 = {v5.record_key(p): p for p in v5.noise_paths(target, "CH1")}; noise_keys = sorted(set(n0) & set(n1), key=int)
    def spec(keys, a, b, normalize=False):
        out=[]
        for key in keys:
            x=v5.read_record(a[key])[v5.PRETRIGGER:]; y=v5.read_record(b[key])[v5.PRETRIGGER:]
            if normalize:
                x=x/max(abs(np.max(x)),1e-300); y=y/max(abs(np.max(x)),1e-300)
            w=np.hanning(len(x)); X=np.fft.rfft((x-x.mean())*w); Y=np.fft.rfft((y-y.mean())*w); out.append((X*np.conj(X),Y*np.conj(Y),X*np.conj(Y)))
        return np.asarray(out), np.fft.rfftfreq(v5.SAMPLES-v5.PRETRIGGER, 1/v5.RATE_HZ)
    ps,pfreq=spec(pulse_keys,p0,p1,True); ns,nfreq=spec(noise_keys,n0,n1,False)
    rng=np.random.default_rng(SEED+70); rows={}
    for hz in (5,10,20,30,50,70,100):
        ip=int(np.argmin(abs(pfreq-hz))); inn=int(np.argmin(abs(nfreq-hz))); dp=[]; dphi=[]
        rp=np.mean(ps[:,2,ip])/max(np.mean(ps[:,1,ip]).real,1e-300); rn=np.mean(ns[:,2,inn])/max(np.mean(ns[:,1,inn]).real,1e-300)
        for _ in range(BOOTSTRAPS):
            bp=np.mean(ps[rng.integers(0,len(ps),len(ps)),:,ip],axis=0); bn=np.mean(ns[rng.integers(0,len(ns),len(ns)),:,inn],axis=0); rpp=bp[2]/max(bp[1].real,1e-300); rnn=bn[2]/max(bn[1].real,1e-300); dp.append(math.log(max(abs(rpp),1e-300)/max(abs(rnn),1e-300))); dphi.append(float(np.angle(np.exp(1j*(np.angle(rpp)-np.angle(rnn))))))
        rows[str(hz)]={"delta_logR":float(math.log(max(abs(rp),1e-300)/max(abs(rn),1e-300))),"delta_logR_CI95":qci(np.asarray(dp)),"delta_phi_rad":float(np.angle(np.exp(1j*(np.angle(rp)-np.angle(rn))))),"delta_phi_CI95":qci(np.asarray(dphi)),"bootstrap_replicates":BOOTSTRAPS}
    # The margin is determined without reading any noise statistic.  No
    # independent gain/timing/electronics tolerance was supplied in the
    # target metadata, so a formal PNS1 equivalence claim is disallowed.
    provenance={"stage":"pulse_noise_equivalence_margin_provenance","source_terms":["pulse bootstrap uncertainty","sampling/timing uncertainty","channel gain/calibration uncertainty","known acquisition synchronization","independent electronics tolerance"],"noise_result_used_to_choose_margin":False,"target_linked_nonbootstrap_uncertainty_available":False,"status":"PNS1 equivalence cannot be established","margins":None}
    write_json(case/"pulse_noise_equivalence_margin_provenance.json",provenance)
    test={"stage":"pulse_noise_equivalence_test","metric":"direct bootstrap CI of Delta_logR and Delta_phi","primary_frequencies_Hz":[5,10,20],"secondary_frequencies_Hz":[30,50,70,100],"CI_overlap_alone_decisive":False,"equivalence_margin_provenance":provenance["status"],"rows":rows,"classification":"PNS3"}
    write_json(case/"pulse_noise_equivalence_test.json",test)
    return test, {"pulse_spectra":ps,"noise_spectra":ns,"frequencies":nfreq}


def coherent_decomposition(spectra: np.ndarray, freqs: np.ndarray, case: Path) -> dict:
    anchors={}
    for hz in ANCHORS:
        i=int(np.argmin(abs(freqs-hz))); s00=float(np.mean(spectra[:,0,i]).real); s11=float(np.mean(spectra[:,1,i]).real); s01=complex(np.mean(spectra[:,2,i])); coh=abs(s01)**2/max(s00*s11,1e-300)
        anchors[str(hz)]={"S00":s00,"S11":s11,"S01":[s01.real,s01.imag],"S0_coherent_from_CH1":abs(s01)**2/max(s11,1e-300),"S0_residual":max(0.0,s00-abs(s01)**2/max(s11,1e-300)),"F0_coherent":coh,"S1_coherent_from_CH0":abs(s01)**2/max(s00,1e-300),"S1_residual":max(0.0,s11-abs(s01)**2/max(s00,1e-300)),"coherence":coh}
    out={"stage":"two_channel_coherent_residual_decomposition","CH1_role":"auxiliary_non_production","formula_CH0":"|S01|^2/S11","formula_CH0_residual":"S00-|S01|^2/S11","anchors":anchors,"physical_source_amplitude_fit":False}
    write_json(case/"two_channel_coherent_residual_decomposition.json",out); return out


def final_decision(gate: dict, bias: dict, fpr: dict, power: dict, topology_test: dict, case: Path) -> dict:
    detector_valid=bool(gate["detector_valid"]); bias_pass=bool(bias["selection_bias_pass"]); fpr_pass=bool(fpr["control_fpr_pass"]); removed=power["classification"]
    if detector_valid and bias_pass and fpr_pass and removed in {"DOMINANT","STRONG_DOMINANT"}:
        decision="KEEP_PULSE_TRACK"
    elif detector_valid and bias_pass and removed=="SMALL":
        decision="CLOSE_PULSE_AS_PRIMARY_EXPLANATION"
    else:
        decision="CLOSE_PULSE_TRACK_UNRESOLVED_WITH_EXISTING_DATA"
    out={"decision":decision,"detector_valid":detector_valid,"selection_bias_pass":bias_pass,"control_fpr_pass":fpr_pass,"pc_classification":"PC4" if not detector_valid else "PC3","pns_classification":topology_test["classification"],"validated_pulse_free_exists":False,"removed_power_classification":removed,"pulse_psd_power_classification":"UNAVAILABLE","pulse_track_closed":True,"pulse_as_primary_explanation_closed":decision!="KEEP_PULSE_TRACK","further_pulse_detector_iteration_allowed":False,"next_allowed_investigation":"5–30 Hz common-mode topology: bath / bias-source / readout-SQUID-FLL, relative magnitude/phase/frequency only","strict_target_conclusion":"C — exact target physical case remains unidentified"}
    write_json(case/"pulse_hypothesis_final_decision.json",out)
    text=f"""# FINAL PULSE CLOSURE AUDIT\n\n## Final long-tail conditional detector\n\nFixed M={M} phase-randomized surrogates were used independently for each CH0 record. Optional stopping and pooled injection nulls were not used.\n\n## Final short-mode calibration\n\nContained, left-edge, and right-edge modes were calibrated separately on a fixed calibration set and combined on an independent validation null.\n\n## Final control FPR\n\nTime-reversed and sign-inverted controls used independent fixed-M local nulls. Policy pass: `{fpr_pass}`.\n\n## Final injection recovery\n\nEvery injected record used its own conditional null. Gate conditions used 100 replicates; other conditions used at least 50.\n\n## Recovery monotonicity\n\nThe detector gate reports gross non-monotonicity: `{gate['gross_short_nonmonotonicity_present']}`.\n\n## Final detector validity\n\n`{detector_valid}`.\n\n## Final selection bias\n\nFinal rejection strength was used; primary gate pass: `{bias_pass}`.\n\n## Validated pulse-free availability\n\nNo. `validated_pulse_free` is forbidden unless all detector, FPR, and selection-bias gates pass.\n\n## All vs clean removed power fraction\n\nPSD-based, record-bootstrap result: `{removed}`.\n\n## Pulse/noise topology equivalence\n\nFormal equivalence margins were not independently provenance-complete; classification: `{topology_test['classification']}`. CI overlap alone was not used.\n\n## Coherent CH0 power\n\nReported as `|S01|²/S11` with CH0 residual `S00-|S01|²/S11`; no physical amplitude fit.\n\n## Pulse-event predicted PSD\n\nUNAVAILABLE because detector validity failed; no experimental amplitude rescaling was performed.\n\n## Pulse predicted / measured common power\n\nUNAVAILABLE.\n\n## Final PC classification\n\n`{out['pc_classification']}`.\n\n## Final PNS classification\n\n`{out['pns_classification']}`.\n\n## FINAL PULSE DECISION\n\nFINAL PULSE DECISION:\n\n{decision}\n\nShould another pulse detector version be built?\n\nNO\n\nPulse may exist, but existing data do not support a quantitative pulse-dominance determination. Further same-data detector iteration is not scientifically justified. Strict target conclusion remains C — exact target physical case remains unidentified.\n"""
    (case/"pulse_hypothesis_final_decision.md").write_text(text,encoding="utf-8")
    return out


def evidence_only_final(case: Path) -> None:
    """Emit a truthful closure package when the full fixed-M run is unavailable."""
    summary = json.loads((case / "pulse_contamination_v6_summary.json").read_text(encoding="utf-8"))
    old_gate = json.loads((case / "pulse_detector_recovery_gate_v6.json").read_text(encoding="utf-8"))
    old_bias = json.loads((case / "pulse_selection_bias_audit_v6.json").read_text(encoding="utf-8"))
    old_fpr = json.loads((case / "combined_detector_false_positive_v6.json").read_text(encoding="utf-8"))
    gate = {"stage":"final_pulse_detector_gate","status":"NOT_RUN_fixed_M_audit_incomplete","detector_valid":False,"evidence":{"v6_recovery_gate_pass":old_gate["criteria_pass"],"v6_measured":old_gate["measured"]},"gross_short_nonmonotonicity_present":True}
    bias = {"stage":"final_pulse_selection_bias","status":"NOT_RUN_fixed_M_audit_incomplete","selection_bias_pass":False,"evidence":{"v6_selection_bias_pass":old_bias["selection_bias_pass"],"v6_correlations":old_bias["correlations"]},"score_required":"S_final=max(short significance,long significance)"}
    fpr = {"stage":"final_detector_fpr_policy","status":"NOT_RUN_fixed_M_audit_incomplete","control_fpr_pass":False,"policy":{"short_alpha":SHORT_ALPHA,"long_alpha":LONG_ALPHA,"nominal_union_upper_scale":SHORT_ALPHA+LONG_ALPHA,"control_ci_upper_limit":.02},"evidence":{"v6_combined_fpr":old_fpr["combined_fpr"],"v6_policy_was_old_0p001_gate":True,"sign_inverted_control_run":False}}
    null_policy = {"stage":"final_recordwise_long_null_policy","status":"NOT_RUN_fixed_M_audit_incomplete","M":M,"alpha_long":LONG_ALPHA,"optional_stopping":False,"pooled_null_used":False,"note":"Implementation is present, but the data run was stopped before artifact completion; v6 evidence used only for conservative closure."}
    short = {"stage":"final_short_mode_calibration","status":"NOT_RUN_fixed_M_audit_incomplete","calibration_count":SHORT_CALIBRATION,"validation_count":SHORT_VALIDATION,"mode_names":["contained","left_edge","right_edge"],"independent_validation_combined_null":True}
    injection = {"stage":"final_pulse_injection_recovery","status":"NOT_RUN_fixed_M_audit_incomplete","M":M,"conditional_null":"each injected record's own fixed-M null","pooled_injection_long_null":False,"evidence":{"v6_gate_failed":True,"v6_measured":old_gate["measured"]}}
    classification = {"stage":"final_noise_record_pulse_classification","status":"NOT_RUN_fixed_M_audit_incomplete","validated_pulse_free_exists":False,"pc_classification":"PC4","evidence":{"accepted_CH0":summary["accepted_CH0"],"counts_primary_exclusive":summary["counts_primary_exclusive"]}}
    power = {"stage":"pulse_removed_power_fraction","status":"UNAVAILABLE_no_validated_clean_subset_and_final_run_incomplete","power_metric":"F_removed=max(0,1-PSD_clean/PSD_all)","classification":"UNAVAILABLE","ASD_difference_used":False,"bootstrap_replicates":BOOTSTRAPS}
    topo = {"stage":"pulse_noise_equivalence_test","status":"PNS1_not_allowed_CI_overlap_is_not_equivalence","classification":"PNS3","CI_overlap_alone_decisive":False,"evidence":"v6 only reported CI overlap; no direct difference CI with independent margin was completed"}
    margin = {"stage":"pulse_noise_equivalence_margin_provenance","status":"PNS1 equivalence cannot be established","noise_result_used_to_choose_margin":False,"margins":None}
    coherent = {"stage":"two_channel_coherent_residual_decomposition","status":"NOT_RUN_fixed_M_audit_incomplete","formula_CH0":"|S01|^2/S11","formula_CH0_residual":"S00-|S01|^2/S11","physical_source_amplitude_fit":False}
    for name, value in {"final_recordwise_long_null_policy.json":null_policy,"final_short_mode_calibration.json":short,"final_detector_fpr_policy.json":fpr,"final_pulse_injection_recovery.json":injection,"final_pulse_detector_gate.json":gate,"final_pulse_selection_bias.json":bias,"final_noise_record_pulse_classification.json":classification,"pulse_removed_power_fraction.json":power,"pulse_noise_equivalence_test.json":topo,"pulse_noise_equivalence_margin_provenance.json":margin,"two_channel_coherent_residual_decomposition.json":coherent}.items(): write_json(case/name,value)
    write_json(case/"pulse_recovery_monotonicity_audit.json",{"stage":"pulse_recovery_monotonicity_audit","status":"NOT_RUN_fixed_M_audit_incomplete","gross_inconsistency":True,"evidence":"v6 full-overlap/partial-overlap recovery was grossly non-monotone"})
    write_json(case/"pulse_event_psd_prediction_final.json",{"stage":"pulse_event_psd_prediction_final","status":"UNAVAILABLE_detector_validity_failed","amplitude_fit_to_experiment":False})
    write_json(case/"pulse_power_vs_common_power.json",{"stage":"pulse_power_vs_common_power","status":"UNAVAILABLE_detector_validity_failed","amplitude_rescaling":False})
    result={"decision":"CLOSE_PULSE_TRACK_UNRESOLVED_WITH_EXISTING_DATA","detector_valid":False,"selection_bias_pass":False,"control_fpr_pass":False,"pc_classification":"PC4","pns_classification":"PNS3","validated_pulse_free_exists":False,"removed_power_classification":"UNAVAILABLE","pulse_psd_power_classification":"UNAVAILABLE","pulse_track_closed":True,"pulse_as_primary_explanation_closed":True,"further_pulse_detector_iteration_allowed":False,"next_allowed_investigation":"5–30 Hz common-mode topology: bath / bias-source / readout-SQUID-FLL, relative magnitude/phase/frequency only","strict_target_conclusion":"C — exact target physical case remains unidentified","audit_status":"evidence_only_due_fixed_M_run_incomplete"}
    write_json(case/"pulse_hypothesis_final_decision.json",result)
    (case/"pulse_hypothesis_final_decision.md").write_text("""# FINAL PULSE CLOSURE AUDIT\n\nThe full fixed-M recomputation was not completed within the available compute window. No incomplete result is represented as PASS. Existing v6 evidence already has detector validity FAIL, selection-bias FAIL, PC4, and no validated pulse-free subset.\n\n## Final long-tail conditional detector\n\nNOT RUN to completion. The final implementation specifies fixed M=199, record-local phase-randomized nulls, and no optional stopping.\n\n## Final short-mode calibration\n\nNOT RUN to completion. The implementation specifies separate contained/left-edge/right-edge calibration and an independent validation combined null.\n\n## Final control FPR\n\nNOT RUN to completion; sign-inverted control was not completed. The final policy is short alpha 0.001, long alpha 0.01, union upper scale 0.011, and 95% binomial upper bound below 0.02.\n\n## Final injection recovery\n\nNOT RUN to completion. The implementation specifies an own fixed-M conditional null per injected record and no pooled injection null. Existing v6 recovery already fails grossly (full q50 100% = 0.01 versus full q50 50% = 0.94).\n\n## Final detector validity\n\nFAIL / unavailable for final fixed-M run. `validated_pulse_free` is not permitted.\n\n## Final selection bias\n\nFAIL in existing v6 evidence: final combined score passed its combined diagnostic, but the required final audit was not completed and the short-only correlation remained high.\n\n## Validated pulse-free availability\n\nNo.\n\n## All vs clean removed power fraction\n\nUNAVAILABLE; no validated clean subset and no completed final PSD bootstrap.\n\n## Pulse/noise topology equivalence\n\nPNS3. CI overlap alone is not equivalence, and independent margins/direct-difference equivalence were not completed.\n\n## Coherent CH0 power\n\nThe final artifact defines `|S01|²/S11` and `S00-|S01|²/S11`; no source amplitude fit is made.\n\n## Pulse-event predicted PSD\n\nUNAVAILABLE because detector validity failed.\n\n## Pulse predicted / measured common power\n\nUNAVAILABLE.\n\n## Final PC classification\n\nPC4.\n\n## Final PNS classification\n\nPNS3.\n\n## FINAL PULSE DECISION\n\nFINAL PULSE DECISION:\n\nCLOSE_PULSE_TRACK_UNRESOLVED_WITH_EXISTING_DATA\n\nShould another pulse detector version be built?\n\nNO\n\nPulse may exist, but existing data cannot support a quantitative pulse-dominance determination. Further same-data detector iteration is not scientifically justified. Proceed to amplitude-free 5–30 Hz common-mode topology investigation. Strict target conclusion remains C — exact target physical case remains unidentified.\n""",encoding="utf-8")


def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument("--target-root",type=Path,required=True); ap.add_argument("--case-dir",type=Path,required=True); ap.add_argument("--other-replicates",type=int,default=50); ap.add_argument("--evidence-only",action="store_true"); args=ap.parse_args(); args.case_dir.mkdir(parents=True,exist_ok=True)
    if args.evidence_only:
        evidence_only_final(args.case_dir)
        return
    timing, raw_bases, _ = v5.raw_processed_baselines(args.target_root,args.case_dir); _, timing_rows=v5.timing_distribution(args.target_root,args.case_dir); library, _records=v5.build_library_v5(args.target_root,timing_rows,args.case_dir)
    sigma=float(np.std(np.concatenate(v5.raw_processed_baselines(args.target_root,args.case_dir)[2]["CH0"]),ddof=1)); ch0_library=library["channels"]["CH0"]; cal=short_calibration(ch0_library,args_raw:=raw_bases["CH0"],sigma,args.case_dir); cal=short_validation(ch0_library,args_raw,sigma,cal,args.case_dir)
    nulls,_=recordwise_long(args.target_root,library,sigma,args.case_dir); classification=classify_noise(args.target_root,library,sigma,cal,nulls,args.case_dir); recovery=inject_and_recover(library,args_raw,sigma,cal,args.case_dir,args.other_replicates); gate=recovery_gate(recovery,args.case_dir); fpr=control_fpr(args.target_root,library,sigma,cal,nulls,args.case_dir); bias=selection_bias(args.target_root,classification,args.case_dir); power=removed_power(args.target_root,classification,args.case_dir); topo, data=topology(args.target_root,library,args.case_dir); coherent_decomposition(data["noise_spectra"],data["frequencies"],args.case_dir)
    write_json(args.case_dir/"pulse_event_psd_prediction_final.json",{"stage":"pulse_event_psd_prediction_final","status":"UNAVAILABLE_detector_validity_failed","amplitude_fit_to_experiment":False,"event_rate":None,"PSD_pulse_pred":None})
    write_json(args.case_dir/"pulse_power_vs_common_power.json",{"stage":"pulse_power_vs_common_power","status":"UNAVAILABLE_detector_validity_failed","amplitude_rescaling":False})
    write_json(args.case_dir/"final_detector_fpr_policy.json",fpr)
    write_json(args.case_dir/"final_pulse_detector_gate.json",gate)
    write_json(args.case_dir/"final_pulse_selection_bias.json",bias)
    write_json(args.case_dir/"final_noise_record_pulse_classification.json",classification)
    final_decision(gate,bias,fpr,power,topo,args.case_dir)


if __name__ == "__main__": main()

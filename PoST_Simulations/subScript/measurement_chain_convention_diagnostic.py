"""Falsify measurement-chain convention variants with detector physics frozen."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Opt_noise as opt  # noqa: E402

ANALOG_CUTOFFS_HZ = (80_000.0, 100_000.0, 120_000.0)
ANALOG_ORDERS = (2, 4, 6)
ANALOG_NORMS = ("phase", "mag", "delay")
DIGITAL_CUTOFFS_HZ = (8_000.0, 10_000.0, 12_000.0)
DIGITAL_ORDERS = (1, 2, 3, 4)
DIGITAL_PASSES = (1, 2)

BASELINE_CHAIN = {
    "analog_cutoff_Hz": 100_000.0,
    "analog_order": 4,
    "analog_norm": "mag",
    "digital_cutoff_Hz": 10_000.0,
    "digital_order": 2,
    "digital_passes": 2,
}


def fit_args(summary):
    f = summary["fit"]
    return SimpleNamespace(
        fit_min_hz=float(f["min_hz"]),
        fit_max_hz=float(f["max_hz"]),
        fit_weight_start_hz=float(f["weight_start_hz"]),
        high_frequency_weight=float(f["high_frequency_weight"]),
        robust_delta_dex=float(f["robust_delta_dex"]),
        fit_points=int(f["points"]),
        absolute_asd_weight=float(f.get("absolute_asd_weight", 0.0)),
    )


def digital_bessel_magnitude(
    frequency_hz,
    rate_hz,
    cutoff_hz,
    order,
    passes,
):
    frequency = np.asarray(frequency_hz, dtype=float)
    rate = float(rate_hz)
    cutoff = float(cutoff_hz)
    order = int(order)
    passes = int(passes)
    if rate <= 0.0 or cutoff <= 0.0 or cutoff >= rate / 2.0:
        raise ValueError("invalid rate or digital cutoff")
    if order < 1 or passes < 1:
        raise ValueError("order and passes must be positive")
    if np.any(frequency < 0.0) or np.any(frequency > rate / 2.0):
        raise ValueError("frequency outside digital Nyquist interval")

    b, a = signal.bessel(order, cutoff / (rate / 2.0), "low")
    angular = 2.0 * np.pi * frequency / rate
    _, response = signal.freqz(b, a, worN=angular)
    return np.abs(response) ** passes


def intrinsic_context(candidate, frequency):
    trial = dict(candidate)
    opt.apply_post_filter_white_fraction(trial)
    frequency = np.asarray(frequency, dtype=float)
    rate = float(trial["rate"])
    alias_frequency = rate - frequency
    query = np.unique(np.concatenate((frequency, alias_frequency)))
    intrinsic = opt.tes_noise_components(trial, query)
    total = np.asarray(intrinsic["total_ch0"], dtype=float)
    main = np.interp(frequency, query, total)
    alias = np.interp(alias_frequency, query, total)
    same_bin = np.isclose(
        alias_frequency,
        frequency,
        rtol=0.0,
        atol=max(rate, 1.0) * 1e-12,
    )
    return {
        "candidate": trial,
        "rate_Hz": rate,
        "frequency_Hz": frequency,
        "alias_frequency_Hz": alias_frequency,
        "main_intrinsic_asd": main,
        "alias_intrinsic_asd": alias,
        "same_bin": same_bin,
        "post_filter_white_asd_A_rtHz": float(
            trial.get("post_filter_white_asd_A_rtHz", 0.0)
        ),
    }


def chain_model(context, config):
    frequency = context["frequency_Hz"]
    alias_frequency = context["alias_frequency_Hz"]
    main_response = opt.hardware_filter_magnitude(
        frequency,
        cutoff_hz=float(config["analog_cutoff_Hz"]),
        order=int(config["analog_order"]),
        norm=str(config["analog_norm"]),
    )
    alias_response = opt.hardware_filter_magnitude(
        alias_frequency,
        cutoff_hz=float(config["analog_cutoff_Hz"]),
        order=int(config["analog_order"]),
        norm=str(config["analog_norm"]),
    )
    main = context["main_intrinsic_asd"] * main_response
    alias = context["alias_intrinsic_asd"] * alias_response
    alias = np.where(context["same_bin"], 0.0, alias)
    detector = np.sqrt(main**2 + alias**2)
    white = float(context["post_filter_white_asd_A_rtHz"])
    pre_analysis = np.sqrt(detector**2 + white**2)
    digital = digital_bessel_magnitude(
        frequency,
        context["rate_Hz"],
        float(config["digital_cutoff_Hz"]),
        int(config["digital_order"]),
        int(config["digital_passes"]),
    )
    absolute = pre_analysis * digital
    normalized = opt.normalize_at(
        frequency,
        absolute,
        reference_hz=opt.ABSOLUTE_ASD_REFERENCE_HZ,
    )
    return normalized


def band_summary(model, target, frequency, args):
    raw = opt.band_fit_diagnostics(model, target, frequency, args)
    out = {}
    for name, row in raw.items():
        mean = float(row["mean_log10_ratio"])
        out[name] = {
            "mean_log10_ratio": mean,
            "geometric_mean_model_over_measurement": float(10.0**mean),
            "rms_log10_ratio": float(row["rms_log10_ratio"]),
            "max_abs_log10_ratio": float(row["max_abs_log10_ratio"]),
        }
    return out


def config_key(config):
    return (
        f"a{int(config['analog_cutoff_Hz'])}_o{int(config['analog_order'])}_"
        f"{config['analog_norm']}_d{int(config['digital_cutoff_Hz'])}_"
        f"o{int(config['digital_order'])}_p{int(config['digital_passes'])}"
    )


def run(candidate, target, frequency, args):
    context = intrinsic_context(candidate, frequency)
    production_model, _ = opt.deterministic_simulated_spectrum(
        dict(candidate), frequency
    )
    reconstructed_baseline = chain_model(context, BASELINE_CHAIN)
    consistency = float(
        np.max(
            np.abs(reconstructed_baseline - production_model)
            / np.maximum(production_model, np.finfo(float).tiny)
        )
    )
    baseline_score = float(
        opt.fit_score(production_model, target, frequency, args)
    )
    baseline_bands = band_summary(
        production_model, target, frequency, args
    )
    mid_name = "5000-15000_Hz"
    high_name = "40000-100000_Hz"
    tail_name = "100000-200000_Hz"
    bm = float(baseline_bands[mid_name]["mean_log10_ratio"])
    bh = float(baseline_bands[high_name]["mean_log10_ratio"])
    bt = float(baseline_bands[tail_name]["mean_log10_ratio"])

    rows = []
    for analog_cutoff in ANALOG_CUTOFFS_HZ:
        for analog_order in ANALOG_ORDERS:
            for analog_norm in ANALOG_NORMS:
                for digital_cutoff in DIGITAL_CUTOFFS_HZ:
                    for digital_order in DIGITAL_ORDERS:
                        for digital_passes in DIGITAL_PASSES:
                            config = {
                                "analog_cutoff_Hz": float(analog_cutoff),
                                "analog_order": int(analog_order),
                                "analog_norm": str(analog_norm),
                                "digital_cutoff_Hz": float(digital_cutoff),
                                "digital_order": int(digital_order),
                                "digital_passes": int(digital_passes),
                            }
                            model = chain_model(context, config)
                            bands = band_summary(
                                model, target, frequency, args
                            )
                            mid = float(
                                bands[mid_name]["mean_log10_ratio"]
                            )
                            high = float(
                                bands[high_name]["mean_log10_ratio"]
                            )
                            tail = float(
                                bands[tail_name]["mean_log10_ratio"]
                            )
                            score = float(
                                opt.fit_score(
                                    model, target, frequency, args
                                )
                            )
                            rows.append(
                                {
                                    "key": config_key(config),
                                    **config,
                                    "is_production_baseline": bool(
                                        config == BASELINE_CHAIN
                                    ),
                                    "shape_score": score,
                                    "score_ratio_to_baseline": float(
                                        score / baseline_score
                                    ),
                                    "bands": bands,
                                    "correct_direction": bool(
                                        mid > bm and high < bh
                                    ),
                                    "improves_mid_high": bool(
                                        abs(mid) < abs(bm)
                                        and abs(high) < abs(bh)
                                    ),
                                    "strict_improvement": bool(
                                        abs(mid) < abs(bm)
                                        and abs(high) < abs(bh)
                                        and abs(tail) <= abs(bt) + 1e-12
                                    ),
                                }
                            )

    correct = [r for r in rows if r["correct_direction"]]
    improve = [r for r in rows if r["improves_mid_high"]]
    strict = [r for r in rows if r["strict_improvement"]]
    best = min(rows, key=lambda row: row["shape_score"])
    best_correct = min(
        correct,
        key=lambda row: row["shape_score"],
        default=None,
    )
    best_strict = min(
        strict,
        key=lambda row: row["shape_score"],
        default=None,
    )

    by_analog_norm = {}
    for norm in ANALOG_NORMS:
        subset = [r for r in rows if r["analog_norm"] == norm]
        by_analog_norm[norm] = min(
            subset, key=lambda row: row["shape_score"]
        )
    by_digital_passes = {}
    for passes in DIGITAL_PASSES:
        subset = [
            r for r in rows if r["digital_passes"] == passes
        ]
        by_digital_passes[str(passes)] = min(
            subset, key=lambda row: row["shape_score"]
        )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "detector_parameters_held_fixed": True,
        "varied_measurement_chain_only": True,
        "baseline_chain": BASELINE_CHAIN,
        "baseline_shape_score": baseline_score,
        "baseline_bands": baseline_bands,
        "baseline_reconstruction_max_relative_error": consistency,
        "baseline_reconstruction_matches_production": bool(
            consistency <= 1e-10
        ),
        "grid": {
            "analog_cutoffs_Hz": list(ANALOG_CUTOFFS_HZ),
            "analog_orders": list(ANALOG_ORDERS),
            "analog_norms": list(ANALOG_NORMS),
            "digital_cutoffs_Hz": list(DIGITAL_CUTOFFS_HZ),
            "digital_orders": list(DIGITAL_ORDERS),
            "digital_passes": list(DIGITAL_PASSES),
            "total_cases": len(rows),
        },
        "best_global_row": best,
        "best_correct_direction_row": best_correct,
        "best_strict_row": best_strict,
        "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction": bool(
            correct
        ),
        "can_improve_5_15k_and_40_100k_absolute_error": bool(improve),
        "can_improve_5_15k_and_40_100k_without_worsening_100_200k": bool(
            strict
        ),
        "best_by_analog_norm": by_analog_norm,
        "best_by_digital_passes": by_digital_passes,
        "rows": rows,
        "guardrail": (
            "This is a convention/cutoff/order falsification screen, not a "
            "new production fit. A better row must be checked against the "
            "actual acquisition and analysis implementation before adoption."
        ),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--output", type=Path)
    a = p.parse_args()

    summary = json.loads(a.summary.read_text(encoding="utf-8"))
    candidate = dict(summary["best_case_parameters"])
    args = fit_args(summary)
    frequency, target, _ = opt.target_spectrum(args)
    result = run(candidate, target, frequency, args)
    result["source_summary"] = str(a.summary)
    output = a.output or a.summary.with_name(
        "measurement_chain_convention_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(output),
        "baseline_reconstruction_matches_production": result[
            "baseline_reconstruction_matches_production"
        ],
        "total_cases": result["grid"]["total_cases"],
        "best_score": result["best_global_row"]["shape_score"],
        "best_score_ratio": result[
            "best_global_row"
        ]["score_ratio_to_baseline"],
        "correct_direction": result[
            "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction"
        ],
        "strict_improvement": result[
            "can_improve_5_15k_and_40_100k_without_worsening_100_200k"
        ],
        "best_key": result["best_global_row"]["key"],
    }, indent=2))


if __name__ == "__main__":
    main()

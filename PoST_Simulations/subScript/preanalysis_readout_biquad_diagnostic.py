"""Falsify pre-ADC readout-transfer biquads against the fresh pre-analysis target.

Detector parameters and the production 100 kHz analog Bessel are frozen.
Only one or two stable minimum-phase complex-conjugate pole/zero sections are
inserted before ADC alias folding. The 10 kHz digital analysis filter is not
part of either the target or model in this diagnostic.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import differential_evolution, least_squares

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    estimate_one_sided_asd,
)

REFERENCE_HZ = 1_000.0
FAMILIES = (0, 1, 2)
CENTER_MIN_HZ_DEFAULT = 1_000.0
CENTER_MAX_HZ_DEFAULT = 300_000.0
Q_MIN_DEFAULT = 0.55
Q_MAX_DEFAULT = 20.0
RMS_SCREEN_DB_DEFAULT = 1.0
MAX_SCREEN_DB_DEFAULT = 3.0


def fit_args(summary):
    fit = summary["fit"]
    return SimpleNamespace(
        fit_min_hz=float(fit["min_hz"]),
        fit_max_hz=float(fit["max_hz"]),
        fit_weight_start_hz=float(fit["weight_start_hz"]),
        high_frequency_weight=float(fit["high_frequency_weight"]),
        robust_delta_dex=float(fit["robust_delta_dex"]),
        fit_points=int(fit["points"]),
        absolute_asd_weight=float(fit.get("absolute_asd_weight", 0.0)),
    )


def read_record(path: Path) -> np.ndarray:
    return np.frombuffer(path.read_bytes()[4:], dtype=np.float64).copy()


def reconstruct_pre_analysis_target(
    comparison_summary: dict,
    experiment_path: Path,
    fit_frequency: np.ndarray,
):
    acquisition = comparison_summary["acquisition"]
    rate = float(acquisition["rate_Hz"])
    sample = int(acquisition["samples"])
    accepted_indices = [
        int(value) for value in acquisition["accepted_record_indices"]
    ]
    raw_dir = experiment_path / "CH0_noise" / "rawdata"
    paths = sorted(raw_dir.glob("CH0_*.dat"))
    if not paths:
        raise FileNotFoundError(f"no CH0 raw records found in {raw_dir}")
    if not accepted_indices:
        raise ValueError("comparison summary contains no accepted records")
    if max(accepted_indices) >= len(paths):
        raise IndexError("accepted record index exceeds raw record list")

    records = [read_record(paths[index]) for index in accepted_indices]
    pre_asd, count = estimate_one_sided_asd(
        records,
        sample,
        rate,
        cutoff=0.0,
        remove_mean=True,
    )
    if count != len(accepted_indices):
        raise RuntimeError("pre-analysis reconstruction changed accepted mask")
    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    normalized = opt.normalize_at(
        frequency,
        pre_asd,
        reference_hz=REFERENCE_HZ,
    )
    target = np.interp(fit_frequency, frequency, normalized)
    return {
        "target": target,
        "full_frequency_Hz": frequency,
        "pre_analysis_asd": pre_asd,
        "accepted_records": int(count),
        "rate_Hz": rate,
        "sample": sample,
    }


def biquad_magnitude(frequency_hz, pole_hz, pole_q, zero_hz, zero_q):
    """DC-normalized analog pole/zero biquad magnitude.

    The section is
      [(s/wz)^2 + s/(Qz*wz) + 1] /
      [(s/wp)^2 + s/(Qp*wp) + 1].
    Positive Q gives stable/minimum-phase LHP roots; Q > 0.5 restricts this
    diagnostic to complex-conjugate pole and zero pairs.
    """

    frequency = np.asarray(frequency_hz, dtype=float)
    fp = float(pole_hz)
    fz = float(zero_hz)
    qp = float(pole_q)
    qz = float(zero_q)
    if fp <= 0.0 or fz <= 0.0:
        raise ValueError("biquad centers must be positive")
    if qp <= 0.5 or qz <= 0.5:
        raise ValueError("biquad Q must exceed 0.5 for complex-conjugate roots")
    xp = frequency / fp
    xz = frequency / fz
    numerator = np.sqrt((1.0 - xz**2) ** 2 + (xz / qz) ** 2)
    denominator = np.sqrt((1.0 - xp**2) ** 2 + (xp / qp) ** 2)
    return numerator / denominator


def sections_magnitude(frequency_hz, sections):
    frequency = np.asarray(frequency_hz, dtype=float)
    response = np.ones_like(frequency)
    for section in sections:
        response *= biquad_magnitude(
            frequency,
            section["pole_Hz"],
            section["pole_Q"],
            section["zero_Hz"],
            section["zero_Q"],
        )
    return response


def intrinsic_context(candidate, frequency):
    trial = dict(candidate)
    opt.apply_post_filter_white_fraction(trial)
    frequency = np.asarray(frequency, dtype=float)
    rate = float(trial["rate"])
    if np.any(frequency < 0.0) or np.any(frequency > rate / 2.0):
        raise ValueError("fit frequencies must lie below Nyquist")
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


def pre_analysis_model(context, sections):
    frequency = context["frequency_Hz"]
    alias_frequency = context["alias_frequency_Hz"]
    hardware_main = opt.hardware_filter_magnitude(
        frequency,
        cutoff_hz=opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=opt.TARGET_HARDWARE_BESSEL_ORDER,
        norm=opt.TARGET_HARDWARE_BESSEL_NORM,
    )
    hardware_alias = opt.hardware_filter_magnitude(
        alias_frequency,
        cutoff_hz=opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=opt.TARGET_HARDWARE_BESSEL_ORDER,
        norm=opt.TARGET_HARDWARE_BESSEL_NORM,
    )
    extra_main = sections_magnitude(frequency, sections)
    extra_alias = sections_magnitude(alias_frequency, sections)

    main = context["main_intrinsic_asd"] * hardware_main * extra_main
    alias = context["alias_intrinsic_asd"] * hardware_alias * extra_alias
    alias = np.where(context["same_bin"], 0.0, alias)
    folded = np.sqrt(main**2 + alias**2)

    # Preserve the production semantics: this empirical white term is added
    # after hardware/alias folding and therefore is not shaped by the
    # diagnostic pre-ADC biquad.
    white = float(context["post_filter_white_asd_A_rtHz"])
    absolute = np.sqrt(folded**2 + white**2)
    return opt.normalize_at(
        frequency,
        absolute,
        reference_hz=REFERENCE_HZ,
    )


def decode(vector, n_sections):
    values = np.asarray(vector, dtype=float)
    if values.size != 4 * int(n_sections):
        raise ValueError("wrong parameter vector length")
    sections = []
    for index in range(int(n_sections)):
        offset = 4 * index
        sections.append(
            {
                "pole_Hz": float(10.0 ** values[offset]),
                "pole_Q": float(10.0 ** values[offset + 1]),
                "zero_Hz": float(10.0 ** values[offset + 2]),
                "zero_Q": float(10.0 ** values[offset + 3]),
            }
        )
    return sorted(sections, key=lambda row: row["pole_Hz"])


def band_summary(model, target, frequency, args):
    raw = opt.band_fit_diagnostics(model, target, frequency, args)
    result = {}
    for name, row in raw.items():
        mean = float(row["mean_log10_ratio"])
        result[name] = {
            "mean_log10_ratio": mean,
            "mean_model_over_target_dB": float(20.0 * mean),
            "geometric_mean_model_over_target": float(10.0**mean),
            "rms_log10_ratio": float(row["rms_log10_ratio"]),
            "max_abs_log10_ratio": float(row["max_abs_log10_ratio"]),
        }
    return result


def residual_db_metrics(model, target, frequency, args):
    residual = 20.0 * np.log10(model / target)
    return {
        "rms_residual_dB": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual_dB": float(np.max(np.abs(residual))),
        "mean_residual_dB": float(np.mean(residual)),
        "bands": {
            key: {
                "mean_residual_dB": float(
                    20.0 * row["mean_log10_ratio"]
                ),
                "rms_residual_dB": float(
                    20.0 * row["rms_log10_ratio"]
                ),
                "max_abs_residual_dB": float(
                    20.0 * row["max_abs_log10_ratio"]
                ),
            }
            for key, row in opt.band_fit_diagnostics(
                model, target, frequency, args
            ).items()
        },
    }


def fit_family(
    n_sections,
    context,
    target,
    frequency,
    args,
    *,
    center_min_hz,
    center_max_hz,
    q_min,
    q_max,
    seed,
    de_maxiter,
    rms_screen_db,
    max_screen_db,
):
    baseline = pre_analysis_model(context, [])
    if n_sections == 0:
        sections = []
        model = baseline
        success = True
        evaluations = 0
        boundary_hits = []
    else:
        log_f_min = np.log10(float(center_min_hz))
        log_f_max = np.log10(float(center_max_hz))
        log_q_min = np.log10(float(q_min))
        log_q_max = np.log10(float(q_max))
        bounds = []
        for _ in range(int(n_sections)):
            bounds.extend(
                [
                    (log_f_min, log_f_max),
                    (log_q_min, log_q_max),
                    (log_f_min, log_f_max),
                    (log_q_min, log_q_max),
                ]
            )

        def model_from_vector(vector):
            return pre_analysis_model(
                context,
                decode(vector, n_sections),
            )

        def objective(vector):
            return float(
                opt.fit_score(
                    model_from_vector(vector),
                    target,
                    frequency,
                    args,
                )
            )

        de = differential_evolution(
            objective,
            bounds,
            seed=int(seed),
            maxiter=int(de_maxiter),
            popsize=10,
            tol=1e-8,
            polish=False,
            workers=1,
            updating="immediate",
        )
        lower = np.asarray([item[0] for item in bounds], dtype=float)
        upper = np.asarray([item[1] for item in bounds], dtype=float)
        ls = least_squares(
            lambda vector: opt.weighted_residual_vector(
                model_from_vector(vector),
                target,
                frequency,
                args,
            ),
            de.x,
            bounds=(lower, upper),
            max_nfev=2500,
        )
        candidate_vectors = [de.x, ls.x]
        best_vector = min(candidate_vectors, key=objective)
        sections = decode(best_vector, n_sections)
        model = model_from_vector(best_vector)
        success = bool(de.success or ls.success)
        evaluations = int(de.nfev + ls.nfev)
        span = upper - lower
        position = (best_vector - lower) / span
        boundary_hits = [
            {
                "parameter_index": int(index),
                "position_fraction": float(value),
                "within_1pct_of_bound": bool(
                    value <= 0.01 or value >= 0.99
                ),
            }
            for index, value in enumerate(position)
        ]

    score = float(opt.fit_score(model, target, frequency, args))
    metrics = residual_db_metrics(model, target, frequency, args)
    correction = model / baseline
    correction = correction / np.interp(
        REFERENCE_HZ, frequency, correction
    )
    passes = bool(
        metrics["rms_residual_dB"] <= float(rms_screen_db)
        and metrics["max_abs_residual_dB"] <= float(max_screen_db)
    )
    return {
        "family": f"{n_sections}_biquad",
        "n_sections": int(n_sections),
        "n_parameters": int(4 * n_sections),
        "stable_minimum_phase_complex_sections": True,
        "success": success,
        "evaluations": evaluations,
        "sections": sections,
        "shape_score": score,
        "residual_metrics": metrics,
        "passes_screen_tolerance": passes,
        "boundary_positions": boundary_hits,
        "_model": model,
        "_correction": correction,
    }


def run(
    candidate,
    target,
    frequency,
    args,
    *,
    center_min_hz=CENTER_MIN_HZ_DEFAULT,
    center_max_hz=CENTER_MAX_HZ_DEFAULT,
    q_min=Q_MIN_DEFAULT,
    q_max=Q_MAX_DEFAULT,
    seed=20260914,
    de_maxiter=100,
    rms_screen_db=RMS_SCREEN_DB_DEFAULT,
    max_screen_db=MAX_SCREEN_DB_DEFAULT,
):
    context = intrinsic_context(candidate, frequency)
    rows = []
    for index, n_sections in enumerate(FAMILIES):
        rows.append(
            fit_family(
                n_sections,
                context,
                target,
                frequency,
                args,
                center_min_hz=center_min_hz,
                center_max_hz=center_max_hz,
                q_min=q_min,
                q_max=q_max,
                seed=seed + index,
                de_maxiter=de_maxiter,
                rms_screen_db=rms_screen_db,
                max_screen_db=max_screen_db,
            )
        )

    baseline = rows[0]
    baseline_bands = band_summary(
        baseline["_model"], target, frequency, args
    )
    mid_name = "5000-15000_Hz"
    high_name = "40000-100000_Hz"
    tail_name = "100000-200000_Hz"
    bm = float(baseline_bands[mid_name]["mean_log10_ratio"])
    bh = float(baseline_bands[high_name]["mean_log10_ratio"])
    bt = float(baseline_bands[tail_name]["mean_log10_ratio"])

    clean_rows = []
    for row in rows:
        bands = band_summary(row["_model"], target, frequency, args)
        mid = float(bands[mid_name]["mean_log10_ratio"])
        high = float(bands[high_name]["mean_log10_ratio"])
        tail = float(bands[tail_name]["mean_log10_ratio"])
        clean_rows.append(
            {
                key: value
                for key, value in row.items()
                if not key.startswith("_")
            }
            | {
                "score_ratio_to_baseline": float(
                    row["shape_score"] / baseline["shape_score"]
                ),
                "bands": bands,
                "correct_direction": bool(mid > bm and high < bh),
                "improves_mid_high": bool(
                    abs(mid) < abs(bm) and abs(high) < abs(bh)
                ),
                "strict_improvement": bool(
                    abs(mid) < abs(bm)
                    and abs(high) < abs(bh)
                    and abs(tail) <= abs(bt) + 1e-12
                ),
            }
        )

    best = min(clean_rows, key=lambda row: row["shape_score"])
    joint = [row for row in clean_rows if row["improves_mid_high"]]
    strict = [row for row in clean_rows if row["strict_improvement"]]
    passing = [
        row for row in clean_rows if row["passes_screen_tolerance"]
    ]

    best_row_internal = min(rows, key=lambda row: row["shape_score"])
    sample_frequencies = (
        1_000.0,
        5_000.0,
        10_000.0,
        20_000.0,
        40_000.0,
        70_000.0,
        100_000.0,
        150_000.0,
        200_000.0,
    )
    sample = []
    for anchor in sample_frequencies:
        index = int(np.argmin(np.abs(frequency - anchor)))
        sample.append(
            {
                "frequency_Hz": float(frequency[index]),
                "target_pre_analysis_normalized": float(target[index]),
                "baseline_model_normalized": float(
                    baseline["_model"][index]
                ),
                "best_model_normalized": float(
                    best_row_internal["_model"][index]
                ),
                "best_end_to_end_correction": float(
                    best_row_internal["_correction"][index]
                ),
                "best_end_to_end_correction_dB": float(
                    20.0
                    * np.log10(
                        best_row_internal["_correction"][index]
                    )
                ),
            }
        )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "detector_parameters_held_fixed": True,
        "digital_analysis_filter_applied": False,
        "target_semantics": (
            "fresh pre-analysis ASD from the exact accepted raw CH0 mask; "
            "mean removal, Hann, power average, one-sided ASD; no 10 kHz "
            "digital analysis filter"
        ),
        "model_semantics": (
            "frozen intrinsic detector ASD -> production 100 kHz analog "
            "Bessel times diagnostic pre-ADC biquad(s) -> first alias fold "
            "-> production post-filter white term; no digital filter"
        ),
        "biquad_definition": (
            "DC-normalized analog second-order pole/zero sections with "
            "complex-conjugate LHP poles and zeros (Q > 0.5)"
        ),
        "search_box": {
            "center_Hz": [
                float(center_min_hz),
                float(center_max_hz),
            ],
            "Q": [float(q_min), float(q_max)],
            "families": ["0_biquad", "1_biquad", "2_biquad"],
        },
        "screen_tolerance": {
            "rms_residual_dB_max": float(rms_screen_db),
            "max_abs_residual_dB_max": float(max_screen_db),
            "diagnostic_not_physical_prior": True,
        },
        "baseline_shape_score": float(baseline["shape_score"]),
        "baseline_bands": baseline_bands,
        "families": clean_rows,
        "best_global_row": best,
        "best_joint_improvement_row": (
            min(joint, key=lambda row: row["shape_score"])
            if joint
            else None
        ),
        "best_strict_row": (
            min(strict, key=lambda row: row["shape_score"])
            if strict
            else None
        ),
        "smallest_passing_family": (
            min(
                passing,
                key=lambda row: (
                    row["n_parameters"],
                    row["shape_score"],
                ),
            )
            if passing
            else None
        ),
        "any_biquad_improves_mid_high": bool(joint),
        "any_biquad_strict_improvement": bool(strict),
        "any_family_passes_screen_tolerance": bool(passing),
        "curve_sample_best_global": sample,
        "guardrail": (
            "A successful biquad fit shows that a structured pre-ADC "
            "readout/electrical transfer has enough shape leverage. It does "
            "not identify the actual circuit or phase response without "
            "independent hardware evidence."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument(
        "--comparison-summary",
        type=Path,
        required=True,
    )
    parser.add_argument("--experiment-path", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--center-min-hz",
        type=float,
        default=CENTER_MIN_HZ_DEFAULT,
    )
    parser.add_argument(
        "--center-max-hz",
        type=float,
        default=CENTER_MAX_HZ_DEFAULT,
    )
    parser.add_argument("--q-min", type=float, default=Q_MIN_DEFAULT)
    parser.add_argument("--q-max", type=float, default=Q_MAX_DEFAULT)
    parser.add_argument("--de-maxiter", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260914)
    parser.add_argument(
        "--rms-screen-db",
        type=float,
        default=RMS_SCREEN_DB_DEFAULT,
    )
    parser.add_argument(
        "--max-screen-db",
        type=float,
        default=MAX_SCREEN_DB_DEFAULT,
    )
    args = parser.parse_args()

    if args.q_min <= 0.5 or args.q_max <= args.q_min:
        raise ValueError("require 0.5 < q_min < q_max")
    if (
        args.center_min_hz <= 0.0
        or args.center_max_hz <= args.center_min_hz
    ):
        raise ValueError("invalid center-frequency search box")

    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    comparison = json.loads(
        args.comparison_summary.read_text(encoding="utf-8")
    )
    experiment_path = (
        args.experiment_path
        if args.experiment_path is not None
        else Path(comparison["experiment_path"])
    )
    fit = fit_args(summary)
    frequency = np.geomspace(
        fit.fit_min_hz,
        fit.fit_max_hz,
        fit.fit_points,
    )
    target_context = reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        frequency,
    )
    candidate = dict(summary["best_case_parameters"])
    result = run(
        candidate,
        target_context["target"],
        frequency,
        fit,
        center_min_hz=args.center_min_hz,
        center_max_hz=args.center_max_hz,
        q_min=args.q_min,
        q_max=args.q_max,
        seed=args.seed,
        de_maxiter=args.de_maxiter,
        rms_screen_db=args.rms_screen_db,
        max_screen_db=args.max_screen_db,
    )
    result["comparison_summary"] = str(args.comparison_summary)
    result["optimizer_summary"] = str(args.summary)
    result["experiment_path"] = str(experiment_path)
    result["accepted_records"] = target_context["accepted_records"]
    result["acquisition"] = {
        "rate_Hz": target_context["rate_Hz"],
        "samples": target_context["sample"],
    }

    output = args.output or args.summary.with_name(
        "preanalysis_readout_biquad_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "accepted_records": result["accepted_records"],
                "baseline_shape_score": result[
                    "baseline_shape_score"
                ],
                "best_family": result["best_global_row"]["family"],
                "best_shape_score": result[
                    "best_global_row"
                ]["shape_score"],
                "best_score_ratio": result[
                    "best_global_row"
                ]["score_ratio_to_baseline"],
                "any_joint_improvement": result[
                    "any_biquad_improves_mid_high"
                ],
                "any_strict_improvement": result[
                    "any_biquad_strict_improvement"
                ],
                "any_screen_pass": result[
                    "any_family_passes_screen_tolerance"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

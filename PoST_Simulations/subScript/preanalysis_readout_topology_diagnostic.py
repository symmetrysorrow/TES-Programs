"""Compare structured pre-ADC readout-transfer topologies at matched complexity.

The detector model, production 100 kHz analog Bessel, ADC first-alias fold,
and production post-filter white term are frozen.  The target is the fresh
pre-analysis ASD from the exact accepted CH0 record mask, with no 10 kHz
digital analysis filter.

The core comparison uses four free parameters for each topology:
  * complex_biquad_4p: second-order pole/zero pair with Q > 0.5,
  * general_second_order_4p: same form with Q > 0, allowing Q < 0.5,
  * real_2p2z_4p: two real first-order poles and two real first-order zeros.

A six-parameter hybrid (general second order + one real lead/lag) is included
as a secondary complexity-ladder check, not as a matched-DOF winner.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, least_squares

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402

REFERENCE_HZ = 1_000.0
CENTER_MIN_HZ_DEFAULT = 1_000.0
CENTER_MAX_HZ_DEFAULT = 300_000.0
COMPLEX_Q_MIN_DEFAULT = 0.55
GENERAL_Q_MIN_DEFAULT = 0.10
Q_MAX_DEFAULT = 20.0
RMS_SCREEN_DB_DEFAULT = 1.0
MAX_SCREEN_DB_DEFAULT = 3.0

CORE_FAMILIES = (
    "complex_biquad_4p",
    "general_second_order_4p",
    "real_2p2z_4p",
)
ALL_FAMILIES = CORE_FAMILIES + ("hybrid_general_plus_leadlag_6p",)


def second_order_magnitude(
    frequency_hz,
    pole_hz,
    pole_q,
    zero_hz,
    zero_q,
):
    """DC-normalized stable/minimum-phase second-order magnitude for Q > 0."""

    frequency = np.asarray(frequency_hz, dtype=float)
    fp = float(pole_hz)
    fz = float(zero_hz)
    qp = float(pole_q)
    qz = float(zero_q)
    if fp <= 0.0 or fz <= 0.0:
        raise ValueError("second-order center frequencies must be positive")
    if qp <= 0.0 or qz <= 0.0:
        raise ValueError("second-order Q values must be positive")

    xp = frequency / fp
    xz = frequency / fz
    numerator = np.sqrt((1.0 - xz**2) ** 2 + (xz / qz) ** 2)
    denominator = np.sqrt((1.0 - xp**2) ** 2 + (xp / qp) ** 2)
    return numerator / denominator


def real_pole_zero_magnitude(frequency_hz, poles_hz, zeros_hz):
    """DC-normalized stable/minimum-phase first-order real pole/zero network."""

    frequency = np.asarray(frequency_hz, dtype=float)
    poles = np.asarray(poles_hz, dtype=float)
    zeros = np.asarray(zeros_hz, dtype=float)
    if np.any(poles <= 0.0) or np.any(zeros <= 0.0):
        raise ValueError("real pole/zero corners must be positive")

    response = np.ones_like(frequency)
    for corner in zeros:
        response *= np.sqrt(1.0 + (frequency / corner) ** 2)
    for corner in poles:
        response /= np.sqrt(1.0 + (frequency / corner) ** 2)
    return response


def equivalent_real_roots(center_hz, q):
    """Return real-root corner frequencies for a second-order factor when Q<=0.5."""

    center = float(center_hz)
    q = float(q)
    if center <= 0.0 or q <= 0.0:
        raise ValueError("center and Q must be positive")
    if q > 0.5:
        return None
    term = np.sqrt(max(0.0, 1.0 / (4.0 * q**2) - 1.0))
    low = center * (1.0 / (2.0 * q) - term)
    high = center * (1.0 / (2.0 * q) + term)
    return [float(low), float(high)]


def root_regime(q):
    q = float(q)
    if q < 0.5:
        return "overdamped_real_roots"
    if np.isclose(q, 0.5, rtol=0.0, atol=1e-9):
        return "critical"
    return "complex_conjugate"


def transfer_magnitude(frequency_hz, family, parameters):
    frequency = np.asarray(frequency_hz, dtype=float)
    if family in {"complex_biquad_4p", "general_second_order_4p"}:
        return second_order_magnitude(
            frequency,
            parameters["pole_Hz"],
            parameters["pole_Q"],
            parameters["zero_Hz"],
            parameters["zero_Q"],
        )
    if family == "real_2p2z_4p":
        return real_pole_zero_magnitude(
            frequency,
            parameters["poles_Hz"],
            parameters["zeros_Hz"],
        )
    if family == "hybrid_general_plus_leadlag_6p":
        return second_order_magnitude(
            frequency,
            parameters["pole_Hz"],
            parameters["pole_Q"],
            parameters["zero_Hz"],
            parameters["zero_Q"],
        ) * real_pole_zero_magnitude(
            frequency,
            [parameters["leadlag_pole_Hz"]],
            [parameters["leadlag_zero_Hz"]],
        )
    raise ValueError(f"unknown family: {family}")


def family_bounds(
    family,
    center_min_hz,
    center_max_hz,
    complex_q_min,
    general_q_min,
    q_max,
):
    lf0 = np.log10(float(center_min_hz))
    lf1 = np.log10(float(center_max_hz))
    lqc0 = np.log10(float(complex_q_min))
    lqg0 = np.log10(float(general_q_min))
    lq1 = np.log10(float(q_max))

    if family == "complex_biquad_4p":
        return (
            ["pole_Hz", "pole_Q", "zero_Hz", "zero_Q"],
            [(lf0, lf1), (lqc0, lq1), (lf0, lf1), (lqc0, lq1)],
        )
    if family == "general_second_order_4p":
        return (
            ["pole_Hz", "pole_Q", "zero_Hz", "zero_Q"],
            [(lf0, lf1), (lqg0, lq1), (lf0, lf1), (lqg0, lq1)],
        )
    if family == "real_2p2z_4p":
        return (
            ["pole_1_Hz", "pole_2_Hz", "zero_1_Hz", "zero_2_Hz"],
            [(lf0, lf1)] * 4,
        )
    if family == "hybrid_general_plus_leadlag_6p":
        return (
            [
                "pole_Hz",
                "pole_Q",
                "zero_Hz",
                "zero_Q",
                "leadlag_pole_Hz",
                "leadlag_zero_Hz",
            ],
            [
                (lf0, lf1),
                (lqg0, lq1),
                (lf0, lf1),
                (lqg0, lq1),
                (lf0, lf1),
                (lf0, lf1),
            ],
        )
    raise ValueError(f"unknown family: {family}")


def decode_family(family, vector):
    values = 10.0 ** np.asarray(vector, dtype=float)
    if family in {"complex_biquad_4p", "general_second_order_4p"}:
        return {
            "pole_Hz": float(values[0]),
            "pole_Q": float(values[1]),
            "zero_Hz": float(values[2]),
            "zero_Q": float(values[3]),
        }
    if family == "real_2p2z_4p":
        return {
            "poles_Hz": sorted([float(values[0]), float(values[1])]),
            "zeros_Hz": sorted([float(values[2]), float(values[3])]),
        }
    if family == "hybrid_general_plus_leadlag_6p":
        return {
            "pole_Hz": float(values[0]),
            "pole_Q": float(values[1]),
            "zero_Hz": float(values[2]),
            "zero_Q": float(values[3]),
            "leadlag_pole_Hz": float(values[4]),
            "leadlag_zero_Hz": float(values[5]),
        }
    raise ValueError(f"unknown family: {family}")


def describe_parameters(family, parameters):
    result = dict(parameters)
    if family in {
        "complex_biquad_4p",
        "general_second_order_4p",
        "hybrid_general_plus_leadlag_6p",
    }:
        result["pole_root_regime"] = root_regime(parameters["pole_Q"])
        result["zero_root_regime"] = root_regime(parameters["zero_Q"])
        result["pole_equivalent_real_root_Hz"] = equivalent_real_roots(
            parameters["pole_Hz"],
            parameters["pole_Q"],
        )
        result["zero_equivalent_real_root_Hz"] = equivalent_real_roots(
            parameters["zero_Hz"],
            parameters["zero_Q"],
        )
    return result


def pre_analysis_model(context, family, parameters):
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

    extra_main = (
        np.ones_like(frequency)
        if family == "unity"
        else transfer_magnitude(frequency, family, parameters)
    )
    extra_alias = (
        np.ones_like(alias_frequency)
        if family == "unity"
        else transfer_magnitude(alias_frequency, family, parameters)
    )

    main = context["main_intrinsic_asd"] * hardware_main * extra_main
    alias = context["alias_intrinsic_asd"] * hardware_alias * extra_alias
    alias = np.where(context["same_bin"], 0.0, alias)
    folded = np.sqrt(main**2 + alias**2)

    white = float(context["post_filter_white_asd_A_rtHz"])
    absolute = np.sqrt(folded**2 + white**2)
    return opt.normalize_at(
        frequency,
        absolute,
        reference_hz=REFERENCE_HZ,
    )


def named_boundary_positions(names, vector, bounds):
    lower = np.asarray([row[0] for row in bounds], dtype=float)
    upper = np.asarray([row[1] for row in bounds], dtype=float)
    vector = np.asarray(vector, dtype=float)
    position = (vector - lower) / (upper - lower)
    return {
        name: {
            "position_fraction": float(value),
            "within_1pct_of_bound": bool(value <= 0.01 or value >= 0.99),
            "nearest_bound": "lower" if value <= 0.5 else "upper",
        }
        for name, value in zip(names, position)
    }


def fit_family(
    family,
    context,
    target,
    frequency,
    args,
    *,
    center_min_hz,
    center_max_hz,
    complex_q_min,
    general_q_min,
    q_max,
    seed,
    de_maxiter,
    rms_screen_db,
    max_screen_db,
):
    names, bounds = family_bounds(
        family,
        center_min_hz,
        center_max_hz,
        complex_q_min,
        general_q_min,
        q_max,
    )

    def model_from_vector(vector):
        params = decode_family(family, vector)
        return pre_analysis_model(context, family, params)

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
    lower = np.asarray([row[0] for row in bounds], dtype=float)
    upper = np.asarray([row[1] for row in bounds], dtype=float)
    ls = least_squares(
        lambda vector: opt.weighted_residual_vector(
            model_from_vector(vector),
            target,
            frequency,
            args,
        ),
        de.x,
        bounds=(lower, upper),
        max_nfev=3000,
    )
    best_vector = min((de.x, ls.x), key=objective)
    params = decode_family(family, best_vector)
    model = model_from_vector(best_vector)
    score = float(opt.fit_score(model, target, frequency, args))
    metrics = base.residual_db_metrics(model, target, frequency, args)
    passes = bool(
        metrics["rms_residual_dB"] <= float(rms_screen_db)
        and metrics["max_abs_residual_dB"] <= float(max_screen_db)
    )
    return {
        "family": family,
        "n_parameters": int(len(best_vector)),
        "matched_4p_core_family": bool(family in CORE_FAMILIES),
        "stable_minimum_phase": True,
        "success": bool(de.success or ls.success),
        "evaluations": int(de.nfev + ls.nfev),
        "parameters": describe_parameters(family, params),
        "shape_score": score,
        "residual_metrics": metrics,
        "passes_screen_tolerance": passes,
        "boundary_positions": named_boundary_positions(
            names,
            best_vector,
            bounds,
        ),
        "_model": model,
    }


def band_flags(model, target, frequency, args, baseline_bands):
    bands = base.band_summary(model, target, frequency, args)
    mid_name = "5000-15000_Hz"
    high_name = "40000-100000_Hz"
    tail_name = "100000-200000_Hz"

    bm = float(baseline_bands[mid_name]["mean_log10_ratio"])
    bh = float(baseline_bands[high_name]["mean_log10_ratio"])
    bt = float(baseline_bands[tail_name]["mean_log10_ratio"])
    mid = float(bands[mid_name]["mean_log10_ratio"])
    high = float(bands[high_name]["mean_log10_ratio"])
    tail = float(bands[tail_name]["mean_log10_ratio"])
    return bands, {
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


def curve_sample(frequency, target, baseline, rows):
    anchors = (
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
    samples = []
    for anchor in anchors:
        index = int(np.argmin(np.abs(frequency - anchor)))
        row = {
            "frequency_Hz": float(frequency[index]),
            "target_pre_analysis_normalized": float(target[index]),
            "baseline_model_normalized": float(baseline[index]),
        }
        for fitted in rows:
            model = fitted["_model"]
            correction = float(model[index] / baseline[index])
            reference_index = int(
                np.argmin(np.abs(frequency - REFERENCE_HZ))
            )
            reference_correction = float(
                model[reference_index] / baseline[reference_index]
            )
            correction /= reference_correction
            key = fitted["family"]
            row[f"{key}_model_normalized"] = float(model[index])
            row[f"{key}_correction_dB"] = float(
                20.0 * np.log10(correction)
            )
        samples.append(row)
    return samples


def run(
    candidate,
    target,
    frequency,
    args,
    *,
    center_min_hz=CENTER_MIN_HZ_DEFAULT,
    center_max_hz=CENTER_MAX_HZ_DEFAULT,
    complex_q_min=COMPLEX_Q_MIN_DEFAULT,
    general_q_min=GENERAL_Q_MIN_DEFAULT,
    q_max=Q_MAX_DEFAULT,
    seed=20260914,
    de_maxiter=100,
    rms_screen_db=RMS_SCREEN_DB_DEFAULT,
    max_screen_db=MAX_SCREEN_DB_DEFAULT,
):
    context = base.intrinsic_context(candidate, frequency)
    baseline = pre_analysis_model(context, "unity", {})
    baseline_score = float(
        opt.fit_score(baseline, target, frequency, args)
    )
    baseline_metrics = base.residual_db_metrics(
        baseline, target, frequency, args
    )
    baseline_bands = base.band_summary(
        baseline, target, frequency, args
    )

    internal_rows = []
    for index, family in enumerate(ALL_FAMILIES):
        internal_rows.append(
            fit_family(
                family,
                context,
                target,
                frequency,
                args,
                center_min_hz=center_min_hz,
                center_max_hz=center_max_hz,
                complex_q_min=complex_q_min,
                general_q_min=general_q_min,
                q_max=q_max,
                seed=seed + index,
                de_maxiter=de_maxiter,
                rms_screen_db=rms_screen_db,
                max_screen_db=max_screen_db,
            )
        )

    clean_rows = []
    for row in internal_rows:
        bands, flags = band_flags(
            row["_model"],
            target,
            frequency,
            args,
            baseline_bands,
        )
        clean_rows.append(
            {
                key: value
                for key, value in row.items()
                if not key.startswith("_")
            }
            | {
                "score_ratio_to_baseline": float(
                    row["shape_score"] / baseline_score
                ),
                "bands": bands,
                **flags,
            }
        )

    core = [
        row for row in clean_rows if row["matched_4p_core_family"]
    ]
    best_core = min(core, key=lambda row: row["shape_score"])
    best_overall = min(clean_rows, key=lambda row: row["shape_score"])
    complex_row = next(
        row for row in clean_rows
        if row["family"] == "complex_biquad_4p"
    )
    general_row = next(
        row for row in clean_rows
        if row["family"] == "general_second_order_4p"
    )
    real_row = next(
        row for row in clean_rows
        if row["family"] == "real_2p2z_4p"
    )

    gp = general_row["parameters"]
    general_prefers_real_roots = bool(
        gp["pole_Q"] < 0.5 or gp["zero_Q"] < 0.5
    )
    any_pass = any(
        row["passes_screen_tolerance"] for row in clean_rows
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
            "Bessel times diagnostic pre-ADC transfer -> first alias fold -> "
            "production post-filter white term; no digital filter"
        ),
        "comparison_design": {
            "matched_4p_core_families": list(CORE_FAMILIES),
            "secondary_complexity_ladder_family": (
                "hybrid_general_plus_leadlag_6p"
            ),
            "no_free_gain": True,
            "all_transfer_DC_gain": 1.0,
            "all_roots_stable_minimum_phase": True,
        },
        "search_box": {
            "center_Hz": [
                float(center_min_hz),
                float(center_max_hz),
            ],
            "complex_biquad_Q": [
                float(complex_q_min),
                float(q_max),
            ],
            "general_second_order_Q": [
                float(general_q_min),
                float(q_max),
            ],
        },
        "screen_tolerance": {
            "rms_residual_dB_max": float(rms_screen_db),
            "max_abs_residual_dB_max": float(max_screen_db),
            "diagnostic_not_physical_prior": True,
        },
        "baseline": {
            "shape_score": baseline_score,
            "residual_metrics": baseline_metrics,
            "bands": baseline_bands,
        },
        "families": clean_rows,
        "best_matched_4p_family": best_core,
        "best_overall_family": best_overall,
        "topology_evidence": {
            "general_second_order_prefers_Q_below_0p5": (
                general_prefers_real_roots
            ),
            "general_score_ratio_to_complex_4p": float(
                general_row["shape_score"]
                / complex_row["shape_score"]
            ),
            "real_2p2z_score_ratio_to_complex_4p": float(
                real_row["shape_score"]
                / complex_row["shape_score"]
            ),
            "general_outperforms_complex_4p": bool(
                general_row["shape_score"]
                < complex_row["shape_score"]
            ),
            "real_2p2z_outperforms_complex_4p": bool(
                real_row["shape_score"]
                < complex_row["shape_score"]
            ),
            "best_4p_is_nonresonant_or_overdamped": bool(
                best_core["family"]
                in {
                    "general_second_order_4p",
                    "real_2p2z_4p",
                }
            ),
            "any_family_passes_screen_tolerance": bool(any_pass),
        },
        "curve_sample": curve_sample(
            frequency,
            target,
            baseline,
            internal_rows,
        ),
        "guardrail": (
            "This is a topology falsification with frozen detector physics. "
            "A better overdamped/real-root fit supports an RC/RL/readout-"
            "network interpretation but does not identify a unique circuit "
            "without schematics or independent transfer measurements."
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
    parser.add_argument(
        "--complex-q-min",
        type=float,
        default=COMPLEX_Q_MIN_DEFAULT,
    )
    parser.add_argument(
        "--general-q-min",
        type=float,
        default=GENERAL_Q_MIN_DEFAULT,
    )
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

    if not (
        0.0 < args.general_q_min < 0.5
        < args.complex_q_min < args.q_max
    ):
        raise ValueError(
            "require 0 < general_q_min < 0.5 < complex_q_min < q_max"
        )
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
    fit = base.fit_args(summary)
    frequency = np.geomspace(
        fit.fit_min_hz,
        fit.fit_max_hz,
        fit.fit_points,
    )
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        frequency,
    )
    result = run(
        dict(summary["best_case_parameters"]),
        target_context["target"],
        frequency,
        fit,
        center_min_hz=args.center_min_hz,
        center_max_hz=args.center_max_hz,
        complex_q_min=args.complex_q_min,
        general_q_min=args.general_q_min,
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
        "preanalysis_readout_topology_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    evidence = result["topology_evidence"]
    print(
        json.dumps(
            {
                "output": str(output),
                "accepted_records": result["accepted_records"],
                "baseline_shape_score": result["baseline"][
                    "shape_score"
                ],
                "best_matched_4p_family": result[
                    "best_matched_4p_family"
                ]["family"],
                "best_matched_4p_score": result[
                    "best_matched_4p_family"
                ]["shape_score"],
                "best_overall_family": result[
                    "best_overall_family"
                ]["family"],
                "best_overall_score": result[
                    "best_overall_family"
                ]["shape_score"],
                "general_prefers_q_below_0p5": evidence[
                    "general_second_order_prefers_Q_below_0p5"
                ],
                "general_score_ratio_to_complex": evidence[
                    "general_score_ratio_to_complex_4p"
                ],
                "real_score_ratio_to_complex": evidence[
                    "real_2p2z_score_ratio_to_complex_4p"
                ],
                "any_screen_pass": evidence[
                    "any_family_passes_screen_tolerance"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

"""Profile the hybrid readout lead/lag pole and first-alias identifiability.

This diagnostic freezes detector physics and the production measurement chain,
then fixes the hybrid lead/lag pole on a wide grid while refitting only the
remaining five transfer parameters.  It also evaluates the exact infinite-pole
limit, where the real lead/lag denominator disappears and only its zero remains.

The target is the fresh pre-analysis ASD from the exact accepted CH0 record
mask.  No 10 kHz digital analysis filter is applied.
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
from subScript import preanalysis_readout_topology_diagnostic as topology  # noqa: E402

REFERENCE_HZ = 1_000.0
CENTER_MIN_HZ_DEFAULT = 1_000.0
CENTER_MAX_HZ_DEFAULT = 300_000.0
GENERAL_Q_MIN_DEFAULT = 0.10
Q_MAX_DEFAULT = 20.0
FIXED_POLE_GRID_HZ_DEFAULT = (
    150_000.0,
    250_000.0,
    300_000.0,
    500_000.0,
    1_000_000.0,
    2_000_000.0,
    5_000_000.0,
)
RMS_SCREEN_DB_DEFAULT = 1.0
MAX_SCREEN_DB_DEFAULT = 3.0
ANCHOR_FREQUENCIES_HZ = (
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


def leadlag_magnitude(
    frequency_hz,
    zero_hz,
    pole_hz,
):
    """DC-normalized real lead/lag magnitude; pole_hz=None means pole -> infinity."""

    frequency = np.asarray(frequency_hz, dtype=float)
    zero = float(zero_hz)
    if zero <= 0.0:
        raise ValueError("lead/lag zero must be positive")
    numerator = np.sqrt(1.0 + (frequency / zero) ** 2)
    if pole_hz is None:
        return numerator
    pole = float(pole_hz)
    if pole <= 0.0:
        raise ValueError("lead/lag pole must be positive")
    return numerator / np.sqrt(1.0 + (frequency / pole) ** 2)


def hybrid_transfer_magnitude(
    frequency_hz,
    parameters,
    fixed_leadlag_pole_hz,
):
    return topology.second_order_magnitude(
        frequency_hz,
        parameters["pole_Hz"],
        parameters["pole_Q"],
        parameters["zero_Hz"],
        parameters["zero_Q"],
    ) * leadlag_magnitude(
        frequency_hz,
        parameters["leadlag_zero_Hz"],
        fixed_leadlag_pole_hz,
    )


def decode(vector):
    values = 10.0 ** np.asarray(vector, dtype=float)
    if values.size != 5:
        raise ValueError("hybrid pole-profile vector must have five parameters")
    return {
        "pole_Hz": float(values[0]),
        "pole_Q": float(values[1]),
        "zero_Hz": float(values[2]),
        "zero_Q": float(values[3]),
        "leadlag_zero_Hz": float(values[4]),
    }


def describe_parameters(parameters, fixed_leadlag_pole_hz):
    result = dict(parameters)
    result["leadlag_pole_Hz"] = (
        None
        if fixed_leadlag_pole_hz is None
        else float(fixed_leadlag_pole_hz)
    )
    result["leadlag_pole_is_infinite"] = bool(
        fixed_leadlag_pole_hz is None
    )
    result["pole_root_regime"] = topology.root_regime(
        parameters["pole_Q"]
    )
    result["zero_root_regime"] = topology.root_regime(
        parameters["zero_Q"]
    )
    result["pole_equivalent_real_root_Hz"] = (
        topology.equivalent_real_roots(
            parameters["pole_Hz"],
            parameters["pole_Q"],
        )
    )
    result["zero_equivalent_real_root_Hz"] = (
        topology.equivalent_real_roots(
            parameters["zero_Hz"],
            parameters["zero_Q"],
        )
    )
    return result


def shaped_components(
    context,
    parameters,
    fixed_leadlag_pole_hz,
):
    """Return absolute main, first-alias, and white ASD before normalization."""

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
    extra_main = hybrid_transfer_magnitude(
        frequency,
        parameters,
        fixed_leadlag_pole_hz,
    )
    extra_alias = hybrid_transfer_magnitude(
        alias_frequency,
        parameters,
        fixed_leadlag_pole_hz,
    )
    main = context["main_intrinsic_asd"] * hardware_main * extra_main
    alias = (
        context["alias_intrinsic_asd"]
        * hardware_alias
        * extra_alias
    )
    alias = np.where(context["same_bin"], 0.0, alias)
    white = np.full_like(
        frequency,
        float(context["post_filter_white_asd_A_rtHz"]),
    )
    return {
        "main_asd": main,
        "alias_asd": alias,
        "white_asd": white,
    }


def pre_analysis_model(
    context,
    parameters,
    fixed_leadlag_pole_hz,
):
    components = shaped_components(
        context,
        parameters,
        fixed_leadlag_pole_hz,
    )
    absolute = np.sqrt(
        components["main_asd"] ** 2
        + components["alias_asd"] ** 2
        + components["white_asd"] ** 2
    )
    return opt.normalize_at(
        context["frequency_Hz"],
        absolute,
        reference_hz=REFERENCE_HZ,
    )


def parameter_bounds(
    center_min_hz,
    center_max_hz,
    general_q_min,
    q_max,
):
    lf0 = np.log10(float(center_min_hz))
    lf1 = np.log10(float(center_max_hz))
    lq0 = np.log10(float(general_q_min))
    lq1 = np.log10(float(q_max))
    names = (
        "pole_Hz",
        "pole_Q",
        "zero_Hz",
        "zero_Q",
        "leadlag_zero_Hz",
    )
    bounds = (
        (lf0, lf1),
        (lq0, lq1),
        (lf0, lf1),
        (lq0, lq1),
        (lf0, lf1),
    )
    return names, bounds


def fit_fixed_pole(
    fixed_pole_hz,
    context,
    target,
    frequency,
    args,
    *,
    center_min_hz,
    center_max_hz,
    general_q_min,
    q_max,
    seed,
    de_maxiter,
    rms_screen_db,
    max_screen_db,
):
    names, bounds = parameter_bounds(
        center_min_hz,
        center_max_hz,
        general_q_min,
        q_max,
    )

    def model_from_vector(vector):
        return pre_analysis_model(
            context,
            decode(vector),
            fixed_pole_hz,
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
        max_nfev=3000,
    )
    best_vector = min((de.x, ls.x), key=objective)
    parameters = decode(best_vector)
    model = model_from_vector(best_vector)
    score = float(opt.fit_score(model, target, frequency, args))
    metrics = base.residual_db_metrics(
        model,
        target,
        frequency,
        args,
    )
    bands = base.band_summary(
        model,
        target,
        frequency,
        args,
    )
    positions = topology.named_boundary_positions(
        names,
        best_vector,
        bounds,
    )
    passes = bool(
        metrics["rms_residual_dB"] <= float(rms_screen_db)
        and metrics["max_abs_residual_dB"] <= float(max_screen_db)
    )
    return {
        "fixed_leadlag_pole_Hz": (
            None if fixed_pole_hz is None else float(fixed_pole_hz)
        ),
        "fixed_leadlag_pole_label": (
            "infinity"
            if fixed_pole_hz is None
            else f"{float(fixed_pole_hz):g}"
        ),
        "leadlag_pole_is_infinite": bool(fixed_pole_hz is None),
        "n_refit_parameters": 5,
        "success": bool(de.success or ls.success),
        "evaluations": int(de.nfev + ls.nfev),
        "parameters": describe_parameters(
            parameters,
            fixed_pole_hz,
        ),
        "shape_score": score,
        "residual_metrics": metrics,
        "bands": bands,
        "passes_screen_tolerance": passes,
        "boundary_positions": positions,
        "_model": model,
        "_parameters_raw": parameters,
    }


def baseline_components(context):
    """Components with the production 100 kHz Bessel and no extra transfer."""

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
    main = context["main_intrinsic_asd"] * hardware_main
    alias = context["alias_intrinsic_asd"] * hardware_alias
    alias = np.where(context["same_bin"], 0.0, alias)
    white = np.full_like(
        frequency,
        float(context["post_filter_white_asd_A_rtHz"]),
    )
    return {
        "main_asd": main,
        "alias_asd": alias,
        "white_asd": white,
    }


def component_fraction_sample(
    frequency,
    components,
):
    rows = []
    main_psd = np.asarray(components["main_asd"], dtype=float) ** 2
    alias_psd = np.asarray(components["alias_asd"], dtype=float) ** 2
    white_psd = np.asarray(components["white_asd"], dtype=float) ** 2
    total_psd = main_psd + alias_psd + white_psd
    tiny = np.finfo(float).tiny

    for anchor in ANCHOR_FREQUENCIES_HZ:
        index = int(np.argmin(np.abs(frequency - anchor)))
        total = max(float(total_psd[index]), tiny)
        rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "main_psd_fraction": float(main_psd[index] / total),
                "first_alias_psd_fraction": float(
                    alias_psd[index] / total
                ),
                "post_filter_white_psd_fraction": float(
                    white_psd[index] / total
                ),
                "alias_over_main_asd": float(
                    np.sqrt(
                        alias_psd[index]
                        / max(float(main_psd[index]), tiny)
                    )
                ),
            }
        )
    return rows


def profile_trend(rows):
    finite = [
        row for row in rows
        if not row["leadlag_pole_is_infinite"]
    ]
    finite = sorted(
        finite,
        key=lambda row: row["fixed_leadlag_pole_Hz"],
    )
    scores = np.asarray(
        [row["shape_score"] for row in finite],
        dtype=float,
    )
    poles = np.asarray(
        [row["fixed_leadlag_pole_Hz"] for row in finite],
        dtype=float,
    )
    best_finite = min(finite, key=lambda row: row["shape_score"])
    infinite = next(
        row for row in rows if row["leadlag_pole_is_infinite"]
    )

    if len(scores) >= 2:
        high_mask = poles >= 300_000.0
        high_scores = scores[high_mask]
        monotonic_high = bool(
            len(high_scores) >= 2
            and np.all(np.diff(high_scores) <= 1e-10)
        )
    else:
        monotonic_high = False

    return {
        "best_finite_pole_Hz": float(
            best_finite["fixed_leadlag_pole_Hz"]
        ),
        "best_finite_shape_score": float(
            best_finite["shape_score"]
        ),
        "infinite_pole_shape_score": float(infinite["shape_score"]),
        "infinite_over_best_finite_score_ratio": float(
            infinite["shape_score"] / best_finite["shape_score"]
        ),
        "score_nonincreasing_from_300k_to_5M": monotonic_high,
        "infinite_is_best_profile_point": bool(
            infinite["shape_score"]
            <= min(row["shape_score"] for row in rows) + 1e-12
        ),
        "finite_minimum_below_5MHz": bool(
            best_finite["fixed_leadlag_pole_Hz"]
            < max(poles)
            and best_finite["shape_score"]
            < finite[-1]["shape_score"] - 1e-8
        ),
        "interpretation_flags": {
            "supports_pole_to_infinity_nonidentifiability": bool(
                monotonic_high
                and infinite["shape_score"]
                <= finite[-1]["shape_score"] * 1.01
            ),
            "supports_finite_high_frequency_corner": bool(
                best_finite["fixed_leadlag_pole_Hz"]
                < max(poles)
                and best_finite["shape_score"]
                < finite[-1]["shape_score"] * 0.98
            ),
        },
    }


def run(
    candidate,
    target,
    frequency,
    args,
    *,
    pole_grid_hz=FIXED_POLE_GRID_HZ_DEFAULT,
    center_min_hz=CENTER_MIN_HZ_DEFAULT,
    center_max_hz=CENTER_MAX_HZ_DEFAULT,
    general_q_min=GENERAL_Q_MIN_DEFAULT,
    q_max=Q_MAX_DEFAULT,
    seed=20260915,
    de_maxiter=100,
    rms_screen_db=RMS_SCREEN_DB_DEFAULT,
    max_screen_db=MAX_SCREEN_DB_DEFAULT,
):
    context = base.intrinsic_context(candidate, frequency)

    base_components = baseline_components(context)
    baseline_absolute = np.sqrt(
        base_components["main_asd"] ** 2
        + base_components["alias_asd"] ** 2
        + base_components["white_asd"] ** 2
    )
    baseline_model = opt.normalize_at(
        frequency,
        baseline_absolute,
        reference_hz=REFERENCE_HZ,
    )
    baseline_score = float(
        opt.fit_score(
            baseline_model,
            target,
            frequency,
            args,
        )
    )

    rows_internal = []
    pole_values = [float(value) for value in pole_grid_hz]
    for index, fixed_pole in enumerate(pole_values):
        rows_internal.append(
            fit_fixed_pole(
                fixed_pole,
                context,
                target,
                frequency,
                args,
                center_min_hz=center_min_hz,
                center_max_hz=center_max_hz,
                general_q_min=general_q_min,
                q_max=q_max,
                seed=seed + index,
                de_maxiter=de_maxiter,
                rms_screen_db=rms_screen_db,
                max_screen_db=max_screen_db,
            )
        )
    rows_internal.append(
        fit_fixed_pole(
            None,
            context,
            target,
            frequency,
            args,
            center_min_hz=center_min_hz,
            center_max_hz=center_max_hz,
            general_q_min=general_q_min,
            q_max=q_max,
            seed=seed + len(pole_values),
            de_maxiter=de_maxiter,
            rms_screen_db=rms_screen_db,
            max_screen_db=max_screen_db,
        )
    )

    for row in rows_internal:
        row["score_ratio_to_baseline"] = float(
            row["shape_score"] / baseline_score
        )

    best = min(rows_internal, key=lambda row: row["shape_score"])
    best_components = shaped_components(
        context,
        best["_parameters_raw"],
        best["fixed_leadlag_pole_Hz"],
    )

    clean_rows = [
        {
            key: value
            for key, value in row.items()
            if not key.startswith("_")
        }
        for row in rows_internal
    ]
    best_clean = next(
        row for row in clean_rows
        if row["fixed_leadlag_pole_label"]
        == best["fixed_leadlag_pole_label"]
    )

    trend = profile_trend(clean_rows)
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
            "Bessel -> hybrid general-second-order times fixed-pole real "
            "lead/lag -> first alias fold -> production post-filter white"
        ),
        "profile_definition": (
            "Fix the hybrid lead/lag pole at each grid point, refit only "
            "second-order pole/Q, second-order zero/Q, and lead/lag zero. "
            "The infinity point is the exact limit with no lead/lag "
            "denominator factor."
        ),
        "fixed_leadlag_pole_grid_Hz": pole_values,
        "includes_exact_infinite_pole_limit": True,
        "refit_parameters": [
            "pole_Hz",
            "pole_Q",
            "zero_Hz",
            "zero_Q",
            "leadlag_zero_Hz",
        ],
        "search_box": {
            "center_Hz": [
                float(center_min_hz),
                float(center_max_hz),
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
            "bands": base.band_summary(
                baseline_model,
                target,
                frequency,
                args,
            ),
            "alias_fraction_sample": component_fraction_sample(
                frequency,
                base_components,
            ),
        },
        "profile_rows": clean_rows,
        "best_profile_row": best_clean,
        "profile_trend": trend,
        "best_profile_alias_fraction_sample": (
            component_fraction_sample(
                frequency,
                best_components,
            )
        ),
        "guardrail": (
            "A pole that runs toward infinity is a non-identifying shape "
            "limit, not evidence for a multi-MHz physical hardware corner. "
            "Alias fractions are model bookkeeping for this first-alias "
            "approximation and do not replace a measured hardware transfer."
        ),
    }


def parse_pole_grid(text):
    values = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        value = float(token)
        if value <= 0.0:
            raise ValueError("all pole-grid values must be positive")
        values.append(value)
    if not values:
        raise ValueError("pole grid must contain at least one value")
    return tuple(values)


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
        "--pole-grid-hz",
        default=",".join(
            f"{value:g}" for value in FIXED_POLE_GRID_HZ_DEFAULT
        ),
    )
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
        "--general-q-min",
        type=float,
        default=GENERAL_Q_MIN_DEFAULT,
    )
    parser.add_argument("--q-max", type=float, default=Q_MAX_DEFAULT)
    parser.add_argument("--de-maxiter", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260915)
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

    pole_grid = parse_pole_grid(args.pole_grid_hz)
    if not (0.0 < args.general_q_min < args.q_max):
        raise ValueError("require 0 < general_q_min < q_max")
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
        pole_grid_hz=pole_grid,
        center_min_hz=args.center_min_hz,
        center_max_hz=args.center_max_hz,
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
        "preanalysis_hybrid_pole_profile_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    trend = result["profile_trend"]
    print(
        json.dumps(
            {
                "output": str(output),
                "accepted_records": result["accepted_records"],
                "baseline_shape_score": result["baseline"][
                    "shape_score"
                ],
                "best_pole_label": result["best_profile_row"][
                    "fixed_leadlag_pole_label"
                ],
                "best_shape_score": result["best_profile_row"][
                    "shape_score"
                ],
                "best_finite_pole_Hz": trend[
                    "best_finite_pole_Hz"
                ],
                "infinite_pole_shape_score": trend[
                    "infinite_pole_shape_score"
                ],
                "supports_pole_to_infinity": trend[
                    "interpretation_flags"
                ]["supports_pole_to_infinity_nonidentifiability"],
                "supports_finite_corner": trend[
                    "interpretation_flags"
                ]["supports_finite_high_frequency_corner"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

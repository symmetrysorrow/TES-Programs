"""Fit a nested effective-numerator complexity ladder on the 1-40 kHz repeat target.

The detector model and profiled post-filter white floor are frozen.  Each family
uses the same second-order pole denominator and increases only the numerator
magnitude-squared complexity:

  order 0: P(x) = 1
  order 1: P(x) = 1 + v^2 x^2
  order 2: P(x) = (1 + u x^2)^2 + v^2 x^2
  order 3: P(x) = [(1 + u x^2)^2 + v^2 x^2] (1 + w^2 x^2)

with x = f / f_ref.

Order 3 is algebraically equivalent to the current infinity-lead-pole hybrid
numerator, but the output reports polynomial coefficients rather than physical
zero/lead-zero corner interpretations.  Orders 0-2 remove numerator degrees of
freedom in a strictly nested way.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, least_squares

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "readout_effective_numerator_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_effective_numerator_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import preanalysis_hybrid_pole_profile_diagnostic as profile  # noqa: E402
from subScript import readout_lowmid_holdout_diagnostic as holdout  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import readout_white_floor_separation_diagnostic as white_sep  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


def effective_numerator_squared(frequency_hz, order, latent, scale_hz):
    frequency = np.asarray(frequency_hz, dtype=float)
    x2 = (frequency / float(scale_hz)) ** 2
    order = int(order)
    if order == 0:
        return np.ones_like(frequency)
    if order == 1:
        v = float(latent["v"])
        return 1.0 + v**2 * x2
    if order in {2, 3}:
        u = float(latent["u"])
        v = float(latent["v"])
        result = (1.0 + u * x2) ** 2 + v**2 * x2
        if order == 3:
            w = float(latent["w"])
            result = result * (1.0 + w**2 * x2)
        return result
    raise ValueError(f"unsupported effective numerator order: {order}")


def polynomial_coefficients(order, latent):
    order = int(order)
    if order == 0:
        return {"c2": 0.0, "c4": 0.0, "c6": 0.0}
    if order == 1:
        v = float(latent["v"])
        return {"c2": v**2, "c4": 0.0, "c6": 0.0}
    u = float(latent["u"])
    v = float(latent["v"])
    base_c2 = 2.0 * u + v**2
    base_c4 = u**2
    if order == 2:
        return {
            "c2": float(base_c2),
            "c4": float(base_c4),
            "c6": 0.0,
        }
    if order == 3:
        w2 = float(latent["w"]) ** 2
        return {
            "c2": float(base_c2 + w2),
            "c4": float(base_c4 + base_c2 * w2),
            "c6": float(base_c4 * w2),
        }
    raise ValueError(f"unsupported effective numerator order: {order}")


def effective_transfer_magnitude(
    frequency_hz,
    pole_hz,
    pole_q,
    order,
    latent,
    scale_hz,
):
    frequency = np.asarray(frequency_hz, dtype=float)
    pole_hz = float(pole_hz)
    pole_q = float(pole_q)
    if pole_hz <= 0.0 or pole_q <= 0.0:
        raise ValueError("pole frequency and Q must be positive")
    numerator_squared = effective_numerator_squared(
        frequency,
        order,
        latent,
        scale_hz,
    )
    if np.any(~np.isfinite(numerator_squared)) or np.any(
        numerator_squared <= 0.0
    ):
        raise ValueError("effective numerator magnitude squared must be positive")
    xp = frequency / pole_hz
    denominator_squared = (
        (1.0 - xp**2) ** 2 + (xp / pole_q) ** 2
    )
    return np.sqrt(numerator_squared / denominator_squared)


def pre_analysis_model(context, parameters, order, scale_hz):
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
    main_transfer = effective_transfer_magnitude(
        frequency,
        parameters["pole_Hz"],
        parameters["pole_Q"],
        order,
        parameters.get("latent", {}),
        scale_hz,
    )
    alias_transfer = effective_transfer_magnitude(
        alias_frequency,
        parameters["pole_Hz"],
        parameters["pole_Q"],
        order,
        parameters.get("latent", {}),
        scale_hz,
    )
    main = (
        context["main_intrinsic_asd"]
        * hardware_main
        * main_transfer
    )
    alias = (
        context["alias_intrinsic_asd"]
        * hardware_alias
        * alias_transfer
    )
    alias = np.where(context["same_bin"], 0.0, alias)
    white = np.full_like(
        frequency,
        float(context["post_filter_white_asd_A_rtHz"]),
    )
    absolute = np.sqrt(main**2 + alias**2 + white**2)
    return opt.normalize_at(
        frequency,
        absolute,
        reference_hz=profile.REFERENCE_HZ,
    )


def hybrid_to_effective(parameters, scale_hz):
    scale = float(scale_hz)
    zero_hz = float(parameters["zero_Hz"])
    zero_q = float(parameters["zero_Q"])
    lead_zero_hz = float(parameters["leadlag_zero_Hz"])
    if zero_hz <= 0.0 or zero_q <= 0.0 or lead_zero_hz <= 0.0:
        raise ValueError("hybrid zero parameters must be positive")
    r = (scale / zero_hz) ** 2
    return {
        "pole_Hz": float(parameters["pole_Hz"]),
        "pole_Q": float(parameters["pole_Q"]),
        "latent": {
            "u": float(-r),
            "v": float(np.sqrt(r) / zero_q),
            "w": float(scale / lead_zero_hz),
        },
    }


def vector_spec(order, optimizer):
    order = int(order)
    names = ["log10_pole_Hz", "log10_pole_Q"]
    bounds = [
        (
            np.log10(float(optimizer["pole_min_Hz"])),
            np.log10(float(optimizer["pole_max_Hz"])),
        ),
        (
            np.log10(float(optimizer["pole_Q_min"])),
            np.log10(float(optimizer["pole_Q_max"])),
        ),
    ]
    if order == 1:
        names += ["log10_v"]
        bounds += [
            (
                np.log10(float(optimizer["v_min"])),
                np.log10(float(optimizer["v_max"])),
            )
        ]
    elif order in {2, 3}:
        names += ["u", "log10_v"]
        bounds += [
            (
                float(optimizer["u_min"]),
                float(optimizer["u_max"]),
            ),
            (
                np.log10(float(optimizer["v_min"])),
                np.log10(float(optimizer["v_max"])),
            ),
        ]
        if order == 3:
            names += ["log10_w"]
            bounds += [
                (
                    np.log10(float(optimizer["w_min"])),
                    np.log10(float(optimizer["w_max"])),
                )
            ]
    elif order != 0:
        raise ValueError(f"unsupported order: {order}")
    return tuple(names), tuple(bounds)


def decode_vector(order, vector):
    vector = np.asarray(vector, dtype=float)
    order = int(order)
    expected = 2 + order
    if vector.size != expected:
        raise ValueError(
            f"order {order} expects {expected} fit parameters, "
            f"received {vector.size}"
        )
    result = {
        "pole_Hz": float(10.0 ** vector[0]),
        "pole_Q": float(10.0 ** vector[1]),
        "latent": {},
    }
    if order == 1:
        result["latent"]["v"] = float(10.0 ** vector[2])
    elif order in {2, 3}:
        result["latent"]["u"] = float(vector[2])
        result["latent"]["v"] = float(10.0 ** vector[3])
        if order == 3:
            result["latent"]["w"] = float(10.0 ** vector[4])
    return result


def warm_vector(order, hybrid_parameters, scale_hz):
    mapped = hybrid_to_effective(hybrid_parameters, scale_hz)
    pole_hz = mapped["pole_Hz"]
    pole_q = mapped["pole_Q"]
    latent = mapped["latent"]
    coeff = polynomial_coefficients(3, latent)
    values = [np.log10(pole_hz), np.log10(pole_q)]
    order = int(order)
    if order == 1:
        v = np.sqrt(max(float(coeff["c2"]), 1.0e-12))
        values.append(np.log10(v))
    elif order in {2, 3}:
        values += [float(latent["u"]), np.log10(float(latent["v"]))]
        if order == 3:
            values.append(np.log10(float(latent["w"])))
    return np.asarray(values, dtype=float)


def vector_within_bounds(vector, bounds):
    vector = np.asarray(vector, dtype=float)
    return bool(
        vector.size == len(bounds)
        and all(
            lower <= value <= upper
            for value, (lower, upper) in zip(vector, bounds)
        )
    )


def fit_order(
    order,
    context,
    target,
    frequency,
    args,
    *,
    scale_hz,
    optimizer,
    shared_parameters,
    free_parameters,
    seed,
):
    names, bounds = vector_spec(order, optimizer)

    def model_from_vector(vector):
        return pre_analysis_model(
            context,
            decode_vector(order, vector),
            order,
            scale_hz,
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
        maxiter=int(optimizer["DE_maxiter"]),
        popsize=10,
        tol=1e-8,
        polish=False,
        workers=1,
        updating="immediate",
    )
    candidates = [("de", np.asarray(de.x, dtype=float))]
    evaluations = int(de.nfev)
    success = bool(de.success)

    lower = np.asarray([item[0] for item in bounds], dtype=float)
    upper = np.asarray([item[1] for item in bounds], dtype=float)

    def residual_vector(vector):
        return opt.weighted_residual_vector(
            model_from_vector(vector),
            target,
            frequency,
            args,
        )

    ls = least_squares(
        residual_vector,
        de.x,
        bounds=(lower, upper),
        max_nfev=3000,
    )
    evaluations += int(ls.nfev)
    success = bool(success or ls.success)
    candidates.append(
        ("de_least_squares", np.asarray(ls.x, dtype=float))
    )

    for label, parameters in (
        ("shared", shared_parameters),
        ("free_reference", free_parameters),
    ):
        warm = warm_vector(order, parameters, scale_hz)
        if not vector_within_bounds(warm, bounds):
            continue
        candidates.append((f"{label}_exact", warm.copy()))
        warm_ls = least_squares(
            residual_vector,
            warm,
            bounds=(lower, upper),
            max_nfev=3000,
        )
        evaluations += int(warm_ls.nfev)
        success = bool(success or warm_ls.success)
        candidates.append(
            (
                f"{label}_least_squares",
                np.asarray(warm_ls.x, dtype=float),
            )
        )

    best_source, best_vector = min(
        candidates,
        key=lambda item: objective(item[1]),
    )
    parameters = decode_vector(order, best_vector)
    model = model_from_vector(best_vector)
    score = objective(best_vector)
    residual = base.residual_db_metrics(
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
    return {
        "order": int(order),
        "n_total_fit_parameters": int(2 + int(order)),
        "n_numerator_parameters": int(order),
        "vector_names": list(names),
        "success": bool(success),
        "evaluations": int(evaluations),
        "best_candidate_source": best_source,
        "parameters": {
            "pole_Hz": parameters["pole_Hz"],
            "pole_Q": parameters["pole_Q"],
            "latent": parameters["latent"],
            "polynomial_coefficients_dimensionless": (
                polynomial_coefficients(
                    order,
                    parameters["latent"],
                )
            ),
        },
        "shape_score": float(score),
        "residual_metrics": residual,
        "bands": bands,
        "_model": model,
        "_parameters_raw": parameters,
    }


def annotate_complexity(rows, free_score, free_rms, screen):
    result = []
    for row in rows:
        clean = {
            key: value
            for key, value in row.items()
            if not key.startswith("_")
        }
        score_ratio = float(row["shape_score"] / free_score)
        rms_delta = float(
            row["residual_metrics"]["rms_residual_dB"]
            - free_rms
        )
        clean["score_ratio_to_free_hybrid"] = score_ratio
        clean["rms_delta_from_free_hybrid_dB"] = rms_delta
        clean["near_free_hybrid"] = bool(
            score_ratio <= float(screen["near_free_score_ratio"])
            and rms_delta
            <= float(screen["near_free_rms_delta_dB"])
        )
        result.append(clean)
    return result


def run(config, config_path: Path):
    manifest_path = ident.resolve_config_path(
        config["manifest"],
        config_path,
    )
    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    cases = shared.normalize_manifest(manifest, manifest_path)
    repeat_case = ident.case_by_label(
        cases,
        config["repeat_case_label"],
    )
    summary = json.loads(
        repeat_case["summary"].read_text(encoding="utf-8")
    )
    comparison, experiment_path, comparison_source = (
        ident.comparison_for_case(repeat_case)
    )

    full_args = base.fit_args(summary)
    full_frequency = np.geomspace(
        full_args.fit_min_hz,
        full_args.fit_max_hz,
        full_args.fit_points,
    )
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        full_frequency,
    )
    full_target = target_context["target"]
    candidate = dict(summary["best_case_parameters"])
    base_context = base.intrinsic_context(
        candidate,
        full_frequency,
    )

    reference_path = ident.resolve_config_path(
        config["reference_transfer"],
        config_path,
    )
    reference_payload = json.loads(
        reference_path.read_text(encoding="utf-8")
    )
    reference = shared.load_shared_transfer(reference_payload)

    free_snapshot_path = ident.resolve_config_path(
        config["lowmid_free_snapshot"],
        config_path,
    )
    free_snapshot = ident.load_free_snapshot(
        free_snapshot_path
    )
    full_context = white_sep.white_scaled_context(
        base_context,
        free_snapshot["white_scale"],
    )
    if not np.isclose(
        float(full_context["post_filter_white_asd_A_rtHz"]),
        free_snapshot["white_asd_A_rtHz"],
        rtol=5.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "tracked low/mid snapshot white ASD does not reproduce"
        )

    fit_cfg = config["fit_region_Hz"]
    fit_min = float(fit_cfg["min"])
    fit_max = float(fit_cfg["max_exclusive"])
    fit_mask = (
        (full_frequency >= fit_min)
        & (full_frequency < fit_max)
    )
    hold_cfg = config["holdout_region_Hz"]
    hold_min = float(hold_cfg["min"])
    hold_max = float(hold_cfg["max"])
    hold_mask = (
        (full_frequency >= hold_min)
        & (full_frequency <= hold_max)
    )
    if np.any(fit_mask & hold_mask):
        raise RuntimeError("fit and holdout masks overlap")

    fit_frequency = full_frequency[fit_mask]
    fit_target = full_target[fit_mask]
    fit_context = holdout.subset_context(
        full_context,
        fit_mask,
    )
    fit_args = ident.band_args(
        full_args,
        fit_min,
        fit_max,
        np.count_nonzero(fit_mask),
    )
    hold_frequency = full_frequency[hold_mask]
    hold_target = full_target[hold_mask]
    hold_args = ident.band_args(
        full_args,
        hold_min,
        hold_max,
        np.count_nonzero(hold_mask),
    )

    scale_hz = float(
        config["effective_numerator"]["reference_scale_Hz"]
    )
    free_hybrid_fit_model = profile.pre_analysis_model(
        fit_context,
        free_snapshot["parameters"],
        None,
    )
    free_score = float(
        opt.fit_score(
            free_hybrid_fit_model,
            fit_target,
            fit_frequency,
            fit_args,
        )
    )
    free_residual = base.residual_db_metrics(
        free_hybrid_fit_model,
        fit_target,
        fit_frequency,
        fit_args,
    )
    if not np.isclose(
        free_score,
        free_snapshot["shape_score"],
        rtol=2e-6,
        atol=1e-12,
    ):
        raise ValueError("free hybrid reference score does not reproduce")

    mapped_free = hybrid_to_effective(
        free_snapshot["parameters"],
        scale_hz,
    )
    mapped_model = pre_analysis_model(
        fit_context,
        mapped_free,
        3,
        scale_hz,
    )
    equivalence_max_abs = float(
        np.max(np.abs(mapped_model - free_hybrid_fit_model))
    )
    equivalence_rms = float(
        np.sqrt(
            np.mean(
                (mapped_model - free_hybrid_fit_model) ** 2
            )
        )
    )
    if equivalence_max_abs > 1.0e-10:
        raise RuntimeError(
            "order-3 effective numerator failed algebraic hybrid equivalence"
        )

    optimizer = config["optimizer"]
    rows = []
    raw_rows = []
    for index, order in enumerate(
        config["effective_numerator"]["orders"]
    ):
        result = fit_order(
            int(order),
            fit_context,
            fit_target,
            fit_frequency,
            fit_args,
            scale_hz=scale_hz,
            optimizer=optimizer,
            shared_parameters=reference["parameters"],
            free_parameters=free_snapshot["parameters"],
            seed=int(optimizer["seed"]) + index,
        )
        raw_rows.append(result)

        full_model = pre_analysis_model(
            full_context,
            result["_parameters_raw"],
            int(order),
            scale_hz,
        )
        hold_metrics = holdout.model_metrics(
            full_model[hold_mask],
            hold_target,
            hold_frequency,
            hold_args,
        )
        full_metrics = holdout.model_metrics(
            full_model,
            full_target,
            full_frequency,
            full_args,
        )
        result["holdout_metrics"] = hold_metrics
        result["full_1_200k_metrics"] = full_metrics

    rows = annotate_complexity(
        raw_rows,
        free_score,
        free_residual["rms_residual_dB"],
        config["minimal_order_screen"],
    )
    near_orders = [
        int(row["order"])
        for row in rows
        if row["near_free_hybrid"]
    ]
    minimal_order = min(near_orders) if near_orders else None

    incremental = []
    for previous, current in zip(rows[:-1], rows[1:]):
        incremental.append(
            {
                "from_order": int(previous["order"]),
                "to_order": int(current["order"]),
                "shape_score_ratio_new_over_old": float(
                    current["shape_score"]
                    / previous["shape_score"]
                ),
                "rms_improvement_dB": float(
                    previous["residual_metrics"][
                        "rms_residual_dB"
                    ]
                    - current["residual_metrics"][
                        "rms_residual_dB"
                    ]
                ),
                "holdout_score_ratio_new_over_old": float(
                    current["holdout_metrics"]["shape_score"]
                    / previous["holdout_metrics"]["shape_score"]
                ),
            }
        )

    order3 = next(
        row for row in rows
        if int(row["order"]) == 3
    )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "fit_region": {
            "min_Hz": fit_min,
            "max_Hz_exclusive": fit_max,
            "points": int(np.count_nonzero(fit_mask)),
        },
        "holdout_region": {
            "min_Hz": hold_min,
            "max_Hz": hold_max,
            "points": int(np.count_nonzero(hold_mask)),
            "optimizer_received_holdout_points": False,
        },
        "repeat_case": {
            "label": repeat_case["label"],
            "comparison_source": comparison_source,
            "experiment_path": str(experiment_path),
            "accepted_records": int(
                target_context["accepted_records"]
            ),
        },
        "white_floor_fixed": {
            "scale": free_snapshot["white_scale"],
            "white_asd_A_rtHz": free_snapshot[
                "white_asd_A_rtHz"
            ],
        },
        "effective_numerator_definition": {
            "reference_scale_Hz": scale_hz,
            "x_definition": "x = f / reference_scale_Hz",
            "polynomial_definition": (
                "P(x)=1+c2*x^2+c4*x^4+c6*x^6; "
                "transfer magnitude numerator=sqrt(P)"
            ),
            "orders": config["effective_numerator"]["semantics"],
        },
        "free_hybrid_reference": {
            "parameters": free_snapshot["parameters"],
            "shape_score": free_score,
            "rms_residual_dB": free_residual[
                "rms_residual_dB"
            ],
            "mapped_order3_latent": mapped_free["latent"],
            "mapped_order3_polynomial_coefficients_dimensionless": (
                polynomial_coefficients(
                    3,
                    mapped_free["latent"],
                )
            ),
        },
        "order3_algebraic_equivalence_check": {
            "max_abs_normalized_model_difference": (
                equivalence_max_abs
            ),
            "rms_normalized_model_difference": equivalence_rms,
            "passes_1e-10": bool(
                equivalence_max_abs <= 1.0e-10
            ),
        },
        "complexity_ladder": rows,
        "incremental_complexity_gain": incremental,
        "minimal_order_screen": config["minimal_order_screen"],
        "minimal_order_near_free_hybrid": minimal_order,
        "numerator_dof_required_to_match_free": minimal_order,
        "interpretation_flags": {
            "pole_only_sufficient": bool(
                minimal_order == 0
            ),
            "one_numerator_degree_sufficient": bool(
                minimal_order is not None
                and minimal_order <= 1
            ),
            "two_numerator_degrees_sufficient": bool(
                minimal_order is not None
                and minimal_order <= 2
            ),
            "three_numerator_degrees_required": bool(
                minimal_order == 3
            ),
            "order3_reproduces_free_hybrid": bool(
                order3["near_free_hybrid"]
                and equivalence_max_abs <= 1.0e-10
            ),
            "individual_physical_zero_corners_identified": False,
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "manifest": str(manifest_path),
            "reference_transfer": str(reference_path),
            "lowmid_free_snapshot": str(free_snapshot_path),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    config = json.loads(
        args.config.read_text(encoding="utf-8")
    )
    result = run(config, args.config)
    output = args.output or DEFAULT_OUTPUT
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(output),
                "minimal_order_near_free_hybrid": result[
                    "minimal_order_near_free_hybrid"
                ],
                "complexity_ladder": [
                    {
                        "order": row["order"],
                        "n_total_fit_parameters": row[
                            "n_total_fit_parameters"
                        ],
                        "shape_score": row["shape_score"],
                        "rms_residual_dB": row[
                            "residual_metrics"
                        ]["rms_residual_dB"],
                        "score_ratio_to_free_hybrid": row[
                            "score_ratio_to_free_hybrid"
                        ],
                        "rms_delta_from_free_hybrid_dB": row[
                            "rms_delta_from_free_hybrid_dB"
                        ],
                        "near_free_hybrid": row[
                            "near_free_hybrid"
                        ],
                        "holdout_shape_score": row[
                            "holdout_metrics"
                        ]["shape_score"],
                    }
                    for row in result["complexity_ladder"]
                ],
                "flags": result["interpretation_flags"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

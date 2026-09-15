"""Compete residual order-2 readout degrees of freedom after detector/DC nuisance.

The 2024-12-05 repeat is fit on 1-40 kHz after fixing R_TES to the best
stability-aware profile point (1.25 times the inherited resistance).  For every
readout family, alpha, beta, C_tes, L, and T_bath are re-optimized.  The
post-filter white ASD remains fixed.

The reference readout is the 2024-12-06 order-2 effective transfer

    P(x) = 1 + c2 x^2 + c4 x^4,  x = f / 40 kHz,

with a second-order pole denominator.  Selected canonical readout coordinates
are released in nested/branch families to identify the minimum residual
magnitude-shape freedom required after the tested detector/DC confounds.

The 40-200 kHz region is a strict holdout.
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
DEFAULT_CONFIG = CONFIG_DIR / "readout_residual_dof_competition_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_residual_dof_competition_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import readout_detector_state_competition_diagnostic as competition  # noqa: E402
from subScript import readout_effective_numerator_diagnostic as effective  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import readout_lowmid_holdout_diagnostic as holdout  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


READOUT_PARAMETER_NAMES = ("pole_Hz", "pole_Q", "c2", "c4")
LOG_READOUT_PARAMETERS = {"pole_Hz", "pole_Q"}
MAIN_NESTED_LADDER = (
    "fixed_reference_readout",
    "pole_frequency_only",
    "pole_section",
    "pole_section_plus_c2",
    "full_order2",
)


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def load_baseline_snapshot(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    baseline = payload["best_profiled_detector_dc_baseline"]
    reference = payload["fixed_reference_order2_readout"]
    local = payload["repeat_local_order2_readout"]
    return {
        "R_TES_Ohm": float(baseline["R_TES_Ohm"]),
        "R_ratio_to_inherited": float(
            baseline["R_ratio_to_inherited"]
        ),
        "detector_candidate": {
            key: float(value)
            for key, value in baseline["candidate"].items()
        },
        "baseline_shape_score": float(baseline["shape_score"]),
        "baseline_rms_dB": float(baseline["rms_residual_dB"]),
        "baseline_holdout_shape_score": float(
            baseline["holdout_shape_score"]
        ),
        "baseline_holdout_rms_dB": float(
            baseline["holdout_rms_residual_dB"]
        ),
        "repeat_local_shape_score": float(
            payload["repeat_local_readout_comparator"]["shape_score"]
        ),
        "repeat_local_rms_dB": float(
            payload["repeat_local_readout_comparator"][
                "rms_residual_dB"
            ]
        ),
        "reference_readout": {
            name: float(reference[name])
            for name in READOUT_PARAMETER_NAMES
        },
        "local_readout": {
            name: float(local[name])
            for name in READOUT_PARAMETER_NAMES
        },
        "reference_scale_Hz": float(
            reference["reference_scale_Hz"]
        ),
        "white_asd_A_rtHz": float(
            payload["fixed_repeat_white_asd_A_rtHz"]
        ),
        "provenance": payload.get("provenance", {}),
    }


def canonical_to_effective(readout):
    pole_hz = float(readout["pole_Hz"])
    pole_q = float(readout["pole_Q"])
    c2 = float(readout["c2"])
    c4 = float(readout["c4"])
    if pole_hz <= 0.0 or pole_q <= 0.0:
        raise ValueError("pole_Hz and pole_Q must be positive")
    if c2 < 0.0 or c4 < 0.0:
        raise ValueError("canonical c2/c4 must be non-negative")
    sqrt_c4 = float(np.sqrt(c4))
    u = -sqrt_c4
    v2 = c2 + 2.0 * sqrt_c4
    if v2 < 0.0:
        raise ValueError("canonical coefficients are not representable")
    return {
        "pole_Hz": pole_hz,
        "pole_Q": pole_q,
        "latent": {
            "u": u,
            "v": float(np.sqrt(v2)),
        },
    }


def canonical_polynomial(frequency_hz, readout, scale_hz):
    frequency = np.asarray(frequency_hz, dtype=float)
    x2 = (frequency / float(scale_hz)) ** 2
    return (
        1.0
        + float(readout["c2"]) * x2
        + float(readout["c4"]) * x2**2
    )


def model_for_candidate(
    detector_candidate,
    readout,
    frequency,
    scale_hz,
    white_asd,
):
    point = opt.tes_operating_point(detector_candidate)
    if not point.get("valid") or not point.get("stable"):
        return None, point
    context = competition.fixed_white_context(
        detector_candidate,
        frequency,
        white_asd,
    )
    model = effective.pre_analysis_model(
        context,
        canonical_to_effective(readout),
        2,
        scale_hz,
    )
    if np.any(~np.isfinite(model)) or np.any(model <= 0.0):
        return None, point
    return model, point


def encode_readout_value(name, value):
    value = float(value)
    if name in LOG_READOUT_PARAMETERS:
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
        return float(np.log10(value))
    return value


def decode_readout_value(name, value):
    value = float(value)
    if name in LOG_READOUT_PARAMETERS:
        return float(10.0**value)
    return value


def joint_vector_bounds(
    detector_names,
    detector_bounds,
    readout_names,
    readout_bounds,
):
    bounds = []
    for name in detector_names:
        lower, upper = detector_bounds[name]
        bounds.append(
            (
                competition.encode_value(name, lower),
                competition.encode_value(name, upper),
            )
        )
    for name in readout_names:
        lower, upper = readout_bounds[name]
        bounds.append(
            (
                encode_readout_value(name, lower),
                encode_readout_value(name, upper),
            )
        )
    return tuple(bounds)


def encode_joint(
    detector_names,
    detector_candidate,
    readout_names,
    readout,
):
    values = [
        competition.encode_value(
            name,
            detector_candidate[name],
        )
        for name in detector_names
    ]
    values += [
        encode_readout_value(name, readout[name])
        for name in readout_names
    ]
    return np.asarray(values, dtype=float)


def decode_joint(
    vector,
    detector_names,
    baseline_detector,
    readout_names,
    baseline_readout,
):
    vector = np.asarray(vector, dtype=float)
    detector = dict(baseline_detector)
    readout = dict(baseline_readout)
    offset = 0
    for name in detector_names:
        detector[name] = competition.decode_value(
            name,
            vector[offset],
        )
        offset += 1
    for name in readout_names:
        readout[name] = decode_readout_value(
            name,
            vector[offset],
        )
        offset += 1
    if offset != vector.size:
        raise ValueError("joint vector length mismatch")
    return detector, readout


def boundary_hits(values, names, bounds, encoder):
    result = {}
    for name in names:
        value = float(values[name])
        lower, upper = bounds[name]
        encoded_value = encoder(name, value)
        encoded_lower = encoder(name, lower)
        encoded_upper = encoder(name, upper)
        span = max(
            abs(encoded_upper - encoded_lower),
            abs(encoded_upper),
            1.0e-12,
        )
        tolerance = span * 1.0e-5
        result[name] = {
            "value": value,
            "lower": float(lower),
            "upper": float(upper),
            "at_lower": bool(
                abs(encoded_value - encoded_lower)
                <= tolerance
            ),
            "at_upper": bool(
                abs(encoded_value - encoded_upper)
                <= tolerance
            ),
        }
    return result


def readout_distance(reference, local, fitted):
    rows = {}
    for name in READOUT_PARAMETER_NAMES:
        ref = float(reference[name])
        loc = float(local[name])
        value = float(fitted[name])
        rows[name] = {
            "reference_value": ref,
            "local_comparator_value": loc,
            "fitted_value": value,
            "fitted_over_reference": (
                float(value / ref)
                if ref != 0.0
                else None
            ),
            "fitted_over_local": (
                float(value / loc)
                if loc != 0.0
                else None
            ),
            "abs_log10_distance_to_reference": (
                float(abs(np.log10(value / ref)))
                if value > 0.0 and ref > 0.0
                else None
            ),
            "abs_log10_distance_to_local": (
                float(abs(np.log10(value / loc)))
                if value > 0.0 and loc > 0.0
                else None
            ),
        }
    return rows


def fit_family(
    *,
    family,
    baseline_detector,
    reference_readout,
    local_readout,
    detector_names,
    detector_bounds,
    readout_bounds,
    frequency,
    target,
    args,
    scale_hz,
    white_asd,
    optimizer_cfg,
    seed,
    warm_solutions=(),
):
    readout_names = tuple(family["readout_parameters"])
    unknown = sorted(
        set(readout_names) - set(READOUT_PARAMETER_NAMES)
    )
    if unknown:
        raise ValueError(
            f"unknown readout parameters for {family['name']}: {unknown}"
        )

    bounds = joint_vector_bounds(
        detector_names,
        detector_bounds,
        readout_names,
        readout_bounds,
    )
    lower = np.asarray(
        [item[0] for item in bounds],
        dtype=float,
    )
    upper = np.asarray(
        [item[1] for item in bounds],
        dtype=float,
    )
    penalty = float(
        optimizer_cfg["instability_penalty"]
    )

    reference_model, reference_point = model_for_candidate(
        baseline_detector,
        reference_readout,
        frequency,
        scale_hz,
        white_asd,
    )
    if reference_model is None:
        raise ValueError(
            "tracked detector/DC baseline is invalid or unstable"
        )
    residual_template = opt.weighted_residual_vector(
        reference_model,
        target,
        frequency,
        args,
    )

    evaluation_count = 0
    instability_rejects = 0

    def evaluate(vector):
        nonlocal evaluation_count, instability_rejects
        evaluation_count += 1
        detector, readout = decode_joint(
            vector,
            detector_names,
            baseline_detector,
            readout_names,
            reference_readout,
        )
        try:
            model, point = model_for_candidate(
                detector,
                readout,
                frequency,
                scale_hz,
                white_asd,
            )
        except Exception:
            model, point = None, {}
        if model is None:
            instability_rejects += 1
            return None, detector, readout, point
        return model, detector, readout, point

    def objective(vector):
        model, _, _, _ = evaluate(vector)
        if model is None:
            return penalty
        return float(
            opt.fit_score(
                model,
                target,
                frequency,
                args,
            )
        )

    def residual(vector):
        model, _, _, _ = evaluate(vector)
        if model is None:
            return np.full_like(
                residual_template,
                np.sqrt(penalty),
                dtype=float,
            )
        return opt.weighted_residual_vector(
            model,
            target,
            frequency,
            args,
        )

    de = differential_evolution(
        objective,
        bounds,
        seed=int(seed),
        maxiter=int(optimizer_cfg["DE_maxiter"]),
        popsize=10,
        tol=1.0e-8,
        polish=False,
        workers=1,
        updating="immediate",
    )
    candidates = [
        ("de", np.asarray(de.x, dtype=float))
    ]
    if objective(de.x) < penalty:
        ls = least_squares(
            residual,
            de.x,
            bounds=(lower, upper),
            max_nfev=int(
                optimizer_cfg[
                    "least_squares_max_nfev"
                ]
            ),
            x_scale="jac",
        )
        candidates.append(
            (
                "de_least_squares",
                np.asarray(ls.x, dtype=float),
            )
        )

    warm_sources = [
        (
            "tracked_detector_dc_reference_readout",
            baseline_detector,
            reference_readout,
        ),
        (
            "tracked_detector_dc_local_readout",
            baseline_detector,
            local_readout,
        ),
        *list(warm_solutions),
    ]
    for label, detector_warm, readout_warm in warm_sources:
        try:
            vector = encode_joint(
                detector_names,
                detector_warm,
                readout_names,
                readout_warm,
            )
        except Exception:
            continue
        if np.any(vector < lower) or np.any(vector > upper):
            continue
        candidates.append(
            (f"{label}_exact", vector.copy())
        )
        if objective(vector) < penalty:
            ls = least_squares(
                residual,
                vector,
                bounds=(lower, upper),
                max_nfev=int(
                    optimizer_cfg[
                        "least_squares_max_nfev"
                    ]
                ),
                x_scale="jac",
            )
            candidates.append(
                (
                    f"{label}_least_squares",
                    np.asarray(ls.x, dtype=float),
                )
            )

    scored = []
    for label, vector in candidates:
        score = objective(vector)
        if score < penalty:
            scored.append((score, label, vector))
    if not scored:
        raise RuntimeError(
            f"no stable solution for readout family {family['name']}"
        )

    score, source, vector = min(
        scored,
        key=lambda item: item[0],
    )
    model, detector, readout, point = evaluate(vector)
    if model is None:
        raise RuntimeError("selected solution became invalid")
    metrics = holdout.model_metrics(
        model,
        target,
        frequency,
        args,
    )

    return {
        "name": family["name"],
        "readout_parameters_varied": list(readout_names),
        "n_readout_parameters": int(len(readout_names)),
        "detector_parameters_varied": list(detector_names),
        "n_detector_parameters": int(len(detector_names)),
        "n_total_free_parameters": int(
            len(detector_names) + len(readout_names)
        ),
        "best_candidate_source": source,
        "detector_candidate": {
            name: float(detector[name])
            for name in detector_names
        },
        "readout": {
            name: float(readout[name])
            for name in READOUT_PARAMETER_NAMES
        },
        "readout_distance": readout_distance(
            reference_readout,
            local_readout,
            readout,
        ),
        "detector_boundary_hits": boundary_hits(
            detector,
            detector_names,
            detector_bounds,
            competition.encode_value,
        ),
        "readout_boundary_hits": boundary_hits(
            readout,
            readout_names,
            readout_bounds,
            encode_readout_value,
        ),
        "operating_point": {
            "valid": bool(point.get("valid")),
            "stable": bool(point.get("stable")),
            "current_A": (
                float(point["current_A"])
                if point.get("current_A") is not None
                else None
            ),
            "joule_power_W": (
                float(point["joule_power_W"])
                if point.get("joule_power_W") is not None
                else None
            ),
        },
        "shape_score": float(score),
        "residual_metrics": metrics["residual_metrics"],
        "bands": metrics["bands"],
        "optimizer": {
            "evaluations_including_candidate_scoring": int(
                evaluation_count
            ),
            "instability_rejects": int(
                instability_rejects
            ),
            "de_success": bool(de.success),
        },
        "_detector_full": detector,
        "_readout_full": readout,
        "_model": model,
    }


def clean_result(row):
    return {
        key: value
        for key, value in row.items()
        if not key.startswith("_")
    }


def run(config, config_path: Path):
    manifest_path = resolve_config_path(
        config["manifest"],
        config_path,
    )
    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    cases = shared.normalize_manifest(
        manifest,
        manifest_path,
    )
    repeat_case = ident.case_by_label(
        cases,
        config["repeat_case_label"],
    )
    summary = json.loads(
        repeat_case["summary"].read_text(
            encoding="utf-8"
        )
    )
    comparison, experiment_path, comparison_source = (
        ident.comparison_for_case(repeat_case)
    )

    snapshot_path = resolve_config_path(
        config["residual_baseline_snapshot"],
        config_path,
    )
    snapshot = load_baseline_snapshot(
        snapshot_path
    )

    inherited_candidate = dict(
        summary["best_case_parameters"]
    )
    baseline_detector = dict(inherited_candidate)
    baseline_detector["R"] = float(
        snapshot["R_TES_Ohm"]
    )
    baseline_detector.update(
        snapshot["detector_candidate"]
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
        raise RuntimeError(
            "fit and holdout masks overlap"
        )

    fit_frequency = full_frequency[fit_mask]
    fit_target = full_target[fit_mask]
    fit_args = competition.band_args(
        full_args,
        fit_min,
        fit_max,
        np.count_nonzero(fit_mask),
    )
    hold_frequency = full_frequency[hold_mask]
    hold_target = full_target[hold_mask]
    hold_args = competition.band_args(
        full_args,
        hold_min,
        hold_max,
        np.count_nonzero(hold_mask),
    )

    reference_readout = snapshot[
        "reference_readout"
    ]
    local_readout = snapshot["local_readout"]
    scale_hz = float(
        config["readout_parameterization"][
            "reference_scale_Hz"
        ]
    )
    if not np.isclose(
        scale_hz,
        snapshot["reference_scale_Hz"],
        rtol=0.0,
        atol=1.0e-12,
    ):
        raise ValueError(
            "readout reference scale does not match tracked snapshot"
        )
    white_asd = float(
        snapshot["white_asd_A_rtHz"]
    )

    detector_names = tuple(
        config["detector_nuisance"]["parameters"]
    )
    if "R" in detector_names:
        raise ValueError(
            "R_TES must remain fixed in residual readout competition"
        )
    detector_bounds = competition.nuisance_bounds(
        inherited_candidate,
        {
            "bounds": config["detector_bounds"],
        },
    )
    readout_bounds = {
        name: tuple(
            float(value)
            for value in config[
                "readout_bounds"
            ][name]
        )
        for name in READOUT_PARAMETER_NAMES
    }

    families = config[
        "residual_readout_families"
    ]
    names = [row["name"] for row in families]
    if len(names) != len(set(names)):
        raise ValueError(
            "residual readout family names must be unique"
        )
    missing_nested = [
        name
        for name in MAIN_NESTED_LADDER
        if name not in names
    ]
    if missing_nested:
        raise ValueError(
            f"missing main nested families: {missing_nested}"
        )

    optimizer_cfg = config["optimizer"]
    rows = []
    warm_solutions = []

    for index, family in enumerate(families):
        row = fit_family(
            family=family,
            baseline_detector=baseline_detector,
            reference_readout=reference_readout,
            local_readout=local_readout,
            detector_names=detector_names,
            detector_bounds=detector_bounds,
            readout_bounds=readout_bounds,
            frequency=fit_frequency,
            target=fit_target,
            args=fit_args,
            scale_hz=scale_hz,
            white_asd=white_asd,
            optimizer_cfg=optimizer_cfg,
            seed=int(optimizer_cfg["seed"])
            + 100 * index,
            warm_solutions=warm_solutions,
        )

        full_model, full_point = model_for_candidate(
            row["_detector_full"],
            row["_readout_full"],
            full_frequency,
            scale_hz,
            white_asd,
        )
        if full_model is None:
            raise RuntimeError(
                f"full-band reconstruction failed for {family['name']}"
            )
        row["holdout_metrics"] = holdout.model_metrics(
            full_model[hold_mask],
            hold_target,
            hold_frequency,
            hold_args,
        )
        row["full_1_200k_metrics"] = holdout.model_metrics(
            full_model,
            full_target,
            full_frequency,
            full_args,
        )
        row["score_ratio_to_repeat_local_readout_comparator"] = float(
            row["shape_score"]
            / snapshot["repeat_local_shape_score"]
        )
        row["rms_delta_to_repeat_local_readout_comparator_dB"] = float(
            row["residual_metrics"]["rms_residual_dB"]
            - snapshot["repeat_local_rms_dB"]
        )
        screen = config["comparison_screen"]
        row["near_repeat_local_readout_comparator"] = bool(
            row[
                "score_ratio_to_repeat_local_readout_comparator"
            ]
            <= float(
                screen["near_repeat_local_score_ratio"]
            )
            and row[
                "rms_delta_to_repeat_local_readout_comparator_dB"
            ]
            <= float(
                screen[
                    "near_repeat_local_rms_delta_dB"
                ]
            )
        )
        rows.append(row)
        warm_solutions.append(
            (
                family["name"],
                row["_detector_full"],
                row["_readout_full"],
            )
        )

    fixed = next(
        row
        for row in rows
        if row["name"] == "fixed_reference_readout"
    )
    anchor_ratio = float(
        fixed["shape_score"]
        / snapshot["baseline_shape_score"]
    )
    if anchor_ratio > 1.0005:
        raise RuntimeError(
            "fixed-reference detector/DC regression anchor did not reproduce"
        )

    nested_rows = [
        next(
            row
            for row in rows
            if row["name"] == name
        )
        for name in MAIN_NESTED_LADDER
    ]
    incremental = []
    for previous, current in zip(
        nested_rows[:-1],
        nested_rows[1:],
    ):
        incremental.append(
            {
                "from_family": previous["name"],
                "to_family": current["name"],
                "added_readout_parameters": [
                    name
                    for name in current[
                        "readout_parameters_varied"
                    ]
                    if name not in previous[
                        "readout_parameters_varied"
                    ]
                ],
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
                    current["holdout_metrics"][
                        "shape_score"
                    ]
                    / previous["holdout_metrics"][
                        "shape_score"
                    ]
                ),
                "material_score_improvement": bool(
                    current["shape_score"]
                    / previous["shape_score"]
                    <= float(
                        config["comparison_screen"][
                            "material_score_improvement_ratio_to_previous_nested"
                        ]
                    )
                ),
            }
        )

    nested_near = [
        row for row in nested_rows
        if row[
            "near_repeat_local_readout_comparator"
        ]
    ]
    minimal_nested = (
        nested_near[0]["name"]
        if nested_near
        else None
    )
    any_near = sorted(
        [
            row for row in rows
            if row[
                "near_repeat_local_readout_comparator"
            ]
        ],
        key=lambda row: (
            row["n_readout_parameters"],
            row["shape_score"],
        ),
    )
    minimal_any = (
        any_near[0]["name"]
        if any_near
        else None
    )

    clean_rows = [
        clean_result(row)
        for row in rows
    ]

    by_name = {
        row["name"]: row
        for row in rows
    }
    pole_only_ok = bool(
        by_name["pole_frequency_only"][
            "near_repeat_local_readout_comparator"
        ]
    )
    pole_section_ok = bool(
        by_name["pole_section"][
            "near_repeat_local_readout_comparator"
        ]
    )
    c2_only_ok = bool(
        by_name["c2_only"][
            "near_repeat_local_readout_comparator"
        ]
    )
    pole_c2_ok = bool(
        by_name["pole_frequency_plus_c2"][
            "near_repeat_local_readout_comparator"
        ]
    )
    pole_section_c2_ok = bool(
        by_name["pole_section_plus_c2"][
            "near_repeat_local_readout_comparator"
        ]
    )
    full_ok = bool(
        by_name["full_order2"][
            "near_repeat_local_readout_comparator"
        ]
    )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "repeat_case": {
            "label": repeat_case["label"],
            "comparison_source": comparison_source,
            "experiment_path": str(experiment_path),
            "accepted_records": int(
                target_context["accepted_records"]
            ),
        },
        "fit_region": {
            "min_Hz": fit_min,
            "max_Hz_exclusive": fit_max,
            "points": int(
                np.count_nonzero(fit_mask)
            ),
        },
        "holdout_region": {
            "min_Hz": hold_min,
            "max_Hz": hold_max,
            "points": int(
                np.count_nonzero(hold_mask)
            ),
            "optimizer_received_holdout_points": False,
        },
        "fixed_detector_dc_baseline": {
            "R_TES_Ohm": float(
                baseline_detector["R"]
            ),
            "R_ratio_to_inherited": snapshot[
                "R_ratio_to_inherited"
            ],
            "tracked_nuisance_candidate": snapshot[
                "detector_candidate"
            ],
            "all_families_reoptimize_detector_nuisance": True,
            "detector_parameters": list(
                detector_names
            ),
        },
        "fixed_white_floor": {
            "white_asd_A_rtHz": white_asd,
            "recomputed_for_trials": False,
        },
        "readout_reference": {
            "order": 2,
            "reference_scale_Hz": scale_hz,
            "parameters": reference_readout,
        },
        "repeat_local_readout_comparator": {
            "parameters": local_readout,
            "shape_score": snapshot[
                "repeat_local_shape_score"
            ],
            "rms_residual_dB": snapshot[
                "repeat_local_rms_dB"
            ],
        },
        "detector_dc_regression_anchor": {
            "tracked_shape_score": snapshot[
                "baseline_shape_score"
            ],
            "recomputed_fixed_reference_shape_score": (
                fixed["shape_score"]
            ),
            "recomputed_over_tracked_score_ratio": (
                anchor_ratio
            ),
        },
        "residual_readout_fits": clean_rows,
        "main_nested_ladder": {
            "families": list(
                MAIN_NESTED_LADDER
            ),
            "incremental_gain": incremental,
            "minimal_family_near_repeat_local": (
                minimal_nested
            ),
        },
        "branch_diagnostics": {
            "c2_only": clean_result(
                by_name["c2_only"]
            ),
            "pole_frequency_plus_c2": clean_result(
                by_name[
                    "pole_frequency_plus_c2"
                ]
            ),
            "minimal_any_family_near_repeat_local": (
                minimal_any
            ),
        },
        "interpretation_flags": {
            "pole_frequency_only_sufficient": pole_only_ok,
            "pole_section_sufficient": pole_section_ok,
            "c2_only_sufficient": c2_only_ok,
            "pole_frequency_plus_c2_sufficient": pole_c2_ok,
            "pole_section_plus_c2_sufficient": (
                pole_section_c2_ok
            ),
            "full_order2_reaches_repeat_local": full_ok,
            "c4_required_after_pole_section_plus_c2": bool(
                full_ok and not pole_section_c2_ok
            ),
            "some_residual_readout_freedom_reaches_repeat_local": bool(
                any_near
            ),
            "physical_electronics_component_identified": False,
            "physical_detector_state_identified": False,
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "manifest": str(manifest_path),
            "residual_baseline_snapshot": str(
                snapshot_path
            ),
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
        json.dumps(
            result,
            indent=2,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(output),
                "detector_dc_regression_anchor": result[
                    "detector_dc_regression_anchor"
                ],
                "minimal_nested_family_near_local": result[
                    "main_nested_ladder"
                ]["minimal_family_near_repeat_local"],
                "minimal_any_family_near_local": result[
                    "branch_diagnostics"
                ]["minimal_any_family_near_repeat_local"],
                "fits": [
                    {
                        "name": row["name"],
                        "n_readout_parameters": row[
                            "n_readout_parameters"
                        ],
                        "shape_score": row[
                            "shape_score"
                        ],
                        "rms_residual_dB": row[
                            "residual_metrics"
                        ]["rms_residual_dB"],
                        "score_ratio_to_local": row[
                            "score_ratio_to_repeat_local_readout_comparator"
                        ],
                        "near_local": row[
                            "near_repeat_local_readout_comparator"
                        ],
                        "holdout_shape_score": row[
                            "holdout_metrics"
                        ]["shape_score"],
                    }
                    for row in result[
                        "residual_readout_fits"
                    ]
                ],
                "flags": result[
                    "interpretation_flags"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

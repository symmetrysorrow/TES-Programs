"""Compete limited detector-state nuisance freedom against a fixed reference readout shape.

The 2024-12-05 repeat target is reconstructed from its tracked accepted-record
mask.  The post-filter white ASD is frozen to the previously profiled repeat
value.  The effective order-2 readout transfer is frozen to the 2024-12-06
reference solution.

Nested detector-state nuisance families are then fit only on 1-40 kHz:

  1) alpha, beta
  2) alpha, beta, C_tes, L
  3) alpha, beta, C_tes, L, T_bath

The 40-200 kHz region is a strict holdout.  The repeat-local order-2 readout fit
with the inherited detector snapshot is retained as the competing comparator.

This is a confound test, not a physical detector-state measurement.
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
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "readout_detector_state_competition_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_detector_state_competition_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import readout_effective_numerator_diagnostic as effective  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import readout_lowmid_holdout_diagnostic as holdout  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


LOG_PARAMETERS = {"C_tes", "L"}


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def load_day_snapshot(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    reference = payload["reference_order2"]
    repeat = payload["repeat_full_order2"]

    def transfer(row):
        return {
            "pole_Hz": float(row["pole_Hz"]),
            "pole_Q": float(row["pole_Q"]),
            "latent": {
                "u": float(row["latent"]["u"]),
                "v": float(row["latent"]["v"]),
            },
        }

    return {
        "reference_transfer": transfer(reference),
        "repeat_transfer": transfer(repeat),
        "reference_canonical": {
            "c2": float(reference["canonical"]["c2"]),
            "c4": float(reference["canonical"]["c4"]),
        },
        "repeat_canonical": {
            "c2": float(repeat["canonical"]["c2"]),
            "c4": float(repeat["canonical"]["c4"]),
        },
        "repeat_white_asd_A_rtHz": float(
            repeat["white_asd_A_rtHz"]
        ),
        "repeat_local_expected_score": float(
            repeat["lowmid_shape_score"]
        ),
        "repeat_local_expected_rms_dB": float(
            repeat["lowmid_rms_residual_dB"]
        ),
        "provenance": payload.get("provenance", {}),
    }


def band_args(original, min_hz, max_hz, points):
    values = dict(vars(original))
    values["fit_min_hz"] = float(min_hz)
    values["fit_max_hz"] = float(max_hz)
    values["fit_points"] = int(points)
    return SimpleNamespace(**values)


def nuisance_bounds(candidate, config):
    bounds_cfg = config["bounds"]
    alpha_cfg = bounds_cfg["alpha"]
    baseline_alpha = float(candidate["alpha"])
    bounds = {
        "alpha": (
            max(
                baseline_alpha
                * float(alpha_cfg["minimum_fraction_of_baseline"]),
                1.0e-12,
            ),
            float(alpha_cfg["maximum"]),
        ),
        "beta": (
            float(bounds_cfg["beta"]["min"]),
            float(bounds_cfg["beta"]["max"]),
        ),
        "C_tes": (
            float(opt.C_TES_FIT_MIN_J_PER_K),
            float(opt.C_TES_FIT_MAX_J_PER_K),
        ),
        "L": (
            float(opt.L_FIT_MIN_H),
            float(opt.L_FIT_MAX_H),
        ),
        "T_bath": (
            float(candidate["T_bath"])
            - float(bounds_cfg["T_bath"]["half_width_K"]),
            float(candidate["T_bath"])
            + float(bounds_cfg["T_bath"]["half_width_K"]),
        ),
    }
    if bounds["alpha"][1] < bounds["alpha"][0]:
        bounds["alpha"] = (
            bounds["alpha"][0],
            bounds["alpha"][0],
        )
    return bounds


def encode_value(name, value):
    value = float(value)
    if name in LOG_PARAMETERS:
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
        return float(np.log10(value))
    return value


def decode_value(name, value):
    value = float(value)
    if name in LOG_PARAMETERS:
        return float(10.0**value)
    return value


def vector_bounds(parameter_names, physical_bounds):
    result = []
    for name in parameter_names:
        lower, upper = physical_bounds[name]
        result.append(
            (
                encode_value(name, lower),
                encode_value(name, upper),
            )
        )
    return tuple(result)


def encode_candidate(parameter_names, candidate):
    return np.asarray(
        [
            encode_value(name, candidate[name])
            for name in parameter_names
        ],
        dtype=float,
    )


def decode_candidate(
    vector,
    parameter_names,
    baseline_candidate,
):
    trial = dict(baseline_candidate)
    for name, encoded in zip(parameter_names, vector):
        trial[name] = decode_value(name, encoded)
    return trial


def fixed_white_context(candidate, frequency, white_asd):
    context = base.intrinsic_context(
        candidate,
        frequency,
    )
    context = dict(context)
    context["post_filter_white_asd_A_rtHz"] = float(
        white_asd
    )
    return context


def model_for_candidate(
    candidate,
    frequency,
    fixed_transfer,
    scale_hz,
    white_asd,
):
    point = opt.tes_operating_point(candidate)
    if not point.get("valid") or not point.get("stable"):
        return None, point
    context = fixed_white_context(
        candidate,
        frequency,
        white_asd,
    )
    model = effective.pre_analysis_model(
        context,
        fixed_transfer,
        2,
        scale_hz,
    )
    if np.any(~np.isfinite(model)) or np.any(model <= 0.0):
        return None, point
    return model, point


def metrics(model, target, frequency, args):
    return holdout.model_metrics(
        model,
        target,
        frequency,
        args,
    )


def parameter_boundary_hits(
    candidate,
    parameter_names,
    physical_bounds,
):
    hits = {}
    for name in parameter_names:
        value = float(candidate[name])
        lower, upper = physical_bounds[name]
        span = max(abs(upper - lower), abs(upper), 1.0e-30)
        tolerance = span * 1.0e-5
        hits[name] = {
            "value": value,
            "lower": float(lower),
            "upper": float(upper),
            "at_lower": bool(
                abs(value - lower) <= tolerance
            ),
            "at_upper": bool(
                abs(value - upper) <= tolerance
            ),
        }
    return hits


def fit_nuisance_family(
    *,
    name,
    parameter_names,
    baseline_candidate,
    fixed_transfer,
    frequency,
    target,
    args,
    scale_hz,
    white_asd,
    physical_bounds,
    optimizer_cfg,
    seed,
    warm_candidates=(),
):
    parameter_names = tuple(parameter_names)
    bounds = vector_bounds(
        parameter_names,
        physical_bounds,
    )
    penalty = float(
        optimizer_cfg["instability_penalty"]
    )

    baseline_model, baseline_point = model_for_candidate(
        baseline_candidate,
        frequency,
        fixed_transfer,
        scale_hz,
        white_asd,
    )
    if baseline_model is None:
        raise ValueError(
            "baseline detector candidate is invalid or unstable"
        )
    residual_template = opt.weighted_residual_vector(
        baseline_model,
        target,
        frequency,
        args,
    )

    evaluation_count = 0
    instability_rejects = 0

    def evaluate(vector):
        nonlocal evaluation_count, instability_rejects
        evaluation_count += 1
        trial = decode_candidate(
            vector,
            parameter_names,
            baseline_candidate,
        )
        try:
            model, point = model_for_candidate(
                trial,
                frequency,
                fixed_transfer,
                scale_hz,
                white_asd,
            )
        except Exception:
            model, point = None, {}
        if model is None:
            instability_rejects += 1
            return None, trial, point
        return model, trial, point

    def objective(vector):
        model, _, _ = evaluate(vector)
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
        model, _, _ = evaluate(vector)
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

    lower = np.asarray(
        [item[0] for item in bounds],
        dtype=float,
    )
    upper = np.asarray(
        [item[1] for item in bounds],
        dtype=float,
    )

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
        )
        candidates.append(
            (
                "de_least_squares",
                np.asarray(ls.x, dtype=float),
            )
        )

    warm_sources = [
        ("baseline", baseline_candidate),
        *list(warm_candidates),
    ]
    for label, warm_candidate in warm_sources:
        try:
            vector = encode_candidate(
                parameter_names,
                warm_candidate,
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
        if score >= penalty:
            continue
        scored.append((score, label, vector))
    if not scored:
        raise RuntimeError(
            f"no stable candidate found for nuisance family {name}"
        )

    score, source, vector = min(
        scored,
        key=lambda item: item[0],
    )
    model, trial, point = evaluate(vector)
    if model is None:
        raise RuntimeError(
            "selected nuisance candidate became invalid"
        )

    result_metrics = metrics(
        model,
        target,
        frequency,
        args,
    )
    baseline_point = opt.tes_operating_point(
        baseline_candidate
    )

    return {
        "name": name,
        "parameters_varied": list(parameter_names),
        "n_free_parameters": len(parameter_names),
        "best_candidate_source": source,
        "candidate": {
            key: float(trial[key])
            for key in parameter_names
        },
        "ratios_to_inherited_detector": {
            key: float(
                trial[key] / baseline_candidate[key]
            )
            if float(baseline_candidate[key]) != 0.0
            else None
            for key in parameter_names
        },
        "boundary_hits": parameter_boundary_hits(
            trial,
            parameter_names,
            physical_bounds,
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
            "current_ratio_to_inherited": (
                float(
                    point["current_A"]
                    / baseline_point["current_A"]
                )
                if point.get("current_A") is not None
                and baseline_point.get("current_A")
                else None
            ),
            "joule_power_ratio_to_inherited": (
                float(
                    point["joule_power_W"]
                    / baseline_point["joule_power_W"]
                )
                if point.get("joule_power_W") is not None
                and baseline_point.get("joule_power_W")
                else None
            ),
        },
        "shape_score": float(score),
        "residual_metrics": result_metrics[
            "residual_metrics"
        ],
        "bands": result_metrics["bands"],
        "optimizer": {
            "evaluations_including_candidate_scoring": int(
                evaluation_count
            ),
            "instability_rejects": int(
                instability_rejects
            ),
            "de_success": bool(de.success),
        },
        "_candidate_full": trial,
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
        config["order2_day_snapshot"],
        config_path,
    )
    snapshot = load_day_snapshot(snapshot_path)

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
    target = target_context["target"]
    inherited_candidate = dict(
        summary["best_case_parameters"]
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
        raise RuntimeError(
            "fit and holdout masks overlap"
        )

    fit_frequency = full_frequency[fit_mask]
    fit_target = target[fit_mask]
    fit_args = band_args(
        full_args,
        fit_min,
        fit_max,
        np.count_nonzero(fit_mask),
    )
    hold_frequency = full_frequency[hold_mask]
    hold_target = target[hold_mask]
    hold_args = band_args(
        full_args,
        hold_min,
        hold_max,
        np.count_nonzero(hold_mask),
    )

    scale_hz = float(
        config["fixed_readout"]["reference_scale_Hz"]
    )
    white_asd = float(
        snapshot["repeat_white_asd_A_rtHz"]
    )
    reference_transfer = snapshot[
        "reference_transfer"
    ]
    repeat_transfer = snapshot[
        "repeat_transfer"
    ]

    inherited_fit_context = fixed_white_context(
        inherited_candidate,
        fit_frequency,
        white_asd,
    )
    fixed_reference_fit_model = (
        effective.pre_analysis_model(
            inherited_fit_context,
            reference_transfer,
            2,
            scale_hz,
        )
    )
    repeat_local_fit_model = (
        effective.pre_analysis_model(
            inherited_fit_context,
            repeat_transfer,
            2,
            scale_hz,
        )
    )
    fixed_reference_fit_metrics = metrics(
        fixed_reference_fit_model,
        fit_target,
        fit_frequency,
        fit_args,
    )
    repeat_local_fit_metrics = metrics(
        repeat_local_fit_model,
        fit_target,
        fit_frequency,
        fit_args,
    )

    if not np.isclose(
        repeat_local_fit_metrics["shape_score"],
        snapshot["repeat_local_expected_score"],
        rtol=2.0e-4,
        atol=1.0e-12,
    ):
        raise ValueError(
            "tracked repeat-local order-2 comparator does not reproduce"
        )

    inherited_full_context = fixed_white_context(
        inherited_candidate,
        full_frequency,
        white_asd,
    )
    fixed_reference_full_model = (
        effective.pre_analysis_model(
            inherited_full_context,
            reference_transfer,
            2,
            scale_hz,
        )
    )
    repeat_local_full_model = (
        effective.pre_analysis_model(
            inherited_full_context,
            repeat_transfer,
            2,
            scale_hz,
        )
    )

    fixed_reference_holdout = metrics(
        fixed_reference_full_model[hold_mask],
        hold_target,
        hold_frequency,
        hold_args,
    )
    repeat_local_holdout = metrics(
        repeat_local_full_model[hold_mask],
        hold_target,
        hold_frequency,
        hold_args,
    )

    physical_bounds = nuisance_bounds(
        inherited_candidate,
        config,
    )
    optimizer_cfg = config["optimizer"]
    family_rows = []
    previous_best = []

    for index, family in enumerate(
        config["nuisance_families"]
    ):
        row = fit_nuisance_family(
            name=family["name"],
            parameter_names=family["parameters"],
            baseline_candidate=inherited_candidate,
            fixed_transfer=reference_transfer,
            frequency=fit_frequency,
            target=fit_target,
            args=fit_args,
            scale_hz=scale_hz,
            white_asd=white_asd,
            physical_bounds=physical_bounds,
            optimizer_cfg=optimizer_cfg,
            seed=int(optimizer_cfg["seed"])
            + index,
            warm_candidates=previous_best,
        )

        best_full_context = fixed_white_context(
            row["_candidate_full"],
            full_frequency,
            white_asd,
        )
        best_full_model = effective.pre_analysis_model(
            best_full_context,
            reference_transfer,
            2,
            scale_hz,
        )
        row["holdout_metrics"] = metrics(
            best_full_model[hold_mask],
            hold_target,
            hold_frequency,
            hold_args,
        )
        row["full_1_200k_metrics"] = metrics(
            best_full_model,
            target,
            full_frequency,
            full_args,
        )
        row["score_ratio_to_fixed_reference_inherited_detector"] = float(
            row["shape_score"]
            / fixed_reference_fit_metrics["shape_score"]
        )
        row["score_ratio_to_repeat_local_readout_comparator"] = float(
            row["shape_score"]
            / repeat_local_fit_metrics["shape_score"]
        )
        row["rms_delta_to_repeat_local_readout_comparator_dB"] = float(
            row["residual_metrics"]["rms_residual_dB"]
            - repeat_local_fit_metrics[
                "residual_metrics"
            ]["rms_residual_dB"]
        )
        family_rows.append(row)
        previous_best.append(
            (
                family["name"],
                row["_candidate_full"],
            )
        )

    screen = config["comparison_screen"]
    for row in family_rows:
        row["material_improvement_vs_fixed_reference"] = bool(
            row[
                "score_ratio_to_fixed_reference_inherited_detector"
            ]
            <= float(
                screen[
                    "material_improvement_vs_fixed_reference_ratio"
                ]
            )
        )
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

    best_nuisance = min(
        family_rows,
        key=lambda row: row["shape_score"],
    )
    detector_can_compete = bool(
        best_nuisance[
            "material_improvement_vs_fixed_reference"
        ]
        and best_nuisance[
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
        "fixed_readout": {
            "order": 2,
            "reference_scale_Hz": scale_hz,
            "source": (
                "2024-12-06 reference order-2 "
                "effective numerator snapshot"
            ),
            "parameters": reference_transfer,
            "canonical": snapshot[
                "reference_canonical"
            ],
        },
        "fixed_repeat_white_floor": {
            "white_asd_A_rtHz": white_asd,
            "source": (
                "tracked 2024-12-05 profiled "
                "post-filter white floor"
            ),
            "recomputed_for_detector_trials": False,
        },
        "inherited_detector_baseline": {
            "nuisance_parameters": {
                name: float(
                    inherited_candidate[name]
                )
                for name in (
                    "alpha",
                    "beta",
                    "C_tes",
                    "L",
                    "T_bath",
                )
            },
            "fixed_reference_readout_metrics": (
                fixed_reference_fit_metrics
            ),
            "fixed_reference_readout_holdout": (
                fixed_reference_holdout
            ),
        },
        "repeat_local_readout_comparator": {
            "detector_state": "inherited/frozen",
            "readout_parameters": repeat_transfer,
            "canonical": snapshot[
                "repeat_canonical"
            ],
            "fit_metrics": repeat_local_fit_metrics,
            "holdout_metrics": repeat_local_holdout,
        },
        "nuisance_bounds": {
            key: {
                "lower": float(value[0]),
                "upper": float(value[1]),
                "log_optimized": bool(
                    key in LOG_PARAMETERS
                ),
            }
            for key, value in physical_bounds.items()
        },
        "detector_state_nuisance_fits": [
            clean_result(row)
            for row in family_rows
        ],
        "best_detector_state_nuisance_family": (
            best_nuisance["name"]
        ),
        "comparison_screen": screen,
        "interpretation_flags": {
            "detector_state_nuisance_materially_improves_fixed_reference_readout": bool(
                best_nuisance[
                    "material_improvement_vs_fixed_reference"
                ]
            ),
            "detector_state_nuisance_reaches_near_repeat_local_readout_fit": bool(
                best_nuisance[
                    "near_repeat_local_readout_comparator"
                ]
            ),
            "limited_detector_state_nuisance_can_compete_with_day_specific_readout_shape": (
                detector_can_compete
            ),
            "day_specific_readout_shape_still_required_within_tested_nuisance_set": bool(
                not detector_can_compete
            ),
            "physical_electronics_drift_identified": False,
            "physical_detector_state_shift_identified": False,
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "manifest": str(manifest_path),
            "order2_day_snapshot": str(
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
                "fixed_reference_score": result[
                    "inherited_detector_baseline"
                ][
                    "fixed_reference_readout_metrics"
                ]["shape_score"],
                "repeat_local_score": result[
                    "repeat_local_readout_comparator"
                ]["fit_metrics"]["shape_score"],
                "best_nuisance_family": result[
                    "best_detector_state_nuisance_family"
                ],
                "nuisance_fits": [
                    {
                        "name": row["name"],
                        "shape_score": row[
                            "shape_score"
                        ],
                        "rms_residual_dB": row[
                            "residual_metrics"
                        ]["rms_residual_dB"],
                        "score_ratio_to_fixed_reference": row[
                            "score_ratio_to_fixed_reference_inherited_detector"
                        ],
                        "score_ratio_to_repeat_local": row[
                            "score_ratio_to_repeat_local_readout_comparator"
                        ],
                        "near_repeat_local": row[
                            "near_repeat_local_readout_comparator"
                        ],
                    }
                    for row in result[
                        "detector_state_nuisance_fits"
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

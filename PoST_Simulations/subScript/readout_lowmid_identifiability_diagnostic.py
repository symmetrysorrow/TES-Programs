"""Profile low/mid hybrid-transfer identifiability after white-floor separation.

The 2024-12-05 repeat target is restricted to 1-40 kHz.  The post-filter white
floor is frozen to the preceding profiled value.  This diagnostic profiles
selected transfer parameters by fixing one value at a time and refitting all
remaining hybrid parameters.  It also fits reduced topologies where selected
parameters are fixed to the shared-reference values.

The purpose is identifiability, not another physical transfer fit.
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
DEFAULT_CONFIG = CONFIG_DIR / "readout_lowmid_identifiability_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_lowmid_identifiability_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import preanalysis_hybrid_pole_profile_diagnostic as profile  # noqa: E402
from subScript import readout_lowmid_holdout_diagnostic as holdout  # noqa: E402
from subScript import readout_white_floor_separation_diagnostic as white_sep  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


PARAMETER_NAMES = tuple(shared.TRANSFER_PARAMETER_NAMES)


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def case_by_label(cases, label):
    matches = [case for case in cases if case["label"] == label]
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one manifest case named {label!r}; "
            f"found {len(matches)}"
        )
    return matches[0]


def comparison_for_case(case):
    if case["comparison_summary"] is not None:
        comparison = json.loads(
            case["comparison_summary"].read_text(encoding="utf-8")
        )
        source = "tracked_comparison_summary"
    else:
        comparison = shared.build_comparison_from_spec(
            case["comparison_spec"]
        )
        source = "tracked_comparison_spec_recomputed_mask"
    experiment_path = (
        case["experiment_path"]
        if case["experiment_path"] is not None
        else Path(comparison["experiment_path"])
    )
    return comparison, experiment_path, source


def load_free_snapshot(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    transfer = payload["transfer"]
    return {
        "parameters": {
            name: float(transfer["parameters"][name])
            for name in PARAMETER_NAMES
        },
        "white_scale": float(payload["white_floor"]["scale"]),
        "white_asd_A_rtHz": float(
            payload["white_floor"]["white_asd_A_rtHz"]
        ),
        "shape_score": float(
            payload["provenance"]["free_fit_shape_score"]
        ),
        "rms_residual_dB": float(
            payload["provenance"]["free_fit_rms_residual_dB"]
        ),
        "max_abs_residual_dB": float(
            payload["provenance"]["free_fit_max_abs_residual_dB"]
        ),
    }


def band_args(original, min_hz, max_hz, points):
    values = dict(vars(original))
    values["fit_min_hz"] = float(min_hz)
    values["fit_max_hz"] = float(max_hz)
    values["fit_points"] = int(points)
    return SimpleNamespace(**values)


def decode_free(vector, free_names, fixed_parameters):
    result = {
        name: float(value)
        for name, value in fixed_parameters.items()
    }
    values = 10.0 ** np.asarray(vector, dtype=float)
    for name, value in zip(free_names, values):
        result[name] = float(value)
    if set(result) != set(PARAMETER_NAMES):
        missing = sorted(set(PARAMETER_NAMES) - set(result))
        raise ValueError(f"incomplete parameter vector: missing {missing}")
    return result


def fit_with_fixed_parameters(
    context,
    target,
    frequency,
    args,
    *,
    fixed_parameters,
    shared_parameters,
    free_reference_parameters,
    center_min_hz,
    center_max_hz,
    general_q_min,
    q_max,
    seed,
    de_maxiter,
):
    fixed_parameters = {
        name: float(value)
        for name, value in fixed_parameters.items()
    }
    unknown = sorted(set(fixed_parameters) - set(PARAMETER_NAMES))
    if unknown:
        raise ValueError(f"unknown fixed parameters: {unknown}")

    all_names, all_bounds = profile.parameter_bounds(
        center_min_hz,
        center_max_hz,
        general_q_min,
        q_max,
    )
    bound_map = dict(zip(all_names, all_bounds))
    free_names = tuple(
        name for name in PARAMETER_NAMES
        if name not in fixed_parameters
    )
    bounds = tuple(bound_map[name] for name in free_names)

    for name, value in fixed_parameters.items():
        lower, upper = bound_map[name]
        log_value = np.log10(value)
        if log_value < lower - 1e-12 or log_value > upper + 1e-12:
            raise ValueError(
                f"fixed {name}={value} lies outside configured bounds"
            )

    def model_from_vector(vector):
        parameters = decode_free(
            vector,
            free_names,
            fixed_parameters,
        )
        return profile.pre_analysis_model(
            context,
            parameters,
            None,
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

    candidate_vectors = []
    evaluations = 0
    success = False

    if free_names:
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
        evaluations += int(de.nfev)
        success = bool(de.success)
        candidate_vectors.append(("de", np.asarray(de.x, dtype=float)))

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
        candidate_vectors.append(
            ("de_least_squares", np.asarray(ls.x, dtype=float))
        )

        warm_sources = (
            ("shared", shared_parameters),
            ("free_reference", free_reference_parameters),
        )
        for label, parameters in warm_sources:
            warm = np.log10(
                np.asarray(
                    [float(parameters[name]) for name in free_names],
                    dtype=float,
                )
            )
            if np.any(warm < lower) or np.any(warm > upper):
                continue
            candidate_vectors.append(
                (f"{label}_exact", warm.copy())
            )
            warm_ls = least_squares(
                residual_vector,
                warm,
                bounds=(lower, upper),
                max_nfev=3000,
            )
            evaluations += int(warm_ls.nfev)
            success = bool(success or warm_ls.success)
            candidate_vectors.append(
                (
                    f"{label}_least_squares",
                    np.asarray(warm_ls.x, dtype=float),
                )
            )
    else:
        candidate_vectors.append(
            ("all_parameters_fixed", np.asarray([], dtype=float))
        )
        success = True

    best_source, best_vector = min(
        candidate_vectors,
        key=lambda item: objective(item[1]),
    )
    parameters = decode_free(
        best_vector,
        free_names,
        fixed_parameters,
    )
    model = model_from_vector(best_vector)
    score = float(
        opt.fit_score(model, target, frequency, args)
    )
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
        "fixed_parameters": fixed_parameters,
        "free_parameters": list(free_names),
        "n_free_parameters": int(len(free_names)),
        "success": bool(success),
        "evaluations": int(evaluations),
        "best_candidate_source": best_source,
        "parameters": parameters,
        "described_parameters": profile.describe_parameters(
            parameters,
            None,
        ),
        "shape_score": score,
        "residual_metrics": residual,
        "bands": bands,
        "_model": model,
    }


def annotate_profile(rows, free_score, free_rms, screen):
    rms_delta_limit = float(screen["near_free_rms_delta_dB"])
    score_ratio_limit = float(screen["near_free_score_ratio"])
    annotated = []
    for row in rows:
        score_ratio = float(row["shape_score"] / free_score)
        rms_delta = float(
            row["residual_metrics"]["rms_residual_dB"] - free_rms
        )
        near_free = bool(
            score_ratio <= score_ratio_limit
            and rms_delta <= rms_delta_limit
        )
        clean = {
            key: value
            for key, value in row.items()
            if not key.startswith("_")
        }
        clean["score_ratio_to_free"] = score_ratio
        clean["rms_delta_from_free_dB"] = rms_delta
        clean["near_free"] = near_free
        annotated.append(clean)
    return annotated


def summarize_profile(parameter_name, rows):
    values = np.asarray(
        [float(row["fixed_parameters"][parameter_name]) for row in rows],
        dtype=float,
    )
    near = np.asarray(
        [bool(row["near_free"]) for row in rows],
        dtype=bool,
    )
    scores = np.asarray(
        [float(row["shape_score"]) for row in rows],
        dtype=float,
    )
    best_index = int(np.argmin(scores))
    near_values = values[near]
    span_ratio = (
        float(np.max(near_values) / np.min(near_values))
        if near_values.size >= 2
        else None
    )
    broad = bool(
        near_values.size >= 3
        and span_ratio is not None
        and span_ratio >= 2.0
    )
    return {
        "parameter": parameter_name,
        "n_profile_points": int(len(rows)),
        "n_near_free_points": int(np.count_nonzero(near)),
        "fraction_near_free": float(np.mean(near)),
        "near_free_min_value": (
            float(np.min(near_values))
            if near_values.size
            else None
        ),
        "near_free_max_value": (
            float(np.max(near_values))
            if near_values.size
            else None
        ),
        "near_free_span_ratio": span_ratio,
        "best_profile_value": float(values[best_index]),
        "best_profile_shape_score": float(scores[best_index]),
        "broad_nonidentifiability_over_tested_range": broad,
    }


def run(config, config_path: Path):
    manifest_path = resolve_config_path(
        config["manifest"],
        config_path,
    )
    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    cases = shared.normalize_manifest(manifest, manifest_path)
    repeat_case = case_by_label(
        cases,
        config["repeat_case_label"],
    )
    summary = json.loads(
        repeat_case["summary"].read_text(encoding="utf-8")
    )
    comparison, experiment_path, comparison_source = (
        comparison_for_case(repeat_case)
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
    candidate = dict(summary["best_case_parameters"])
    base_context = base.intrinsic_context(
        candidate,
        full_frequency,
    )

    reference_path = resolve_config_path(
        config["reference_transfer"],
        config_path,
    )
    reference_payload = json.loads(
        reference_path.read_text(encoding="utf-8")
    )
    reference = shared.load_shared_transfer(reference_payload)

    free_snapshot_path = resolve_config_path(
        config["lowmid_free_snapshot"],
        config_path,
    )
    free_snapshot = load_free_snapshot(free_snapshot_path)
    context = white_sep.white_scaled_context(
        base_context,
        free_snapshot["white_scale"],
    )
    if not np.isclose(
        float(context["post_filter_white_asd_A_rtHz"]),
        free_snapshot["white_asd_A_rtHz"],
        rtol=5.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "tracked low/mid snapshot white ASD does not match "
            "detector-snapshot baseline"
        )

    fit_region = config["fit_region_Hz"]
    fit_min = float(fit_region["min"])
    fit_max = float(fit_region["max_exclusive"])
    mask = (
        (full_frequency >= fit_min)
        & (full_frequency < fit_max)
    )
    frequency = full_frequency[mask]
    target = target_context["target"][mask]
    fit_context = holdout.subset_context(context, mask)
    fit_args = band_args(
        full_args,
        fit_min,
        fit_max,
        np.count_nonzero(mask),
    )

    free_model = profile.pre_analysis_model(
        fit_context,
        free_snapshot["parameters"],
        None,
    )
    free_score = float(
        opt.fit_score(
            free_model,
            target,
            frequency,
            fit_args,
        )
    )
    free_residual = base.residual_db_metrics(
        free_model,
        target,
        frequency,
        fit_args,
    )
    if not np.isclose(
        free_score,
        free_snapshot["shape_score"],
        rtol=2e-6,
        atol=1e-12,
    ):
        raise ValueError(
            "low/mid free snapshot score does not reproduce"
        )

    optimizer = config["optimizer"]
    screen = config["identifiability_screen"]
    profiles = {}
    profile_summaries = {}

    seed_base = int(optimizer["seed"])
    profile_index = 0
    for parameter_name, values in config["fixed_profiles"].items():
        rows = []
        for value in values:
            result = fit_with_fixed_parameters(
                fit_context,
                target,
                frequency,
                fit_args,
                fixed_parameters={
                    parameter_name: float(value)
                },
                shared_parameters=reference["parameters"],
                free_reference_parameters=free_snapshot[
                    "parameters"
                ],
                center_min_hz=float(
                    optimizer["center_min_Hz"]
                ),
                center_max_hz=float(
                    optimizer["center_max_Hz"]
                ),
                general_q_min=float(
                    optimizer["general_Q_min"]
                ),
                q_max=float(optimizer["Q_max"]),
                seed=seed_base + profile_index,
                de_maxiter=int(optimizer["DE_maxiter"]),
            )
            rows.append(result)
            profile_index += 1
        annotated = annotate_profile(
            rows,
            free_score,
            free_residual["rms_residual_dB"],
            screen,
        )
        profiles[parameter_name] = annotated
        profile_summaries[parameter_name] = summarize_profile(
            parameter_name,
            annotated,
        )

    reduced = []
    for index, spec in enumerate(config["reduced_topologies"]):
        fixed = {
            name: float(reference["parameters"][name])
            for name in spec["fixed_parameters"]
        }
        result = fit_with_fixed_parameters(
            fit_context,
            target,
            frequency,
            fit_args,
            fixed_parameters=fixed,
            shared_parameters=reference["parameters"],
            free_reference_parameters=free_snapshot["parameters"],
            center_min_hz=float(optimizer["center_min_Hz"]),
            center_max_hz=float(optimizer["center_max_Hz"]),
            general_q_min=float(optimizer["general_Q_min"]),
            q_max=float(optimizer["Q_max"]),
            seed=seed_base + 1000 + index,
            de_maxiter=int(optimizer["DE_maxiter"]),
        )
        clean = {
            key: value
            for key, value in result.items()
            if not key.startswith("_")
        }
        clean["name"] = spec["name"]
        clean["score_ratio_to_free"] = float(
            result["shape_score"] / free_score
        )
        clean["rms_delta_from_free_dB"] = float(
            result["residual_metrics"]["rms_residual_dB"]
            - free_residual["rms_residual_dB"]
        )
        clean["near_free"] = bool(
            clean["score_ratio_to_free"]
            <= float(screen["near_free_score_ratio"])
            and clean["rms_delta_from_free_dB"]
            <= float(screen["near_free_rms_delta_dB"])
        )
        reduced.append(clean)

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "fit_region": {
            "min_Hz": fit_min,
            "max_Hz_exclusive": fit_max,
            "points": int(np.count_nonzero(mask)),
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
        "free_low_mid_reference": {
            "parameters": free_snapshot["parameters"],
            "tracked_shape_score": free_snapshot["shape_score"],
            "recomputed_shape_score": free_score,
            "tracked_rms_residual_dB": free_snapshot[
                "rms_residual_dB"
            ],
            "recomputed_rms_residual_dB": free_residual[
                "rms_residual_dB"
            ],
        },
        "identifiability_screen": screen,
        "fixed_parameter_profiles": profiles,
        "profile_summaries": profile_summaries,
        "reduced_topology_fits": reduced,
        "interpretation_flags": {
            "zero_frequency_broadly_nonidentified": bool(
                profile_summaries["zero_Hz"][
                    "broad_nonidentifiability_over_tested_range"
                ]
            ),
            "lead_zero_broadly_nonidentified": bool(
                profile_summaries["leadlag_zero_Hz"][
                    "broad_nonidentifiability_over_tested_range"
                ]
            ),
            "pole_frequency_broadly_nonidentified": bool(
                profile_summaries["pole_Hz"][
                    "broad_nonidentifiability_over_tested_range"
                ]
            ),
            "shared_zero_section_reduced_fit_near_free": bool(
                next(
                    row["near_free"]
                    for row in reduced
                    if row["name"]
                    == "shared_zero_section_fixed"
                )
            ),
            "shared_lead_zero_reduced_fit_near_free": bool(
                next(
                    row["near_free"]
                    for row in reduced
                    if row["name"]
                    == "shared_lead_zero_fixed"
                )
            ),
            "shared_zero_and_lead_reduced_fit_near_free": bool(
                next(
                    row["near_free"]
                    for row in reduced
                    if row["name"]
                    == "shared_zero_section_and_lead_zero_fixed"
                )
            ),
            "physical_readout_component_identified": False,
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
                "profile_summaries": result["profile_summaries"],
                "reduced_topologies": [
                    {
                        "name": row["name"],
                        "shape_score": row["shape_score"],
                        "score_ratio_to_free": row[
                            "score_ratio_to_free"
                        ],
                        "rms_delta_from_free_dB": row[
                            "rms_delta_from_free_dB"
                        ],
                        "near_free": row["near_free"],
                    }
                    for row in result["reduced_topology_fits"]
                ],
                "flags": result["interpretation_flags"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

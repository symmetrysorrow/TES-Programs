"""Separate repeat-case high-frequency white-floor leverage from readout-transfer drift.

Stage 0 evaluates the 2024-12-05 repeat with the reference shared transfer and
its inherited production post-filter white floor.

Stage 1 keeps every detector and shared-transfer parameter fixed and profiles
only the post-filter white ASD amplitude.

Stage 2 freezes the Stage-1 white floor, then refits the five parameters of the
same infinity-pole hybrid transfer with the exact shared transfer retained as a
warm-start candidate.

This ordering tests whether the previously observed repeat-local transfer drift
was partly compensating for a case-specific high-frequency white floor.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "readout_white_floor_separation_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_white_floor_separation_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import preanalysis_hybrid_pole_profile_diagnostic as profile  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


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


def load_previous_local_transfer(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    transfer = payload["transfer"]
    parameters = {
        name: float(transfer["parameters"][name])
        for name in shared.TRANSFER_PARAMETER_NAMES
    }
    fixed_pole = transfer.get("fixed_leadlag_pole_Hz")
    if transfer.get("leadlag_pole_is_infinite", fixed_pole is None):
        fixed_pole = None
    elif fixed_pole is None:
        raise ValueError(
            "finite previous local lead/lag pole requires a value"
        )
    else:
        fixed_pole = float(fixed_pole)
    return {
        "parameters": parameters,
        "fixed_leadlag_pole_Hz": fixed_pole,
        "provenance": payload.get("provenance", {}),
    }


def white_scaled_context(context, scale):
    scale = float(scale)
    if scale < 0.0 or not np.isfinite(scale):
        raise ValueError("white scale must be finite and non-negative")
    result = dict(context)
    baseline = float(context["post_filter_white_asd_A_rtHz"])
    result["post_filter_white_asd_A_rtHz"] = baseline * scale
    return result


def model_summary(model, target, frequency, fit_args, *, rms_screen, max_screen):
    score = float(
        opt.fit_score(
            model,
            target,
            frequency,
            fit_args,
        )
    )
    residual = base.residual_db_metrics(
        model,
        target,
        frequency,
        fit_args,
    )
    bands = base.band_summary(
        model,
        target,
        frequency,
        fit_args,
    )
    return {
        "shape_score": score,
        "residual_metrics": residual,
        "bands": bands,
        "passes_screen_tolerance": bool(
            residual["rms_residual_dB"] <= float(rms_screen)
            and residual["max_abs_residual_dB"] <= float(max_screen)
        ),
    }


def shared_model_for_scale(
    context,
    scale,
    shared_transfer,
):
    trial_context = white_scaled_context(context, scale)
    return profile.pre_analysis_model(
        trial_context,
        shared_transfer["parameters"],
        shared_transfer["fixed_leadlag_pole_Hz"],
    )


def profile_white_floor(
    context,
    target,
    frequency,
    fit_args,
    shared_transfer,
    profile_cfg,
):
    baseline_white_asd = float(
        context["post_filter_white_asd_A_rtHz"]
    )
    if baseline_white_asd <= 0.0:
        raise ValueError(
            "baseline post-filter white ASD must be positive "
            "for multiplicative profiling"
        )

    min_scale = float(profile_cfg["min_positive_scale"])
    max_scale = float(profile_cfg["max_scale"])
    points = int(profile_cfg["grid_points"])
    if min_scale <= 0.0 or max_scale <= min_scale or points < 2:
        raise ValueError("invalid white-floor profile grid")

    def score_for_scale(scale):
        model = shared_model_for_scale(
            context,
            scale,
            shared_transfer,
        )
        return float(
            opt.fit_score(
                model,
                target,
                frequency,
                fit_args,
            )
        )

    scales = list(
        np.geomspace(min_scale, max_scale, points)
    )
    if bool(profile_cfg.get("include_exact_zero", True)):
        scales.append(0.0)
    if bool(profile_cfg.get("include_baseline_scale_one", True)):
        scales.append(1.0)
    scales = sorted(set(float(value) for value in scales))

    profile_rows = [
        {
            "white_scale": scale,
            "white_asd_A_rtHz": baseline_white_asd * scale,
            "shape_score": score_for_scale(scale),
        }
        for scale in scales
    ]

    optimized = minimize_scalar(
        lambda log10_scale: score_for_scale(
            10.0 ** float(log10_scale)
        ),
        bounds=(np.log10(min_scale), np.log10(max_scale)),
        method="bounded",
        options={"xatol": 1.0e-7, "maxiter": 300},
    )
    optimized_scale = float(10.0 ** optimized.x)
    candidates = list(profile_rows)
    candidates.append(
        {
            "white_scale": optimized_scale,
            "white_asd_A_rtHz": (
                baseline_white_asd * optimized_scale
            ),
            "shape_score": float(optimized.fun),
            "source": "bounded_scalar_optimization",
        }
    )
    best = min(candidates, key=lambda row: row["shape_score"])
    baseline_row = min(
        candidates,
        key=lambda row: abs(float(row["white_scale"]) - 1.0),
    )
    return {
        "baseline_white_asd_A_rtHz": baseline_white_asd,
        "grid": profile_rows,
        "optimizer_success": bool(optimized.success),
        "optimizer_message": str(optimized.message),
        "best": {
            **best,
            "score_ratio_to_scale_one": float(
                best["shape_score"]
                / baseline_row["shape_score"]
            ),
        },
        "scale_one_shape_score": float(
            baseline_row["shape_score"]
        ),
    }


def parameter_drift(reference, local):
    return shared.parameter_drift(reference, local)


def anchor_rows(
    frequency,
    target,
    stage0_model,
    stage1_model,
    stage2_model,
):
    rows = []
    for anchor in shared.ANCHOR_FREQUENCIES_HZ:
        index = int(np.argmin(np.abs(frequency - anchor)))
        rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "target_pre_analysis_normalized": float(
                    target[index]
                ),
                "stage0_shared_baseline_white": float(
                    stage0_model[index]
                ),
                "stage1_shared_profiled_white": float(
                    stage1_model[index]
                ),
                "stage2_profiled_white_local_transfer": float(
                    stage2_model[index]
                ),
            }
        )
    return rows


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
    fit_args = base.fit_args(summary)
    frequency = np.geomspace(
        fit_args.fit_min_hz,
        fit_args.fit_max_hz,
        fit_args.fit_points,
    )
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        frequency,
    )
    candidate = dict(summary["best_case_parameters"])
    context = base.intrinsic_context(
        candidate,
        frequency,
    )

    reference_transfer_path = resolve_config_path(
        config["reference_transfer"],
        config_path,
    )
    reference_payload = json.loads(
        reference_transfer_path.read_text(encoding="utf-8")
    )
    reference_transfer = shared.load_shared_transfer(
        reference_payload
    )

    previous_local_path = resolve_config_path(
        config["previous_repeat_local_transfer"],
        config_path,
    )
    previous_local = load_previous_local_transfer(
        previous_local_path
    )

    refit_cfg = config["transfer_refit"]
    rms_screen = float(refit_cfg["RMS_screen_dB"])
    max_screen = float(refit_cfg["max_screen_dB"])

    stage0_model = shared_model_for_scale(
        context,
        1.0,
        reference_transfer,
    )
    stage0 = model_summary(
        stage0_model,
        target_context["target"],
        frequency,
        fit_args,
        rms_screen=rms_screen,
        max_screen=max_screen,
    )

    white_profile = profile_white_floor(
        context,
        target_context["target"],
        frequency,
        fit_args,
        reference_transfer,
        config["white_floor_profile"],
    )
    best_scale = float(
        white_profile["best"]["white_scale"]
    )
    profiled_context = white_scaled_context(
        context,
        best_scale,
    )
    stage1_model = profile.pre_analysis_model(
        profiled_context,
        reference_transfer["parameters"],
        reference_transfer["fixed_leadlag_pole_Hz"],
    )
    stage1 = model_summary(
        stage1_model,
        target_context["target"],
        frequency,
        fit_args,
        rms_screen=rms_screen,
        max_screen=max_screen,
    )

    local = profile.fit_fixed_pole(
        reference_transfer["fixed_leadlag_pole_Hz"],
        profiled_context,
        target_context["target"],
        frequency,
        fit_args,
        center_min_hz=float(refit_cfg["center_min_Hz"]),
        center_max_hz=float(refit_cfg["center_max_Hz"]),
        general_q_min=float(refit_cfg["general_Q_min"]),
        q_max=float(refit_cfg["Q_max"]),
        seed=int(refit_cfg["seed"]),
        de_maxiter=int(refit_cfg["DE_maxiter"]),
        rms_screen_db=rms_screen,
        max_screen_db=max_screen,
        initial_parameters=reference_transfer["parameters"],
    )
    stage2_model = local["_model"]
    stage2 = {
        key: value
        for key, value in local.items()
        if not key.startswith("_")
    }

    score_tolerance = max(
        1.0e-12,
        abs(stage1["shape_score"]) * 1.0e-9,
    )
    if (
        float(stage2["shape_score"])
        > float(stage1["shape_score"]) + score_tolerance
    ):
        raise RuntimeError(
            "warm-started transfer refit is worse than the "
            "shared-transfer candidate at the same white floor"
        )

    previous_drift = parameter_drift(
        reference_transfer["parameters"],
        previous_local["parameters"],
    )
    new_drift = parameter_drift(
        reference_transfer["parameters"],
        local["_parameters_raw"],
    )

    shared_over_local_after_white = float(
        stage1["shape_score"] / stage2["shape_score"]
    )
    white_score_ratio = float(
        stage1["shape_score"] / stage0["shape_score"]
    )
    drift_ratio = (
        float(
            new_drift["rms_log10_parameter_ratio"]
            / previous_drift["rms_log10_parameter_ratio"]
        )
        if previous_drift["rms_log10_parameter_ratio"] > 0.0
        else None
    )

    tolerance = float(
        refit_cfg["shared_local_score_tolerance"]
    )
    white_material = bool(white_score_ratio <= 0.80)
    shared_close = bool(
        shared_over_local_after_white <= tolerance
    )
    drift_reduced = bool(
        new_drift["rms_log10_parameter_ratio"]
        < previous_drift["rms_log10_parameter_ratio"]
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
        "separation_order": [
            "stage0: shared transfer + inherited production white floor",
            "stage1: shared transfer fixed; profile only post-filter white ASD",
            "stage2: freeze profiled white ASD; refit only the five hybrid transfer parameters",
        ],
        "white_floor_semantics": (
            "diagnostic post-filter ASD amplitude added after analog "
            "filter/alias folding; absolute CH0 calibration remains unresolved"
        ),
        "stage0_shared_transfer_baseline_white": {
            **stage0,
            "white_scale": 1.0,
            "white_asd_A_rtHz": float(
                context["post_filter_white_asd_A_rtHz"]
            ),
        },
        "white_floor_profile": white_profile,
        "stage1_shared_transfer_profiled_white": {
            **stage1,
            "white_scale": best_scale,
            "white_asd_A_rtHz": float(
                profiled_context[
                    "post_filter_white_asd_A_rtHz"
                ]
            ),
            "score_ratio_to_stage0": white_score_ratio,
        },
        "stage2_profiled_white_transfer_refit": {
            **stage2,
            "shared_over_local_score_ratio": (
                shared_over_local_after_white
            ),
            "shared_within_local_score_tolerance": bool(
                shared_close
            ),
        },
        "transfer_drift_comparison": {
            "previous_local_snapshot": {
                "path": str(previous_local_path),
                "parameters": previous_local["parameters"],
                "drift_from_shared": previous_drift,
            },
            "after_white_profile": {
                "parameters": local["_parameters_raw"],
                "drift_from_shared": new_drift,
            },
            "rms_log10_drift_after_over_before": drift_ratio,
            "drift_reduced_after_white_profile": drift_reduced,
        },
        "anchor_models": anchor_rows(
            frequency,
            target_context["target"],
            stage0_model,
            stage1_model,
            stage2_model,
        ),
        "interpretation_flags": {
            "white_floor_materially_improves_shared_fit": white_material,
            "shared_transfer_close_to_local_after_white_profile": shared_close,
            "transfer_parameter_drift_reduced_after_white_profile": drift_reduced,
            "supports_white_floor_absorbing_previous_transfer_drift": bool(
                white_material
                and shared_close
                and drift_reduced
            ),
            "physical_white_source_identified": False,
            "electronics_transfer_drift_identified": False,
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "manifest": str(manifest_path),
            "reference_transfer": str(reference_transfer_path),
            "previous_repeat_local_transfer": str(
                previous_local_path
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
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(output),
                "accepted_records": result["repeat_case"][
                    "accepted_records"
                ],
                "best_white_scale": result[
                    "white_floor_profile"
                ]["best"]["white_scale"],
                "stage0_score": result[
                    "stage0_shared_transfer_baseline_white"
                ]["shape_score"],
                "stage1_score": result[
                    "stage1_shared_transfer_profiled_white"
                ]["shape_score"],
                "stage2_score": result[
                    "stage2_profiled_white_transfer_refit"
                ]["shape_score"],
                "shared_over_local_after_white": result[
                    "stage2_profiled_white_transfer_refit"
                ]["shared_over_local_score_ratio"],
                "drift_reduced_after_white": result[
                    "transfer_drift_comparison"
                ]["drift_reduced_after_white_profile"],
                "supports_white_floor_absorption": result[
                    "interpretation_flags"
                ][
                    "supports_white_floor_absorbing_previous_transfer_drift"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

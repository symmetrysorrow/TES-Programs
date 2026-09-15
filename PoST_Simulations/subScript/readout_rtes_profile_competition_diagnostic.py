"""Profile R_TES while refitting nested detector-state nuisance families.

The 2024-12-05 repeat target is fit on 1-40 kHz with:
  * the 2024-12-06 order-2 effective readout transfer fixed,
  * the profiled 2024-12-05 post-filter white ASD fixed,
  * R_TES fixed externally at each profile point,
  * nested detector-state nuisance parameters refit internally.

R_TES is never a simultaneous continuous optimizer coordinate.  This preserves
profile identifiability and directly tests the remaining DC operating-point
confound left by the detector-state competition diagnostic.

The 40-200 kHz region is a strict holdout.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "readout_rtes_profile_competition_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_rtes_profile_competition_diagnostic.json"
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


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def load_detector_state_snapshot(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    best = payload["best_detector_state_family"]
    return {
        "inherited_R_TES_Ohm": float(
            payload["provenance"]["inherited_R_TES_Ohm"]
        ),
        "fixed_reference_shape_score": float(
            payload["fixed_reference_readout_baseline"]["shape_score"]
        ),
        "fixed_reference_rms_dB": float(
            payload["fixed_reference_readout_baseline"]["rms_residual_dB"]
        ),
        "repeat_local_shape_score": float(
            payload["repeat_local_readout_comparator"]["shape_score"]
        ),
        "repeat_local_rms_dB": float(
            payload["repeat_local_readout_comparator"]["rms_residual_dB"]
        ),
        "best_family_name": str(best["name"]),
        "best_candidate": {
            key: float(value)
            for key, value in best["candidate"].items()
        },
        "best_shape_score": float(best["shape_score"]),
        "best_rms_dB": float(best["rms_residual_dB"]),
    }


def fit_args_for_band(original, min_hz, max_hz, points):
    return competition.band_args(
        original,
        min_hz,
        max_hz,
        points,
    )


def full_model(
    candidate,
    frequency,
    transfer,
    scale_hz,
    white_asd,
):
    context = competition.fixed_white_context(
        candidate,
        frequency,
        white_asd,
    )
    return effective.pre_analysis_model(
        context,
        transfer,
        2,
        scale_hz,
    )


def clean_row(row):
    return {
        key: value
        for key, value in row.items()
        if not key.startswith("_")
    }


def profile_point(
    *,
    ratio,
    inherited_candidate,
    families,
    reference_transfer,
    full_frequency,
    fit_frequency,
    fit_target,
    fit_args,
    hold_frequency,
    hold_target,
    hold_args,
    full_target,
    full_args,
    hold_mask,
    scale_hz,
    white_asd,
    physical_bounds,
    optimizer_cfg,
    seed,
    warm_candidates,
    repeat_local_score,
    repeat_local_rms,
    screen,
):
    profiled_candidate = dict(inherited_candidate)
    baseline_r = float(inherited_candidate["R"])
    profiled_r = baseline_r * float(ratio)
    profiled_candidate["R"] = profiled_r

    point = opt.tes_operating_point(profiled_candidate)
    base = {
        "R_ratio_to_inherited": float(ratio),
        "R_TES_Ohm": float(profiled_r),
        "R_TES_mOhm": float(profiled_r * 1.0e3),
        "baseline_operating_point": {
            "valid": bool(point.get("valid")),
            "stable": bool(point.get("stable")),
            "reason": point.get("reason"),
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
    }
    baseline_metrics = None
    baseline_model_error = None
    if point.get("valid") and point.get("stable"):
        try:
            baseline_model = full_model(
                profiled_candidate,
                fit_frequency,
                reference_transfer,
                scale_hz,
                white_asd,
            )
            baseline_metrics = holdout.model_metrics(
                baseline_model,
                fit_target,
                fit_frequency,
                fit_args,
            )
        except Exception as exc:
            baseline_model_error = str(exc)

    family_rows = []
    previous = list(warm_candidates)

    for family_index, family in enumerate(families):
        try:
            row = competition.fit_nuisance_family(
                name=family["name"],
                parameter_names=family["parameters"],
                baseline_candidate=profiled_candidate,
                fixed_transfer=reference_transfer,
                frequency=fit_frequency,
                target=fit_target,
                args=fit_args,
                scale_hz=scale_hz,
                white_asd=white_asd,
                physical_bounds=physical_bounds,
                optimizer_cfg=optimizer_cfg,
                seed=int(seed) + family_index,
                warm_candidates=previous,
                allow_unstable_baseline=True,
            )
        except Exception as exc:
            family_rows.append(
                {
                    "name": family["name"],
                    "parameters_varied": list(family["parameters"]),
                    "status": "fit_error",
                    "error": str(exc),
                }
            )
            continue

        best_full_model = full_model(
            row["_candidate_full"],
            full_frequency,
            reference_transfer,
            scale_hz,
            white_asd,
        )
        row["holdout_metrics"] = holdout.model_metrics(
            best_full_model[hold_mask],
            hold_target,
            hold_frequency,
            hold_args,
        )
        row["full_1_200k_metrics"] = holdout.model_metrics(
            best_full_model,
            full_target,
            full_frequency,
            full_args,
        )
        row["score_ratio_to_repeat_local_readout_comparator"] = float(
            row["shape_score"] / repeat_local_score
        )
        row["rms_delta_to_repeat_local_readout_comparator_dB"] = float(
            row["residual_metrics"]["rms_residual_dB"]
            - repeat_local_rms
        )
        row["near_repeat_local_readout_comparator"] = bool(
            row["score_ratio_to_repeat_local_readout_comparator"]
            <= float(screen["near_repeat_local_score_ratio"])
            and row["rms_delta_to_repeat_local_readout_comparator_dB"]
            <= float(screen["near_repeat_local_rms_delta_dB"])
        )
        row["status"] = "evaluated"
        family_rows.append(clean_row(row))
        previous.append(
            (
                family["name"],
                row["_candidate_full"],
            )
        )

    evaluated = [
        row
        for row in family_rows
        if row.get("status") == "evaluated"
    ]
    best_family = (
        min(
            evaluated,
            key=lambda row: row["shape_score"],
        )
        if evaluated
        else None
    )
    base.update(
        {
            "status": (
                "evaluated"
                if evaluated
                else "no_stable_solution_at_fixed_R"
            ),
            "baseline_model_error": baseline_model_error,
            "inherited_nuisance_state_was_stable": bool(
                point.get("valid") and point.get("stable")
            ),
            "fixed_reference_readout_inherited_nuisance_metrics": (
                baseline_metrics
            ),
            "family_fits": family_rows,
            "best_family": best_family,
        }
    )
    return base


def run(config, config_path: Path):
    if int(config["fixed_readout"]["effective_numerator_order"]) != 2:
        raise ValueError("R_TES profile diagnostic requires order-2 fixed readout")
    for family in config["nuisance_families"]:
        if "R" in family["parameters"]:
            raise ValueError("R_TES must remain an outer fixed profile coordinate")

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

    day_snapshot_path = resolve_config_path(
        config["order2_day_snapshot"],
        config_path,
    )
    day_snapshot = competition.load_day_snapshot(
        day_snapshot_path
    )
    detector_snapshot_path = resolve_config_path(
        config["detector_state_snapshot"],
        config_path,
    )
    detector_snapshot = load_detector_state_snapshot(
        detector_snapshot_path
    )

    inherited_candidate = dict(
        summary["best_case_parameters"]
    )
    inherited_r = float(inherited_candidate["R"])
    if not np.isclose(
        inherited_r,
        detector_snapshot["inherited_R_TES_Ohm"],
        rtol=1.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "detector-state snapshot R_TES does not match repeat detector snapshot"
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
        raise RuntimeError("fit and holdout masks overlap")

    fit_frequency = full_frequency[fit_mask]
    fit_target = full_target[fit_mask]
    fit_args = fit_args_for_band(
        full_args,
        fit_min,
        fit_max,
        np.count_nonzero(fit_mask),
    )
    hold_frequency = full_frequency[hold_mask]
    hold_target = full_target[hold_mask]
    hold_args = fit_args_for_band(
        full_args,
        hold_min,
        hold_max,
        np.count_nonzero(hold_mask),
    )

    reference_transfer = day_snapshot[
        "reference_transfer"
    ]
    repeat_transfer = day_snapshot[
        "repeat_transfer"
    ]
    scale_hz = float(
        config["fixed_readout"]["reference_scale_Hz"]
    )
    white_asd = float(
        day_snapshot["repeat_white_asd_A_rtHz"]
    )

    repeat_local_model = full_model(
        inherited_candidate,
        fit_frequency,
        repeat_transfer,
        scale_hz,
        white_asd,
    )
    repeat_local_metrics = holdout.model_metrics(
        repeat_local_model,
        fit_target,
        fit_frequency,
        fit_args,
    )
    if not np.isclose(
        repeat_local_metrics["shape_score"],
        detector_snapshot["repeat_local_shape_score"],
        rtol=2.0e-4,
        atol=1.0e-12,
    ):
        raise ValueError(
            "repeat-local order-2 comparator does not reproduce"
        )

    physical_bounds = competition.nuisance_bounds(
        inherited_candidate,
        {
            "bounds": config["bounds"],
        },
    )
    optimizer_cfg = config["optimizer"]
    families = config["nuisance_families"]
    profile_ratios = [
        float(value)
        for value in config["R_TES_profile"]["ratios"]
    ]
    if any(value <= 0.0 for value in profile_ratios):
        raise ValueError("R_TES profile ratios must be positive")
    if len(set(profile_ratios)) != len(profile_ratios):
        raise ValueError("R_TES profile ratios must be unique")
    if 1.0 not in profile_ratios:
        raise ValueError("R_TES profile must include ratio 1.0")

    inherited_point = opt.tes_operating_point(
        inherited_candidate
    )
    if not inherited_point.get("valid") or not inherited_point.get("stable"):
        raise ValueError(
            "inherited detector operating point is invalid or unstable"
        )

    warm_anchor = dict(inherited_candidate)
    warm_anchor.update(
        detector_snapshot["best_candidate"]
    )
    rows = []
    for index, ratio in enumerate(profile_ratios):
        rows.append(
            profile_point(
                ratio=ratio,
                inherited_candidate=inherited_candidate,
                families=families,
                reference_transfer=reference_transfer,
                full_frequency=full_frequency,
                fit_frequency=fit_frequency,
                fit_target=fit_target,
                fit_args=fit_args,
                hold_frequency=hold_frequency,
                hold_target=hold_target,
                hold_args=hold_args,
                full_target=full_target,
                full_args=full_args,
                hold_mask=hold_mask,
                scale_hz=scale_hz,
                white_asd=white_asd,
                physical_bounds=physical_bounds,
                optimizer_cfg=optimizer_cfg,
                seed=int(optimizer_cfg["seed"])
                + 100 * index,
                warm_candidates=[
                    (
                        "tracked_R1_detector_state_best",
                        warm_anchor,
                    )
                ],
                repeat_local_score=repeat_local_metrics[
                    "shape_score"
                ],
                repeat_local_rms=repeat_local_metrics[
                    "residual_metrics"
                ]["rms_residual_dB"],
                screen=config["comparison_screen"],
            )
        )

    evaluated = [
        row
        for row in rows
        if row.get("status") == "evaluated"
        and row.get("best_family") is not None
    ]
    if not evaluated:
        raise RuntimeError(
            "no R_TES profile point produced a stable nuisance fit"
        )
    best = min(
        evaluated,
        key=lambda row: row["best_family"]["shape_score"],
    )
    r1 = next(
        row
        for row in evaluated
        if np.isclose(
            row["R_ratio_to_inherited"],
            1.0,
            rtol=0.0,
            atol=1.0e-12,
        )
    )
    tracked_anchor_ratio = float(
        r1["best_family"]["shape_score"]
        / detector_snapshot["best_shape_score"]
    )
    if tracked_anchor_ratio > 1.0005:
        raise RuntimeError(
            "R/R0=1 profile failed to reproduce the tracked detector-state optimum"
        )

    screen = config["comparison_screen"]
    modest_min = float(
        screen["modest_R_shift_ratio_min"]
    )
    modest_max = float(
        screen["modest_R_shift_ratio_max"]
    )
    modest = [
        row
        for row in evaluated
        if modest_min
        <= row["R_ratio_to_inherited"]
        <= modest_max
    ]
    best_modest = min(
        modest,
        key=lambda row: row["best_family"]["shape_score"],
    )

    for row in evaluated:
        row["best_family"]["score_ratio_to_R1_best_detector_nuisance"] = float(
            row["best_family"]["shape_score"]
            / r1["best_family"]["shape_score"]
        )
        row["best_family"]["materially_better_than_R1_detector_nuisance"] = bool(
            row["best_family"][
                "score_ratio_to_R1_best_detector_nuisance"
            ]
            <= float(
                screen[
                    "materially_better_than_R1_detector_nuisance_score_ratio"
                ]
            )
        )

    stable_ratios = [
        float(row["R_ratio_to_inherited"])
        for row in evaluated
    ]
    unresolved = [
        row for row in rows
        if row.get("status") != "evaluated"
    ]
    any_near_local = any(
        row["best_family"][
            "near_repeat_local_readout_comparator"
        ]
        for row in evaluated
    )
    modest_near_local = any(
        row["best_family"][
            "near_repeat_local_readout_comparator"
        ]
        for row in modest
    )
    best_material = bool(
        best["best_family"][
            "materially_better_than_R1_detector_nuisance"
        ]
    )
    best_at_profile_edge = bool(
        np.isclose(
            best["R_ratio_to_inherited"],
            min(profile_ratios),
        )
        or np.isclose(
            best["R_ratio_to_inherited"],
            max(profile_ratios),
        )
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
            "points": int(np.count_nonzero(fit_mask)),
        },
        "holdout_region": {
            "min_Hz": hold_min,
            "max_Hz": hold_max,
            "points": int(np.count_nonzero(hold_mask)),
            "optimizer_received_holdout_points": False,
        },
        "fixed_readout": {
            "order": 2,
            "reference_scale_Hz": scale_hz,
            "parameters": reference_transfer,
            "canonical": day_snapshot[
                "reference_canonical"
            ],
        },
        "fixed_repeat_white_floor": {
            "white_asd_A_rtHz": white_asd,
            "recomputed_for_profile_trials": False,
        },
        "repeat_local_readout_comparator": {
            "shape_score": repeat_local_metrics[
                "shape_score"
            ],
            "rms_residual_dB": repeat_local_metrics[
                "residual_metrics"
            ]["rms_residual_dB"],
        },
        "inherited_detector": {
            "R_TES_Ohm": inherited_r,
            "R_TES_mOhm": inherited_r * 1.0e3,
            "current_A": float(
                inherited_point["current_A"]
            ),
            "joule_power_W": float(
                inherited_point["joule_power_W"]
            ),
        },
        "R_TES_profile": {
            "ratios": profile_ratios,
            "rows": rows,
            "stability_summary": {
                "n_profile_points": int(len(rows)),
                "n_points_with_stable_fitted_solution": int(len(evaluated)),
                "n_points_without_stable_fitted_solution": int(len(unresolved)),
                "stable_solution_R_ratio_min": float(min(stable_ratios)),
                "stable_solution_R_ratio_max": float(max(stable_ratios)),
                "low_R_below_inherited_stable_solution_found": bool(
                    any(value < 1.0 for value in stable_ratios)
                ),
                "unresolved_R_ratios": [
                    float(row["R_ratio_to_inherited"])
                    for row in unresolved
                ],
            },
            "R1_regression_anchor": {
                "tracked_best_shape_score": detector_snapshot[
                    "best_shape_score"
                ],
                "recomputed_profile_best_shape_score": r1[
                    "best_family"
                ]["shape_score"],
                "recomputed_over_tracked_score_ratio": (
                    tracked_anchor_ratio
                ),
            },
        },
        "best_overall_profile_point": {
            "R_ratio_to_inherited": best[
                "R_ratio_to_inherited"
            ],
            "R_TES_Ohm": best["R_TES_Ohm"],
            "R_TES_mOhm": best["R_TES_mOhm"],
            "best_family": best["best_family"],
        },
        "best_modest_R_shift_profile_point": {
            "R_ratio_to_inherited": best_modest[
                "R_ratio_to_inherited"
            ],
            "R_TES_Ohm": best_modest["R_TES_Ohm"],
            "R_TES_mOhm": best_modest["R_TES_mOhm"],
            "best_family": best_modest["best_family"],
        },
        "comparison_screen": screen,
        "interpretation_flags": {
            "some_R_profile_point_reaches_near_repeat_local_readout_fit": bool(
                any_near_local
            ),
            "modest_R_shift_reaches_near_repeat_local_readout_fit": bool(
                modest_near_local
            ),
            "R_profile_materially_improves_over_R1_detector_nuisance": (
                best_material
            ),
            "best_R_profile_point_is_at_tested_edge": (
                best_at_profile_edge
            ),
            "low_R_stable_nuisance_solution_found": bool(
                any(value < 1.0 for value in stable_ratios)
            ),
            "all_tested_R_points_have_stable_nuisance_solution": bool(
                len(unresolved) == 0
            ),
            "profiled_R_TES_can_compete_with_day_specific_readout_shape": bool(
                any_near_local
            ),
            "day_specific_readout_shape_still_required_within_profiled_R_and_tested_nuisance_set": bool(
                not any_near_local
            ),
            "physical_R_TES_shift_identified": False,
            "physical_electronics_drift_identified": False,
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "manifest": str(manifest_path),
            "order2_day_snapshot": str(
                day_snapshot_path
            ),
            "detector_state_snapshot": str(
                detector_snapshot_path
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
                "R1_regression_anchor": result[
                    "R_TES_profile"
                ]["R1_regression_anchor"],
                "best_overall": {
                    "R_ratio": result[
                        "best_overall_profile_point"
                    ]["R_ratio_to_inherited"],
                    "R_TES_mOhm": result[
                        "best_overall_profile_point"
                    ]["R_TES_mOhm"],
                    "family": result[
                        "best_overall_profile_point"
                    ]["best_family"]["name"],
                    "shape_score": result[
                        "best_overall_profile_point"
                    ]["best_family"]["shape_score"],
                    "rms_residual_dB": result[
                        "best_overall_profile_point"
                    ]["best_family"][
                        "residual_metrics"
                    ]["rms_residual_dB"],
                    "near_repeat_local": result[
                        "best_overall_profile_point"
                    ]["best_family"][
                        "near_repeat_local_readout_comparator"
                    ],
                },
                "best_modest_R_shift": {
                    "R_ratio": result[
                        "best_modest_R_shift_profile_point"
                    ]["R_ratio_to_inherited"],
                    "shape_score": result[
                        "best_modest_R_shift_profile_point"
                    ]["best_family"]["shape_score"],
                    "near_repeat_local": result[
                        "best_modest_R_shift_profile_point"
                    ]["best_family"][
                        "near_repeat_local_readout_comparator"
                    ],
                },
                "flags": result[
                    "interpretation_flags"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

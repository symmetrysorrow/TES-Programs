"""Fit repeat-case readout transfer on 1-40 kHz and validate 40-200 kHz as holdout.

The post-filter white floor is frozen to the Git-tracked best value from the
preceding white-floor separation diagnostic.  Only 1-40 kHz samples are passed
to the hybrid-transfer optimizer.  The fitted parameters are then evaluated on
40-200 kHz without any further adjustment.

This tests whether the full-band repeat fit's high-frequency parameter shifts
(particularly the ~75 kHz zero) are already demanded by low/mid-frequency data
or instead arise from fitting the high-frequency region itself.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "readout_lowmid_holdout_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_lowmid_holdout_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import preanalysis_hybrid_pole_profile_diagnostic as profile  # noqa: E402
from subScript import readout_white_floor_separation_diagnostic as white_sep  # noqa: E402
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


def load_full_band_snapshot(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    white = payload["white_floor"]
    transfer = payload["full_band_transfer_fit"]
    parameters = {
        name: float(transfer["parameters"][name])
        for name in shared.TRANSFER_PARAMETER_NAMES
    }
    fixed_pole = transfer.get("fixed_leadlag_pole_Hz")
    if transfer.get("leadlag_pole_is_infinite", fixed_pole is None):
        fixed_pole = None
    elif fixed_pole is None:
        raise ValueError(
            "finite full-band lead/lag pole requires a value"
        )
    else:
        fixed_pole = float(fixed_pole)
    return {
        "white_scale": float(white["best_white_scale"]),
        "white_asd_A_rtHz": float(white["best_white_asd_A_rtHz"]),
        "parameters": parameters,
        "fixed_leadlag_pole_Hz": fixed_pole,
        "expected_full_band_shape_score": float(
            payload["provenance"][
                "stage2_full_band_local_shape_score"
            ]
        ),
        "provenance": payload.get("provenance", {}),
    }


def subset_context(context, mask):
    mask = np.asarray(mask, dtype=bool)
    full_frequency = np.asarray(
        context["frequency_Hz"],
        dtype=float,
    )
    if mask.shape != full_frequency.shape:
        raise ValueError("context subset mask shape mismatch")
    result = {}
    for key, value in context.items():
        if isinstance(value, np.ndarray) and value.shape == mask.shape:
            result[key] = value[mask]
        else:
            result[key] = value
    return result


def band_args(original, min_hz, max_hz, points):
    values = dict(vars(original))
    values["fit_min_hz"] = float(min_hz)
    values["fit_max_hz"] = float(max_hz)
    values["fit_points"] = int(points)
    return SimpleNamespace(**values)


def model_metrics(model, target, frequency, args):
    return {
        "shape_score": float(
            opt.fit_score(model, target, frequency, args)
        ),
        "residual_metrics": base.residual_db_metrics(
            model,
            target,
            frequency,
            args,
        ),
        "bands": base.band_summary(
            model,
            target,
            frequency,
            args,
        ),
    }


def drift_toward_shared(reference, full_band, low_mid):
    rows = {}
    closer_count = 0
    for name in shared.TRANSFER_PARAMETER_NAMES:
        ref = float(reference[name])
        full = float(full_band[name])
        low = float(low_mid[name])
        full_distance = abs(float(np.log10(full / ref)))
        low_distance = abs(float(np.log10(low / ref)))
        closer = bool(low_distance < full_distance)
        closer_count += int(closer)
        rows[name] = {
            "shared_value": ref,
            "full_band_value": full,
            "low_mid_value": low,
            "full_band_abs_log10_distance_from_shared": full_distance,
            "low_mid_abs_log10_distance_from_shared": low_distance,
            "low_mid_is_closer_to_shared": closer,
        }
    return {
        "parameters": rows,
        "n_parameters_closer_to_shared": int(closer_count),
        "fraction_parameters_closer_to_shared": float(
            closer_count / len(shared.TRANSFER_PARAMETER_NAMES)
        ),
    }


def anchor_rows(
    frequency,
    target,
    shared_model,
    full_band_model,
    low_mid_model,
):
    rows = []
    for anchor in shared.ANCHOR_FREQUENCIES_HZ:
        index = int(np.argmin(np.abs(frequency - anchor)))
        rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "region": (
                    "fit"
                    if frequency[index] <= 40000.0
                    else "holdout"
                ),
                "target_pre_analysis_normalized": float(
                    target[index]
                ),
                "shared_transfer_model": float(
                    shared_model[index]
                ),
                "full_band_local_model": float(
                    full_band_model[index]
                ),
                "low_mid_fit_model": float(
                    low_mid_model[index]
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
    frequency = np.geomspace(
        full_args.fit_min_hz,
        full_args.fit_max_hz,
        full_args.fit_points,
    )
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        frequency,
    )
    target = target_context["target"]
    candidate = dict(summary["best_case_parameters"])
    base_context = base.intrinsic_context(
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

    snapshot_path = resolve_config_path(
        config["white_profiled_full_band_snapshot"],
        config_path,
    )
    snapshot = load_full_band_snapshot(snapshot_path)
    context = white_sep.white_scaled_context(
        base_context,
        snapshot["white_scale"],
    )
    computed_white = float(
        context["post_filter_white_asd_A_rtHz"]
    )
    if not np.isclose(
        computed_white,
        snapshot["white_asd_A_rtHz"],
        rtol=5.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "white-profile snapshot does not match detector-snapshot "
            "baseline white ASD"
        )

    shared_full_model = profile.pre_analysis_model(
        context,
        reference_transfer["parameters"],
        reference_transfer["fixed_leadlag_pole_Hz"],
    )
    full_band_model = profile.pre_analysis_model(
        context,
        snapshot["parameters"],
        snapshot["fixed_leadlag_pole_Hz"],
    )

    fit_cfg = config["low_mid_fit"]
    fit_min = float(fit_cfg["min_Hz"])
    fit_max = float(fit_cfg["max_Hz"])
    fit_mask = (
        (frequency >= fit_min)
        & (frequency < fit_max)
    )
    if np.count_nonzero(fit_mask) < 10:
        raise ValueError("low/mid fit mask contains too few points")
    fit_frequency = frequency[fit_mask]
    fit_target = target[fit_mask]
    fit_context = subset_context(context, fit_mask)
    fit_args = band_args(
        full_args,
        fit_min,
        fit_max,
        np.count_nonzero(fit_mask),
    )

    screen_cfg = config["screen"]
    low_mid_fit = profile.fit_fixed_pole(
        reference_transfer["fixed_leadlag_pole_Hz"],
        fit_context,
        fit_target,
        fit_frequency,
        fit_args,
        center_min_hz=float(fit_cfg["center_min_Hz"]),
        center_max_hz=float(fit_cfg["center_max_Hz"]),
        general_q_min=float(fit_cfg["general_Q_min"]),
        q_max=float(fit_cfg["Q_max"]),
        seed=int(fit_cfg["seed"]),
        de_maxiter=int(fit_cfg["DE_maxiter"]),
        rms_screen_db=float(screen_cfg["RMS_dB"]),
        max_screen_db=float(screen_cfg["max_abs_dB"]),
        initial_parameters=reference_transfer["parameters"],
    )
    low_mid_parameters = low_mid_fit["_parameters_raw"]
    low_mid_full_model = profile.pre_analysis_model(
        context,
        low_mid_parameters,
        reference_transfer["fixed_leadlag_pole_Hz"],
    )

    low_mid_shared_metrics = model_metrics(
        shared_full_model[fit_mask],
        fit_target,
        fit_frequency,
        fit_args,
    )
    low_mid_full_band_metrics = model_metrics(
        full_band_model[fit_mask],
        fit_target,
        fit_frequency,
        fit_args,
    )
    low_mid_fitted_metrics = model_metrics(
        low_mid_full_model[fit_mask],
        fit_target,
        fit_frequency,
        fit_args,
    )

    hold_cfg = config["high_frequency_holdout"]
    hold_min = float(hold_cfg["min_Hz"])
    hold_max = float(hold_cfg["max_Hz"])
    hold_mask = (
        (frequency >= hold_min)
        & (frequency <= hold_max)
    )
    if np.any(fit_mask & hold_mask):
        raise RuntimeError(
            "fit and holdout masks overlap; strict holdout violated"
        )
    if np.count_nonzero(hold_mask) < 10:
        raise ValueError("holdout mask contains too few points")
    hold_frequency = frequency[hold_mask]
    hold_target = target[hold_mask]
    hold_args = band_args(
        full_args,
        hold_min,
        hold_max,
        np.count_nonzero(hold_mask),
    )
    hold_shared = model_metrics(
        shared_full_model[hold_mask],
        hold_target,
        hold_frequency,
        hold_args,
    )
    hold_full_band = model_metrics(
        full_band_model[hold_mask],
        hold_target,
        hold_frequency,
        hold_args,
    )
    hold_low_mid = model_metrics(
        low_mid_full_model[hold_mask],
        hold_target,
        hold_frequency,
        hold_args,
    )

    full_band_computed = model_metrics(
        full_band_model,
        target,
        frequency,
        full_args,
    )
    low_mid_full_metrics = model_metrics(
        low_mid_full_model,
        target,
        frequency,
        full_args,
    )

    shared_drift_full = shared.parameter_drift(
        reference_transfer["parameters"],
        snapshot["parameters"],
    )
    shared_drift_low_mid = shared.parameter_drift(
        reference_transfer["parameters"],
        low_mid_parameters,
    )
    movement = drift_toward_shared(
        reference_transfer["parameters"],
        snapshot["parameters"],
        low_mid_parameters,
    )

    fit_score_ratio_to_shared = float(
        low_mid_fitted_metrics["shape_score"]
        / low_mid_shared_metrics["shape_score"]
    )
    holdout_score_ratio_to_shared = float(
        hold_low_mid["shape_score"]
        / hold_shared["shape_score"]
    )
    holdout_score_ratio_to_full_band = float(
        hold_low_mid["shape_score"]
        / hold_full_band["shape_score"]
    )
    drift_ratio = float(
        shared_drift_low_mid["rms_log10_parameter_ratio"]
        / shared_drift_full["rms_log10_parameter_ratio"]
    )

    tolerance = float(
        hold_cfg["relative_score_tolerance_to_shared"]
    )
    lowmid_needs_local = bool(fit_score_ratio_to_shared <= 0.80)
    holdout_not_worse_than_shared = bool(
        holdout_score_ratio_to_shared <= tolerance
    )
    zero_returns = bool(
        movement["parameters"]["zero_Hz"][
            "low_mid_is_closer_to_shared"
        ]
    )
    overall_returns = bool(drift_ratio < 1.0)

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "strict_holdout_semantics": True,
        "repeat_case": {
            "label": repeat_case["label"],
            "comparison_source": comparison_source,
            "experiment_path": str(experiment_path),
            "accepted_records": int(
                target_context["accepted_records"]
            ),
        },
        "white_floor_fixed": {
            "scale": snapshot["white_scale"],
            "white_asd_A_rtHz": snapshot["white_asd_A_rtHz"],
            "source_snapshot": str(snapshot_path),
        },
        "fit_region": {
            "min_Hz": fit_min,
            "max_Hz": fit_max,
            "points": int(np.count_nonzero(fit_mask)),
            "optimizer_received_holdout_points": False,
        },
        "holdout_region": {
            "min_Hz": hold_min,
            "max_Hz": hold_max,
            "points": int(np.count_nonzero(hold_mask)),
            "evaluated_only_after_fit": True,
        },
        "reference_shared_transfer": reference_transfer["parameters"],
        "full_band_white_profiled_snapshot": {
            "parameters": snapshot["parameters"],
            "expected_shape_score": snapshot[
                "expected_full_band_shape_score"
            ],
            "recomputed_shape_score": full_band_computed[
                "shape_score"
            ],
            "drift_from_shared": shared_drift_full,
        },
        "low_mid_only_transfer_fit": {
            **{
                key: value
                for key, value in low_mid_fit.items()
                if not key.startswith("_")
            },
            "parameters_raw": low_mid_parameters,
            "drift_from_shared": shared_drift_low_mid,
            "fit_score_ratio_to_shared": fit_score_ratio_to_shared,
        },
        "low_mid_region_metrics": {
            "shared_transfer": low_mid_shared_metrics,
            "full_band_local_transfer": low_mid_full_band_metrics,
            "low_mid_fit_transfer": low_mid_fitted_metrics,
        },
        "high_frequency_holdout_metrics": {
            "shared_transfer": hold_shared,
            "full_band_local_transfer": hold_full_band,
            "low_mid_fit_transfer": hold_low_mid,
            "low_mid_over_shared_score_ratio": (
                holdout_score_ratio_to_shared
            ),
            "low_mid_over_full_band_score_ratio": (
                holdout_score_ratio_to_full_band
            ),
        },
        "full_1_200k_metrics_for_low_mid_fit": (
            low_mid_full_metrics
        ),
        "parameter_return_toward_shared": {
            **movement,
            "rms_log10_drift_low_mid_over_full_band": drift_ratio,
        },
        "anchors": anchor_rows(
            frequency,
            target,
            shared_full_model,
            full_band_model,
            low_mid_full_model,
        ),
        "interpretation_flags": {
            "low_mid_data_require_case_specific_transfer": (
                lowmid_needs_local
            ),
            "low_mid_fit_holdout_not_materially_worse_than_shared": (
                holdout_not_worse_than_shared
            ),
            "zero_frequency_returns_toward_shared_when_high_is_held_out": (
                zero_returns
            ),
            "overall_transfer_drift_reduced_when_high_is_held_out": (
                overall_returns
            ),
            "supports_full_band_zero_drift_being_high_frequency_nuisance": bool(
                zero_returns and overall_returns
            ),
            "physical_readout_component_identified": False,
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "manifest": str(manifest_path),
            "reference_transfer": str(reference_transfer_path),
            "white_profiled_full_band_snapshot": str(
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
                "low_mid_parameters": result[
                    "low_mid_only_transfer_fit"
                ]["parameters_raw"],
                "low_mid_score_ratio_to_shared": result[
                    "low_mid_only_transfer_fit"
                ]["fit_score_ratio_to_shared"],
                "holdout_low_mid_over_shared_score_ratio": result[
                    "high_frequency_holdout_metrics"
                ]["low_mid_over_shared_score_ratio"],
                "drift_low_mid_over_full_band": result[
                    "parameter_return_toward_shared"
                ]["rms_log10_drift_low_mid_over_full_band"],
                "zero_returns_toward_shared": result[
                    "interpretation_flags"
                ][
                    "zero_frequency_returns_toward_shared_when_high_is_held_out"
                ],
                "supports_high_frequency_zero_nuisance": result[
                    "interpretation_flags"
                ][
                    "supports_full_band_zero_drift_being_high_frequency_nuisance"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

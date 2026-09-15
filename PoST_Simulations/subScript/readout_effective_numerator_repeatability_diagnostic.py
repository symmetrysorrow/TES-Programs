"""Order-2 effective-numerator block stability and adjacent-day repeatability.

The 2024-12-05 repeat is split into deterministic non-overlapping blocks with
the same accepted-record count as the 2024-12-06 reference.  Each block and
each full-day target is fit independently with the four-parameter order-2
effective-numerator model on 1-40 kHz:

    P(x) = (1 + u x^2)^2 + v^2 x^2
         = 1 + c2 x^2 + c4 x^4

Only the canonical shape parameters pole_Hz, pole_Q, c2, and c4 are used for
cross-block/day comparisons.  Latent u/v are optimization coordinates only.

The 40-200 kHz region remains a strict holdout for every fit.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = (
    CONFIG_DIR / "readout_effective_numerator_repeatability_config.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_effective_numerator_repeatability_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import readout_day_to_day_drift_diagnostic as day_drift  # noqa: E402
from subScript import readout_effective_numerator_diagnostic as effective  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import readout_lowmid_holdout_diagnostic as holdout  # noqa: E402
from subScript import readout_white_floor_separation_diagnostic as white_sep  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


CANONICAL_PARAMETER_NAMES = (
    "pole_Hz",
    "pole_Q",
    "c2",
    "c4",
)


def comparison_with_indices(comparison, accepted_indices):
    result = copy.deepcopy(comparison)
    acquisition = result["acquisition"]
    acquisition["accepted_record_indices"] = [
        int(value) for value in accepted_indices
    ]
    return result


def canonical_parameters(fit_result):
    parameters = fit_result["parameters"]
    coefficients = parameters[
        "polynomial_coefficients_dimensionless"
    ]
    return {
        "pole_Hz": float(parameters["pole_Hz"]),
        "pole_Q": float(parameters["pole_Q"]),
        "c2": float(coefficients["c2"]),
        "c4": float(coefficients["c4"]),
    }


def normalized_transfer(frequency, fit_result, scale_hz):
    raw = fit_result["_parameters_raw"]
    magnitude = effective.effective_transfer_magnitude(
        frequency,
        raw["pole_Hz"],
        raw["pole_Q"],
        2,
        raw["latent"],
        scale_hz,
    )
    return opt.normalize_at(
        frequency,
        magnitude,
        reference_hz=1000.0,
    )


def fit_dataset(
    *,
    label,
    summary,
    comparison,
    experiment_path,
    white_scale,
    accepted_indices,
    config,
    reference_parameters,
    free_reference_parameters,
    seed,
):
    full_args = base.fit_args(summary)
    frequency = np.geomspace(
        full_args.fit_min_hz,
        full_args.fit_max_hz,
        full_args.fit_points,
    )
    selected_comparison = (
        comparison
        if accepted_indices is None
        else comparison_with_indices(
            comparison,
            accepted_indices,
        )
    )
    target_context = base.reconstruct_pre_analysis_target(
        selected_comparison,
        experiment_path,
        frequency,
    )
    candidate = dict(summary["best_case_parameters"])
    context = base.intrinsic_context(
        candidate,
        frequency,
    )
    context = white_sep.white_scaled_context(
        context,
        white_scale,
    )

    fit_cfg = config["fit_region_Hz"]
    fit_min = float(fit_cfg["min"])
    fit_max = float(fit_cfg["max_exclusive"])
    fit_mask = (
        (frequency >= fit_min)
        & (frequency < fit_max)
    )
    hold_cfg = config["holdout_region_Hz"]
    hold_min = float(hold_cfg["min"])
    hold_max = float(hold_cfg["max"])
    hold_mask = (
        (frequency >= hold_min)
        & (frequency <= hold_max)
    )
    if np.any(fit_mask & hold_mask):
        raise RuntimeError(
            "fit and holdout masks overlap"
        )

    fit_frequency = frequency[fit_mask]
    fit_target = target_context["target"][fit_mask]
    fit_context = holdout.subset_context(
        context,
        fit_mask,
    )
    fit_args = ident.band_args(
        full_args,
        fit_min,
        fit_max,
        np.count_nonzero(fit_mask),
    )

    optimizer = config["optimizer"]
    scale_hz = float(config["reference_scale_Hz"])
    fit_result = effective.fit_order(
        2,
        fit_context,
        fit_target,
        fit_frequency,
        fit_args,
        scale_hz=scale_hz,
        optimizer=optimizer,
        shared_parameters=reference_parameters,
        free_parameters=free_reference_parameters,
        seed=int(seed),
    )

    full_model = effective.pre_analysis_model(
        context,
        fit_result["_parameters_raw"],
        2,
        scale_hz,
    )
    hold_frequency = frequency[hold_mask]
    hold_target = target_context["target"][hold_mask]
    hold_args = ident.band_args(
        full_args,
        hold_min,
        hold_max,
        np.count_nonzero(hold_mask),
    )
    fit_metrics = holdout.model_metrics(
        full_model[fit_mask],
        fit_target,
        fit_frequency,
        fit_args,
    )
    holdout_metrics = holdout.model_metrics(
        full_model[hold_mask],
        hold_target,
        hold_frequency,
        hold_args,
    )
    full_metrics = holdout.model_metrics(
        full_model,
        target_context["target"],
        frequency,
        full_args,
    )

    return {
        "label": label,
        "accepted_records": int(
            target_context["accepted_records"]
        ),
        "white_scale": float(white_scale),
        "white_asd_A_rtHz": float(
            context["post_filter_white_asd_A_rtHz"]
        ),
        "canonical_parameters": canonical_parameters(
            fit_result
        ),
        "optimization_parameters": {
            "u": float(
                fit_result["_parameters_raw"]["latent"]["u"]
            ),
            "v": float(
                fit_result["_parameters_raw"]["latent"]["v"]
            ),
        },
        "best_candidate_source": fit_result[
            "best_candidate_source"
        ],
        "evaluations": int(fit_result["evaluations"]),
        "fit_metrics": fit_metrics,
        "holdout_metrics": holdout_metrics,
        "full_1_200k_metrics": full_metrics,
        "_frequency_Hz": frequency,
        "_fit_mask": fit_mask,
        "_holdout_mask": hold_mask,
        "_transfer_normalized": normalized_transfer(
            frequency,
            fit_result,
            scale_hz,
        ),
        "_accepted_indices": (
            None
            if accepted_indices is None
            else [int(value) for value in accepted_indices]
        ),
    }


def clean_fit(result):
    return {
        key: value
        for key, value in result.items()
        if not key.startswith("_")
    }


def shape_difference_db(
    numerator,
    denominator,
    mask,
):
    residual_db = 20.0 * np.log10(
        np.asarray(numerator, dtype=float)
        / np.asarray(denominator, dtype=float)
    )
    values = residual_db[np.asarray(mask, dtype=bool)]
    return {
        "rms_dB": float(
            np.sqrt(np.mean(values**2))
        ),
        "mean_dB": float(np.mean(values)),
        "max_abs_dB": float(
            np.max(np.abs(values))
        ),
    }


def parameter_repeatability(
    reference_parameters,
    repeat_full_parameters,
    block_parameters,
    sigma_threshold,
):
    rows = {}
    outside_count = 0
    over_sigma = 0

    for name in CANONICAL_PARAMETER_NAMES:
        reference_value = float(reference_parameters[name])
        repeat_full_value = float(
            repeat_full_parameters[name]
        )
        values = np.asarray(
            [
                float(parameters[name])
                for parameters in block_parameters
            ],
            dtype=float,
        )
        median = float(np.median(values))
        mean = float(np.mean(values))
        std = float(
            np.std(values, ddof=1)
            if values.size > 1
            else 0.0
        )
        minimum = float(np.min(values))
        maximum = float(np.max(values))
        outside = bool(
            reference_value < minimum
            or reference_value > maximum
        )
        outside_count += int(outside)

        row = {
            "reference_value": reference_value,
            "repeat_full_value": repeat_full_value,
            "repeat_block_values": [
                float(value) for value in values
            ],
            "repeat_block_median": median,
            "repeat_block_mean": mean,
            "repeat_block_sample_std": std,
            "repeat_block_min": minimum,
            "repeat_block_max": maximum,
            "reference_outside_repeat_block_range": outside,
            "reference_over_repeat_full_ratio": (
                float(reference_value / repeat_full_value)
                if repeat_full_value != 0.0
                else None
            ),
        }

        positive = bool(
            reference_value > 0.0
            and repeat_full_value > 0.0
            and np.all(values > 0.0)
        )
        if positive:
            logs = np.log10(values)
            log_median = float(np.median(logs))
            log_std = float(
                np.std(logs, ddof=1)
                if logs.size > 1
                else 0.0
            )
            reference_log = float(
                np.log10(reference_value)
            )
            repeat_full_log = float(
                np.log10(repeat_full_value)
            )
            sigma = (
                abs(reference_log - log_median)
                / log_std
                if log_std > 0.0
                else None
            )
            exceeds_sigma = bool(
                sigma is not None
                and sigma >= float(sigma_threshold)
            )
            over_sigma += int(exceeds_sigma)
            row.update(
                {
                    "repeat_block_log10_median": (
                        log_median
                    ),
                    "repeat_block_log10_sample_std": (
                        log_std
                    ),
                    "reference_minus_repeat_full_log10": (
                        reference_log - repeat_full_log
                    ),
                    "reference_minus_block_median_log10": (
                        reference_log - log_median
                    ),
                    "reference_shift_over_repeat_block_log10_std": (
                        None
                        if sigma is None
                        else float(sigma)
                    ),
                    "reference_shift_exceeds_sigma_threshold": (
                        exceeds_sigma
                    ),
                }
            )
        else:
            row.update(
                {
                    "repeat_block_log10_median": None,
                    "repeat_block_log10_sample_std": None,
                    "reference_minus_repeat_full_log10": None,
                    "reference_minus_block_median_log10": None,
                    "reference_shift_over_repeat_block_log10_std": None,
                    "reference_shift_exceeds_sigma_threshold": False,
                }
            )
        rows[name] = row

    return {
        "parameters": rows,
        "n_parameters": len(CANONICAL_PARAMETER_NAMES),
        "n_reference_outside_repeat_block_range": int(
            outside_count
        ),
        "n_reference_shift_exceeds_sigma_threshold": int(
            over_sigma
        ),
    }


def run(config, config_path: Path):
    manifest_path = ident.resolve_config_path(
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
    reference_case = ident.case_by_label(
        cases,
        config["reference_case_label"],
    )
    repeat_case = ident.case_by_label(
        cases,
        config["repeat_case_label"],
    )

    reference_summary = json.loads(
        reference_case["summary"].read_text(
            encoding="utf-8"
        )
    )
    repeat_summary = json.loads(
        repeat_case["summary"].read_text(
            encoding="utf-8"
        )
    )
    (
        reference_comparison,
        reference_path,
        reference_source,
    ) = ident.comparison_for_case(reference_case)
    (
        repeat_comparison,
        repeat_path,
        repeat_source,
    ) = ident.comparison_for_case(repeat_case)

    reference_transfer_path = ident.resolve_config_path(
        config["reference_transfer"],
        config_path,
    )
    reference_transfer_payload = json.loads(
        reference_transfer_path.read_text(
            encoding="utf-8"
        )
    )
    reference_transfer = shared.load_shared_transfer(
        reference_transfer_payload
    )

    free_snapshot_path = ident.resolve_config_path(
        config["lowmid_free_snapshot"],
        config_path,
    )
    free_snapshot = ident.load_free_snapshot(
        free_snapshot_path
    )

    reference_indices = [
        int(value)
        for value in reference_comparison[
            "acquisition"
        ]["accepted_record_indices"]
    ]
    repeat_indices = [
        int(value)
        for value in repeat_comparison[
            "acquisition"
        ]["accepted_record_indices"]
    ]

    block_cfg = config["block_size"]
    if block_cfg["mode"] == (
        "reference_accepted_record_count"
    ):
        block_size = len(reference_indices)
    else:
        block_size = int(
            block_cfg["fallback_records"]
        )
    blocks, remainder = day_drift.make_blocks(
        repeat_indices,
        block_size,
    )
    if len(blocks) < 2:
        raise ValueError(
            "repeat run does not contain at least "
            "two full comparison blocks"
        )

    reference_white_scale = float(
        config["white_floor"][
            "reference_scale_from_detector_snapshot"
        ]
    )
    repeat_white_scale = float(
        free_snapshot["white_scale"]
    )
    seed = int(config["optimizer"]["seed"])

    reference_fit = fit_dataset(
        label=reference_case["label"],
        summary=reference_summary,
        comparison=reference_comparison,
        experiment_path=reference_path,
        white_scale=reference_white_scale,
        accepted_indices=None,
        config=config,
        reference_parameters=reference_transfer[
            "parameters"
        ],
        free_reference_parameters=free_snapshot[
            "parameters"
        ],
        seed=seed,
    )
    repeat_full_fit = fit_dataset(
        label=repeat_case["label"] + ":full",
        summary=repeat_summary,
        comparison=repeat_comparison,
        experiment_path=repeat_path,
        white_scale=repeat_white_scale,
        accepted_indices=None,
        config=config,
        reference_parameters=reference_transfer[
            "parameters"
        ],
        free_reference_parameters=free_snapshot[
            "parameters"
        ],
        seed=seed + 1,
    )

    np.testing.assert_allclose(
        reference_fit["_frequency_Hz"],
        repeat_full_fit["_frequency_Hz"],
        rtol=0.0,
        atol=1.0e-12,
    )

    block_fits = []
    for index, indices in enumerate(blocks):
        block_fits.append(
            fit_dataset(
                label=(
                    repeat_case["label"]
                    + f":block_{index + 1}"
                ),
                summary=repeat_summary,
                comparison=repeat_comparison,
                experiment_path=repeat_path,
                white_scale=repeat_white_scale,
                accepted_indices=indices,
                config=config,
                reference_parameters=reference_transfer[
                    "parameters"
                ],
                free_reference_parameters=free_snapshot[
                    "parameters"
                ],
                seed=seed + 100 + index,
            )
        )

    frequency = reference_fit["_frequency_Hz"]
    fit_mask = reference_fit["_fit_mask"]
    holdout_mask = reference_fit["_holdout_mask"]
    block_fit_shape_rows = []
    block_hold_shape_rows = []
    for block in block_fits:
        block_fit_shape_rows.append(
            shape_difference_db(
                block["_transfer_normalized"],
                repeat_full_fit[
                    "_transfer_normalized"
                ],
                fit_mask,
            )
        )
        block_hold_shape_rows.append(
            shape_difference_db(
                block["_transfer_normalized"],
                repeat_full_fit[
                    "_transfer_normalized"
                ],
                holdout_mask,
            )
        )

    day_fit_shape = shape_difference_db(
        reference_fit["_transfer_normalized"],
        repeat_full_fit["_transfer_normalized"],
        fit_mask,
    )
    day_hold_shape = shape_difference_db(
        reference_fit["_transfer_normalized"],
        repeat_full_fit["_transfer_normalized"],
        holdout_mask,
    )
    max_block_fit_rms = float(
        max(
            row["rms_dB"]
            for row in block_fit_shape_rows
        )
    )
    max_block_hold_rms = float(
        max(
            row["rms_dB"]
            for row in block_hold_shape_rows
        )
    )
    fit_shape_ratio = (
        float(
            day_fit_shape["rms_dB"]
            / max_block_fit_rms
        )
        if max_block_fit_rms > 0.0
        else None
    )
    hold_shape_ratio = (
        float(
            day_hold_shape["rms_dB"]
            / max_block_hold_rms
        )
        if max_block_hold_rms > 0.0
        else None
    )

    screen = config["repeatability_screen"]
    parameter_stats = parameter_repeatability(
        reference_fit["canonical_parameters"],
        repeat_full_fit["canonical_parameters"],
        [
            row["canonical_parameters"]
            for row in block_fits
        ],
        screen["parameter_shift_sigma_threshold"],
    )
    day_shape_exceeds = bool(
        fit_shape_ratio is not None
        and fit_shape_ratio
        > float(
            screen[
                "day_shape_rms_over_block_max_threshold"
            ]
        )
    )
    multiple_parameter_shifts = bool(
        parameter_stats[
            "n_reference_shift_exceeds_sigma_threshold"
        ]
        >= 2
    )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "model": {
            "order": 2,
            "n_total_fit_parameters": 4,
            "canonical_parameters": list(
                CANONICAL_PARAMETER_NAMES
            ),
            "reference_scale_Hz": float(
                config["reference_scale_Hz"]
            ),
            "numerator": (
                "P(x)=1+c2*x^2+c4*x^4, "
                "x=f/reference_scale_Hz"
            ),
            "latent_u_v_are_not_physical_comparison_parameters": True,
        },
        "fit_region": config["fit_region_Hz"],
        "holdout_region": {
            **config["holdout_region_Hz"],
            "optimizer_received_holdout_points": False,
        },
        "white_floor_policy": {
            "reference_scale": reference_white_scale,
            "repeat_scale": repeat_white_scale,
            "reference_white_from_detector_snapshot": True,
            "repeat_white_from_lowmid_free_snapshot": True,
            "white_floor_refit_per_block": False,
        },
        "reference": {
            "comparison_source": reference_source,
            "experiment_path": str(reference_path),
            "fit": clean_fit(reference_fit),
        },
        "repeat_full": {
            "comparison_source": repeat_source,
            "experiment_path": str(repeat_path),
            "fit": clean_fit(repeat_full_fit),
        },
        "repeat_blocks": {
            "block_size_records": int(block_size),
            "n_full_nonoverlapping_blocks": int(
                len(blocks)
            ),
            "unused_remainder_records": int(remainder),
            "fits": [
                {
                    **clean_fit(row),
                    "accepted_index_first": int(
                        row["_accepted_indices"][0]
                    ),
                    "accepted_index_last": int(
                        row["_accepted_indices"][-1]
                    ),
                }
                for row in block_fits
            ],
        },
        "parameter_repeatability": parameter_stats,
        "transfer_shape_repeatability": {
            "reference_vs_repeat_full_fit_region": (
                day_fit_shape
            ),
            "reference_vs_repeat_full_holdout_region": (
                day_hold_shape
            ),
            "repeat_block_vs_repeat_full_fit_region": (
                block_fit_shape_rows
            ),
            "repeat_block_vs_repeat_full_holdout_region": (
                block_hold_shape_rows
            ),
            "max_repeat_block_fit_region_rms_dB": (
                max_block_fit_rms
            ),
            "max_repeat_block_holdout_region_rms_dB": (
                max_block_hold_rms
            ),
            "day_fit_region_rms_over_block_max_ratio": (
                fit_shape_ratio
            ),
            "day_holdout_region_rms_over_block_max_ratio": (
                hold_shape_ratio
            ),
        },
        "interpretation_flags": {
            "day_order2_transfer_shape_exceeds_within_repeat_block_variation": (
                day_shape_exceeds
            ),
            "at_least_two_canonical_parameters_shift_beyond_block_sigma_screen": (
                multiple_parameter_shifts
            ),
            "day_specific_order2_shape_shift_screen_passes": bool(
                day_shape_exceeds
                and multiple_parameter_shifts
            ),
            "physical_electronics_drift_identified": False,
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "manifest": str(manifest_path),
            "reference_transfer": str(
                reference_transfer_path
            ),
            "lowmid_free_snapshot": str(
                free_snapshot_path
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
        json.dumps(result, indent=2, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(output),
                "reference_parameters": result[
                    "reference"
                ]["fit"]["canonical_parameters"],
                "repeat_full_parameters": result[
                    "repeat_full"
                ]["fit"]["canonical_parameters"],
                "n_repeat_blocks": result[
                    "repeat_blocks"
                ]["n_full_nonoverlapping_blocks"],
                "day_fit_shape_rms_over_block_max": result[
                    "transfer_shape_repeatability"
                ][
                    "day_fit_region_rms_over_block_max_ratio"
                ],
                "n_parameters_outside_block_range": result[
                    "parameter_repeatability"
                ][
                    "n_reference_outside_repeat_block_range"
                ],
                "n_parameters_over_sigma_screen": result[
                    "parameter_repeatability"
                ][
                    "n_reference_shift_exceeds_sigma_threshold"
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

"""Cross-validate one shared pre-ADC readout transfer across multiple datasets.

A Git-tracked reference transfer and Git-tracked case manifest are used by
default.  For every manifest case, the case-specific detector snapshot is held
fixed.  The experimental target is reconstructed from that case's exact
accepted CH0 raw-record mask with no 10 kHz digital analysis filter.

Command-line overrides remain available, but the standard run requires no
external JSON files outside the repository.

The shared transfer is never refit per case.  Optionally, the same topology is
also refit locally for each case to quantify how close the shared transfer is to
the case-specific best shape.
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
DEFAULT_MANIFEST = CONFIG_DIR / "shared_readout_cross_dataset_manifest.json"
DEFAULT_REFERENCE_PROFILE = CONFIG_DIR / "shared_readout_reference_transfer.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "shared_readout_cross_dataset_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import preanalysis_hybrid_pole_profile_diagnostic as profile  # noqa: E402
from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    accepted_noise_indices,
)

REFERENCE_HZ = 1_000.0
RMS_SCREEN_DB_DEFAULT = 1.0
MAX_SCREEN_DB_DEFAULT = 3.0
SHARED_LOCAL_SCORE_TOLERANCE_DEFAULT = 1.25
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
TRANSFER_PARAMETER_NAMES = (
    "pole_Hz",
    "pole_Q",
    "zero_Hz",
    "zero_Q",
    "leadlag_zero_Hz",
)


def resolve_manifest_path(value, manifest_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (manifest_path.parent / path).resolve()


def load_shared_transfer(reference_profile: dict):
    row = reference_profile.get("best_profile_row")
    if not isinstance(row, dict):
        raise ValueError("reference profile has no best_profile_row")
    parameters = row.get("parameters")
    if not isinstance(parameters, dict):
        raise ValueError("best_profile_row has no parameters")

    shared = {}
    for name in TRANSFER_PARAMETER_NAMES:
        if name not in parameters:
            raise ValueError(f"reference transfer missing {name}")
        value = float(parameters[name])
        if value <= 0.0:
            raise ValueError(f"reference transfer {name} must be positive")
        shared[name] = value

    pole_is_infinite = bool(
        parameters.get(
            "leadlag_pole_is_infinite",
            row.get("leadlag_pole_is_infinite", False),
        )
    )
    pole_value = parameters.get(
        "leadlag_pole_Hz",
        row.get("fixed_leadlag_pole_Hz"),
    )
    if pole_is_infinite:
        fixed_pole_hz = None
    else:
        if pole_value is None:
            raise ValueError(
                "finite reference lead/lag pole requires leadlag_pole_Hz"
            )
        fixed_pole_hz = float(pole_value)
        if fixed_pole_hz <= 0.0:
            raise ValueError("reference lead/lag pole must be positive")

    return {
        "parameters": shared,
        "fixed_leadlag_pole_Hz": fixed_pole_hz,
        "leadlag_pole_is_infinite": bool(fixed_pole_hz is None),
        "reference_shape_score": (
            float(row["shape_score"])
            if row.get("shape_score") is not None
            else None
        ),
    }


def normalize_manifest(manifest: dict, manifest_path: Path):
    raw_cases = manifest.get("cases")
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ValueError("manifest must contain a non-empty cases list")

    seen = set()
    cases = []
    for index, raw in enumerate(raw_cases):
        if not isinstance(raw, dict):
            raise ValueError(f"manifest case {index} must be an object")
        label = str(raw.get("label", f"case_{index}")).strip()
        if not label:
            raise ValueError(f"manifest case {index} has an empty label")
        if label in seen:
            raise ValueError(f"duplicate case label: {label}")
        seen.add(label)

        if "summary" not in raw:
            raise ValueError(f"case {label} requires summary")
        comparison_summary = raw.get("comparison_summary")
        comparison_spec = raw.get("comparison_spec")
        if (comparison_summary is None) == (comparison_spec is None):
            raise ValueError(
                f"case {label} requires exactly one of "
                "comparison_summary or comparison_spec"
            )
        role = str(raw.get("role", "validation")).strip().lower()
        allowed_roles = {"reference", "validation", "repeat_validation"}
        if role not in allowed_roles:
            raise ValueError(
                f"case {label} role must be reference, validation, "
                "or repeat_validation"
            )

        row = {
            "label": label,
            "role": role,
            "summary": resolve_manifest_path(
                raw["summary"], manifest_path
            ),
            "comparison_summary": (
                resolve_manifest_path(comparison_summary, manifest_path)
                if comparison_summary is not None
                else None
            ),
            "comparison_spec": (
                resolve_manifest_path(comparison_spec, manifest_path)
                if comparison_spec is not None
                else None
            ),
            "allow_missing_raw_data": bool(
                raw.get("allow_missing_raw_data", False)
            ),
        }
        if raw.get("experiment_path") is not None:
            row["experiment_path"] = Path(raw["experiment_path"])
        else:
            row["experiment_path"] = None
        cases.append(row)
    return cases


def build_comparison_from_spec(spec_path: Path):
    """Build the minimal comparison-summary contract from a tracked spec.

    This reproduces the production noise-record acceptance mask from raw
    records, so repeat-validation datasets do not require an untracked
    generated comparison_summary.json.
    """

    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    experiment_path = Path(spec["experiment_path"])
    acquisition = spec["acquisition"]
    rate = float(acquisition["rate_Hz"])
    samples = int(acquisition["samples"])
    cutoff = float(acquisition["cutoff_Hz"])
    channel = str(acquisition.get("channel", "CH0"))
    acceptance = spec.get("acceptance", {})
    max_peak_to_peak = float(
        acceptance.get("max_peak_to_peak_raw_units", 0.04)
    )
    remove_mean = bool(acceptance.get("remove_mean", True))
    apply_raw_range = bool(
        acceptance.get("apply_to_raw", True)
    )
    apply_processed_range = bool(
        acceptance.get("apply_to_processed", True)
    )

    raw_dir = experiment_path / f"{channel}_noise" / "rawdata"
    paths = sorted(raw_dir.glob(f"{channel}_*.dat"))
    if not paths:
        raise FileNotFoundError(
            f"no {channel} raw records found in {raw_dir}"
        )

    def records():
        for path in paths:
            yield base.read_record(path)

    def range_ok(values):
        values = np.asarray(values, dtype=float)
        return float(np.max(values) - np.min(values)) <= max_peak_to_peak

    accepted_indices = accepted_noise_indices(
        records(),
        samples,
        rate,
        cutoff=cutoff,
        remove_mean=remove_mean,
        accept_raw=range_ok if apply_raw_range else None,
        accept_processed=(
            range_ok if apply_processed_range else None
        ),
    )
    if not accepted_indices:
        raise ValueError(
            f"comparison spec {spec_path} accepted no records"
        )

    return {
        "status": "generated_in_memory_from_git_tracked_spec",
        "experiment_path": str(experiment_path),
        "comparison_spec": str(spec_path),
        "acquisition": {
            "rate_Hz": rate,
            "samples": samples,
            "cutoff_Hz": cutoff,
            "channel": channel,
            "accepted_records": int(len(accepted_indices)),
            "accepted_record_indices": [
                int(value) for value in accepted_indices
            ],
        },
        "acceptance": {
            "max_peak_to_peak_raw_units": max_peak_to_peak,
            "remove_mean": remove_mean,
            "apply_to_raw": apply_raw_range,
            "apply_to_processed": apply_processed_range,
            "selection_function": "accepted_noise_indices",
        },
    }


def normalized_model_from_components(frequency, components):
    absolute = np.sqrt(
        np.asarray(components["main_asd"], dtype=float) ** 2
        + np.asarray(components["alias_asd"], dtype=float) ** 2
        + np.asarray(components["white_asd"], dtype=float) ** 2
    )
    return opt.normalize_at(
        frequency,
        absolute,
        reference_hz=REFERENCE_HZ,
    )


def band_change_summary(baseline_bands, shared_bands):
    result = {}
    for name in baseline_bands:
        before = float(
            baseline_bands[name]["mean_model_over_target_dB"]
        )
        after = float(
            shared_bands[name]["mean_model_over_target_dB"]
        )
        result[name] = {
            "baseline_mean_model_over_target_dB": before,
            "shared_mean_model_over_target_dB": after,
            "absolute_error_change_dB": float(
                abs(after) - abs(before)
            ),
            "improves_absolute_mean_error": bool(
                abs(after) < abs(before)
            ),
        }
    return result


def parameter_drift(shared, local):
    rows = {}
    log_values = []
    for name in TRANSFER_PARAMETER_NAMES:
        reference = float(shared[name])
        fitted = float(local[name])
        ratio = fitted / reference
        log10_ratio = float(np.log10(ratio))
        log_values.append(log10_ratio)
        rows[name] = {
            "shared_value": reference,
            "local_value": fitted,
            "local_over_shared_ratio": float(ratio),
            "log10_ratio": log10_ratio,
        }
    return {
        "parameters": rows,
        "rms_log10_parameter_ratio": float(
            np.sqrt(np.mean(np.asarray(log_values) ** 2))
        ),
        "max_abs_log10_parameter_ratio": float(
            np.max(np.abs(log_values))
        ),
    }


def anchor_sample(
    frequency,
    target,
    baseline_model,
    shared_model,
    local_model=None,
):
    rows = []
    for anchor in ANCHOR_FREQUENCIES_HZ:
        index = int(np.argmin(np.abs(frequency - anchor)))
        row = {
            "frequency_Hz": float(frequency[index]),
            "target_pre_analysis_normalized": float(target[index]),
            "baseline_model_normalized": float(
                baseline_model[index]
            ),
            "shared_model_normalized": float(shared_model[index]),
            "shared_over_baseline_correction_dB": float(
                20.0
                * np.log10(
                    shared_model[index] / baseline_model[index]
                )
            ),
        }
        if local_model is not None:
            row["local_best_model_normalized"] = float(
                local_model[index]
            )
        rows.append(row)
    return rows


def evaluation_flags(
    baseline_score,
    shared_score,
    baseline_bands,
    shared_bands,
    local_score,
    shared_local_tolerance,
):
    mid = "5000-15000_Hz"
    high = "40000-100000_Hz"
    tail = "100000-200000_Hz"

    def mean_abs(bands, name):
        return abs(float(bands[name]["mean_model_over_target_dB"]))

    score_improves = shared_score < baseline_score
    mid_improves = mean_abs(shared_bands, mid) < mean_abs(
        baseline_bands, mid
    )
    high_improves = mean_abs(shared_bands, high) < mean_abs(
        baseline_bands, high
    )
    tail_not_worse = mean_abs(shared_bands, tail) <= mean_abs(
        baseline_bands, tail
    ) + 1e-12

    shared_local_ratio = None
    close_to_local = None
    if local_score is not None:
        shared_local_ratio = float(shared_score / local_score)
        close_to_local = bool(
            shared_local_ratio <= float(shared_local_tolerance)
        )

    return {
        "shared_improves_shape_score": bool(score_improves),
        "shared_improves_5_15k": bool(mid_improves),
        "shared_improves_40_100k": bool(high_improves),
        "shared_does_not_worsen_100_200k_mean": bool(
            tail_not_worse
        ),
        "shared_strict_broad_improvement": bool(
            score_improves
            and mid_improves
            and high_improves
            and tail_not_worse
        ),
        "shared_over_local_best_score_ratio": shared_local_ratio,
        "shared_within_local_score_tolerance": close_to_local,
    }


def evaluate_case(
    case,
    shared_transfer,
    *,
    local_refit,
    shared_local_tolerance,
    seed,
    de_maxiter,
    rms_screen_db,
    max_screen_db,
):
    summary = json.loads(
        case["summary"].read_text(encoding="utf-8")
    )
    if case["comparison_summary"] is not None:
        comparison = json.loads(
            case["comparison_summary"].read_text(encoding="utf-8")
        )
        comparison_source = "tracked_comparison_summary"
    else:
        comparison = build_comparison_from_spec(
            case["comparison_spec"]
        )
        comparison_source = "tracked_comparison_spec_recomputed_mask"
    fit_args = base.fit_args(summary)
    frequency = np.geomspace(
        fit_args.fit_min_hz,
        fit_args.fit_max_hz,
        fit_args.fit_points,
    )

    experiment_path = (
        case["experiment_path"]
        if case["experiment_path"] is not None
        else Path(comparison["experiment_path"])
    )
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        frequency,
    )
    candidate = dict(summary["best_case_parameters"])
    context = base.intrinsic_context(candidate, frequency)

    acquisition_rate = float(target_context["rate_Hz"])
    model_rate = float(context["rate_Hz"])
    if not np.isclose(
        acquisition_rate,
        model_rate,
        rtol=0.0,
        atol=max(acquisition_rate, 1.0) * 1e-12,
    ):
        raise ValueError(
            f"{case['label']}: acquisition rate {acquisition_rate} "
            f"does not match detector summary rate {model_rate}"
        )

    baseline_components = profile.baseline_components(context)
    baseline_model = normalized_model_from_components(
        frequency,
        baseline_components,
    )
    shared_model = profile.pre_analysis_model(
        context,
        shared_transfer["parameters"],
        shared_transfer["fixed_leadlag_pole_Hz"],
    )

    baseline_score = float(
        opt.fit_score(
            baseline_model,
            target_context["target"],
            frequency,
            fit_args,
        )
    )
    shared_score = float(
        opt.fit_score(
            shared_model,
            target_context["target"],
            frequency,
            fit_args,
        )
    )
    baseline_metrics = base.residual_db_metrics(
        baseline_model,
        target_context["target"],
        frequency,
        fit_args,
    )
    shared_metrics = base.residual_db_metrics(
        shared_model,
        target_context["target"],
        frequency,
        fit_args,
    )
    baseline_bands = base.band_summary(
        baseline_model,
        target_context["target"],
        frequency,
        fit_args,
    )
    shared_bands = base.band_summary(
        shared_model,
        target_context["target"],
        frequency,
        fit_args,
    )

    shared_components = profile.shaped_components(
        context,
        shared_transfer["parameters"],
        shared_transfer["fixed_leadlag_pole_Hz"],
    )
    shared_screen = bool(
        shared_metrics["rms_residual_dB"] <= float(rms_screen_db)
        and shared_metrics["max_abs_residual_dB"]
        <= float(max_screen_db)
    )

    local_clean = None
    local_model = None
    local_score = None
    drift = None
    if local_refit:
        local = profile.fit_fixed_pole(
            shared_transfer["fixed_leadlag_pole_Hz"],
            context,
            target_context["target"],
            frequency,
            fit_args,
            center_min_hz=1_000.0,
            center_max_hz=300_000.0,
            general_q_min=0.10,
            q_max=20.0,
            seed=seed,
            de_maxiter=de_maxiter,
            rms_screen_db=rms_screen_db,
            max_screen_db=max_screen_db,
            initial_parameters=shared_transfer["parameters"],
        )
        local_model = local["_model"]
        local_score = float(local["shape_score"])
        local_score_tolerance = max(
            1.0e-12,
            abs(shared_score) * 1.0e-9,
        )
        if local_score > shared_score + local_score_tolerance:
            raise RuntimeError(
                f"{case['label']}: warm-started local fit "
                f"({local_score}) is worse than the included shared "
                f"candidate ({shared_score})"
            )
        local_clean = {
            key: value
            for key, value in local.items()
            if not key.startswith("_")
        }
        drift = parameter_drift(
            shared_transfer["parameters"],
            local["_parameters_raw"],
        )
        local_clean["shared_candidate_score"] = shared_score
        local_clean["local_not_worse_than_shared"] = bool(
            local_score <= shared_score + local_score_tolerance
        )
        local_clean["local_minus_shared_score"] = float(
            local_score - shared_score
        )

    flags = evaluation_flags(
        baseline_score,
        shared_score,
        baseline_bands,
        shared_bands,
        local_score,
        shared_local_tolerance,
    )

    return {
        "status": "evaluated",
        "label": case["label"],
        "role": case["role"],
        "summary": str(case["summary"]),
        "comparison_summary": (
            str(case["comparison_summary"])
            if case["comparison_summary"] is not None
            else None
        ),
        "comparison_spec": (
            str(case["comparison_spec"])
            if case["comparison_spec"] is not None
            else None
        ),
        "comparison_source": comparison_source,
        "experiment_path": str(experiment_path),
        "accepted_records": int(target_context["accepted_records"]),
        "acquisition": {
            "rate_Hz": acquisition_rate,
            "samples": int(target_context["sample"]),
        },
        "baseline": {
            "shape_score": baseline_score,
            "residual_metrics": baseline_metrics,
            "bands": baseline_bands,
        },
        "shared_transfer": {
            "shape_score": shared_score,
            "score_ratio_to_baseline": float(
                shared_score / baseline_score
            ),
            "residual_metrics": shared_metrics,
            "bands": shared_bands,
            "passes_screen_tolerance": shared_screen,
            "alias_fraction_sample": (
                profile.component_fraction_sample(
                    frequency,
                    shared_components,
                )
            ),
        },
        "band_change": band_change_summary(
            baseline_bands,
            shared_bands,
        ),
        "local_refit": local_clean,
        "local_parameter_drift_from_shared": drift,
        "flags": flags,
        "curve_sample": anchor_sample(
            frequency,
            target_context["target"],
            baseline_model,
            shared_model,
            local_model,
        ),
    }


def geometric_mean(values):
    values = np.asarray(values, dtype=float)
    if np.any(values <= 0.0):
        raise ValueError("geometric mean requires positive values")
    return float(np.exp(np.mean(np.log(values))))


def aggregate_results(rows, shared_local_tolerance):
    validation = [row for row in rows if row["role"] == "validation"]
    reference = [row for row in rows if row["role"] == "reference"]
    score_ratios = [
        row["shared_transfer"]["score_ratio_to_baseline"]
        for row in rows
    ]
    validation_ratios = [
        row["shared_transfer"]["score_ratio_to_baseline"]
        for row in validation
    ]
    local_ratios = [
        row["flags"]["shared_over_local_best_score_ratio"]
        for row in rows
        if row["flags"]["shared_over_local_best_score_ratio"]
        is not None
    ]

    all_score_improve = all(
        row["flags"]["shared_improves_shape_score"] for row in rows
    )
    all_validation_score_improve = bool(validation) and all(
        row["flags"]["shared_improves_shape_score"]
        for row in validation
    )
    all_validation_mid_high = bool(validation) and all(
        row["flags"]["shared_improves_5_15k"]
        and row["flags"]["shared_improves_40_100k"]
        for row in validation
    )
    all_validation_close_local = bool(validation) and all(
        row["flags"]["shared_within_local_score_tolerance"] is True
        for row in validation
    )

    return {
        "n_cases": int(len(rows)),
        "n_reference_cases": int(len(reference)),
        "n_validation_cases": int(len(validation)),
        "has_validation_cases": bool(validation),
        "all_cases_shared_improve_score": bool(all_score_improve),
        "all_validation_cases_shared_improve_score": bool(
            all_validation_score_improve
        ),
        "all_validation_cases_improve_mid_and_high": bool(
            all_validation_mid_high
        ),
        "all_validation_cases_within_local_score_tolerance": bool(
            all_validation_close_local
        ),
        "shared_local_score_tolerance": float(
            shared_local_tolerance
        ),
        "geometric_mean_shared_score_ratio_to_baseline": (
            geometric_mean(score_ratios)
        ),
        "validation_geometric_mean_shared_score_ratio_to_baseline": (
            geometric_mean(validation_ratios)
            if validation_ratios
            else None
        ),
        "median_shared_over_local_best_score_ratio": (
            float(np.median(local_ratios))
            if local_ratios
            else None
        ),
        "max_shared_over_local_best_score_ratio": (
            float(np.max(local_ratios))
            if local_ratios
            else None
        ),
        "interpretation_flags": {
            "supports_shared_readout_transfer_across_validation_cases": bool(
                validation
                and all_validation_score_improve
                and all_validation_mid_high
                and (
                    not local_ratios
                    or all_validation_close_local
                )
            ),
            "validation_dataset_required_for_cross_validation_claim": bool(
                not validation
            ),
        },
    }


def run(
    manifest,
    manifest_path,
    reference_profile,
    *,
    local_refit=True,
    shared_local_tolerance=SHARED_LOCAL_SCORE_TOLERANCE_DEFAULT,
    seed=20260915,
    de_maxiter=100,
    rms_screen_db=RMS_SCREEN_DB_DEFAULT,
    max_screen_db=MAX_SCREEN_DB_DEFAULT,
):
    shared = load_shared_transfer(reference_profile)
    cases = normalize_manifest(manifest, manifest_path)
    rows = []
    for index, case in enumerate(cases):
        rows.append(
            evaluate_case(
                case,
                shared,
                local_refit=local_refit,
                shared_local_tolerance=shared_local_tolerance,
                seed=seed + index,
                de_maxiter=de_maxiter,
                rms_screen_db=rms_screen_db,
                max_screen_db=max_screen_db,
            )
        )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "shared_readout_transfer_held_fixed_across_cases": True,
        "detector_parameters_are_case_specific_and_frozen": True,
        "digital_analysis_filter_applied": False,
        "target_semantics": (
            "for every case: fresh pre-analysis ASD from that case's exact "
            "accepted raw CH0 mask; mean removal, Hann, power average, "
            "one-sided ASD; no 10 kHz digital analysis filter"
        ),
        "model_semantics": (
            "case-specific frozen intrinsic detector ASD -> production "
            "100 kHz analog Bessel -> one shared pre-ADC hybrid transfer -> "
            "first alias fold -> case-specific production post-filter white"
        ),
        "shared_transfer": {
            **shared,
            "parameters": {
                **shared["parameters"],
                "leadlag_pole_Hz": shared[
                    "fixed_leadlag_pole_Hz"
                ],
                "leadlag_pole_is_infinite": shared[
                    "leadlag_pole_is_infinite"
                ],
            },
        },
        "local_refit_enabled": bool(local_refit),
        "local_refit_purpose": (
            "diagnostic comparator only: quantify how close the fixed "
            "shared transfer is to each case-specific best transfer; local "
            "refits are not used to score shared-transfer generalization"
        ),
        "screen_tolerance": {
            "rms_residual_dB_max": float(rms_screen_db),
            "max_abs_residual_dB_max": float(max_screen_db),
            "diagnostic_not_physical_prior": True,
        },
        "shared_local_score_tolerance": float(
            shared_local_tolerance
        ),
        "cases": rows,
        "aggregate": aggregate_results(
            rows,
            shared_local_tolerance,
        ),
        "guardrail": (
            "A shared-transfer improvement on the reference case alone is "
            "not cross-validation. A readout-origin claim requires one or "
            "more independent validation cases acquired through the same "
            "readout chain. Large local-transfer drift across cases weakens "
            "the interpretation even if every case can be fit individually."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=(
            "Git-tracked case manifest; defaults to "
            "PoST_Simulations/config/shared_readout_cross_dataset_manifest.json"
        ),
    )
    parser.add_argument(
        "--reference-profile",
        type=Path,
        default=DEFAULT_REFERENCE_PROFILE,
        help=(
            "Git-tracked shared transfer; defaults to "
            "PoST_Simulations/config/shared_readout_reference_transfer.json"
        ),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-local-refit", action="store_true")
    parser.add_argument(
        "--shared-local-score-tolerance",
        type=float,
        default=SHARED_LOCAL_SCORE_TOLERANCE_DEFAULT,
    )
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

    if args.shared_local_score_tolerance < 1.0:
        raise ValueError(
            "shared-local score tolerance must be >= 1"
        )

    manifest = json.loads(
        args.manifest.read_text(encoding="utf-8")
    )
    reference_profile = json.loads(
        args.reference_profile.read_text(encoding="utf-8")
    )
    result = run(
        manifest,
        args.manifest,
        reference_profile,
        local_refit=not args.skip_local_refit,
        shared_local_tolerance=args.shared_local_score_tolerance,
        seed=args.seed,
        de_maxiter=args.de_maxiter,
        rms_screen_db=args.rms_screen_db,
        max_screen_db=args.max_screen_db,
    )
    result["manifest"] = str(args.manifest)
    result["reference_profile"] = str(args.reference_profile)

    output = args.output or DEFAULT_OUTPUT
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    aggregate = result["aggregate"]
    print(
        json.dumps(
            {
                "output": str(output),
                "n_cases": aggregate["n_cases"],
                "n_validation_cases": aggregate[
                    "n_validation_cases"
                ],
                "shared_transfer": result["shared_transfer"][
                    "parameters"
                ],
                "validation_geometric_mean_score_ratio": aggregate[
                    "validation_geometric_mean_shared_score_ratio_to_baseline"
                ],
                "median_shared_over_local_best_score_ratio": aggregate[
                    "median_shared_over_local_best_score_ratio"
                ],
                "supports_shared_readout_transfer": aggregate[
                    "interpretation_flags"
                ][
                    "supports_shared_readout_transfer_across_validation_cases"
                ],
                "needs_validation_dataset": aggregate[
                    "interpretation_flags"
                ][
                    "validation_dataset_required_for_cross_validation_claim"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

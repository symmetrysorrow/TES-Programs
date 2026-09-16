"""Fit the smooth repeat-case continuum after removing only narrow spectral lines.

The line mask is determined once from the experimental native 5 Hz spectrum
before any TES/readout model is evaluated.  Only narrow peaks above 40 kHz are
eligible for repair.  Their intervals are replaced by a local median floor in
log-ASD; broad shoulders/humps are deliberately preserved.

The repaired continuum is then fit over 1--200 kHz using the existing
stability-aware detector/DC baseline.  The current minimal readout family
(pole + Q + c2, with c4 fixed to the tracked reference value) is compared with
full order-2 after profiling the post-filter white floor.  Added c4 is accepted
only when it gives a material continuum improvement and does not run to a
bound.

This is a diagnostic/effective measurement-shape fit.  It does not identify a
physical electronics component and does not modify production noise data.
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
DEFAULT_CONFIG = CONFIG_DIR / "line_robust_continuum_fit_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "line_robust_continuum_fit_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "line_robust_continuum_fit_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import high_frequency_line_noise_diagnostic as line_diag  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import readout_detector_state_competition_diagnostic as competition  # noqa: E402
from subScript import readout_lowmid_holdout_diagnostic as holdout  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as residual  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def repair_narrow_lines(
    native_frequency,
    native_asd,
    config,
):
    frequency = np.asarray(native_frequency, dtype=float)
    asd = np.asarray(native_asd, dtype=float)
    normalized = opt.normalize_at(
        frequency,
        asd,
        reference_hz=base.REFERENCE_HZ,
    )
    detect_min = float(config["detect_min_Hz"])
    detect_max = float(config["detect_max_Hz"])
    detect_mask = (
        (frequency >= detect_min)
        & (frequency <= detect_max)
    )
    detected, _baseline, _excess = line_diag.detect_lines(
        frequency[detect_mask],
        normalized[detect_mask],
        baseline_width_hz=float(config["baseline_width_Hz"]),
        min_excess_db=float(config["minimum_excess_dB"]),
        min_prominence_db=float(config["minimum_prominence_dB"]),
        min_spacing_hz=float(config["minimum_spacing_Hz"]),
        max_peaks=int(config["max_detected_peaks"]),
    )

    bin_hz = float(np.median(np.diff(frequency)))
    kernel = line_diag.odd_kernel_bins(
        float(config["baseline_width_Hz"]),
        bin_hz,
    )
    baseline_db = line_diag.local_baseline_db(
        normalized,
        kernel,
    )
    repaired = np.asarray(normalized, dtype=float).copy()
    repair_mask = np.zeros_like(frequency, dtype=bool)
    accepted_lines = []
    rejected_broad = []

    max_width = float(config["max_line_width_Hz"])
    minimum_half_width = float(
        config["minimum_repair_half_width_Hz"]
    )
    expansion = float(config["width_expansion_factor"])

    for row in detected:
        width = float(row["half_prominence_width_Hz"])
        if width > max_width:
            rejected_broad.append(
                {
                    **row,
                    "reason": "width_exceeds_narrow_line_limit",
                }
            )
            continue
        center = float(row["frequency_Hz"])
        half_width = max(
            minimum_half_width,
            0.5 * width * expansion,
        )
        local = np.abs(frequency - center) <= half_width
        local &= detect_mask
        repair_mask |= local
        accepted_lines.append(
            {
                **row,
                "repair_half_width_Hz": float(half_width),
                "repair_min_Hz": float(center - half_width),
                "repair_max_Hz": float(center + half_width),
                "repaired_native_bins": int(np.count_nonzero(local)),
            }
        )

    repaired[repair_mask] = 10.0 ** (
        baseline_db[repair_mask] / 20.0
    )
    return {
        "raw_normalized": normalized,
        "continuum_normalized": repaired,
        "repair_mask": repair_mask,
        "accepted_lines": accepted_lines,
        "rejected_broad_features": rejected_broad,
        "native_bin_Hz": bin_hz,
        "repaired_fraction_full_grid": float(
            np.mean(repair_mask)
        ),
        "repaired_fraction_detection_band": float(
            np.mean(repair_mask[detect_mask])
        ),
    }


def profile_white_scale(
    detector,
    readout,
    frequency,
    target,
    fit_args,
    scale_hz,
    baseline_white_asd,
    profile_cfg,
):
    min_scale = float(profile_cfg["min_scale"])
    max_scale = float(profile_cfg["max_scale"])
    if not 0.0 < min_scale < max_scale:
        raise ValueError("invalid white profile scale bounds")

    def objective(log10_scale):
        scale = float(10.0 ** float(log10_scale))
        model, _point = residual.model_for_candidate(
            detector,
            readout,
            frequency,
            scale_hz,
            float(baseline_white_asd) * scale,
        )
        if model is None:
            return 1.0e6
        return float(
            opt.fit_score(
                model,
                target,
                frequency,
                fit_args,
            )
        )

    optimized = minimize_scalar(
        objective,
        bounds=(np.log10(min_scale), np.log10(max_scale)),
        method="bounded",
        options={"xatol": 1.0e-6, "maxiter": 200},
    )
    candidates = [
        {
            "scale": 1.0,
            "shape_score": objective(0.0),
            "source": "baseline_scale_one",
        },
        {
            "scale": min_scale,
            "shape_score": objective(np.log10(min_scale)),
            "source": "lower_bound",
        },
        {
            "scale": max_scale,
            "shape_score": objective(np.log10(max_scale)),
            "source": "upper_bound",
        },
        {
            "scale": float(10.0 ** optimized.x),
            "shape_score": float(optimized.fun),
            "source": "bounded_scalar_optimization",
        },
    ]
    best = min(candidates, key=lambda row: row["shape_score"])
    scale = float(best["scale"])
    span = np.log10(max_scale) - np.log10(min_scale)
    tol = max(span, 1.0e-12) * 1.0e-5
    at_lower = abs(np.log10(scale) - np.log10(min_scale)) <= tol
    at_upper = abs(np.log10(scale) - np.log10(max_scale)) <= tol
    return {
        "baseline_white_asd_A_rtHz": float(baseline_white_asd),
        "best_scale": scale,
        "best_white_asd_A_rtHz": float(baseline_white_asd) * scale,
        "shape_score_with_detector_readout_frozen": float(
            best["shape_score"]
        ),
        "best_source": best["source"],
        "optimizer_success": bool(optimized.success),
        "optimizer_message": str(optimized.message),
        "at_lower_bound": bool(at_lower),
        "at_upper_bound": bool(at_upper),
        "candidates": candidates,
    }


def metrics_in_region(
    model,
    target,
    frequency,
    full_args,
    min_hz,
    max_hz,
):
    mask = (
        (frequency >= float(min_hz))
        & (frequency <= float(max_hz))
    )
    args = competition.band_args(
        full_args,
        float(min_hz),
        float(max_hz),
        int(np.count_nonzero(mask)),
    )
    return holdout.model_metrics(
        np.asarray(model)[mask],
        np.asarray(target)[mask],
        np.asarray(frequency)[mask],
        args,
    )


def fit_with_profiled_white(
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
    fit_args,
    scale_hz,
    baseline_white_asd,
    optimizer_cfg,
    profile_cfg,
    seed,
    warm_solutions=(),
):
    current_white = float(baseline_white_asd)
    history = []
    warm = list(warm_solutions)
    latest = None

    outer_iterations = int(profile_cfg.get("outer_iterations", 1))
    for outer in range(max(1, outer_iterations)):
        fit_family = dict(family)
        fit_family["name"] = (
            f"{family['name']}:white_iteration_{outer}"
        )
        row = residual.fit_family(
            family=fit_family,
            baseline_detector=baseline_detector,
            reference_readout=reference_readout,
            local_readout=local_readout,
            detector_names=detector_names,
            detector_bounds=detector_bounds,
            readout_bounds=readout_bounds,
            frequency=frequency,
            target=target,
            args=fit_args,
            scale_hz=scale_hz,
            white_asd=current_white,
            optimizer_cfg=optimizer_cfg,
            seed=int(seed) + outer * 17,
            warm_solutions=warm,
        )
        white_profile = profile_white_scale(
            row["_detector_full"],
            row["_readout_full"],
            frequency,
            target,
            fit_args,
            scale_hz,
            baseline_white_asd,
            profile_cfg,
        )
        history.append(
            {
                "iteration": outer,
                "fit_with_white_asd_A_rtHz": current_white,
                "shape_score_before_white_profile": float(
                    row["shape_score"]
                ),
                "rms_before_white_profile_dB": float(
                    row["residual_metrics"]["rms_residual_dB"]
                ),
                "white_profile": white_profile,
            }
        )
        current_white = float(
            white_profile["best_white_asd_A_rtHz"]
        )
        warm = [
            (
                f"{family['name']}:outer_{outer}",
                row["_detector_full"],
                row["_readout_full"],
            )
        ]
        latest = row

    final_family = dict(family)
    final_family["name"] = family["name"]
    final = residual.fit_family(
        family=final_family,
        baseline_detector=baseline_detector,
        reference_readout=reference_readout,
        local_readout=local_readout,
        detector_names=detector_names,
        detector_bounds=detector_bounds,
        readout_bounds=readout_bounds,
        frequency=frequency,
        target=target,
        args=fit_args,
        scale_hz=scale_hz,
        white_asd=current_white,
        optimizer_cfg=optimizer_cfg,
        seed=int(seed) + 997,
        warm_solutions=warm,
    )
    final["profiled_white_asd_A_rtHz"] = current_white
    final["profiled_white_scale_to_tracked"] = float(
        current_white / float(baseline_white_asd)
    )
    final["white_profile_history"] = history
    final["_initial_outer_solution"] = latest
    return final


def added_c4_is_at_bound(row):
    hit = row.get("readout_boundary_hits", {}).get("c4")
    if hit is None:
        return False
    return bool(hit["at_lower"] or hit["at_upper"])


def choose_recommended(pole_row, full_row, screen):
    ratio = float(
        full_row["shape_score"] / pole_row["shape_score"]
    )
    rms_improvement = float(
        pole_row["residual_metrics"]["rms_residual_dB"]
        - full_row["residual_metrics"]["rms_residual_dB"]
    )
    c4_bound = added_c4_is_at_bound(full_row)
    material = bool(
        ratio <= float(
            screen["full_order2_material_score_ratio"]
        )
        and rms_improvement
        >= float(
            screen["full_order2_min_rms_improvement_dB"]
        )
    )
    if bool(screen.get("reject_if_added_c4_at_bound", True)):
        material = bool(material and not c4_bound)
    selected = full_row if material else pole_row
    return {
        "selected_family": selected["name"],
        "full_order2_score_ratio_to_pole_section_plus_c2": ratio,
        "full_order2_rms_improvement_dB": rms_improvement,
        "full_order2_c4_at_bound": c4_bound,
        "full_order2_material_improvement": material,
        "rule": screen,
    }


def clean_row(row):
    return {
        key: value
        for key, value in row.items()
        if not key.startswith("_")
    }


def run(
    config,
    config_path: Path,
    experiment_path_override=None,
    readout_bounds_override=None,
):
    base_config_path = resolve_config_path(
        config["base_residual_config"],
        config_path,
    )
    base_config = json.loads(
        base_config_path.read_text(encoding="utf-8")
    )

    manifest_path = residual.resolve_config_path(
        base_config["manifest"],
        base_config_path,
    )
    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    cases = shared.normalize_manifest(manifest, manifest_path)
    repeat_case = ident.case_by_label(
        cases,
        base_config["repeat_case_label"],
    )
    summary = json.loads(
        repeat_case["summary"].read_text(encoding="utf-8")
    )
    comparison, experiment_path, comparison_source = (
        ident.comparison_for_case(repeat_case)
    )
    if experiment_path_override is not None:
        experiment_path = experiment_path_override

    snapshot_path = residual.resolve_config_path(
        base_config["residual_baseline_snapshot"],
        base_config_path,
    )
    snapshot = residual.load_baseline_snapshot(snapshot_path)

    inherited_candidate = dict(summary["best_case_parameters"])
    baseline_detector = dict(inherited_candidate)
    baseline_detector["R"] = float(snapshot["R_TES_Ohm"])
    baseline_detector.update(snapshot["detector_candidate"])

    full_args = base.fit_args(summary)
    fit_min = float(config["fit_region_Hz"]["min"])
    fit_max = float(config["fit_region_Hz"]["max"])
    full_frequency = np.geomspace(
        fit_min,
        fit_max,
        int(full_args.fit_points),
    )
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        full_frequency,
    )

    repair = repair_narrow_lines(
        target_context["full_frequency_Hz"],
        target_context["pre_analysis_asd"],
        config["line_repair"],
    )
    raw_target = np.interp(
        full_frequency,
        target_context["full_frequency_Hz"],
        repair["raw_normalized"],
    )
    continuum_target = np.interp(
        full_frequency,
        target_context["full_frequency_Hz"],
        repair["continuum_normalized"],
    )
    fit_args = competition.band_args(
        full_args,
        fit_min,
        fit_max,
        len(full_frequency),
    )

    detector_names = tuple(
        base_config["detector_nuisance"]["parameters"]
    )
    detector_bounds = competition.nuisance_bounds(
        inherited_candidate,
        {"bounds": base_config["detector_bounds"]},
    )
    readout_bounds = {
        name: tuple(
            float(value)
            for value in base_config["readout_bounds"][name]
        )
        for name in residual.READOUT_PARAMETER_NAMES
    }
    if readout_bounds_override is not None:
        unknown = sorted(
            set(readout_bounds_override)
            - set(residual.READOUT_PARAMETER_NAMES)
        )
        if unknown:
            raise ValueError(
                f"unknown readout bound overrides: {unknown}"
            )
        for name, values in readout_bounds_override.items():
            if len(values) != 2:
                raise ValueError(
                    f"readout bound override for {name} must have two values"
                )
            lower, upper = (float(values[0]), float(values[1]))
            if upper < lower:
                raise ValueError(
                    f"readout bound override for {name} is reversed"
                )
            readout_bounds[name] = (lower, upper)
    reference_readout = snapshot["reference_readout"]
    local_readout = snapshot["local_readout"]
    scale_hz = float(
        base_config["readout_parameterization"][
            "reference_scale_Hz"
        ]
    )
    tracked_white = float(snapshot["white_asd_A_rtHz"])
    optimizer_cfg = base_config["optimizer"]

    family_by_name = {
        row["name"]: row
        for row in config["families"]
    }
    required = {"pole_section_plus_c2", "full_order2"}
    if not required.issubset(family_by_name):
        raise ValueError(
            f"continuum config must contain families {sorted(required)}"
        )

    # A direct fixed-white reconstruction of the current minimal family is a
    # useful comparator: it isolates the gain from the new line repair and
    # white-floor profiling without adding a readout degree of freedom.
    fixed_current = residual.fit_family(
        family={
            "name": "pole_section_plus_c2_fixed_white",
            "readout_parameters": family_by_name[
                "pole_section_plus_c2"
            ]["readout_parameters"],
        },
        baseline_detector=baseline_detector,
        reference_readout=reference_readout,
        local_readout=local_readout,
        detector_names=detector_names,
        detector_bounds=detector_bounds,
        readout_bounds=readout_bounds,
        frequency=full_frequency,
        target=continuum_target,
        args=fit_args,
        scale_hz=scale_hz,
        white_asd=tracked_white,
        optimizer_cfg=optimizer_cfg,
        seed=int(optimizer_cfg["seed"]) + 3000,
    )

    pole = fit_with_profiled_white(
        family=family_by_name["pole_section_plus_c2"],
        baseline_detector=baseline_detector,
        reference_readout=reference_readout,
        local_readout=local_readout,
        detector_names=detector_names,
        detector_bounds=detector_bounds,
        readout_bounds=readout_bounds,
        frequency=full_frequency,
        target=continuum_target,
        fit_args=fit_args,
        scale_hz=scale_hz,
        baseline_white_asd=tracked_white,
        optimizer_cfg=optimizer_cfg,
        profile_cfg=config["white_profile"],
        seed=int(optimizer_cfg["seed"]) + 4000,
        warm_solutions=[
            (
                "fixed_white_continuum_fit",
                fixed_current["_detector_full"],
                fixed_current["_readout_full"],
            )
        ],
    )

    full = fit_with_profiled_white(
        family=family_by_name["full_order2"],
        baseline_detector=baseline_detector,
        reference_readout=reference_readout,
        local_readout=local_readout,
        detector_names=detector_names,
        detector_bounds=detector_bounds,
        readout_bounds=readout_bounds,
        frequency=full_frequency,
        target=continuum_target,
        fit_args=fit_args,
        scale_hz=scale_hz,
        baseline_white_asd=tracked_white,
        optimizer_cfg=optimizer_cfg,
        profile_cfg=config["white_profile"],
        seed=int(optimizer_cfg["seed"]) + 5000,
        warm_solutions=[
            (
                "profiled_pole_section_plus_c2",
                pole["_detector_full"],
                pole["_readout_full"],
            )
        ],
    )

    selection = choose_recommended(
        pole,
        full,
        config["selection_screen"],
    )
    selected = (
        full
        if selection["selected_family"] == full["name"]
        else pole
    )

    for row in (fixed_current, pole, full):
        white = float(
            row.get(
                "profiled_white_asd_A_rtHz",
                tracked_white,
            )
        )
        model, point = residual.model_for_candidate(
            row["_detector_full"],
            row["_readout_full"],
            full_frequency,
            scale_hz,
            white,
        )
        if model is None:
            raise RuntimeError(
                f"full-band model reconstruction failed for {row['name']}: {point}"
            )
        row["_model_full"] = model
        row["continuum_metrics_full"] = holdout.model_metrics(
            model,
            continuum_target,
            full_frequency,
            fit_args,
        )
        row["raw_experiment_metrics_full"] = holdout.model_metrics(
            model,
            raw_target,
            full_frequency,
            fit_args,
        )
        row["continuum_metrics_1_40k"] = metrics_in_region(
            model,
            continuum_target,
            full_frequency,
            full_args,
            1000.0,
            40000.0,
        )
        row["continuum_metrics_40_200k"] = metrics_in_region(
            model,
            continuum_target,
            full_frequency,
            full_args,
            40000.0,
            200000.0,
        )

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "no_notch_applied_to_production_data": True,
        "mask_selected_before_model_fit": True,
        "repeat_case": {
            "label": repeat_case["label"],
            "comparison_source": comparison_source,
            "experiment_path": str(experiment_path),
            "accepted_records": int(
                target_context["accepted_records"]
            ),
        },
        "line_repair": {
            "definition": (
                "Detect native-grid narrow peaks from experiment only; "
                "replace only their fixed intervals by the local median "
                "log-ASD floor before interpolation to the optimizer grid."
            ),
            "config": config["line_repair"],
            "native_fft_bin_Hz": repair["native_bin_Hz"],
            "accepted_narrow_lines": repair["accepted_lines"],
            "rejected_broad_features": repair[
                "rejected_broad_features"
            ],
            "repaired_fraction_full_native_grid": repair[
                "repaired_fraction_full_grid"
            ],
            "repaired_fraction_detection_band": repair[
                "repaired_fraction_detection_band"
            ],
        },
        "fit": {
            "min_Hz": fit_min,
            "max_Hz": fit_max,
            "points": len(full_frequency),
            "target": "line-repaired smooth continuum",
            "broad_features_preserved": True,
            "tracked_R_TES_Ohm": float(
                baseline_detector["R"]
            ),
            "detector_parameters_reoptimized": list(
                detector_names
            ),
            "tracked_white_asd_A_rtHz": tracked_white,
        },
        "fits": {
            "pole_section_plus_c2_fixed_white": clean_row(
                fixed_current
            ),
            "pole_section_plus_c2": clean_row(pole),
            "full_order2": clean_row(full),
        },
        "selection": selection,
        "recommended": {
            "family": selected["name"],
            "detector_candidate": {
                key: float(value)
                for key, value in selected[
                    "_detector_full"
                ].items()
                if key in detector_names
            },
            "R_TES_Ohm": float(
                selected["_detector_full"]["R"]
            ),
            "readout": {
                key: float(value)
                for key, value in selected[
                    "_readout_full"
                ].items()
            },
            "white_asd_A_rtHz": float(
                selected.get(
                    "profiled_white_asd_A_rtHz",
                    tracked_white,
                )
            ),
            "continuum_rms_dB": float(
                selected["continuum_metrics_full"][
                    "residual_metrics"
                ]["rms_residual_dB"]
            ),
            "continuum_1_40k_rms_dB": float(
                selected["continuum_metrics_1_40k"][
                    "residual_metrics"
                ]["rms_residual_dB"]
            ),
            "continuum_40_200k_rms_dB": float(
                selected["continuum_metrics_40_200k"][
                    "residual_metrics"
                ]["rms_residual_dB"]
            ),
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "base_residual_config": str(base_config_path),
            "manifest": str(manifest_path),
            "residual_baseline_snapshot": str(snapshot_path),
        },
        "_plot": {
            "frequency_Hz": full_frequency.tolist(),
            "raw_target": raw_target.tolist(),
            "continuum_target": continuum_target.tolist(),
            "fixed_white_model": fixed_current[
                "_model_full"
            ].tolist(),
            "pole_model": pole["_model_full"].tolist(),
            "full_model": full["_model_full"].tolist(),
            "recommended_family": selection[
                "selected_family"
            ],
        },
    }
    return result


def make_plot(result, output: Path, show=False):
    import matplotlib.pyplot as plt

    plot = result["_plot"]
    frequency = np.asarray(plot["frequency_Hz"], dtype=float)
    raw = np.asarray(plot["raw_target"], dtype=float)
    continuum = np.asarray(plot["continuum_target"], dtype=float)
    models = {
        "fixed white pole+Q+c2": np.asarray(
            plot["fixed_white_model"], dtype=float
        ),
        "profiled white pole+Q+c2": np.asarray(
            plot["pole_model"], dtype=float
        ),
        "profiled white full order-2": np.asarray(
            plot["full_model"], dtype=float
        ),
    }

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(11, 9),
        sharex=True,
        gridspec_kw={"height_ratios": [1.6, 1.0, 0.8]},
    )
    top, middle, bottom = axes

    top.plot(
        frequency,
        raw,
        linewidth=1.2,
        alpha=0.55,
        label="Experiment (raw target)",
    )
    top.plot(
        frequency,
        continuum,
        linewidth=2.0,
        label="Line-repaired continuum target",
    )
    for label, model in models.items():
        top.plot(
            frequency,
            model,
            linewidth=1.35,
            label=label,
        )
    for row in result["line_repair"]["accepted_narrow_lines"]:
        top.axvspan(
            row["repair_min_Hz"],
            row["repair_max_Hz"],
            alpha=0.08,
        )
    top.set_xscale("log")
    top.set_yscale("log")
    top.set_ylabel("Normalized ASD")
    top.set_title("Narrow-line-robust continuum fit")
    top.legend(frameon=False, fontsize=8)
    top.grid(True, which="both", alpha=0.25)

    for label, model in models.items():
        residual_db = 20.0 * np.log10(model / continuum)
        middle.plot(
            frequency,
            residual_db,
            linewidth=1.2,
            label=label,
        )
    middle.axhline(0.0, linewidth=0.9)
    middle.axhline(1.0, linestyle="--", linewidth=0.7)
    middle.axhline(-1.0, linestyle="--", linewidth=0.7)
    middle.set_xscale("log")
    middle.set_ylabel("Model / continuum [dB]")
    middle.grid(True, which="both", alpha=0.25)
    middle.legend(frameon=False, fontsize=8)

    repair_delta = 20.0 * np.log10(
        np.maximum(raw, np.finfo(float).tiny)
        / np.maximum(continuum, np.finfo(float).tiny)
    )
    bottom.plot(frequency, repair_delta, linewidth=1.0)
    bottom.axhline(0.0, linewidth=0.8)
    bottom.set_xscale("log")
    bottom.set_xlabel("Frequency [Hz]")
    bottom.set_ylabel("Raw / floor [dB]")
    bottom.set_title(
        "Only fixed narrow-line intervals are replaced"
    )
    bottom.grid(True, which="both", alpha=0.25)

    fig.suptitle(
        (
            f"{result['repeat_case']['label']} | "
            f"recommended: {result['recommended']['family']} | "
            f"continuum RMS "
            f"{result['recommended']['continuum_rms_dB']:.3f} dB"
        ),
        fontsize=10,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def cleaned_result(result):
    return {
        key: value
        for key, value in result.items()
        if key != "_plot"
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
    )
    parser.add_argument("--experiment-path", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    config = json.loads(
        args.config.read_text(encoding="utf-8")
    )
    result = run(
        config,
        args.config,
        experiment_path_override=args.experiment_path,
    )
    make_plot(result, args.figure, show=args.show)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            cleaned_result(result),
            indent=2,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "accepted_records": result["repeat_case"][
                    "accepted_records"
                ],
                "repaired_lines": [
                    {
                        "frequency_Hz": row["frequency_Hz"],
                        "width_Hz": row[
                            "half_prominence_width_Hz"
                        ],
                        "excess_dB": row[
                            "excess_over_local_baseline_dB"
                        ],
                    }
                    for row in result["line_repair"][
                        "accepted_narrow_lines"
                    ]
                ],
                "selection": result["selection"],
                "recommended": result["recommended"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

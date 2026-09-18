"""Cross-day repeatability test for the residual XXF-1 c2 shape freedom.

The documented Magnicon XXF-1 Connector Box LPF is represented as a second-order
Bessel section, c4 is fixed to zero, and the residual c2 numerator term is fit
independently on the 2024-12-05 repeat and 2024-12-06 reference acquisitions.

The diagnostic asks two deliberately narrow questions:

1. Is c2 materially required on both adjacent-day same-condition datasets?
2. If so, does the locally fitted Magnicon+c2 readout shape cross-apply to the
   other day after only detector nuisance parameters and the white floor are
   re-optimized?

For c4=0 and c2>0, the effective numerator magnitude is exactly

    sqrt(1 + c2 * (f / f_ref)^2),

which is algebraically identical in magnitude to one first-order zero with

    f_zero = f_ref / sqrt(c2).

That equivalent zero frequency is reported only as a compact coordinate for
repeatability.  It is not a claim that a physical zero has been identified.

Adjacent-day same-condition repeatability cannot distinguish a stable readout
feature from a detector-model deficiency that repeats at the same operating
point.  Independent operating-point validation or direct transfer measurement
is still required for a physical-origin claim.
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
DEFAULT_CONFIG = (
    CONFIG_DIR / "magnicon_xxf1_lpf_c2_cross_day_repeatability_config.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import line_robust_continuum_fit_diagnostic as continuum  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_necessity_diagnostic as c2diag  # noqa: E402
from subScript import magnicon_xxf1_lpf_c4zero_cutoff_profile_diagnostic as c4prof  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import readout_detector_state_competition_diagnostic as competition  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as residual  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def equivalent_first_order_zero_hz(c2: float, reference_scale_hz: float):
    """Magnitude-equivalent first-order zero for sqrt(1+c2*(f/f_ref)^2)."""
    c2 = float(c2)
    reference_scale_hz = float(reference_scale_hz)
    if reference_scale_hz <= 0.0:
        raise ValueError("reference_scale_hz must be positive")
    if c2 < 0.0:
        raise ValueError("c2 must be non-negative")
    if c2 == 0.0:
        return None
    return float(reference_scale_hz / np.sqrt(c2))


def symmetric_fractional_difference(a: float, b: float):
    a = float(a)
    b = float(b)
    if a <= 0.0 or b <= 0.0:
        return None
    return float(abs(a - b) / (0.5 * (a + b)))


def _problem_for_case(
    magnicon_config: dict,
    magnicon_config_path: Path,
    case_label: str,
    experiment_path_override=None,
):
    """Build the same continuum problem as xxf1._problem for a chosen case.

    The residual baseline snapshot is intentionally shared between the two
    same-condition days.  This preserves the detector/DC convention of the
    current Magnicon diagnostics while allowing detector nuisance parameters
    to refit locally on each target.
    """
    continuum_config_path = continuum.resolve_config_path(
        magnicon_config["base_continuum_config"],
        magnicon_config_path,
    )
    continuum_config = json.loads(
        continuum_config_path.read_text(encoding="utf-8")
    )
    residual_config_path = continuum.resolve_config_path(
        continuum_config["base_residual_config"],
        continuum_config_path,
    )
    residual_config = json.loads(
        residual_config_path.read_text(encoding="utf-8")
    )

    manifest_path = residual.resolve_config_path(
        residual_config["manifest"],
        residual_config_path,
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = shared.normalize_manifest(manifest, manifest_path)
    case = ident.case_by_label(cases, case_label)

    summary = json.loads(case["summary"].read_text(encoding="utf-8"))
    comparison, experiment_path, comparison_source = ident.comparison_for_case(case)
    if experiment_path_override is not None:
        experiment_path = Path(experiment_path_override)

    snapshot_path = residual.resolve_config_path(
        residual_config["residual_baseline_snapshot"],
        residual_config_path,
    )
    snapshot = residual.load_baseline_snapshot(snapshot_path)

    inherited_candidate = dict(summary["best_case_parameters"])
    baseline_detector = dict(inherited_candidate)
    baseline_detector["R"] = float(snapshot["R_TES_Ohm"])
    baseline_detector.update(snapshot["detector_candidate"])

    full_args = base.fit_args(summary)
    fit_min = float(continuum_config["fit_region_Hz"]["min"])
    fit_max = float(continuum_config["fit_region_Hz"]["max"])
    frequency = np.geomspace(fit_min, fit_max, int(full_args.fit_points))
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        frequency,
    )
    repair = continuum.repair_narrow_lines(
        target_context["full_frequency_Hz"],
        target_context["pre_analysis_asd"],
        continuum_config["line_repair"],
    )
    raw_target = np.interp(
        frequency,
        target_context["full_frequency_Hz"],
        repair["raw_normalized"],
    )
    target = np.interp(
        frequency,
        target_context["full_frequency_Hz"],
        repair["continuum_normalized"],
    )
    fit_args = competition.band_args(
        full_args,
        fit_min,
        fit_max,
        len(frequency),
    )

    detector_names = tuple(
        residual_config["detector_nuisance"]["parameters"]
    )
    detector_bounds = competition.nuisance_bounds(
        inherited_candidate,
        {"bounds": residual_config["detector_bounds"]},
    )
    readout_bounds = {
        name: tuple(
            float(value)
            for value in residual_config["readout_bounds"][name]
        )
        for name in residual.READOUT_PARAMETER_NAMES
    }
    readout_bounds["c2"] = (
        0.0,
        float(magnicon_config["readout_profile"]["c2_upper_bound"]),
    )

    return {
        "continuum_config_path": continuum_config_path,
        "residual_config_path": residual_config_path,
        "manifest_path": manifest_path,
        "snapshot_path": snapshot_path,
        "continuum_config": continuum_config,
        "residual_config": residual_config,
        "snapshot": snapshot,
        "repeat_case": case,
        "comparison_source": comparison_source,
        "experiment_path": experiment_path,
        "accepted_records": int(target_context["accepted_records"]),
        "baseline_detector": baseline_detector,
        "detector_names": detector_names,
        "detector_bounds": detector_bounds,
        "readout_bounds": readout_bounds,
        "frequency": frequency,
        "target": target,
        "raw_target": raw_target,
        "full_args": full_args,
        "fit_args": fit_args,
        "scale_hz": float(
            residual_config["readout_parameterization"]["reference_scale_Hz"]
        ),
        "tracked_white": float(snapshot["white_asd_A_rtHz"]),
        "optimizer_cfg": residual_config["optimizer"],
        "repair": repair,
    }


def _fit_local_case(
    *,
    problem: dict,
    c2_config: dict,
    magnicon_config: dict,
    seed_offset: int,
):
    bounds = dict(problem["readout_bounds"])
    bounds["c2"] = (0.0, float(c2_config["c2_upper_bound"]))
    snapshot = problem["snapshot"]

    current_reference, free_zero_reference = c4prof._fit_references(
        problem,
        snapshot,
        bounds,
    )
    nominal_cutoff = float(
        magnicon_config["magnicon_filter"]["nominal_cutoff_Hz"]
    )
    tolerance = float(
        magnicon_config["magnicon_filter"]["cutoff_tolerance_fraction"]
    )
    norm = str(c2_config["primary_normalization"])
    branch = c2diag._fit_branch(
        problem=problem,
        norm=norm,
        nominal_cutoff=nominal_cutoff,
        tolerance=tolerance,
        bounds=bounds,
        current_reference=current_reference,
        free_zero_reference=free_zero_reference,
        screen=c2_config["materiality_screen"],
        seed_offset=seed_offset,
    )
    return {
        "bounds": bounds,
        "current_reference": current_reference,
        "free_zero_reference": free_zero_reference,
        "branch": branch,
    }


def _cross_apply(
    *,
    source_local: dict,
    target_problem: dict,
    target_local: dict,
    seed_offset: int,
):
    source_readout = {
        key: float(value)
        for key, value in source_local["branch"]["profiled_c2_free"][
            "_readout_full"
        ].items()
    }
    row = xxf1._fit_variant(
        problem=target_problem,
        name="cross_applied_magnicon_c2_readout",
        readout_parameters=(),
        reference_readout=source_readout,
        local_readout=source_readout,
        readout_bounds=target_local["bounds"],
        seed_offset=seed_offset,
        warm_solutions=[
            (
                "target_local_detector_with_source_readout",
                target_local["branch"]["profiled_c2_free"]["_detector_full"],
                source_readout,
            )
        ],
    )
    return row


def _cross_apply_pass(comparison: dict, screen: dict) -> bool:
    return bool(
        float(comparison["continuum_rms_delta_dB"])
        <= float(screen["max_cross_applied_rms_degradation_dB"])
        and float(comparison["shape_score_ratio"])
        <= float(screen["max_cross_applied_shape_score_ratio"])
    )


def _c2_required(branch: dict) -> bool:
    return bool(
        branch["profiled_c2_material_improvement"]
        and not branch["profiled_c2_zero_passes_current_screen"]
        and branch["profiled_c2_free_passes_current_screen"]
    )


def _local_summary(
    problem: dict,
    local: dict,
):
    branch = local["branch"]
    free = branch["profiled_c2_free"]
    zero = branch["profiled_c2_zero"]
    c2 = float(free["readout"]["c2"])
    zero_hz = equivalent_first_order_zero_hz(c2, problem["scale_hz"])
    return {
        "label": problem["repeat_case"]["label"],
        "role": problem["repeat_case"].get("role"),
        "comparison_source": problem["comparison_source"],
        "experiment_path": str(problem["experiment_path"]),
        "accepted_records": int(problem["accepted_records"]),
        "n_repaired_lines": int(len(problem["repair"]["accepted_lines"])),
        "readout_reference_scale_Hz": float(problem["scale_hz"]),
        "profiled_cutoff_Hz": float(
            branch["profiled_c2_free_cutoff_Hz"]
        ),
        "profiled_cutoff_boundary_state": branch[
            "profiled_c2_free_cutoff_boundary_state"
        ],
        "profiled_c2": c2,
        "equivalent_first_order_zero_Hz": zero_hz,
        "equivalent_zero_semantics": (
            "algebraic magnitude coordinate only: "
            "sqrt(1+c2*(f/f_ref)^2) = sqrt(1+(f/f_zero)^2)"
        ),
        "c2_required_by_nested_screen": _c2_required(branch),
        "profiled_c2_gain": branch["profiled_c2_gain"],
        "c2_zero_passes_current_screen": bool(
            branch["profiled_c2_zero_passes_current_screen"]
        ),
        "c2_free_passes_current_screen": bool(
            branch["profiled_c2_free_passes_current_screen"]
        ),
        "c2_zero_fit": xxf1._clean(zero),
        "c2_free_fit": xxf1._clean(free),
    }


def _classify(
    *,
    reference_summary: dict,
    repeat_summary: dict,
    cross_reference_from_repeat: dict,
    cross_repeat_from_reference: dict,
    screen: dict,
):
    both_required = bool(
        reference_summary["c2_required_by_nested_screen"]
        and repeat_summary["c2_required_by_nested_screen"]
    )

    ref_zero = reference_summary["equivalent_first_order_zero_Hz"]
    rep_zero = repeat_summary["equivalent_first_order_zero_Hz"]
    zero_fractional_difference = (
        symmetric_fractional_difference(ref_zero, rep_zero)
        if ref_zero is not None and rep_zero is not None
        else None
    )
    zero_repeatable = bool(
        zero_fractional_difference is not None
        and zero_fractional_difference
        <= float(screen["max_equivalent_zero_fractional_difference"])
    )

    ref_cross_pass = _cross_apply_pass(
        cross_reference_from_repeat,
        screen,
    )
    rep_cross_pass = _cross_apply_pass(
        cross_repeat_from_reference,
        screen,
    )
    both_cross_pass = bool(ref_cross_pass and rep_cross_pass)

    if both_required and zero_repeatable and both_cross_pass:
        classification = (
            "c2_residual_repeatable_across_adjacent_days_"
            "supports_stable_readout_shape_candidate"
        )
    elif both_required and both_cross_pass:
        classification = (
            "c2_residual_cross_applies_but_local_coordinate_drift_detected"
        )
    elif both_required:
        classification = "c2_required_but_not_cross_day_repeatable"
    else:
        classification = "c2_requirement_not_reproduced_on_both_days"

    return {
        "classification": classification,
        "c2_required_on_both_days": both_required,
        "equivalent_zero_fractional_difference": zero_fractional_difference,
        "equivalent_zero_repeatable_within_screen": zero_repeatable,
        "repeat_readout_cross_applied_to_reference_passes": ref_cross_pass,
        "reference_readout_cross_applied_to_repeat_passes": rep_cross_pass,
        "bidirectional_cross_application_passes": both_cross_pass,
        "cross_day_screen": screen,
    }


def run(config: dict, config_path: Path):
    c2_config_path = resolve_config_path(
        config["base_c2_config"],
        config_path,
    )
    c2_config = json.loads(
        c2_config_path.read_text(encoding="utf-8")
    )
    profile_config_path = c2diag.resolve_config_path(
        c2_config["base_profile_config"],
        c2_config_path,
    )
    profile_config = json.loads(
        profile_config_path.read_text(encoding="utf-8")
    )
    magnicon_config_path = c4prof.resolve_config_path(
        profile_config["base_magnicon_config"],
        profile_config_path,
    )
    magnicon_config = json.loads(
        magnicon_config_path.read_text(encoding="utf-8")
    )

    requested_norm = str(config["normalization"])
    base_primary = str(c2_config["primary_normalization"])
    if requested_norm != base_primary:
        c2_config = dict(c2_config)
        c2_config["primary_normalization"] = requested_norm

    reference_problem = _problem_for_case(
        magnicon_config,
        magnicon_config_path,
        str(config["reference_case_label"]),
    )
    repeat_problem = _problem_for_case(
        magnicon_config,
        magnicon_config_path,
        str(config["repeat_case_label"]),
    )

    if not np.allclose(
        reference_problem["frequency"],
        repeat_problem["frequency"],
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError("reference and repeat fit grids differ")
    if not np.isclose(
        reference_problem["scale_hz"],
        repeat_problem["scale_hz"],
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError("reference and repeat readout scales differ")

    reference_local = _fit_local_case(
        problem=reference_problem,
        c2_config=c2_config,
        magnicon_config=magnicon_config,
        seed_offset=9100,
    )
    repeat_local = _fit_local_case(
        problem=repeat_problem,
        c2_config=c2_config,
        magnicon_config=magnicon_config,
        seed_offset=9700,
    )

    reference_from_repeat = _cross_apply(
        source_local=repeat_local,
        target_problem=reference_problem,
        target_local=reference_local,
        seed_offset=10300,
    )
    repeat_from_reference = _cross_apply(
        source_local=reference_local,
        target_problem=repeat_problem,
        target_local=repeat_local,
        seed_offset=10400,
    )

    ref_local_best = reference_local["branch"]["profiled_c2_free"]
    rep_local_best = repeat_local["branch"]["profiled_c2_free"]
    comp_ref_from_repeat = c4prof._comparison_to(
        reference_from_repeat,
        ref_local_best,
    )
    comp_rep_from_reference = c4prof._comparison_to(
        repeat_from_reference,
        rep_local_best,
    )

    reference_summary = _local_summary(
        reference_problem,
        reference_local,
    )
    repeat_summary = _local_summary(
        repeat_problem,
        repeat_local,
    )
    interpretation = _classify(
        reference_summary=reference_summary,
        repeat_summary=repeat_summary,
        cross_reference_from_repeat=comp_ref_from_repeat,
        cross_repeat_from_reference=comp_rep_from_reference,
        screen=config["cross_day_screen"],
    )

    c2_ref = float(reference_summary["profiled_c2"])
    c2_rep = float(repeat_summary["profiled_c2"])
    c2_ratio = (
        float(c2_rep / c2_ref)
        if c2_ref > 0.0
        else None
    )

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "measurement_time_filter_state_known": False,
        "tested_hypothesis": (
            "Conditional on the XXF-1 10 kHz second-order Bessel LPF being ON, "
            "test whether the residual c2 shape term is reproducible across "
            "the 2024-12-05 repeat and 2024-12-06 reference acquisitions."
        ),
        "normalization": requested_norm,
        "manual_reference": magnicon_config["manual_reference"],
        "repeatability_semantics": {
            "detector_nuisance": (
                "refit independently per day using the same residual-model "
                "bounds and the same fixed R_TES/DC baseline convention"
            ),
            "white_floor": "profiled independently per day",
            "line_repair": (
                "narrow experimental lines detected and repaired independently "
                "per day before continuum fitting"
            ),
            "readout_local_fit": (
                "Magnicon Bessel cutoff constrained to manual +/-2.5%; "
                "Bessel Q fixed by normalization; c4=0; c2 free"
            ),
            "cross_application": (
                "source-day cutoff/Q/c2/c4 frozen on target day; only target "
                "detector nuisance and white floor are re-optimized"
            ),
        },
        "reference_day": reference_summary,
        "repeat_day": repeat_summary,
        "cross_day_coordinate_comparison": {
            "reference_c2": c2_ref,
            "repeat_c2": c2_rep,
            "repeat_over_reference_c2_ratio": c2_ratio,
            "reference_equivalent_zero_Hz": reference_summary[
                "equivalent_first_order_zero_Hz"
            ],
            "repeat_equivalent_zero_Hz": repeat_summary[
                "equivalent_first_order_zero_Hz"
            ],
            "equivalent_zero_fractional_difference": interpretation[
                "equivalent_zero_fractional_difference"
            ],
        },
        "cross_application": {
            "repeat_readout_on_reference": {
                "source_readout": {
                    key: float(value)
                    for key, value in repeat_local["branch"][
                        "profiled_c2_free"
                    ]["readout"].items()
                },
                "target_refit": xxf1._clean(reference_from_repeat),
                "comparison_to_reference_local_best": comp_ref_from_repeat,
                "passes_screen": interpretation[
                    "repeat_readout_cross_applied_to_reference_passes"
                ],
            },
            "reference_readout_on_repeat": {
                "source_readout": {
                    key: float(value)
                    for key, value in reference_local["branch"][
                        "profiled_c2_free"
                    ]["readout"].items()
                },
                "target_refit": xxf1._clean(repeat_from_reference),
                "comparison_to_repeat_local_best": comp_rep_from_reference,
                "passes_screen": interpretation[
                    "reference_readout_cross_applied_to_repeat_passes"
                ],
            },
        },
        "interpretation": {
            **interpretation,
            "origin_guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_c2_config": str(c2_config_path),
            "base_profile_config": str(profile_config_path),
            "base_magnicon_config": str(magnicon_config_path),
            "manifest": str(reference_problem["manifest_path"]),
            "residual_baseline_snapshot": str(
                reference_problem["snapshot_path"]
            ),
        },
        "_plot": {
            "frequency_Hz": reference_problem["frequency"].tolist(),
            "reference_target": reference_problem["target"].tolist(),
            "repeat_target": repeat_problem["target"].tolist(),
            "reference_local": ref_local_best["_model_full"].tolist(),
            "repeat_local": rep_local_best["_model_full"].tolist(),
            "reference_from_repeat": reference_from_repeat[
                "_model_full"
            ].tolist(),
            "repeat_from_reference": repeat_from_reference[
                "_model_full"
            ].tolist(),
        },
    }
    return result


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    f = np.asarray(p["frequency_Hz"], dtype=float)
    ref_target = np.asarray(p["reference_target"], dtype=float)
    rep_target = np.asarray(p["repeat_target"], dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.0), sharex=True)
    top, middle, bottom = axes

    top.loglog(f, ref_target, label="12/06 continuum target")
    top.loglog(f, np.asarray(p["reference_local"]), label="12/06 local readout")
    top.loglog(f, rep_target, label="12/05 continuum target")
    top.loglog(f, np.asarray(p["repeat_local"]), label="12/05 local readout")
    top.set_ylabel("Normalized ASD")
    top.legend(frameon=False, fontsize=8)
    top.grid(True, which="both", alpha=0.2)

    curves = [
        ("12/06 local", p["reference_local"], ref_target),
        ("12/06 with 12/05 readout", p["reference_from_repeat"], ref_target),
        ("12/05 local", p["repeat_local"], rep_target),
        ("12/05 with 12/06 readout", p["repeat_from_reference"], rep_target),
    ]
    for label, model, target in curves:
        residual_db = 20.0 * np.log10(
            np.asarray(model, dtype=float) / target
        )
        middle.semilogx(f, residual_db, label=label)
    middle.axhline(0.0, linewidth=1.0)
    middle.set_ylabel("Model / continuum [dB]")
    middle.legend(frameon=False, fontsize=8)
    middle.grid(True, which="both", alpha=0.2)

    labels = ["12/06", "12/05"]
    c2_values = [
        result["reference_day"]["profiled_c2"],
        result["repeat_day"]["profiled_c2"],
    ]
    zero_values = [
        result["reference_day"]["equivalent_first_order_zero_Hz"] / 1000.0,
        result["repeat_day"]["equivalent_first_order_zero_Hz"] / 1000.0,
    ]
    x = np.arange(2, dtype=float)
    bottom.plot(x, c2_values, marker="o", label="c2")
    bottom.plot(
        x,
        zero_values,
        marker="o",
        label="equiv. zero [kHz]",
    )
    bottom.set_xticks(x, labels)
    bottom.set_ylabel("Cross-day coordinates")
    bottom.set_xlabel("Acquisition")
    bottom.grid(True, alpha=0.2)
    bottom.legend(frameon=False)

    fig.suptitle(result["interpretation"]["classification"], fontsize=10)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def cleaned_result(result: dict) -> dict:
    return {key: value for key, value in result.items() if key != "_plot"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = run(config, args.config)
    make_plot(result, args.figure, show=args.show)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(cleaned_result(result), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "classification": result["interpretation"]["classification"],
                "reference_c2": result["reference_day"]["profiled_c2"],
                "repeat_c2": result["repeat_day"]["profiled_c2"],
                "reference_equivalent_zero_Hz": result["reference_day"][
                    "equivalent_first_order_zero_Hz"
                ],
                "repeat_equivalent_zero_Hz": result["repeat_day"][
                    "equivalent_first_order_zero_Hz"
                ],
                "bidirectional_cross_application_passes": result[
                    "interpretation"
                ]["bidirectional_cross_application_passes"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

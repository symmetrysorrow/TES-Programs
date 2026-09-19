"""Identify which TES nuisance freedoms compensate the residual c2 term.

The preceding diagnostics established two useful facts:
  * the fitted Magnicon Bessel denominator is essentially identical across the
    adjacent-day datasets, while
  * the phenomenological c2 numerator is not directly repeatable.

This diagnostic therefore stops treating c2 as a putative fixed readout
component and asks a narrower identifiability question.  For each day:

  1. fit the current local Magnicon+Bessel+c2 solution;
  2. freeze that day's Magnicon pole/Q and profiled white floor;
  3. set c2=0 and keep every TES nuisance frozen;
  4. repeat c2=0 while releasing alpha, beta, C_tes, L, or T_bath one at a time;
  5. repeat once with all five TES nuisance parameters released.

The recovery of fit quality from steps 4--5 measures compensation capacity, not
physical causation.  A parameter that recovers the c2 penalty may simply be
correlated with a missing readout or detector-model term.
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
DEFAULT_CONFIG = CONFIG_DIR / "magnicon_xxf1_c2_detector_compensation_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_c2_detector_compensation_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_c2_detector_compensation_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_necessity_diagnostic as c2diag  # noqa: E402
from subScript import magnicon_xxf1_lpf_c4zero_cutoff_profile_diagnostic as c4prof  # noqa: E402
from subScript import readout_lowmid_holdout_diagnostic as holdout  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as residual  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def _rms_from_metrics(metrics: dict) -> float:
    return float(metrics["residual_metrics"]["rms_residual_dB"])


def _rms_from_fit(row: dict) -> float:
    if "continuum_metrics_full" in row:
        return _rms_from_metrics(row["continuum_metrics_full"])
    return float(row["residual_metrics"]["rms_residual_dB"])


def recovery_summary(
    *,
    c2_free_rms: float,
    frozen_c2_zero_rms: float,
    all_detector_c2_zero_rms: float,
    single_rows: dict,
):
    c2_free_rms = float(c2_free_rms)
    frozen = float(frozen_c2_zero_rms)
    all_free = float(all_detector_c2_zero_rms)
    total_recovery = float(frozen - all_free)
    removal_penalty = float(frozen - c2_free_rms)
    remaining_penalty = float(all_free - c2_free_rms)

    rows = {}
    for name, rms in single_rows.items():
        rms = float(rms)
        recovered = float(frozen - rms)
        fraction = (
            float(recovered / total_recovery)
            if total_recovery > 1.0e-12
            else None
        )
        rows[str(name)] = {
            "rms_residual_dB": rms,
            "rms_recovery_from_frozen_c2_zero_dB": recovered,
            "fraction_of_all_detector_recovery": fraction,
            "remaining_rms_penalty_to_c2_free_dB": float(
                rms - c2_free_rms
            ),
        }

    top_name = None
    top_recovery = None
    if rows:
        top_name = max(
            rows,
            key=lambda name: rows[name][
                "rms_recovery_from_frozen_c2_zero_dB"
            ],
        )
        top_recovery = rows[top_name]

    return {
        "c2_free_rms_dB": c2_free_rms,
        "frozen_c2_zero_rms_dB": frozen,
        "all_detector_c2_zero_rms_dB": all_free,
        "c2_removal_penalty_with_detector_frozen_dB": removal_penalty,
        "all_detector_rms_recovery_dB": total_recovery,
        "remaining_rms_penalty_after_all_detector_refit_dB": (
            remaining_penalty
        ),
        "single_parameter_recovery": rows,
        "top_single_parameter": top_name,
        "top_single": top_recovery,
    }


def classify_recovery(summary: dict, screen: dict) -> dict:
    total = float(summary["all_detector_rms_recovery_dB"])
    top = summary.get("top_single")
    fraction = (
        None
        if top is None
        else top["fraction_of_all_detector_recovery"]
    )
    material = bool(
        total
        >= float(screen["min_material_all_detector_rms_recovery_dB"])
    )

    if not material:
        classification = "detector_nuisance_does_not_materially_compensate_c2"
    elif fraction is not None and fraction >= float(
        screen["single_parameter_strong_fraction_of_all_recovery"]
    ):
        classification = "one_parameter_recovers_most_detector_compensation"
    elif fraction is not None and fraction >= float(
        screen["single_parameter_mixed_fraction_of_all_recovery"]
    ):
        classification = "mixed_compensation_with_leading_parameter"
    else:
        classification = "distributed_or_correlated_detector_compensation"

    return {
        "classification": classification,
        "all_detector_recovery_is_material": material,
        "top_single_parameter": summary.get("top_single_parameter"),
        "top_single_fraction_of_all_detector_recovery": fraction,
        "screen": screen,
    }


def parameter_shift_summary(
    c2_free_detector: dict,
    c2_zero_detector: dict,
    detector_names,
    detector_bounds: dict,
    boundary_hits: dict,
):
    rows = {}
    for name in detector_names:
        before = float(c2_free_detector[name])
        after = float(c2_zero_detector[name])
        lower, upper = (
            float(value) for value in detector_bounds[name]
        )
        span = float(upper - lower)
        ratio = (
            float(after / before)
            if before != 0.0
            else None
        )
        rows[name] = {
            "c2_free_value": before,
            "c2_zero_all_detector_value": after,
            "delta": float(after - before),
            "ratio_c2_zero_over_c2_free": ratio,
            "delta_as_fraction_of_allowed_span": (
                float((after - before) / span)
                if span > 0.0
                else None
            ),
            "c2_zero_boundary_hit": boundary_hits.get(name),
        }
    return rows


def _fixed_model_row(
    *,
    problem: dict,
    detector: dict,
    readout: dict,
    white_asd: float,
    name: str,
):
    model, point = residual.model_for_candidate(
        detector,
        readout,
        problem["frequency"],
        problem["scale_hz"],
        white_asd,
    )
    if model is None:
        raise RuntimeError(
            f"fixed model failed for {name}: {point}"
        )
    metrics = holdout.model_metrics(
        model,
        problem["target"],
        problem["frequency"],
        problem["fit_args"],
    )
    score = float(
        opt.fit_score(
            model,
            problem["target"],
            problem["frequency"],
            problem["fit_args"],
        )
    )
    return {
        "name": name,
        "shape_score": score,
        "residual_metrics": metrics["residual_metrics"],
        "bands": metrics["bands"],
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
        "_model": model,
    }


def _fit_detector_subset(
    *,
    problem: dict,
    baseline_detector: dict,
    zero_readout: dict,
    fixed_white_asd: float,
    detector_names,
    readout_bounds: dict,
    warm_detectors,
    seed_offset: int,
):
    names = tuple(detector_names)
    if not names:
        raise ValueError("detector subset must not be empty")

    warm_solutions = [
        (
            label,
            detector,
            zero_readout,
        )
        for label, detector in warm_detectors
    ]
    row = residual.fit_family(
        family={
            "name": (
                "magnicon_fixed_c2_zero_detector_"
                + "_".join(names)
            ),
            "readout_parameters": [],
        },
        baseline_detector=baseline_detector,
        reference_readout=zero_readout,
        local_readout=zero_readout,
        detector_names=names,
        detector_bounds=problem["detector_bounds"],
        readout_bounds=readout_bounds,
        frequency=problem["frequency"],
        target=problem["target"],
        args=problem["fit_args"],
        scale_hz=problem["scale_hz"],
        white_asd=fixed_white_asd,
        optimizer_cfg=problem["optimizer_cfg"],
        seed=int(problem["optimizer_cfg"]["seed"]) + int(seed_offset),
        warm_solutions=warm_solutions,
    )
    return row


def _clean_fit(row: dict):
    return {
        key: value
        for key, value in row.items()
        if not key.startswith("_")
    }


def _load_stack(config: dict, config_path: Path):
    base_cross_path = resolve_config_path(
        config["base_cross_day_config"],
        config_path,
    )
    base_cross = json.loads(
        base_cross_path.read_text(encoding="utf-8")
    )
    c2_config_path = crossday.resolve_config_path(
        base_cross["base_c2_config"],
        base_cross_path,
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
    return {
        "base_cross_path": base_cross_path,
        "base_cross": base_cross,
        "c2_config_path": c2_config_path,
        "c2_config": c2_config,
        "profile_config_path": profile_config_path,
        "magnicon_config_path": magnicon_config_path,
        "magnicon_config": magnicon_config,
    }


def _run_case(
    *,
    problem: dict,
    c2_config: dict,
    magnicon_config: dict,
    detector_names,
    screen: dict,
    seed_base: int,
):
    local = crossday._fit_local_case(
        problem=problem,
        c2_config=c2_config,
        magnicon_config=magnicon_config,
        seed_offset=seed_base,
    )
    branch = local["branch"]
    free = branch["profiled_c2_free"]
    prior_zero = branch["profiled_c2_zero"]

    c2_free_detector = dict(free["_detector_full"])
    c2_free_readout = dict(free["_readout_full"])
    fixed_white = float(free["profiled_white_asd_A_rtHz"])

    zero_readout = dict(c2_free_readout)
    zero_readout["c2"] = 0.0
    zero_readout["c4"] = 0.0

    frozen = _fixed_model_row(
        problem=problem,
        detector=c2_free_detector,
        readout=zero_readout,
        white_asd=fixed_white,
        name="c2_zero_detector_frozen_white_frozen",
    )

    warm_detectors = [
        ("c2_free_detector", c2_free_detector),
        ("prior_profiled_c2_zero_detector", prior_zero["_detector_full"]),
    ]

    single = {}
    for index, name in enumerate(detector_names):
        single[name] = _fit_detector_subset(
            problem=problem,
            baseline_detector=c2_free_detector,
            zero_readout=zero_readout,
            fixed_white_asd=fixed_white,
            detector_names=(name,),
            readout_bounds=local["bounds"],
            warm_detectors=warm_detectors,
            seed_offset=seed_base + 100 + 50 * index,
        )

    all_free = _fit_detector_subset(
        problem=problem,
        baseline_detector=c2_free_detector,
        zero_readout=zero_readout,
        fixed_white_asd=fixed_white,
        detector_names=tuple(detector_names),
        readout_bounds=local["bounds"],
        warm_detectors=warm_detectors,
        seed_offset=seed_base + 500,
    )

    summary = recovery_summary(
        c2_free_rms=_rms_from_fit(free),
        frozen_c2_zero_rms=float(
            frozen["residual_metrics"]["rms_residual_dB"]
        ),
        all_detector_c2_zero_rms=float(
            all_free["residual_metrics"]["rms_residual_dB"]
        ),
        single_rows={
            name: row["residual_metrics"]["rms_residual_dB"]
            for name, row in single.items()
        },
    )
    interpretation = classify_recovery(summary, screen)

    shifts = parameter_shift_summary(
        c2_free_detector,
        all_free["_detector_full"],
        detector_names,
        problem["detector_bounds"],
        all_free.get("detector_boundary_hits", {}),
    )

    return {
        "case": {
            "label": problem["repeat_case"]["label"],
            "role": problem["repeat_case"].get("role"),
            "comparison_source": problem["comparison_source"],
            "experiment_path": str(problem["experiment_path"]),
            "accepted_records": int(problem["accepted_records"]),
            "n_repaired_lines": int(
                len(problem["repair"]["accepted_lines"])
            ),
        },
        "fixed_readout_for_attribution": {
            "semantics": (
                "Magnicon pole/Q frozen to this day's local c2-free fit; "
                "c2 then forced to zero; c4=0"
            ),
            "pole_Hz": float(zero_readout["pole_Hz"]),
            "pole_Q": float(zero_readout["pole_Q"]),
            "c2": 0.0,
            "c4": 0.0,
        },
        "fixed_white_asd_A_rtHz": fixed_white,
        "c2_free_reference": {
            "readout": {
                key: float(value)
                for key, value in c2_free_readout.items()
            },
            "detector": {
                name: float(c2_free_detector[name])
                for name in detector_names
            },
            "white_asd_A_rtHz": fixed_white,
            "continuum_rms_dB": _rms_from_fit(free),
            "shape_score": float(free["shape_score"]),
        },
        "c2_zero_detector_frozen": _clean_fit(frozen),
        "single_parameter_refits": {
            name: _clean_fit(row)
            for name, row in single.items()
        },
        "all_detector_refit": _clean_fit(all_free),
        "prior_c2_zero_with_white_profiled": {
            "continuum_rms_dB": _rms_from_fit(prior_zero),
            "shape_score": float(prior_zero["shape_score"]),
            "white_asd_A_rtHz": float(
                prior_zero["profiled_white_asd_A_rtHz"]
            ),
        },
        "recovery": summary,
        "detector_shift_c2_free_to_c2_zero_all_refit": shifts,
        "interpretation": interpretation,
        "_plot": {
            "frequency_Hz": problem["frequency"].tolist(),
            "target": problem["target"].tolist(),
            "c2_free": free["_model_full"].tolist(),
            "c2_zero_frozen": frozen["_model"].tolist(),
            "c2_zero_all_detector": all_free["_model"].tolist(),
        },
    }


def run(config: dict, config_path: Path):
    stack = _load_stack(config, config_path)
    base_cross = stack["base_cross"]
    c2_config = dict(stack["c2_config"])
    c2_config["primary_normalization"] = str(
        base_cross["normalization"]
    )
    magnicon_config = stack["magnicon_config"]

    detector_names = tuple(config["detector_parameters"])
    expected = {"alpha", "beta", "C_tes", "L", "T_bath"}
    unknown = sorted(set(detector_names) - expected)
    if unknown:
        raise ValueError(f"unknown detector parameters: {unknown}")

    reference_problem = crossday._problem_for_case(
        magnicon_config,
        stack["magnicon_config_path"],
        str(base_cross["reference_case_label"]),
    )
    repeat_problem = crossday._problem_for_case(
        magnicon_config,
        stack["magnicon_config_path"],
        str(base_cross["repeat_case_label"]),
    )

    reference = _run_case(
        problem=reference_problem,
        c2_config=c2_config,
        magnicon_config=magnicon_config,
        detector_names=detector_names,
        screen=config["attribution_screen"],
        seed_base=12100,
    )
    repeat = _run_case(
        problem=repeat_problem,
        c2_config=c2_config,
        magnicon_config=magnicon_config,
        detector_names=detector_names,
        screen=config["attribution_screen"],
        seed_base=13100,
    )

    ref_top = reference["interpretation"]["top_single_parameter"]
    rep_top = repeat["interpretation"]["top_single_parameter"]

    return {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "tested_question": (
            "When phenomenological c2 is removed while Magnicon pole/Q and "
            "white are frozen, which TES nuisance freedoms can recover the "
            "lost continuum fit quality?"
        ),
        "attribution_semantics": {
            "readout": (
                "per day: freeze Magnicon pole/Q to local c2-free solution, "
                "then force c2=c4=0"
            ),
            "white": (
                "per day: freeze white ASD to the local c2-free profiled value"
            ),
            "detector": (
                "start from local c2-free detector solution; release each "
                "nuisance individually and finally all configured nuisances"
            ),
            "causal_claim_allowed": False,
        },
        "reference_day": {
            key: value
            for key, value in reference.items()
            if key != "_plot"
        },
        "repeat_day": {
            key: value
            for key, value in repeat.items()
            if key != "_plot"
        },
        "cross_day_summary": {
            "reference_top_single_parameter": ref_top,
            "repeat_top_single_parameter": rep_top,
            "same_top_single_parameter": bool(
                ref_top is not None and ref_top == rep_top
            ),
            "reference_classification": reference[
                "interpretation"
            ]["classification"],
            "repeat_classification": repeat[
                "interpretation"
            ]["classification"],
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "base_cross_day_config": str(stack["base_cross_path"]),
            "base_c2_config": str(stack["c2_config_path"]),
            "base_profile_config": str(stack["profile_config_path"]),
            "base_magnicon_config": str(
                stack["magnicon_config_path"]
            ),
        },
        "_plot": {
            "reference": reference["_plot"],
            "repeat": repeat["_plot"],
        },
    }


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10.0, 8.0))
    top, bottom = axes

    for label, key in [
        ("12/06", "reference_day"),
        ("12/05", "repeat_day"),
    ]:
        recovery = result[key]["recovery"]
        names = list(recovery["single_parameter_recovery"])
        values = [
            recovery["single_parameter_recovery"][name][
                "rms_recovery_from_frozen_c2_zero_dB"
            ]
            for name in names
        ]
        x = np.arange(len(names), dtype=float)
        top.plot(x, values, marker="o", label=label)
    top.set_xticks(x, names)
    top.set_ylabel("RMS recovery [dB]")
    top.set_title("Single TES-nuisance recovery after forcing c2=0")
    top.grid(True, alpha=0.2)
    top.legend(frameon=False)

    for label, plot_key in [
        ("12/06", "reference"),
        ("12/05", "repeat"),
    ]:
        p = result["_plot"][plot_key]
        f = np.asarray(p["frequency_Hz"], dtype=float)
        target = np.asarray(p["target"], dtype=float)
        for curve_label, curve_key in [
            ("c2 free", "c2_free"),
            ("c2=0 frozen TES", "c2_zero_frozen"),
            ("c2=0 all TES free", "c2_zero_all_detector"),
        ]:
            residual_db = 20.0 * np.log10(
                np.asarray(p[curve_key], dtype=float) / target
            )
            bottom.semilogx(
                f,
                residual_db,
                label=f"{label} {curve_label}",
            )
    bottom.axhline(0.0, linewidth=1.0)
    bottom.set_xlabel("Frequency [Hz]")
    bottom.set_ylabel("Model / continuum [dB]")
    bottom.grid(True, which="both", alpha=0.2)
    bottom.legend(frameon=False, fontsize=8)

    fig.suptitle("TES compensation of removed c2", fontsize=10)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def cleaned_result(result: dict):
    return {
        key: value
        for key, value in result.items()
        if key != "_plot"
    }


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
        json.dumps(cleaned_result(result), indent=2, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "reference": result["reference_day"][
                    "interpretation"
                ],
                "repeat": result["repeat_day"]["interpretation"],
                "cross_day_summary": result["cross_day_summary"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

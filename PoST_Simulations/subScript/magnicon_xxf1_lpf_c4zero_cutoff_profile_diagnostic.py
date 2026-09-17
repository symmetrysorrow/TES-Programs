"""Profile the documented XXF-1 10 kHz LPF with c4 fixed to zero.

This follows the first Magnicon LPF diagnostic, where the strongest physically
interesting branch was the documented second-order Bessel filter with c4 removed.
Here the filter cutoff is allowed to move only inside the manual tolerance
(10 kHz +/-2.5%), while Q remains the exact Bessel value for the chosen
normalization and c2 is the only residual smooth numerator freedom.

The measurement-time filter state remains unknown.  This diagnostic therefore
checks consistency with the LPF-ON hypothesis; it does not establish that the
Connector Box filter was actually enabled during the acquisition.
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
DEFAULT_CONFIG = CONFIG_DIR / "magnicon_xxf1_lpf_c4zero_cutoff_profile_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_c4zero_cutoff_profile_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_c4zero_cutoff_profile_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def _cutoff_boundary_state(value, lower, upper, relative_tolerance):
    value = float(value)
    lower = float(lower)
    upper = float(upper)
    span = max(upper - lower, abs(upper), 1.0)
    tolerance = span * float(relative_tolerance)
    at_lower = bool(abs(value - lower) <= tolerance)
    at_upper = bool(abs(value - upper) <= tolerance)
    return {
        "value_Hz": value,
        "lower_Hz": lower,
        "upper_Hz": upper,
        "at_lower": at_lower,
        "at_upper": at_upper,
        "interior": bool(not at_lower and not at_upper),
    }


def _comparison_to(row, reference):
    reference_score = float(reference["shape_score"])
    return {
        "shape_score_ratio": float(row["shape_score"] / reference_score),
        "continuum_rms_delta_dB": float(
            xxf1._rms(row) - xxf1._rms(reference)
        ),
        "continuum_1_40k_rms_delta_dB": float(
            xxf1._rms(row, "continuum_metrics_1_40k")
            - xxf1._rms(reference, "continuum_metrics_1_40k")
        ),
        "continuum_40_200k_rms_delta_dB": float(
            xxf1._rms(row, "continuum_metrics_40_200k")
            - xxf1._rms(reference, "continuum_metrics_40_200k")
        ),
    }


def _passes_current_reference(row, current_reference, screen):
    comp = _comparison_to(row, current_reference)
    return bool(
        comp["continuum_rms_delta_dB"]
        <= float(screen["max_rms_degradation_to_current_free_pole_dB"])
        and comp["shape_score_ratio"]
        <= float(screen["max_shape_score_ratio_to_current_free_pole"])
    )


def _passes_c4zero_reference(row, c4zero_reference, screen):
    comp = _comparison_to(row, c4zero_reference)
    return bool(
        comp["continuum_rms_delta_dB"]
        <= float(screen["max_rms_degradation_to_free_pole_c4zero_dB"])
    )


def _reference_summary(row):
    return {
        "readout": {
            key: float(value)
            for key, value in row["readout"].items()
        },
        "detector_candidate": {
            key: float(value)
            for key, value in row["detector_candidate"].items()
        },
        "white_asd_A_rtHz": float(row["profiled_white_asd_A_rtHz"]),
        "continuum_rms_dB": xxf1._rms(row),
        "continuum_1_40k_rms_dB": xxf1._rms(
            row, "continuum_metrics_1_40k"
        ),
        "continuum_40_200k_rms_dB": xxf1._rms(
            row, "continuum_metrics_40_200k"
        ),
        "shape_score": float(row["shape_score"]),
        "detector_boundary_hits": row.get("detector_boundary_hits", {}),
        "readout_boundary_hits": row.get("readout_boundary_hits", {}),
    }


def _fit_references(problem, snapshot, bounds):
    current = xxf1._fit_variant(
        problem=problem,
        name="current_free_pole_plus_c2_legacy_c4",
        readout_parameters=("pole_Hz", "pole_Q", "c2"),
        reference_readout=dict(snapshot["reference_readout"]),
        local_readout=dict(snapshot["local_readout"]),
        readout_bounds=bounds,
        seed_offset=7100,
    )

    reference_zero = dict(snapshot["reference_readout"])
    reference_zero["c4"] = 0.0
    local_zero = dict(snapshot["local_readout"])
    local_zero["c4"] = 0.0
    free_zero = xxf1._fit_variant(
        problem=problem,
        name="free_pole_plus_c2_c4_zero",
        readout_parameters=("pole_Hz", "pole_Q", "c2"),
        reference_readout=reference_zero,
        local_readout=local_zero,
        readout_bounds=bounds,
        seed_offset=7150,
        warm_solutions=[
            (
                "current_free_pole",
                current["_detector_full"],
                {
                    **current["_readout_full"],
                    "c4": 0.0,
                },
            )
        ],
    )
    return current, free_zero


def _fit_magnicon_branch(
    *,
    problem,
    norm,
    nominal_cutoff,
    tolerance,
    bounds,
    current_reference,
    free_zero_reference,
    screen,
    seed_offset,
):
    nominal_coordinates = xxf1.second_order_bessel_canonical(
        nominal_cutoff, norm
    )
    low_cutoff = nominal_cutoff * (1.0 - tolerance)
    high_cutoff = nominal_cutoff * (1.0 + tolerance)
    low_coordinates = xxf1.second_order_bessel_canonical(low_cutoff, norm)
    high_coordinates = xxf1.second_order_bessel_canonical(high_cutoff, norm)

    reference = {
        "pole_Hz": float(nominal_coordinates["pole_Hz"]),
        "pole_Q": float(nominal_coordinates["pole_Q"]),
        "c2": 0.0,
        "c4": 0.0,
    }

    nominal = xxf1._fit_variant(
        problem=problem,
        name=f"magnicon_{norm}_nominal_c4zero_plus_c2",
        readout_parameters=("c2",),
        reference_readout=reference,
        local_readout=reference,
        readout_bounds=bounds,
        seed_offset=seed_offset,
        warm_solutions=[
            (
                "free_pole_c4zero",
                free_zero_reference["_detector_full"],
                free_zero_reference["_readout_full"],
            )
        ],
    )

    tolerance_bounds = dict(bounds)
    tolerance_bounds["pole_Hz"] = (
        float(low_coordinates["pole_Hz"]),
        float(high_coordinates["pole_Hz"]),
    )
    profiled = xxf1._fit_variant(
        problem=problem,
        name=f"magnicon_{norm}_tolerance_c4zero_plus_c2",
        readout_parameters=("pole_Hz", "c2"),
        reference_readout=reference,
        local_readout=reference,
        readout_bounds=tolerance_bounds,
        seed_offset=seed_offset + 50,
        warm_solutions=[
            (
                "magnicon_nominal_c4zero",
                nominal["_detector_full"],
                nominal["_readout_full"],
            ),
            (
                "free_pole_c4zero",
                free_zero_reference["_detector_full"],
                free_zero_reference["_readout_full"],
            ),
        ],
    )

    fitted_cutoff = xxf1.equivalent_cutoff_from_pole(
        profiled["readout"]["pole_Hz"],
        nominal_coordinates,
    )
    boundary = _cutoff_boundary_state(
        fitted_cutoff,
        low_cutoff,
        high_cutoff,
        screen["cutoff_bound_relative_tolerance"],
    )

    nominal_comp_current = _comparison_to(nominal, current_reference)
    profiled_comp_current = _comparison_to(profiled, current_reference)
    profiled_comp_zero = _comparison_to(profiled, free_zero_reference)

    return {
        "norm": norm,
        "manual_cutoff_range_Hz": [float(low_cutoff), float(high_cutoff)],
        "nominal_coordinates": nominal_coordinates,
        "fitted_cutoff_Hz": float(fitted_cutoff),
        "fitted_cutoff_fraction_from_nominal": float(
            fitted_cutoff / nominal_cutoff - 1.0
        ),
        "cutoff_boundary_state": boundary,
        "nominal": nominal,
        "profiled": profiled,
        "comparison_nominal_to_current_free_pole": nominal_comp_current,
        "comparison_profiled_to_current_free_pole": profiled_comp_current,
        "comparison_profiled_to_free_pole_c4zero": profiled_comp_zero,
        "profile_improvement_over_nominal_rms_dB": float(
            xxf1._rms(nominal) - xxf1._rms(profiled)
        ),
        "passes_current_free_pole_screen": _passes_current_reference(
            profiled, current_reference, screen
        ),
        "passes_free_pole_c4zero_screen": _passes_c4zero_reference(
            profiled, free_zero_reference, screen
        ),
    }


def _classify(branches, screen):
    ranked = sorted(
        branches.values(),
        key=lambda row: float(row["profiled"]["shape_score"]),
    )
    best = ranked[0]
    pass_current = bool(best["passes_current_free_pole_screen"])
    pass_zero = bool(best["passes_free_pole_c4zero_screen"])
    interior = bool(best["cutoff_boundary_state"]["interior"])

    if pass_current and pass_zero and interior:
        classification = (
            "documented_lpf_c4zero_within_tolerance_replaces_free_pole_supported"
        )
    elif pass_current and pass_zero and not interior:
        classification = (
            "documented_lpf_c4zero_numerically_close_but_cutoff_hits_manual_tolerance"
        )
    elif pass_zero:
        classification = (
            "documented_lpf_c4zero_matches_c4zero_free_pole_but_not_current_reference"
        )
    else:
        classification = "documented_lpf_c4zero_profile_not_close_enough"

    return {
        "classification": classification,
        "best_normalization": best["norm"],
        "best_fitted_cutoff_Hz": float(best["fitted_cutoff_Hz"]),
        "best_cutoff_boundary_state": best["cutoff_boundary_state"],
        "best_c2": float(best["profiled"]["readout"]["c2"]),
        "best_continuum_rms_dB": xxf1._rms(best["profiled"]),
        "best_passes_current_free_pole_screen": pass_current,
        "best_passes_free_pole_c4zero_screen": pass_zero,
        "comparison_screen": screen,
    }


def run(config, config_path: Path, experiment_path_override=None):
    base_config_path = resolve_config_path(
        config["base_magnicon_config"], config_path
    )
    base_config = json.loads(base_config_path.read_text(encoding="utf-8"))
    problem = xxf1._problem(
        base_config,
        base_config_path,
        experiment_path_override=experiment_path_override,
    )
    snapshot = problem["snapshot"]

    bounds = dict(problem["readout_bounds"])
    bounds["c2"] = (0.0, float(config["c2_upper_bound"]))

    current_reference, free_zero_reference = _fit_references(
        problem, snapshot, bounds
    )

    nominal_cutoff = float(
        base_config["magnicon_filter"]["nominal_cutoff_Hz"]
    )
    tolerance = float(
        base_config["magnicon_filter"]["cutoff_tolerance_fraction"]
    )
    screen = config["comparison_screen"]

    branches = {}
    for index, norm in enumerate(config["normalization_candidates"]):
        branches[norm] = _fit_magnicon_branch(
            problem=problem,
            norm=str(norm),
            nominal_cutoff=nominal_cutoff,
            tolerance=tolerance,
            bounds=bounds,
            current_reference=current_reference,
            free_zero_reference=free_zero_reference,
            screen=screen,
            seed_offset=7200 + index * 300,
        )

    interpretation = _classify(branches, screen)

    return {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "measurement_time_filter_state_known": False,
        "tested_hypothesis": (
            "XXF-1 Connector Box 10 kHz second-order Bessel LPF was ON; "
            "c4 is unnecessary once that documented filter is explicit."
        ),
        "manual_reference": base_config["manual_reference"],
        "profile_definition": {
            "c4_fixed": 0.0,
            "residual_readout_parameter": "c2",
            "c2_upper_bound": float(config["c2_upper_bound"]),
            "nominal_cutoff_Hz": nominal_cutoff,
            "manual_cutoff_tolerance_fraction": tolerance,
            "normalization_candidates": list(config["normalization_candidates"]),
            "bessel_Q_fixed_by_normalization": True,
        },
        "repeat_case": {
            "label": problem["repeat_case"]["label"],
            "comparison_source": problem["comparison_source"],
            "experiment_path": str(problem["experiment_path"]),
            "accepted_records": int(problem["accepted_records"]),
        },
        "references": {
            "current_free_pole_plus_c2_legacy_c4": _reference_summary(
                current_reference
            ),
            "free_pole_plus_c2_c4_zero": _reference_summary(
                free_zero_reference
            ),
            "c4_zero_cost_relative_to_current": _comparison_to(
                free_zero_reference, current_reference
            ),
        },
        "branches": {
            norm: {
                "manual_cutoff_range_Hz": row["manual_cutoff_range_Hz"],
                "nominal_coordinates": row["nominal_coordinates"],
                "fitted_cutoff_Hz": row["fitted_cutoff_Hz"],
                "fitted_cutoff_fraction_from_nominal": row[
                    "fitted_cutoff_fraction_from_nominal"
                ],
                "cutoff_boundary_state": row["cutoff_boundary_state"],
                "nominal_fit": xxf1._clean(row["nominal"]),
                "profiled_fit": xxf1._clean(row["profiled"]),
                "comparison_nominal_to_current_free_pole": row[
                    "comparison_nominal_to_current_free_pole"
                ],
                "comparison_profiled_to_current_free_pole": row[
                    "comparison_profiled_to_current_free_pole"
                ],
                "comparison_profiled_to_free_pole_c4zero": row[
                    "comparison_profiled_to_free_pole_c4zero"
                ],
                "profile_improvement_over_nominal_rms_dB": row[
                    "profile_improvement_over_nominal_rms_dB"
                ],
                "passes_current_free_pole_screen": row[
                    "passes_current_free_pole_screen"
                ],
                "passes_free_pole_c4zero_screen": row[
                    "passes_free_pole_c4zero_screen"
                ],
            }
            for norm, row in branches.items()
        },
        "interpretation": {
            **interpretation,
            "guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_magnicon_config": str(base_config_path),
            "base_continuum_config": str(problem["continuum_config_path"]),
            "base_residual_config": str(problem["residual_config_path"]),
            "manifest": str(problem["manifest_path"]),
            "residual_baseline_snapshot": str(problem["snapshot_path"]),
        },
        "_plot": {
            "frequency_Hz": problem["frequency"].tolist(),
            "target": problem["target"].tolist(),
            "current_reference": current_reference["_model_full"].tolist(),
            "free_zero_reference": free_zero_reference["_model_full"].tolist(),
            "profiled_models": {
                norm: row["profiled"]["_model_full"].tolist()
                for norm, row in branches.items()
            },
        },
    }


def make_plot(result, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    frequency = np.asarray(p["frequency_Hz"], dtype=float)
    target = np.asarray(p["target"], dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.0), sharex=True)
    top, middle, bottom = axes

    top.loglog(frequency, target, label="line-repaired continuum")
    top.loglog(
        frequency,
        p["current_reference"],
        label="free pole + legacy c4",
    )
    top.loglog(
        frequency,
        p["free_zero_reference"],
        label="free pole, c4=0",
    )
    for norm, model in p["profiled_models"].items():
        top.loglog(
            frequency,
            np.asarray(model, dtype=float),
            label=f"Magnicon {norm}, c4=0",
        )
    top.set_ylabel("Normalized ASD")
    top.legend(frameon=False, fontsize=8)
    top.grid(True, which="both", alpha=0.2)

    for label, model in [
        ("free pole + legacy c4", p["current_reference"]),
        ("free pole, c4=0", p["free_zero_reference"]),
    ]:
        residual_db = 20.0 * np.log10(
            np.asarray(model, dtype=float) / target
        )
        middle.semilogx(frequency, residual_db, label=label)
    for norm, model in p["profiled_models"].items():
        residual_db = 20.0 * np.log10(
            np.asarray(model, dtype=float) / target
        )
        middle.semilogx(
            frequency,
            residual_db,
            label=f"Magnicon {norm}, c4=0",
        )
    middle.axhline(0.0, linewidth=1.0)
    middle.set_ylabel("Model / continuum [dB]")
    middle.grid(True, which="both", alpha=0.2)

    names = list(result["branches"])
    cutoffs = [
        result["branches"][name]["fitted_cutoff_Hz"] / 1000.0
        for name in names
    ]
    c2_values = [
        result["branches"][name]["profiled_fit"]["readout"]["c2"]
        for name in names
    ]
    x = np.arange(len(names), dtype=float)
    bottom.plot(x, cutoffs, marker="o", label="cutoff [kHz]")
    bottom.plot(x, c2_values, marker="o", label="c2")
    bottom.set_xticks(x, names)
    bottom.set_ylabel("Profile coordinates")
    bottom.set_xlabel("Bessel normalization")
    bottom.grid(True, alpha=0.2)
    bottom.legend(frameon=False)

    fig.suptitle(result["interpretation"]["classification"], fontsize=10)
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
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--experiment-path", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = run(
        config,
        args.config,
        experiment_path_override=args.experiment_path,
    )
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
                "interpretation": result["interpretation"],
                "branches": {
                    norm: {
                        "fitted_cutoff_Hz": row["fitted_cutoff_Hz"],
                        "c2": row["profiled_fit"]["readout"]["c2"],
                        "continuum_rms_dB": row["profiled_fit"][
                            "continuum_metrics_full"
                        ]["residual_metrics"]["rms_residual_dB"],
                    }
                    for norm, row in result["branches"].items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

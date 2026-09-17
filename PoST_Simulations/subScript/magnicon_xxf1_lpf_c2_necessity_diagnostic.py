"""Test whether residual c2 is needed once the documented XXF-1 LPF is explicit.

The preceding diagnostics found that the documented Magnicon XXF-1 Connector Box
10 kHz second-order Bessel LPF, with c4 fixed to zero, can reproduce the smooth
continuum almost as well as the previous free pole/Q effective readout model.
This diagnostic performs the next nested ablation: compare c2=0 against c2 free
at the same Bessel normalization and with the same cutoff freedom.

For each tested Bessel normalization it fits four branches on the identical
line-repaired 1--200 kHz target:

  1. nominal 10 kHz cutoff, c2=0,
  2. nominal 10 kHz cutoff, c2 free,
  3. cutoff profiled only within the manual +/-2.5% tolerance, c2=0,
  4. the same cutoff profile, c2 free.

All branches re-optimize the same detector nuisance parameters and profile the
same post-filter white floor.  c4 is fixed to zero throughout.  A material c2
increment is therefore assessed by a strictly nested comparison with identical
cutoff freedom.

The measurement-time filter state remains unknown.  Results are conditional on
the LPF-ON hypothesis and must not be read as proof that the Connector Box filter
was enabled during the acquisition.
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
DEFAULT_CONFIG = CONFIG_DIR / "magnicon_xxf1_lpf_c2_necessity_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_c2_necessity_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_c2_necessity_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import magnicon_xxf1_lpf_c4zero_cutoff_profile_diagnostic as c4prof  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def _nested_gain(c2_free: dict, c2_zero: dict) -> dict:
    """Return the improvement from releasing c2 at identical cutoff freedom."""
    return {
        "shape_score_ratio_c2_free_to_c2_zero": float(
            c2_free["shape_score"] / float(c2_zero["shape_score"])
        ),
        "continuum_rms_improvement_dB": float(
            xxf1._rms(c2_zero) - xxf1._rms(c2_free)
        ),
        "continuum_1_40k_rms_improvement_dB": float(
            xxf1._rms(c2_zero, "continuum_metrics_1_40k")
            - xxf1._rms(c2_free, "continuum_metrics_1_40k")
        ),
        "continuum_40_200k_rms_improvement_dB": float(
            xxf1._rms(c2_zero, "continuum_metrics_40_200k")
            - xxf1._rms(c2_free, "continuum_metrics_40_200k")
        ),
    }


def _c2_material(gain: dict, screen: dict) -> bool:
    return bool(
        gain["shape_score_ratio_c2_free_to_c2_zero"]
        <= float(screen["max_c2_free_shape_score_ratio_to_c2_zero"])
        and gain["continuum_rms_improvement_dB"]
        >= float(screen["min_c2_rms_improvement_dB"])
    )


def _passes_current(row: dict, current_reference: dict, screen: dict, prefix: str) -> bool:
    comp = c4prof._comparison_to(row, current_reference)
    return bool(
        comp["continuum_rms_delta_dB"]
        <= float(screen[f"max_{prefix}_rms_degradation_to_current_free_pole_dB"])
        and comp["shape_score_ratio"]
        <= float(screen[f"max_{prefix}_shape_score_ratio_to_current_free_pole"])
    )


def _fit_branch(
    *,
    problem: dict,
    norm: str,
    nominal_cutoff: float,
    tolerance: float,
    bounds: dict,
    current_reference: dict,
    free_zero_reference: dict,
    screen: dict,
    seed_offset: int,
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

    nominal_zero = xxf1._fit_variant(
        problem=problem,
        name=f"magnicon_{norm}_nominal_c2_zero",
        readout_parameters=(),
        reference_readout=reference,
        local_readout=reference,
        readout_bounds=bounds,
        seed_offset=seed_offset,
        warm_solutions=[
            (
                "free_pole_c4zero_detector",
                free_zero_reference["_detector_full"],
                reference,
            )
        ],
    )

    nominal_free = xxf1._fit_variant(
        problem=problem,
        name=f"magnicon_{norm}_nominal_c2_free",
        readout_parameters=("c2",),
        reference_readout=reference,
        local_readout=reference,
        readout_bounds=bounds,
        seed_offset=seed_offset + 25,
        warm_solutions=[
            (
                "nominal_c2_zero",
                nominal_zero["_detector_full"],
                nominal_zero["_readout_full"],
            ),
            (
                "free_pole_c4zero_detector",
                free_zero_reference["_detector_full"],
                reference,
            ),
        ],
    )

    tolerance_bounds = dict(bounds)
    tolerance_bounds["pole_Hz"] = (
        float(low_coordinates["pole_Hz"]),
        float(high_coordinates["pole_Hz"]),
    )

    profiled_zero = xxf1._fit_variant(
        problem=problem,
        name=f"magnicon_{norm}_tolerance_c2_zero",
        readout_parameters=("pole_Hz",),
        reference_readout=reference,
        local_readout=reference,
        readout_bounds=tolerance_bounds,
        seed_offset=seed_offset + 50,
        warm_solutions=[
            (
                "nominal_c2_zero",
                nominal_zero["_detector_full"],
                nominal_zero["_readout_full"],
            )
        ],
    )

    profiled_free = xxf1._fit_variant(
        problem=problem,
        name=f"magnicon_{norm}_tolerance_c2_free",
        readout_parameters=("pole_Hz", "c2"),
        reference_readout=reference,
        local_readout=reference,
        readout_bounds=tolerance_bounds,
        seed_offset=seed_offset + 75,
        warm_solutions=[
            (
                "profiled_c2_zero",
                profiled_zero["_detector_full"],
                profiled_zero["_readout_full"],
            ),
            (
                "nominal_c2_free",
                nominal_free["_detector_full"],
                nominal_free["_readout_full"],
            ),
        ],
    )

    cutoff_zero = xxf1.equivalent_cutoff_from_pole(
        profiled_zero["readout"]["pole_Hz"], nominal_coordinates
    )
    cutoff_free = xxf1.equivalent_cutoff_from_pole(
        profiled_free["readout"]["pole_Hz"], nominal_coordinates
    )
    boundary_zero = c4prof._cutoff_boundary_state(
        cutoff_zero,
        low_cutoff,
        high_cutoff,
        screen["cutoff_bound_relative_tolerance"],
    )
    boundary_free = c4prof._cutoff_boundary_state(
        cutoff_free,
        low_cutoff,
        high_cutoff,
        screen["cutoff_bound_relative_tolerance"],
    )

    nominal_gain = _nested_gain(nominal_free, nominal_zero)
    profiled_gain = _nested_gain(profiled_free, profiled_zero)

    return {
        "norm": norm,
        "manual_cutoff_range_Hz": [float(low_cutoff), float(high_cutoff)],
        "nominal_coordinates": nominal_coordinates,
        "nominal_c2_zero": nominal_zero,
        "nominal_c2_free": nominal_free,
        "profiled_c2_zero": profiled_zero,
        "profiled_c2_free": profiled_free,
        "profiled_c2_zero_cutoff_Hz": float(cutoff_zero),
        "profiled_c2_free_cutoff_Hz": float(cutoff_free),
        "profiled_c2_zero_cutoff_boundary_state": boundary_zero,
        "profiled_c2_free_cutoff_boundary_state": boundary_free,
        "cutoff_shift_when_c2_released_Hz": float(cutoff_free - cutoff_zero),
        "nominal_c2_gain": nominal_gain,
        "profiled_c2_gain": profiled_gain,
        "nominal_c2_material_improvement": _c2_material(nominal_gain, screen),
        "profiled_c2_material_improvement": _c2_material(profiled_gain, screen),
        "profiled_c2_zero_to_current": c4prof._comparison_to(
            profiled_zero, current_reference
        ),
        "profiled_c2_free_to_current": c4prof._comparison_to(
            profiled_free, current_reference
        ),
        "profiled_c2_zero_to_free_pole_c4zero": c4prof._comparison_to(
            profiled_zero, free_zero_reference
        ),
        "profiled_c2_free_to_free_pole_c4zero": c4prof._comparison_to(
            profiled_free, free_zero_reference
        ),
        "profiled_c2_zero_passes_current_screen": _passes_current(
            profiled_zero, current_reference, screen, "no_c2"
        ),
        "profiled_c2_free_passes_current_screen": _passes_current(
            profiled_free, current_reference, screen, "c2_free"
        ),
    }


def _classify_primary(branch: dict, screen: dict) -> dict:
    material = bool(branch["profiled_c2_material_improvement"])
    zero_pass = bool(branch["profiled_c2_zero_passes_current_screen"])
    free_pass = bool(branch["profiled_c2_free_passes_current_screen"])
    zero_boundary = branch["profiled_c2_zero_cutoff_boundary_state"]
    free_boundary = branch["profiled_c2_free_cutoff_boundary_state"]

    if zero_pass and not material:
        classification = "c2_not_materially_required_under_documented_lpf"
    elif zero_pass and material:
        classification = "c2_materially_improves_fit_but_is_not_required_for_current_screen"
    elif (not zero_pass) and free_pass and material:
        classification = "c2_required_to_meet_current_fit_screen_under_documented_lpf"
    elif (not zero_pass) and free_pass:
        classification = "c2_helps_meet_current_screen_but_increment_is_below_materiality_threshold"
    else:
        classification = "documented_lpf_branch_not_close_enough_even_with_c2"

    c2_hit = branch["profiled_c2_free"].get("readout_boundary_hits", {}).get("c2", {})
    return {
        "classification": classification,
        "primary_normalization": branch["norm"],
        "profiled_c2_material_improvement": material,
        "profiled_c2_zero_passes_current_screen": zero_pass,
        "profiled_c2_free_passes_current_screen": free_pass,
        "profiled_c2_zero_cutoff_Hz": float(
            branch["profiled_c2_zero_cutoff_Hz"]
        ),
        "profiled_c2_free_cutoff_Hz": float(
            branch["profiled_c2_free_cutoff_Hz"]
        ),
        "profiled_c2_zero_cutoff_boundary_state": zero_boundary,
        "profiled_c2_free_cutoff_boundary_state": free_boundary,
        "profiled_c2": float(branch["profiled_c2_free"]["readout"]["c2"]),
        "profiled_c2_boundary_hit": bool(
            c2_hit.get("at_lower", False) or c2_hit.get("at_upper", False)
        ),
        "profiled_c2_gain": branch["profiled_c2_gain"],
        "materiality_screen": screen,
    }


def _reference_summary(row: dict) -> dict:
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


def run(config: dict, config_path: Path, experiment_path_override=None):
    profile_config_path = resolve_config_path(
        config["base_profile_config"], config_path
    )
    profile_config = json.loads(
        profile_config_path.read_text(encoding="utf-8")
    )
    magnicon_config_path = c4prof.resolve_config_path(
        profile_config["base_magnicon_config"], profile_config_path
    )
    magnicon_config = json.loads(
        magnicon_config_path.read_text(encoding="utf-8")
    )

    problem = xxf1._problem(
        magnicon_config,
        magnicon_config_path,
        experiment_path_override=experiment_path_override,
    )
    snapshot = problem["snapshot"]

    bounds = dict(problem["readout_bounds"])
    bounds["c2"] = (0.0, float(config["c2_upper_bound"]))

    current_reference, free_zero_reference = c4prof._fit_references(
        problem, snapshot, bounds
    )

    nominal_cutoff = float(
        magnicon_config["magnicon_filter"]["nominal_cutoff_Hz"]
    )
    tolerance = float(
        magnicon_config["magnicon_filter"]["cutoff_tolerance_fraction"]
    )
    screen = config["materiality_screen"]

    branches = {}
    for index, norm in enumerate(config["normalization_candidates"]):
        branches[str(norm)] = _fit_branch(
            problem=problem,
            norm=str(norm),
            nominal_cutoff=nominal_cutoff,
            tolerance=tolerance,
            bounds=bounds,
            current_reference=current_reference,
            free_zero_reference=free_zero_reference,
            screen=screen,
            seed_offset=8100 + index * 400,
        )

    primary = str(config["primary_normalization"])
    if primary not in branches:
        raise ValueError(
            f"primary normalization {primary!r} was not evaluated"
        )
    interpretation = _classify_primary(branches[primary], screen)

    best_zero_norm = min(
        branches,
        key=lambda name: float(branches[name]["profiled_c2_zero"]["shape_score"]),
    )
    best_free_norm = min(
        branches,
        key=lambda name: float(branches[name]["profiled_c2_free"]["shape_score"]),
    )

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "measurement_time_filter_state_known": False,
        "tested_hypothesis": (
            "XXF-1 Connector Box 10 kHz second-order Bessel LPF was ON; "
            "test whether residual phenomenological c2 is materially required "
            "after c4 is fixed to zero."
        ),
        "manual_reference": magnicon_config["manual_reference"],
        "profile_definition": {
            "c4_fixed": 0.0,
            "nominal_cutoff_Hz": nominal_cutoff,
            "manual_cutoff_tolerance_fraction": tolerance,
            "primary_normalization": primary,
            "normalization_candidates": list(config["normalization_candidates"]),
            "c2_upper_bound": float(config["c2_upper_bound"]),
            "nested_rule": (
                "Compare c2=0 and c2-free fits with identical Bessel "
                "normalization, cutoff freedom, detector nuisance freedom, "
                "and white-floor profiling."
            ),
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
            "c4_zero_cost_relative_to_current": c4prof._comparison_to(
                free_zero_reference, current_reference
            ),
        },
        "branches": {
            norm: {
                "manual_cutoff_range_Hz": row["manual_cutoff_range_Hz"],
                "nominal_coordinates": row["nominal_coordinates"],
                "profiled_c2_zero_cutoff_Hz": row[
                    "profiled_c2_zero_cutoff_Hz"
                ],
                "profiled_c2_free_cutoff_Hz": row[
                    "profiled_c2_free_cutoff_Hz"
                ],
                "profiled_c2_zero_cutoff_boundary_state": row[
                    "profiled_c2_zero_cutoff_boundary_state"
                ],
                "profiled_c2_free_cutoff_boundary_state": row[
                    "profiled_c2_free_cutoff_boundary_state"
                ],
                "cutoff_shift_when_c2_released_Hz": row[
                    "cutoff_shift_when_c2_released_Hz"
                ],
                "nominal_c2_zero_fit": xxf1._clean(row["nominal_c2_zero"]),
                "nominal_c2_free_fit": xxf1._clean(row["nominal_c2_free"]),
                "profiled_c2_zero_fit": xxf1._clean(row["profiled_c2_zero"]),
                "profiled_c2_free_fit": xxf1._clean(row["profiled_c2_free"]),
                "nominal_c2_gain": row["nominal_c2_gain"],
                "profiled_c2_gain": row["profiled_c2_gain"],
                "nominal_c2_material_improvement": row[
                    "nominal_c2_material_improvement"
                ],
                "profiled_c2_material_improvement": row[
                    "profiled_c2_material_improvement"
                ],
                "profiled_c2_zero_to_current": row[
                    "profiled_c2_zero_to_current"
                ],
                "profiled_c2_free_to_current": row[
                    "profiled_c2_free_to_current"
                ],
                "profiled_c2_zero_to_free_pole_c4zero": row[
                    "profiled_c2_zero_to_free_pole_c4zero"
                ],
                "profiled_c2_free_to_free_pole_c4zero": row[
                    "profiled_c2_free_to_free_pole_c4zero"
                ],
                "profiled_c2_zero_passes_current_screen": row[
                    "profiled_c2_zero_passes_current_screen"
                ],
                "profiled_c2_free_passes_current_screen": row[
                    "profiled_c2_free_passes_current_screen"
                ],
            }
            for norm, row in branches.items()
        },
        "interpretation": {
            **interpretation,
            "best_profiled_c2_zero_normalization": best_zero_norm,
            "best_profiled_c2_free_normalization": best_free_norm,
            "guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_profile_config": str(profile_config_path),
            "base_magnicon_config": str(magnicon_config_path),
            "base_continuum_config": str(problem["continuum_config_path"]),
            "base_residual_config": str(problem["residual_config_path"]),
            "manifest": str(problem["manifest_path"]),
            "residual_baseline_snapshot": str(problem["snapshot_path"]),
        },
        "_plot": {
            "frequency_Hz": problem["frequency"].tolist(),
            "target": problem["target"].tolist(),
            "current_reference": current_reference["_model_full"].tolist(),
            "primary": {
                "profiled_c2_zero": branches[primary]["profiled_c2_zero"][
                    "_model_full"
                ].tolist(),
                "profiled_c2_free": branches[primary]["profiled_c2_free"][
                    "_model_full"
                ].tolist(),
                "nominal_c2_zero": branches[primary]["nominal_c2_zero"][
                    "_model_full"
                ].tolist(),
                "nominal_c2_free": branches[primary]["nominal_c2_free"][
                    "_model_full"
                ].tolist(),
            },
        },
    }
    return result


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    frequency = np.asarray(p["frequency_Hz"], dtype=float)
    target = np.asarray(p["target"], dtype=float)
    models = p["primary"]

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.0), sharex=True)
    top, middle, bottom = axes

    top.loglog(frequency, target, label="line-repaired continuum")
    top.loglog(
        frequency,
        np.asarray(p["current_reference"], dtype=float),
        label="current free pole reference",
    )
    top.loglog(
        frequency,
        np.asarray(models["profiled_c2_zero"], dtype=float),
        label="Magnicon profiled, c2=0",
    )
    top.loglog(
        frequency,
        np.asarray(models["profiled_c2_free"], dtype=float),
        label="Magnicon profiled, c2 free",
    )
    top.set_ylabel("Normalized ASD")
    top.legend(frameon=False, fontsize=8)
    top.grid(True, which="both", alpha=0.2)

    for label, model in [
        ("current free pole", p["current_reference"]),
        ("profiled c2=0", models["profiled_c2_zero"]),
        ("profiled c2 free", models["profiled_c2_free"]),
    ]:
        residual_db = 20.0 * np.log10(
            np.asarray(model, dtype=float) / target
        )
        middle.semilogx(frequency, residual_db, label=label)
    middle.axhline(0.0, linewidth=1.0)
    middle.set_ylabel("Model / continuum [dB]")
    middle.grid(True, which="both", alpha=0.2)
    middle.legend(frameon=False, fontsize=8)

    branch = result["branches"][
        result["interpretation"]["primary_normalization"]
    ]
    labels = [
        "nominal c2=0",
        "nominal c2 free",
        "profiled c2=0",
        "profiled c2 free",
    ]
    rows = [
        branch["nominal_c2_zero_fit"],
        branch["nominal_c2_free_fit"],
        branch["profiled_c2_zero_fit"],
        branch["profiled_c2_free_fit"],
    ]
    rms_values = [
        row["continuum_metrics_full"]["residual_metrics"]["rms_residual_dB"]
        for row in rows
    ]
    x = np.arange(len(labels), dtype=float)
    bottom.plot(x, rms_values, marker="o")
    bottom.set_xticks(x, labels, rotation=15)
    bottom.set_ylabel("Continuum RMS [dB]")
    bottom.grid(True, alpha=0.2)

    fig.suptitle(result["interpretation"]["classification"], fontsize=10)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def cleaned_result(result: dict) -> dict:
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

    primary = result["interpretation"]["primary_normalization"]
    branch = result["branches"][primary]
    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "interpretation": result["interpretation"],
                "primary_branch": {
                    "normalization": primary,
                    "profiled_c2_zero_cutoff_Hz": branch[
                        "profiled_c2_zero_cutoff_Hz"
                    ],
                    "profiled_c2_free_cutoff_Hz": branch[
                        "profiled_c2_free_cutoff_Hz"
                    ],
                    "profiled_c2": branch["profiled_c2_free_fit"]["readout"][
                        "c2"
                    ],
                    "profiled_c2_gain": branch["profiled_c2_gain"],
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

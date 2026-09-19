"""Test residual c2 after anchoring both L and C_tes to physical scales.

The preceding physical-L diagnostic showed that fixing L near its independently
motivated circuit value makes the fitted c2 transfer strongly repeatable across
2024-12-05 and 2024-12-06.  Those same fits pushed C_tes to the configured lower
bound.  This diagnostic therefore removes that remaining escape direction.

C_tes is anchored to Opt_noise.C_TES_MATERIAL_J_PER_K, which is computed from
the repository's TES bilayer geometry, density, and specific heat.  This is
independent of the noise fit, but it is a material estimate rather than direct
calorimetry.  A multiplier grid tests sensitivity to that estimate.

For every C_tes multiplier and each adjacent-day dataset:
  * L is fixed to the external 1e-10 H circuit scale;
  * C_tes is fixed to the selected material-estimate multiple;
  * Magnicon is fixed to nominal 10 kHz phase-normalized 2nd-order Bessel;
  * c4=0;
  * alpha, beta, T_bath and the white floor are refit;
  * c2=0 and c2-free models are compared with identical freedoms.

The fitted c2-free readout transfers are then compared directly across days.
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
    CONFIG_DIR
    / "magnicon_xxf1_physical_L_Ctes_anchor_c2_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_physical_L_Ctes_anchor_c2_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_physical_L_Ctes_anchor_c2_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_necessity_diagnostic as c2diag  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_transfer_shape_stability_diagnostic as transferdiag  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def _multiplier_key(value: float) -> str:
    value = float(value)
    if value.is_integer():
        return f"{int(value)}x"
    return f"{value:g}x"


def _problem_with_fixed_detector_values(
    problem: dict,
    fixed_values: dict,
) -> dict:
    result = dict(problem)
    baseline = dict(problem["baseline_detector"])
    detector_names = list(problem["detector_names"])

    for name, value in fixed_values.items():
        value = float(value)
        if name not in problem["detector_bounds"]:
            raise ValueError(f"unknown detector parameter {name}")
        lower, upper = (
            float(v) for v in problem["detector_bounds"][name]
        )
        if value < lower or value > upper:
            raise ValueError(
                f"fixed {name}={value:g} outside bounds "
                f"[{lower:g}, {upper:g}]"
            )
        baseline[name] = value
        detector_names = [item for item in detector_names if item != name]

    result["baseline_detector"] = baseline
    result["detector_names"] = tuple(detector_names)
    return result


def _nested_summary(c2_free: dict, c2_zero: dict, screen: dict) -> dict:
    gain = c2diag._nested_gain(c2_free, c2_zero)
    return {
        "c2_material_improvement": bool(
            c2diag._c2_material(gain, screen)
        ),
        "gain": gain,
        "c2_free_rms_dB": xxf1._rms(c2_free),
        "c2_zero_rms_dB": xxf1._rms(c2_zero),
        "c2_free_1_40k_rms_dB": xxf1._rms(
            c2_free, "continuum_metrics_1_40k"
        ),
        "c2_zero_1_40k_rms_dB": xxf1._rms(
            c2_zero, "continuum_metrics_1_40k"
        ),
        "c2_free_40_200k_rms_dB": xxf1._rms(
            c2_free, "continuum_metrics_40_200k"
        ),
        "c2_zero_40_200k_rms_dB": xxf1._rms(
            c2_zero, "continuum_metrics_40_200k"
        ),
    }


def _fit_case(
    *,
    problem: dict,
    fixed_L_H: float,
    fixed_C_tes_J_per_K: float,
    readout_reference: dict,
    c2_upper_bound: float,
    materiality_screen: dict,
    seed_offset: int,
):
    fixed_problem = _problem_with_fixed_detector_values(
        problem,
        {
            "L": float(fixed_L_H),
            "C_tes": float(fixed_C_tes_J_per_K),
        },
    )
    bounds = dict(fixed_problem["readout_bounds"])
    bounds["c2"] = (0.0, float(c2_upper_bound))

    zero = xxf1._fit_variant(
        problem=fixed_problem,
        name="magnicon_nominal_phase_physical_L_Ctes_c2_zero",
        readout_parameters=(),
        reference_readout=readout_reference,
        local_readout=readout_reference,
        readout_bounds=bounds,
        seed_offset=seed_offset,
    )
    free = xxf1._fit_variant(
        problem=fixed_problem,
        name="magnicon_nominal_phase_physical_L_Ctes_c2_free",
        readout_parameters=("c2",),
        reference_readout=readout_reference,
        local_readout=readout_reference,
        readout_bounds=bounds,
        seed_offset=seed_offset + 50,
        warm_solutions=[
            (
                "c2_zero_solution",
                zero["_detector_full"],
                zero["_readout_full"],
            )
        ],
    )

    return {
        "fixed_L_H": float(fixed_L_H),
        "fixed_C_tes_J_per_K": float(fixed_C_tes_J_per_K),
        "detector_parameters_varied": list(
            fixed_problem["detector_names"]
        ),
        "c2_zero": zero,
        "c2_free": free,
        "nested": _nested_summary(
            free,
            zero,
            materiality_screen,
        ),
    }


def _equivalent_zero_hz(c2: float, scale_hz: float):
    c2 = float(c2)
    if c2 <= 0.0:
        return None
    return float(scale_hz / np.sqrt(c2))


def _case_summary(row: dict, scale_hz: float) -> dict:
    free = row["c2_free"]
    zero = row["c2_zero"]
    c2 = float(free["readout"]["c2"])
    return {
        "fixed_L_H": float(row["fixed_L_H"]),
        "fixed_C_tes_J_per_K": float(row["fixed_C_tes_J_per_K"]),
        "detector_parameters_varied": list(
            row["detector_parameters_varied"]
        ),
        "c2_free": {
            "c2": c2,
            "equivalent_first_order_zero_Hz": _equivalent_zero_hz(
                c2, scale_hz
            ),
            "equivalent_zero_semantics": (
                "magnitude coordinate only: "
                "sqrt(1+c2*(f/f_ref)^2)=sqrt(1+(f/f_zero)^2)"
            ),
            "white_asd_A_rtHz": float(
                free["profiled_white_asd_A_rtHz"]
            ),
            "detector_candidate": {
                key: float(value)
                for key, value in free["detector_candidate"].items()
            },
            "detector_boundary_hits": free.get(
                "detector_boundary_hits", {}
            ),
            "readout_boundary_hits": free.get(
                "readout_boundary_hits", {}
            ),
            "continuum_rms_dB": xxf1._rms(free),
            "continuum_1_40k_rms_dB": xxf1._rms(
                free, "continuum_metrics_1_40k"
            ),
            "continuum_40_200k_rms_dB": xxf1._rms(
                free, "continuum_metrics_40_200k"
            ),
            "shape_score": float(free["shape_score"]),
        },
        "c2_zero": {
            "white_asd_A_rtHz": float(
                zero["profiled_white_asd_A_rtHz"]
            ),
            "detector_candidate": {
                key: float(value)
                for key, value in zero["detector_candidate"].items()
            },
            "detector_boundary_hits": zero.get(
                "detector_boundary_hits", {}
            ),
            "continuum_rms_dB": xxf1._rms(zero),
            "continuum_1_40k_rms_dB": xxf1._rms(
                zero, "continuum_metrics_1_40k"
            ),
            "continuum_40_200k_rms_dB": xxf1._rms(
                zero, "continuum_metrics_40_200k"
            ),
            "shape_score": float(zero["shape_score"]),
        },
        "nested_c2_test": row["nested"],
    }


def _direct_transfer_difference(
    *,
    frequency,
    scale_hz: float,
    reference_readout: dict,
    repeat_readout: dict,
    normalization_reference_hz: float,
    bands,
):
    reference = transferdiag.readout_components(
        frequency,
        reference_readout,
        scale_hz,
        normalization_reference_hz,
    )
    repeat = transferdiag.readout_components(
        frequency,
        repeat_readout,
        scale_hz,
        normalization_reference_hz,
    )
    difference = transferdiag.difference_db(
        repeat["total"],
        reference["total"],
    )
    return {
        "semantics": (
            "repeat/reference c2-free readout transfer, each normalized at "
            f"{float(normalization_reference_hz):g} Hz"
        ),
        "full_1_200k": transferdiag.global_metrics(difference),
        "bands": transferdiag.band_metrics(
            frequency,
            difference,
            bands,
        ),
        "_difference_dB": difference,
    }


def _transfer_passes(metrics: dict, screen: dict) -> bool:
    full = metrics["full_1_200k"]
    return bool(
        float(full["rms_difference_dB"])
        <= float(screen["max_full_band_rms_difference_dB"])
        and float(full["max_abs_difference_dB"])
        <= float(screen["max_abs_difference_dB"])
    )


def _classify(
    reference_material: bool,
    repeat_material: bool,
    transfer_passes: bool,
) -> str:
    if not reference_material and not repeat_material:
        return "physical_L_Ctes_anchors_remove_material_c2_need_on_both_days"
    if reference_material and repeat_material and transfer_passes:
        return "physical_L_Ctes_anchors_leave_repeatable_material_c2_residual"
    if reference_material and repeat_material:
        return "physical_L_Ctes_anchors_leave_nonrepeatable_material_c2_residual"
    return "physical_L_Ctes_anchors_leave_day_dependent_c2_requirement"


def run(config: dict, config_path: Path):
    stack = compensation._load_stack(config, config_path)
    base_cross = stack["base_cross"]
    c2_config = stack["c2_config"]
    magnicon_config = stack["magnicon_config"]

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
    if not np.allclose(
        reference_problem["frequency"],
        repeat_problem["frequency"],
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError("reference/repeat frequency grids differ")
    if not np.isclose(
        reference_problem["scale_hz"],
        repeat_problem["scale_hz"],
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError("reference/repeat readout scales differ")

    filter_cfg = config["magnicon_filter"]
    canonical = xxf1.second_order_bessel_canonical(
        float(filter_cfg["cutoff_Hz"]),
        str(filter_cfg["normalization"]),
    )
    readout_reference = {
        "pole_Hz": float(canonical["pole_Hz"]),
        "pole_Q": float(canonical["pole_Q"]),
        "c2": 0.0,
        "c4": float(filter_cfg["c4_fixed"]),
    }
    if readout_reference["c4"] != 0.0:
        raise ValueError("this diagnostic requires c4=0")

    fixed_L = float(config["external_L_H"])
    material_ctes = float(opt.C_TES_MATERIAL_J_PER_K)
    ctes_cfg = config["C_tes_anchor"]
    multipliers = [float(v) for v in ctes_cfg["multipliers"]]
    primary_multiplier = float(ctes_cfg["primary_multiplier"])
    if not any(
        np.isclose(v, primary_multiplier, rtol=0.0, atol=1e-12)
        for v in multipliers
    ):
        raise ValueError("primary C_tes multiplier is absent from grid")

    rows = {}
    plot_rows = {}
    signatures = []

    for index, multiplier in enumerate(multipliers):
        key = _multiplier_key(multiplier)
        fixed_ctes = material_ctes * multiplier

        reference = _fit_case(
            problem=reference_problem,
            fixed_L_H=fixed_L,
            fixed_C_tes_J_per_K=fixed_ctes,
            readout_reference=readout_reference,
            c2_upper_bound=float(c2_config["c2_upper_bound"]),
            materiality_screen=c2_config["materiality_screen"],
            seed_offset=19100 + 600 * index,
        )
        repeat = _fit_case(
            problem=repeat_problem,
            fixed_L_H=fixed_L,
            fixed_C_tes_J_per_K=fixed_ctes,
            readout_reference=readout_reference,
            c2_upper_bound=float(c2_config["c2_upper_bound"]),
            materiality_screen=c2_config["materiality_screen"],
            seed_offset=19400 + 600 * index,
        )

        transfer = _direct_transfer_difference(
            frequency=reference_problem["frequency"],
            scale_hz=reference_problem["scale_hz"],
            reference_readout=reference["c2_free"]["_readout_full"],
            repeat_readout=repeat["c2_free"]["_readout_full"],
            normalization_reference_hz=float(
                config["transfer_normalization_reference_Hz"]
            ),
            bands=config["bands_Hz"],
        )
        transfer_pass = _transfer_passes(
            transfer,
            config["transfer_repeatability_screen"],
        )
        ref_material = bool(
            reference["nested"]["c2_material_improvement"]
        )
        rep_material = bool(
            repeat["nested"]["c2_material_improvement"]
        )
        classification = _classify(
            ref_material,
            rep_material,
            transfer_pass,
        )
        signatures.append(
            (ref_material, rep_material, transfer_pass)
        )

        ref_c2 = float(reference["c2_free"]["readout"]["c2"])
        rep_c2 = float(repeat["c2_free"]["readout"]["c2"])
        c2_sym = float(
            2.0 * abs(ref_c2 - rep_c2) / (abs(ref_c2) + abs(rep_c2))
        )

        rows[key] = {
            "C_tes_multiplier": multiplier,
            "fixed_L_H": fixed_L,
            "fixed_C_tes_J_per_K": fixed_ctes,
            "reference_day": _case_summary(
                reference, reference_problem["scale_hz"]
            ),
            "repeat_day": _case_summary(
                repeat, repeat_problem["scale_hz"]
            ),
            "cross_day_c2": {
                "reference_c2": ref_c2,
                "repeat_c2": rep_c2,
                "symmetric_fractional_difference": c2_sym,
            },
            "cross_day_c2_free_transfer": {
                key2: value
                for key2, value in transfer.items()
                if key2 != "_difference_dB"
            },
            "transfer_repeatability_passes": transfer_pass,
            "classification": classification,
        }
        plot_rows[key] = {
            "multiplier": multiplier,
            "reference_c2": ref_c2,
            "repeat_c2": rep_c2,
            "reference_c2_free_rms": xxf1._rms(
                reference["c2_free"]
            ),
            "repeat_c2_free_rms": xxf1._rms(
                repeat["c2_free"]
            ),
            "reference_c2_zero_rms": xxf1._rms(
                reference["c2_zero"]
            ),
            "repeat_c2_zero_rms": xxf1._rms(
                repeat["c2_zero"]
            ),
            "transfer_difference_dB": transfer[
                "_difference_dB"
            ].tolist(),
        }

    primary_key = _multiplier_key(primary_multiplier)
    primary = rows[primary_key]
    sensitivity_consistent = bool(
        all(signature == signatures[0] for signature in signatures)
    )

    return {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "tested_question": (
            "After fixing both L and C_tes to independently motivated "
            "physical scales, does a material and cross-day-repeatable c2 "
            "residual remain?"
        ),
        "physical_anchors": {
            "L_H": fixed_L,
            "C_tes_material_J_per_K": material_ctes,
            "C_tes_material_source": ctes_cfg["source"],
            "C_tes_material_source_semantics": ctes_cfg[
                "source_semantics"
            ],
            "C_tes_sensitivity_multipliers": multipliers,
            "primary_C_tes_multiplier": primary_multiplier,
            "sensitivity_grid_is_not_uncertainty_interval": True,
        },
        "magnicon_fixed_model": {
            "normalization": str(filter_cfg["normalization"]),
            "cutoff_Hz": float(filter_cfg["cutoff_Hz"]),
            "pole_Hz": float(canonical["pole_Hz"]),
            "pole_Q": float(canonical["pole_Q"]),
            "c4": 0.0,
        },
        "fit_semantics": {
            "L": "fixed externally; not optimized",
            "C_tes": "fixed to each material-estimate sensitivity row",
            "remaining_detector_nuisance": (
                "alpha, beta, T_bath optimized independently per day"
            ),
            "white_floor": "profiled independently per day",
            "nested_test": (
                "c2=0 versus c2 free with identical physical anchors, "
                "Magnicon transfer, nuisance freedoms, and white profiling"
            ),
        },
        "sensitivity_rows": rows,
        "primary_result": {
            "key": primary_key,
            **primary,
        },
        "interpretation": {
            "classification": primary["classification"],
            "primary_reference_c2_material": bool(
                primary["reference_day"]["nested_c2_test"][
                    "c2_material_improvement"
                ]
            ),
            "primary_repeat_c2_material": bool(
                primary["repeat_day"]["nested_c2_test"][
                    "c2_material_improvement"
                ]
            ),
            "primary_transfer_repeatability_passes": bool(
                primary["transfer_repeatability_passes"]
            ),
            "primary_c2_symmetric_fractional_difference": float(
                primary["cross_day_c2"][
                    "symmetric_fractional_difference"
                ]
            ),
            "sensitivity_signature_consistent_across_grid": (
                sensitivity_consistent
            ),
            "guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_cross_day_config": str(stack["base_cross_path"]),
            "base_c2_config": str(stack["c2_config_path"]),
            "base_magnicon_config": str(
                stack["magnicon_config_path"]
            ),
        },
        "_plot": {
            "frequency_Hz": reference_problem["frequency"].tolist(),
            "rows": plot_rows,
        },
    }


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    rows = list(p["rows"].items())
    multipliers = np.asarray(
        [row["multiplier"] for _, row in rows],
        dtype=float,
    )

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.0))
    top, middle, bottom = axes

    top.plot(
        multipliers,
        [row["reference_c2"] for _, row in rows],
        marker="o",
        label="12/06 c2",
    )
    top.plot(
        multipliers,
        [row["repeat_c2"] for _, row in rows],
        marker="o",
        label="12/05 c2",
    )
    top.set_xscale("log")
    top.set_ylabel("Fitted c2")
    top.grid(True, which="both", alpha=0.2)
    top.legend(frameon=False)

    middle.plot(
        multipliers,
        [row["reference_c2_free_rms"] for _, row in rows],
        marker="o",
        label="12/06 c2 free",
    )
    middle.plot(
        multipliers,
        [row["reference_c2_zero_rms"] for _, row in rows],
        marker="o",
        label="12/06 c2=0",
    )
    middle.plot(
        multipliers,
        [row["repeat_c2_free_rms"] for _, row in rows],
        marker="o",
        label="12/05 c2 free",
    )
    middle.plot(
        multipliers,
        [row["repeat_c2_zero_rms"] for _, row in rows],
        marker="o",
        label="12/05 c2=0",
    )
    middle.set_xscale("log")
    middle.set_ylabel("Continuum RMS [dB]")
    middle.grid(True, which="both", alpha=0.2)
    middle.legend(frameon=False, fontsize=8)

    frequency = np.asarray(p["frequency_Hz"], dtype=float)
    for key, row in rows:
        bottom.semilogx(
            frequency,
            np.asarray(row["transfer_difference_dB"], dtype=float),
            label=f"Ctes {key}",
        )
    bottom.axhline(0.0, linewidth=1.0)
    bottom.set_xlabel("Frequency [Hz]")
    bottom.set_ylabel("12/05 / 12/06 readout [dB]")
    bottom.grid(True, which="both", alpha=0.2)
    bottom.legend(frameon=False)

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

    primary = result["primary_result"]
    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "classification": result["interpretation"][
                    "classification"
                ],
                "fixed_L_H": primary["fixed_L_H"],
                "fixed_C_tes_J_per_K": primary[
                    "fixed_C_tes_J_per_K"
                ],
                "reference_c2": primary["reference_day"][
                    "c2_free"
                ]["c2"],
                "repeat_c2": primary["repeat_day"][
                    "c2_free"
                ]["c2"],
                "c2_symmetric_fractional_difference": result[
                    "interpretation"
                ]["primary_c2_symmetric_fractional_difference"],
                "reference_c2_material": result["interpretation"][
                    "primary_reference_c2_material"
                ],
                "repeat_c2_material": result["interpretation"][
                    "primary_repeat_c2_material"
                ],
                "transfer_repeatability_passes": result[
                    "interpretation"
                ]["primary_transfer_repeatability_passes"],
                "sensitivity_signature_consistent": result[
                    "interpretation"
                ]["sensitivity_signature_consistent_across_grid"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

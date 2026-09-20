"""Plot c2=0 versus c2-free anchored fits over the full 1-200 kHz domain.

This diagnostic is deliberately visual.  It reruns the same nested fits used by
magnicon_xxf1_physical_L_Ctes_anchor_c2_diagnostic.py for selected fixed C_tes
anchors, then overlays for both 2024-12-06 and 2024-12-05:

  * the line-repaired continuum target used by the smooth fit;
  * the raw experimental ASD (including narrow lines), for provenance;
  * the best nested c2=0 model;
  * the best c2-free model.

The lower panels show model/continuum residuals in dB.  One PNG is produced per
selected C_tes anchor, and the plotted numerical arrays plus fit summaries are
written to JSON.
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
WORK_DIR = ROOT / ".noise_optimization_work_rsh_sweep"

DEFAULT_CONFIG = CONFIG_DIR / "magnicon_xxf1_c2_on_off_fullband_plot_config.json"
DEFAULT_OUTPUT = WORK_DIR / "magnicon_xxf1_c2_on_off_fullband_plot.json"
DEFAULT_FIGURE_DIR = WORK_DIR

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import magnicon_xxf1_physical_L_Ctes_anchor_c2_diagnostic as anchored  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def residual_db(model, target):
    model = np.asarray(model, dtype=float)
    target = np.asarray(target, dtype=float)
    if model.shape != target.shape:
        raise ValueError("model and target shapes differ")
    if np.any(model <= 0.0) or np.any(target <= 0.0):
        raise ValueError("model and target must be positive")
    return 20.0 * np.log10(model / target)


def multiplier_key(value: float) -> str:
    return anchored._multiplier_key(float(value))


def multiplier_slug(value: float) -> str:
    key = multiplier_key(value)
    return key.replace(".", "p").replace("-", "m")


def figure_path(figure_dir: Path, multiplier: float) -> Path:
    return figure_dir / (
        "magnicon_xxf1_c2_on_off_fullband_Ctes_"
        + multiplier_slug(multiplier)
        + ".png"
    )


def display_mask(frequency, display_cfg):
    frequency = np.asarray(frequency, dtype=float)
    low = float(display_cfg["min"])
    high = float(display_cfg["max"])
    if not 0.0 < low < high:
        raise ValueError("invalid display frequency range")
    mask = (frequency >= low) & (frequency <= high)
    if np.count_nonzero(mask) < 2:
        raise ValueError("display range contains fewer than two points")
    return mask


def fit_summary(row: dict) -> dict:
    free = row["c2_free"]
    zero = row["c2_zero"]
    return {
        "c2_free": {
            "c2": float(free["readout"]["c2"]),
            "white_asd_A_rtHz": float(free["profiled_white_asd_A_rtHz"]),
            "continuum_rms_dB": xxf1._rms(free),
            "continuum_1_40k_rms_dB": xxf1._rms(
                free, "continuum_metrics_1_40k"
            ),
            "continuum_40_200k_rms_dB": xxf1._rms(
                free, "continuum_metrics_40_200k"
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
        },
        "c2_zero": {
            "white_asd_A_rtHz": float(zero["profiled_white_asd_A_rtHz"]),
            "continuum_rms_dB": xxf1._rms(zero),
            "continuum_1_40k_rms_dB": xxf1._rms(
                zero, "continuum_metrics_1_40k"
            ),
            "continuum_40_200k_rms_dB": xxf1._rms(
                zero, "continuum_metrics_40_200k"
            ),
            "detector_candidate": {
                key: float(value)
                for key, value in zero["detector_candidate"].items()
            },
            "detector_boundary_hits": zero.get(
                "detector_boundary_hits", {}
            ),
        },
        "nested_c2_test": row["nested"],
    }


def _day_plot_payload(problem: dict, row: dict) -> dict:
    target = np.asarray(problem["target"], dtype=float)
    raw_target = np.asarray(problem["raw_target"], dtype=float)
    zero_model = np.asarray(row["c2_zero"]["_model_full"], dtype=float)
    free_model = np.asarray(row["c2_free"]["_model_full"], dtype=float)
    return {
        "continuum_target": target,
        "raw_experimental_ASD": raw_target,
        "c2_zero_model": zero_model,
        "c2_free_model": free_model,
        "c2_zero_residual_dB": residual_db(zero_model, target),
        "c2_free_residual_dB": residual_db(free_model, target),
    }


def _fit_selected_anchor(
    *,
    multiplier: float,
    index: int,
    fixed_L_H: float,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
    c2_upper_bound: float,
    materiality_screen: dict,
):
    fixed_C_tes = material_C_tes * float(multiplier)
    reference = anchored._fit_case(
        problem=reference_problem,
        fixed_L_H=fixed_L_H,
        fixed_C_tes_J_per_K=fixed_C_tes,
        readout_reference=readout_reference,
        c2_upper_bound=c2_upper_bound,
        materiality_screen=materiality_screen,
        seed_offset=23100 + 600 * index,
    )
    repeat = anchored._fit_case(
        problem=repeat_problem,
        fixed_L_H=fixed_L_H,
        fixed_C_tes_J_per_K=fixed_C_tes,
        readout_reference=readout_reference,
        c2_upper_bound=c2_upper_bound,
        materiality_screen=materiality_screen,
        seed_offset=23400 + 600 * index,
    )
    return fixed_C_tes, reference, repeat


def run(config: dict, config_path: Path):
    base_config_path = resolve_config_path(
        config["base_anchor_config"],
        config_path,
    )
    base_config = json.loads(
        base_config_path.read_text(encoding="utf-8")
    )

    stack = compensation._load_stack(base_config, base_config_path)
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

    filter_cfg = base_config["magnicon_filter"]
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

    fixed_L = float(base_config["external_L_H"])
    material_C_tes = float(opt.C_TES_MATERIAL_J_PER_K)
    selected = [
        float(value)
        for value in config["C_tes_multipliers_to_plot"]
    ]
    if not selected:
        raise ValueError("no C_tes multipliers selected")

    available = {
        float(value)
        for value in base_config["C_tes_anchor"]["multipliers"]
    }
    missing = [
        value for value in selected if value not in available
    ]
    if missing:
        raise ValueError(
            "selected C_tes multipliers are absent from base anchor grid: "
            + ", ".join(f"{value:g}" for value in missing)
        )

    frequency = np.asarray(reference_problem["frequency"], dtype=float)
    mask = display_mask(
        frequency,
        config["frequency_display_Hz"],
    )

    cases = {}
    internal_plot = {}

    for index, multiplier in enumerate(selected):
        fixed_C_tes, reference, repeat = _fit_selected_anchor(
            multiplier=multiplier,
            index=index,
            fixed_L_H=fixed_L,
            material_C_tes=material_C_tes,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
            readout_reference=readout_reference,
            c2_upper_bound=float(c2_config["c2_upper_bound"]),
            materiality_screen=c2_config["materiality_screen"],
        )
        key = multiplier_key(multiplier)
        reference_payload = _day_plot_payload(
            reference_problem, reference
        )
        repeat_payload = _day_plot_payload(
            repeat_problem, repeat
        )

        cases[key] = {
            "C_tes_multiplier": multiplier,
            "fixed_C_tes_J_per_K": fixed_C_tes,
            "fixed_L_H": fixed_L,
            "reference_day_2024_12_06": {
                "fit": fit_summary(reference),
                "arrays": {
                    name: values.tolist()
                    for name, values in reference_payload.items()
                },
            },
            "repeat_day_2024_12_05": {
                "fit": fit_summary(repeat),
                "arrays": {
                    name: values.tolist()
                    for name, values in repeat_payload.items()
                },
            },
        }
        internal_plot[key] = {
            "multiplier": multiplier,
            "fixed_C_tes_J_per_K": fixed_C_tes,
            "reference": reference_payload,
            "repeat": repeat_payload,
        }

    return {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "tested_question": (
            "What do the nested c2=0 and c2-free best fits look like over "
            "the complete 1-200 kHz continuum-fit domain?"
        ),
        "plot_semantics": {
            "frequency_domain_Hz": {
                "min": float(frequency[mask][0]),
                "max": float(frequency[mask][-1]),
            },
            "continuum_target": (
                "line-repaired experimental ASD used by the smooth "
                "continuum objective"
            ),
            "raw_experimental_ASD": (
                "experimental ASD before narrow-line repair; visual "
                "provenance only, not the smooth fit target"
            ),
            "residual_dB": "20*log10(model / line-repaired continuum target)",
            "Magnicon": {
                "normalization": str(filter_cfg["normalization"]),
                "cutoff_Hz": float(filter_cfg["cutoff_Hz"]),
                "pole_Hz": float(canonical["pole_Hz"]),
                "pole_Q": float(canonical["pole_Q"]),
                "c4": 0.0,
            },
        },
        "physical_anchor_context": {
            "fixed_L_H": fixed_L,
            "C_tes_material_J_per_K": material_C_tes,
            "selected_C_tes_multipliers": selected,
        },
        "frequency_Hz": frequency.tolist(),
        "cases": cases,
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "base_anchor_config": str(base_config_path),
        },
        "_display_mask": mask,
        "_plot": internal_plot,
    }


def make_case_plot(
    result: dict,
    key: str,
    output: Path,
    include_raw=True,
    show=False,
):
    import matplotlib.pyplot as plt

    frequency = np.asarray(result["frequency_Hz"], dtype=float)
    mask = np.asarray(result["_display_mask"], dtype=bool)
    p = result["_plot"][key]
    summary = result["cases"][key]

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(13.5, 8.5),
        sharex="col",
        gridspec_kw={"height_ratios": [2.0, 1.0]},
    )

    day_specs = [
        (
            0,
            "2024-12-06 reference",
            p["reference"],
            summary["reference_day_2024_12_06"]["fit"],
        ),
        (
            1,
            "2024-12-05 repeat",
            p["repeat"],
            summary["repeat_day_2024_12_05"]["fit"],
        ),
    ]

    for column, day_label, payload, fit in day_specs:
        top = axes[0, column]
        bottom = axes[1, column]

        if include_raw:
            top.loglog(
                frequency[mask],
                payload["raw_experimental_ASD"][mask],
                linewidth=0.9,
                alpha=0.55,
                label="Raw experimental ASD",
            )
        top.loglog(
            frequency[mask],
            payload["continuum_target"][mask],
            linewidth=1.5,
            label="Line-repaired continuum target",
        )
        top.loglog(
            frequency[mask],
            payload["c2_zero_model"][mask],
            linewidth=1.4,
            label="Best fit: c2=0",
        )
        top.loglog(
            frequency[mask],
            payload["c2_free_model"][mask],
            linewidth=1.4,
            label="Best fit: c2 free",
        )
        top.set_title(
            day_label
            + "\n"
            + (
                f"c2={fit['c2_free']['c2']:.4g}, "
                f"RMS: {fit['c2_zero']['continuum_rms_dB']:.3f} "
                f"→ {fit['c2_free']['continuum_rms_dB']:.3f} dB"
            ),
            fontsize=10,
        )
        top.set_ylabel("Normalized ASD")
        top.grid(True, which="both", alpha=0.2)
        top.legend(frameon=False, fontsize=8)

        bottom.semilogx(
            frequency[mask],
            payload["c2_zero_residual_dB"][mask],
            linewidth=1.2,
            label="c2=0",
        )
        bottom.semilogx(
            frequency[mask],
            payload["c2_free_residual_dB"][mask],
            linewidth=1.2,
            label="c2 free",
        )
        bottom.axhline(0.0, linewidth=1.0)
        bottom.set_xlabel("Frequency [Hz]")
        bottom.set_ylabel("Model / continuum [dB]")
        bottom.grid(True, which="both", alpha=0.2)
        bottom.legend(frameon=False, fontsize=8)

    multiplier = float(summary["C_tes_multiplier"])
    fixed_ctes = float(summary["fixed_C_tes_J_per_K"])
    fixed_l = float(summary["fixed_L_H"])
    fig.suptitle(
        (
            f"c2 on/off full-band fit — Ctes={multiplier:g}× material "
            f"({fixed_ctes:.4e} J/K), L={fixed_l:.3e} H"
        ),
        fontsize=11,
    )
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
        if key not in {"_display_mask", "_plot"}
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=DEFAULT_FIGURE_DIR,
    )
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = run(config, args.config)

    figure_files = {}
    include_raw = bool(
        config.get("include_raw_experimental_ASD", True)
    )
    for key, row in result["cases"].items():
        path = figure_path(
            args.figure_dir,
            float(row["C_tes_multiplier"]),
        )
        make_case_plot(
            result,
            key,
            path,
            include_raw=include_raw,
            show=args.show,
        )
        figure_files[key] = str(path)

    clean = cleaned_result(result)
    clean["figure_files"] = figure_files
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(clean, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    console_cases = {}
    for key, row in clean["cases"].items():
        ref = row["reference_day_2024_12_06"]["fit"]
        rep = row["repeat_day_2024_12_05"]["fit"]
        console_cases[key] = {
            "figure": figure_files[key],
            "reference_c2": ref["c2_free"]["c2"],
            "reference_rms_c2_zero_dB": ref["c2_zero"][
                "continuum_rms_dB"
            ],
            "reference_rms_c2_free_dB": ref["c2_free"][
                "continuum_rms_dB"
            ],
            "repeat_c2": rep["c2_free"]["c2"],
            "repeat_rms_c2_zero_dB": rep["c2_zero"][
                "continuum_rms_dB"
            ],
            "repeat_rms_c2_free_dB": rep["c2_free"][
                "continuum_rms_dB"
            ],
        }

    print(
        json.dumps(
            {
                "output": str(args.output),
                "frequency_domain_Hz": clean["plot_semantics"][
                    "frequency_domain_Hz"
                ],
                "cases": console_cases,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

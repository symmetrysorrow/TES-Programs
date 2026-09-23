"""Magnicon XXF-1 connector-box 10 kHz LPF ON/OFF diagnostic.

This is a deliberately narrow acquisition-state test.  It keeps c2=0 and the
same relaxed thermal-backbone profile in both branches, then asks whether the
measured continuum is better described by:

  * OFF: no Magnicon connector-box 10 kHz attenuation (unity), or
  * ON: the documented second-order phase-normalized Bessel response.

The OFF branch uses diagnostic order=0 in the shared filter-order fitter.  Its
cutoff coordinate is retained but mathematically inert so ON and OFF have
exactly the same nuisance-parameter count.  The external SRS SIM965 model,
alias handling, detector freedoms, and per-day white floors are unchanged.
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

DEFAULT_CONFIG = (
    CONFIG_DIR / "magnicon_xxf1_filter_state_thermal_backbone_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR / "magnicon_xxf1_filter_state_thermal_backbone_diagnostic.json"
)

REFERENCE_HZ = 1000.0
STATE_OFF = "off"
STATE_ON = "on"
ORDER_OFF = 0
ORDER_ON = 2

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_filter_order_thermal_backbone_diagnostic as orderdiag  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_thermal_backbone_pure_zero_diagnostic as thermal  # noqa: E402


def compare_states(off: dict, on: dict, config: dict) -> dict:
    screen = config["materiality_screen"]
    score_ratio = float(off["joint_shape_score"] / on["joint_shape_score"])
    rms_improvement = float(
        on["joint_continuum_rms_dB"] - off["joint_continuum_rms_dB"]
    )
    off_materially_better = bool(
        score_ratio <= float(screen["max_off_over_on_score_ratio"])
        and rms_improvement >= float(screen["min_off_rms_improvement_dB"])
    )
    off_good = bool(
        off["joint_continuum_rms_dB"]
        <= float(screen["target_good_fit_rms_dB"])
    )
    on_good = bool(
        on["joint_continuum_rms_dB"]
        <= float(screen["target_good_fit_rms_dB"])
    )

    if off_materially_better and off_good:
        classification = "filter_off_reaches_good_fit_and_beats_documented_on"
    elif off_materially_better:
        classification = "filter_off_materially_improves_but_does_not_fully_close_residual"
    else:
        classification = "filter_off_does_not_materially_resolve_residual"

    direction = (
        "off_better"
        if rms_improvement > 0.0
        else "on_better"
        if rms_improvement < 0.0
        else "tie"
    )
    return {
        "classification": classification,
        "direction_by_joint_rms": direction,
        "filter_off_materially_better": off_materially_better,
        "filter_off_reaches_target_good_fit": off_good,
        "filter_on_reaches_target_good_fit": on_good,
        "off_over_on_shape_score_ratio": score_ratio,
        "filter_off_rms_improvement_dB": rms_improvement,
        "off_minus_on_joint_rms_dB": float(
            off["joint_continuum_rms_dB"] - on["joint_continuum_rms_dB"]
        ),
        "screen": screen,
    }


def transfer_comparison(off: dict, on: dict, *, norm: str) -> dict:
    frequencies = np.asarray(
        [2e3, 5e3, 10e3, 20e3, 40e3, 100e3, 200e3],
        dtype=float,
    )
    off_mag = orderdiag.normalized_magnicon_magnitude(
        frequencies,
        order=ORDER_OFF,
        cutoff_hz=off["solution"]["magnicon_cutoff_Hz"],
        norm=norm,
    )
    on_mag = orderdiag.normalized_magnicon_magnitude(
        frequencies,
        order=ORDER_ON,
        cutoff_hz=on["solution"]["magnicon_cutoff_Hz"],
        norm=norm,
    )
    ratio_db = 20.0 * np.log10(off_mag / on_mag)
    return {
        "reference_Hz": REFERENCE_HZ,
        "frequency_Hz": frequencies.tolist(),
        "off_over_on_dB": ratio_db.tolist(),
        "off_branch_is_unity": True,
    }


def run(config: dict, config_path: Path) -> dict:
    states = tuple(str(value) for value in config["states"])
    if states != (STATE_OFF, STATE_ON):
        raise ValueError('this diagnostic expects states ["off", "on"]')
    if int(config["magnicon_filter"]["documented_order"]) != ORDER_ON:
        raise ValueError("documented Magnicon order must remain two")

    thermal_path = thermal.resolve_config_path(
        config["base_thermal_backbone_config"],
        config_path,
    )
    thermal_config = json.loads(thermal_path.read_text(encoding="utf-8"))
    separated_path = thermal.resolve_config_path(
        thermal_config["base_separated_substrate_config"],
        thermal_path,
    )
    separated_config = json.loads(separated_path.read_text(encoding="utf-8"))

    stack = compensation._load_stack(separated_config, separated_path)
    base_cross = stack["base_cross"]
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

    material_C_tes = float(opt.C_TES_MATERIAL_J_PER_K)
    common = {
        "config": config,
        "thermal_config": thermal_config,
        "separated_config": separated_config,
        "material_C_tes": material_C_tes,
        "reference_problem": reference_problem,
        "repeat_problem": repeat_problem,
    }
    off = orderdiag.fit_order(order=ORDER_OFF, **common)
    on = orderdiag.fit_order(order=ORDER_ON, **common)

    comparison = compare_states(off, on, config)
    norm = str(config["magnicon_filter"]["normalization"])
    transfer = transfer_comparison(off, on, norm=norm)

    return {
        "diagnostic_only": True,
        "production_default_topology_unchanged": True,
        "tested_question": (
            "With c2 fixed to zero and all relaxed thermal/detector nuisance "
            "freedoms held identical, does removing the Magnicon Connector Box "
            "10 kHz anti-alias filter materially improve the two-day continuum fit?"
        ),
        "manual_constraint": {
            "documented_order_when_on": ORDER_ON,
            "statement": str(config["magnicon_filter"]["manual_statement"]),
            "cutoff_profile_Hz": [
                float(config["magnicon_filter"]["cutoff_Hz"]["min"]),
                float(config["magnicon_filter"]["cutoff_Hz"]["max"]),
            ],
            "normalization_when_on": norm,
        },
        "fit_semantics": {
            "c2": 0.0,
            "same_number_of_free_parameters_per_state": bool(
                off["n_free_parameters"] == on["n_free_parameters"]
            ),
            "off_branch_cutoff_coordinate_inert": True,
            "shared_thermal_parameters": [
                "C_tes",
                "L",
                "C_substrate",
                "G_tes-substrate/G_substrate-bath",
                "G_tes-bath scale to inherited snapshot",
                "Pb absorber thickness",
            ],
            "day_specific": [
                "alpha",
                "beta",
                "T_bath",
                "post-filter white ASD",
            ],
            "external_analog_filter": {
                "model": "SRS SIM965",
                "front_panel_cutoff_Hz": float(
                    opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ
                ),
                "order": int(opt.TARGET_HARDWARE_BESSEL_ORDER),
                "scipy_norm": str(opt.TARGET_HARDWARE_BESSEL_NORM),
            },
        },
        "filter_off": {
            key: value for key, value in off.items() if not key.startswith("_")
        },
        "filter_on": {
            key: value for key, value in on.items() if not key.startswith("_")
        },
        "comparison": comparison,
        "filter_transfer_comparison": transfer,
        "interpretation": {
            "classification": comparison["classification"],
            "lab_priority_if_off_wins": (
                "Verify the historical XXF-1 Connector Box LPF switch/state and "
                "which output path fed the SIM965/DAQ before interpreting c2 as "
                "missing detector physics."
            ),
            "guardrail": str(config["guardrail"]),
        },
        "inputs": {
            "config": str(config_path),
            "base_thermal_backbone_config": str(thermal_path),
            "base_separated_substrate_config": str(separated_path),
            "base_cross_day_config": str(stack["base_cross_path"]),
            "base_magnicon_config": str(stack["magnicon_config_path"]),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = run(config, args.config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "classification": result["comparison"]["classification"],
                "filter_off": {
                    "joint_rms_dB": result["filter_off"][
                        "joint_continuum_rms_dB"
                    ],
                    "joint_shape_score": result["filter_off"][
                        "joint_shape_score"
                    ],
                },
                "filter_on": {
                    "joint_rms_dB": result["filter_on"][
                        "joint_continuum_rms_dB"
                    ],
                    "joint_shape_score": result["filter_on"][
                        "joint_shape_score"
                    ],
                    "cutoff_Hz": result["filter_on"]["solution"][
                        "magnicon_cutoff_Hz"
                    ],
                },
                "comparison": result["comparison"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

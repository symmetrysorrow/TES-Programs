"""Broad effective-C_tes cross-day joint c2 diagnostic.

The earlier shared-L/shared-C_tes test restricted C_tes to 0.25--1 times a
bare-bilayer material estimate.  Historical PoST inputs instead used
C_tes=7.9e-12 J/K, roughly an order of magnitude larger.  In the reduced
five-state model C_tes is the heat capacity of the entire TES thermal state,
so this diagnostic deliberately treats it as an *effective local-node* heat
capacity and profiles a broad absolute range that includes the historical
PoST value.

This module reuses the established nested cross-day fitter:
  * L and C_tes are shared across 2024-12-06 and 2024-12-05;
  * alpha, beta, T_bath, and post-filter white floor remain day-specific;
  * Magnicon is fixed to the nominal 10 kHz phase-normalized second-order
    Bessel response with c4=0;
  * the nested comparison is shared c2=0 versus one shared c2 free.

The historical 7.9e-12 J/K value is a simulation reference.  The repository
does not document a derivation that explicitly assigns it to silicon, so this
script must not be read as evidence for that microscopic decomposition.
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
WORK_DIR = ROOT / ".noise_optimization_work_rsh_sweep"

DEFAULT_CONFIG = (
    CONFIG_DIR
    / "magnicon_xxf1_effective_Ctes_joint_c2_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR
    / "magnicon_xxf1_effective_Ctes_joint_c2_diagnostic.json"
)
DEFAULT_FIGURE = (
    WORK_DIR
    / "magnicon_xxf1_effective_Ctes_joint_c2_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_shared_L_Ctes_joint_c2_diagnostic as shared  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def absolute_bounds_to_material_multipliers(
    lower_J_per_K: float,
    upper_J_per_K: float,
    material_J_per_K: float,
) -> tuple[float, float]:
    lower = float(lower_J_per_K)
    upper = float(upper_J_per_K)
    material = float(material_J_per_K)
    if not 0.0 < lower < upper:
        raise ValueError("effective C_tes bounds must satisfy 0 < min < max")
    if not material > 0.0:
        raise ValueError("material C_tes must be positive")
    return lower / material, upper / material


def log_position(value: float, lower: float, upper: float) -> float:
    value = float(value)
    lower = float(lower)
    upper = float(upper)
    if not 0.0 < lower < upper:
        raise ValueError("invalid positive bounds")
    if not lower <= value <= upper:
        raise ValueError("value lies outside bounds")
    return float(
        (np.log(value) - np.log(lower))
        / (np.log(upper) - np.log(lower))
    )


def build_profile_config(
    config: dict,
    config_path: Path,
) -> tuple[dict, Path, dict]:
    base_path = resolve_config_path(config["base_joint_config"], config_path)
    base = json.loads(base_path.read_text(encoding="utf-8"))
    derived = copy.deepcopy(base)

    effective = config["effective_C_tes_profile_J_per_K"]
    ctes_min = float(effective["min"])
    ctes_max = float(effective["max"])
    historical = float(effective["historical_PoST_reference"])
    material = float(opt.C_TES_MATERIAL_J_PER_K)
    mult_min, mult_max = absolute_bounds_to_material_multipliers(
        ctes_min,
        ctes_max,
        material,
    )
    if not ctes_min <= historical <= ctes_max:
        raise ValueError("historical PoST C_tes reference is outside profile")

    l_cfg = config["shared_L_profile_H"]
    l_min = float(l_cfg["min"])
    l_max = float(l_cfg["max"])
    if not 0.0 < l_min < l_max:
        raise ValueError("invalid L profile bounds")

    derived["purpose"] = config["purpose"]
    derived["shared_physical_profile"]["C_tes_material_multiplier"] = {
        "min": mult_min,
        "max": mult_max,
        "semantics": (
            "derived from the absolute effective-C_tes diagnostic range; "
            "the bare-bilayer material estimate is used only as the fitter's "
            "coordinate scale, not as a physical prior center"
        ),
    }
    derived["shared_physical_profile"]["L_H"] = {
        "min": l_min,
        "max": l_max,
        "semantics": str(l_cfg["semantics"]),
    }
    derived["guardrail"] = str(config["guardrail"])

    # Use an independent deterministic optimizer seed stream and slightly
    # increase global coverage because the C_tes interval is much broader.
    optimizer = derived["optimizer"]
    optimizer["seed_offset"] = int(optimizer["seed_offset"]) + 4000
    optimizer["DE_maxiter"] = max(int(optimizer["DE_maxiter"]), 160)
    optimizer["DE_popsize"] = max(int(optimizer["DE_popsize"]), 12)
    optimizer["least_squares_max_nfev"] = max(
        int(optimizer["least_squares_max_nfev"]),
        6000,
    )

    context = {
        "material_C_tes_J_per_K": material,
        "effective_C_tes_bounds_J_per_K": {
            "min": ctes_min,
            "max": ctes_max,
        },
        "effective_C_tes_bounds_over_bare_material": {
            "min": mult_min,
            "max": mult_max,
        },
        "historical_PoST_C_tes_J_per_K": historical,
        "historical_PoST_over_bare_material": historical / material,
        "historical_PoST_log_position_in_profile": log_position(
            historical,
            ctes_min,
            ctes_max,
        ),
        "shared_L_bounds_H": {
            "min": l_min,
            "max": l_max,
        },
    }
    return derived, base_path, context


def effective_classification(result: dict) -> str:
    nested = result["nested_test"]
    free_solution = result["c2_free"]["solution"]
    zero_solution = result["c2_zero"]["solution"]

    if not bool(nested["shared_c2_material_improvement"]):
        return "broad_effective_Ctes_profile_removes_material_c2_need"

    c2_hit = free_solution.get("shared_c2_boundary_hit")
    if c2_hit and bool(c2_hit.get("at_upper")):
        return (
            "shared_c2_material_but_c2_upper_bound_limited_with_"
            "broad_effective_Ctes_profile"
        )

    for solution in (zero_solution, free_solution):
        hit = solution["shared_physical_boundary_hits"]["C_tes"]
        if bool(hit["at_lower"]) or bool(hit["at_upper"]):
            return (
                "shared_c2_material_but_effective_Ctes_profile_"
                "boundary_active"
            )

    return "shared_c2_material_after_broad_effective_Ctes_profile"


def _reference_comparison(result: dict, context: dict) -> dict:
    historical = float(context["historical_PoST_C_tes_J_per_K"])
    zero = result["c2_zero"]["solution"]
    free = result["c2_free"]["solution"]
    ctes_min = float(context["effective_C_tes_bounds_J_per_K"]["min"])
    ctes_max = float(context["effective_C_tes_bounds_J_per_K"]["max"])

    def row(solution: dict) -> dict:
        ctes = float(solution["shared_C_tes_J_per_K"])
        return {
            "best_shared_C_tes_J_per_K": ctes,
            "best_over_historical_PoST": ctes / historical,
            "best_log_position_in_profile": log_position(
                ctes,
                ctes_min,
                ctes_max,
            ),
            "shared_L_H": float(solution["shared_L_H"]),
            "shared_c2": float(solution["shared_c2"]),
            "C_tes_boundary_hit": solution[
                "shared_physical_boundary_hits"
            ]["C_tes"],
        }

    return {
        "historical_PoST_C_tes_J_per_K": historical,
        "c2_zero": row(zero),
        "c2_free": row(free),
    }


def run(config: dict, config_path: Path) -> dict:
    derived, base_path, context = build_profile_config(
        config,
        config_path,
    )

    result = shared.run(derived, config_path)
    classification = effective_classification(result)

    result["tested_question"] = (
        "If C_tes is treated as the effective heat capacity of the reduced "
        "TES thermal node and profiled across an absolute range that includes "
        "the historical PoST value 7.9e-12 J/K, does a shared residual c2 "
        "remain materially necessary across 2024-12-05/06?"
    )
    result["effective_C_tes_context"] = context
    result["historical_PoST_reference_comparison"] = (
        _reference_comparison(result, context)
    )
    result["interpretation"]["base_joint_classification"] = result[
        "interpretation"
    ]["classification"]
    result["interpretation"]["classification"] = classification
    result["interpretation"]["C_tes_semantics"] = (
        "effective reduced-node heat capacity; may include rapidly "
        "thermalized substrate/membrane material, but the historical PoST "
        "value is not documented as an explicit silicon calculation"
    )
    result["interpretation"]["historical_PoST_C_tes_J_per_K"] = float(
        context["historical_PoST_C_tes_J_per_K"]
    )
    result["interpretation"]["guardrail"] = str(config["guardrail"])
    result["inputs"]["effective_C_tes_config"] = str(config_path)
    result["inputs"]["base_joint_config"] = str(base_path)
    return result


def make_plot(result: dict, output: Path, show=False) -> None:
    # Reuse the established full-band nested-fit visualization.  The title
    # will report the new effective-C_tes classification and the fitted C_tes
    # in bare-material units; absolute values are preserved in the JSON.
    shared.make_plot(result, output, show=show)


def cleaned_result(result: dict) -> dict:
    return shared.cleaned_result(result)


def main() -> None:
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

    zero = result["c2_zero"]["solution"]
    free = result["c2_free"]["solution"]
    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "classification": result["interpretation"][
                    "classification"
                ],
                "historical_PoST_C_tes_J_per_K": result[
                    "effective_C_tes_context"
                ]["historical_PoST_C_tes_J_per_K"],
                "c2_zero": {
                    "shared_C_tes_J_per_K": zero[
                        "shared_C_tes_J_per_K"
                    ],
                    "C_tes_over_historical_PoST": (
                        zero["shared_C_tes_J_per_K"]
                        / result["effective_C_tes_context"][
                            "historical_PoST_C_tes_J_per_K"
                        ]
                    ),
                    "shared_L_H": zero["shared_L_H"],
                    "joint_rms_dB": result["c2_zero"][
                        "joint_continuum_rms_dB"
                    ],
                },
                "c2_free": {
                    "shared_C_tes_J_per_K": free[
                        "shared_C_tes_J_per_K"
                    ],
                    "C_tes_over_historical_PoST": (
                        free["shared_C_tes_J_per_K"]
                        / result["effective_C_tes_context"][
                            "historical_PoST_C_tes_J_per_K"
                        ]
                    ),
                    "shared_L_H": free["shared_L_H"],
                    "shared_c2": free["shared_c2"],
                    "joint_rms_dB": result["c2_free"][
                        "joint_continuum_rms_dB"
                    ],
                },
                "joint_rms_improvement_dB": result["nested_test"][
                    "joint_rms_improvement_dB"
                ],
                "joint_shape_score_ratio": result["nested_test"][
                    "joint_shape_score_ratio_c2_free_to_c2_zero"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

"""Scan TES inductance through the Elmer-configured 12.3 nH scale.

For each fixed-L row this diagnostic uses exactly the relaxed thermal-backbone
profile from magnicon_xxf1_thermal_backbone_pure_zero_diagnostic and compares:

  * c2 = 0
  * c2 free

The L grid spans the historical noise-model scale (0.1 nH), the Elmer project
input (12.3 nH), and deliberately larger values up to 100 nH.  The Elmer value
is treated only as a configured circuit input, not as an independently inferred
inductance.

To reuse the existing optimizer without duplicating its physics, each nominal
fixed L is represented by an extremely narrow relative profile.  The width is
reported explicitly and is small enough to be numerically fixed for this
sensitivity test.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
WORK_DIR = ROOT / ".noise_optimization_work_rsh_sweep"

DEFAULT_CONFIG = (
    CONFIG_DIR
    / "magnicon_xxf1_high_L_thermal_backbone_c2_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR
    / "magnicon_xxf1_high_L_thermal_backbone_c2_diagnostic.json"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import magnicon_xxf1_thermal_backbone_pure_zero_diagnostic as thermal  # noqa: E402


def _resolve(value, anchor: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (anchor.parent / path).resolve()


def _parse_nh_expression(value: str) -> float:
    match = re.fullmatch(
        r"\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*\[nH\]\s*",
        str(value),
    )
    if match is None:
        raise ValueError(f"expected an nH expression, got {value!r}")
    result = float(match.group(1))
    if result <= 0.0:
        raise ValueError("Elmer L expression must be positive")
    return result


def _base_config_at_fixed_L(
    base_config: dict,
    fixed_L_H: float,
    relative_half_width: float,
) -> dict:
    fixed = float(fixed_L_H)
    rel = float(relative_half_width)
    if fixed <= 0.0:
        raise ValueError("fixed L must be positive")
    if not 0.0 < rel < 1.0:
        raise ValueError("relative L half-width must lie in (0, 1)")

    result = copy.deepcopy(base_config)
    bounds = result["shared_physical_profile"]["L_H"]
    bounds["min"] = fixed * (1.0 - rel)
    bounds["max"] = fixed * (1.0 + rel)
    bounds["nominal"] = fixed
    bounds["semantics"] = (
        "Numerically fixed L row for high-L sensitivity diagnostic; "
        f"target={fixed:.12g} H, relative half-width={rel:.3g}."
    )
    return result


def _strip_models(branch: dict) -> dict:
    return {
        key: value for key, value in branch.items() if not key.startswith("_")
    }


def _nested_metrics(zero: dict, c2: dict, screen: dict) -> dict:
    zero_score = float(zero["joint_shape_score"])
    c2_score = float(c2["joint_shape_score"])
    zero_rms = float(zero["joint_continuum_rms_dB"])
    c2_rms = float(c2["joint_continuum_rms_dB"])
    ratio = c2_score / zero_score
    improvement = zero_rms - c2_rms
    material = bool(
        ratio <= float(screen["max_c2_over_zero_shape_score_ratio"])
        and improvement >= float(screen["min_joint_rms_improvement_dB"])
    )
    hit = c2["solution"].get("readout_boundary_hit")
    censored = bool(hit and hit.get("at_upper"))

    c2_value = float(c2["solution"]["shared_c2"])
    scale_hz = float(screen["_scale_hz"])
    equivalent_zero = (
        float(scale_hz / np.sqrt(c2_value)) if c2_value > 0.0 else None
    )
    zero_good = bool(
        zero_rms <= float(screen["target_good_fit_rms_dB"])
    )

    if zero_good and not material:
        classification = "fixed_L_c2_zero_reaches_good_fit_without_material_c2"
    elif not material:
        classification = "fixed_L_removes_material_c2_need_but_fit_not_good"
    elif censored:
        classification = "fixed_L_still_requires_c2_and_c2_is_upper_censored"
    else:
        classification = "fixed_L_still_requires_material_c2"

    return {
        "classification": classification,
        "c2_material_improvement": material,
        "c2_upper_bound_censored": censored,
        "c2_over_zero_shape_score_ratio": float(ratio),
        "joint_rms_improvement_dB": float(improvement),
        "c2_zero_reaches_target_good_fit": zero_good,
        "fitted_c2": c2_value,
        "equivalent_first_order_zero_Hz": equivalent_zero,
    }


def _band_metrics(branch: dict) -> dict:
    solution = branch["solution"]
    return {
        "joint_continuum_rms_dB": float(branch["joint_continuum_rms_dB"]),
        "joint_shape_score": float(branch["joint_shape_score"]),
        "reference_day": {
            key: float(value)
            for key, value in solution["reference_day"]["metrics"].items()
        },
        "repeat_day": {
            key: float(value)
            for key, value in solution["repeat_day"]["metrics"].items()
        },
    }


def _row_summary(
    *,
    nominal_L_nH: float,
    zero: dict,
    c2: dict,
    screen: dict,
) -> dict:
    fitted_L_zero = float(zero["solution"]["shared_L_H"])
    fitted_L_c2 = float(c2["solution"]["shared_L_H"])
    return {
        "nominal_fixed_L_nH": float(nominal_L_nH),
        "nominal_fixed_L_H": float(nominal_L_nH) * 1.0e-9,
        "realized_L_H": {
            "c2_zero": fitted_L_zero,
            "c2_reference": fitted_L_c2,
        },
        "nested_c2_test": _nested_metrics(zero, c2, screen),
        "c2_zero_metrics": _band_metrics(zero),
        "c2_reference_metrics": _band_metrics(c2),
        "c2_zero": _strip_models(zero),
        "c2_reference": _strip_models(c2),
    }


def _best_stable_row(rows: dict) -> tuple[str | None, dict | None]:
    stable = [
        (key, row)
        for key, row in rows.items()
        if row.get("status") == "ok"
    ]
    if not stable:
        return None, None
    return min(
        stable,
        key=lambda item: item[1]["c2_zero_metrics"][
            "joint_continuum_rms_dB"
        ],
    )


def run(config: dict, config_path: Path) -> dict:
    thermal_path = _resolve(
        config["base_thermal_backbone_config"], config_path
    )
    thermal_config = json.loads(thermal_path.read_text(encoding="utf-8"))
    base_path = _resolve(
        thermal_config["base_separated_substrate_config"], thermal_path
    )
    base_config = json.loads(base_path.read_text(encoding="utf-8"))

    elmer_cfg = config["elmer_reference"]
    elmer_path = _resolve(elmer_cfg["project"], config_path)
    elmer_project = json.loads(elmer_path.read_text(encoding="utf-8"))
    elmer_expression = str(
        elmer_project["parameter_expressions"]["L_tes"]
    )
    expected_expression = str(
        elmer_cfg["expected_parameter_expression"]
    )
    if elmer_expression != expected_expression:
        raise ValueError(
            "Elmer L_tes expression changed: "
            f"expected {expected_expression!r}, got {elmer_expression!r}"
        )
    elmer_L_nH = _parse_nh_expression(elmer_expression)

    stack = compensation._load_stack(base_config, base_path)
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
    if not np.isclose(
        reference_problem["scale_hz"],
        repeat_problem["scale_hz"],
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError("reference/repeat readout scales differ")

    filter_cfg = base_config["magnicon_filter"]
    canonical = xxf1.second_order_bessel_canonical(
        float(filter_cfg["cutoff_Hz"]),
        str(filter_cfg["normalization"]),
    )
    readout_reference = {
        "pole_Hz": float(canonical["pole_Hz"]),
        "pole_Q": float(canonical["pole_Q"]),
        "c2": 0.0,
        "c4": 0.0,
    }
    material_C_tes = float(opt.C_TES_MATERIAL_J_PER_K)
    rel_width = float(config["numerical_fixed_L_relative_half_width"])

    screen = dict(config["materiality_screen"])
    screen["_scale_hz"] = float(reference_problem["scale_hz"])

    grid = [float(value) for value in config["L_grid_nH"]]
    if not any(np.isclose(value, elmer_L_nH) for value in grid):
        raise ValueError("L grid must explicitly contain the Elmer L_tes value")

    rows = {}
    previous_zero = None
    for index, L_nH in enumerate(grid):
        key = f"{L_nH:g}_nH"
        fixed_base = _base_config_at_fixed_L(
            base_config,
            L_nH * 1.0e-9,
            rel_width,
        )
        row_config = copy.deepcopy(thermal_config)
        row_config["optimizer"]["seed_offset"] = (
            int(thermal_config["optimizer"]["seed_offset"])
            + 5000
            + 3000 * index
        )

        common = {
            "config": row_config,
            "base_config": fixed_base,
            "material_C_tes": material_C_tes,
            "reference_problem": reference_problem,
            "repeat_problem": repeat_problem,
            "readout_reference": readout_reference,
        }

        try:
            zero = thermal.fit_branch(
                branch=thermal.BRANCH_ZERO,
                warm_from=previous_zero,
                **common,
            )
            c2 = thermal.fit_branch(
                branch=thermal.BRANCH_C2,
                warm_from=zero,
                **common,
            )
        except Exception as exc:
            rows[key] = {
                "status": "fit_failed",
                "nominal_fixed_L_nH": L_nH,
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
            previous_zero = None
            continue

        summary = _row_summary(
            nominal_L_nH=L_nH,
            zero=zero,
            c2=c2,
            screen=screen,
        )
        summary["status"] = "ok"
        rows[key] = summary
        previous_zero = zero

    best_key, best = _best_stable_row(rows)
    elmer_key = f"{elmer_L_nH:g}_nH"
    elmer_row = rows[elmer_key]

    baseline_key = f"{grid[0]:g}_nH"
    baseline = rows.get(baseline_key)
    best_improvement = None
    if (
        best is not None
        and baseline is not None
        and baseline.get("status") == "ok"
    ):
        best_improvement = float(
            baseline["c2_zero_metrics"]["joint_continuum_rms_dB"]
            - best["c2_zero_metrics"]["joint_continuum_rms_dB"]
        )

    if elmer_row.get("status") != "ok":
        elmer_classification = "elmer_L_row_not_stably_fitted"
    else:
        elmer_classification = elmer_row["nested_c2_test"]["classification"]

    any_no_c2 = any(
        row.get("status") == "ok"
        and not row["nested_c2_test"]["c2_material_improvement"]
        for row in rows.values()
    )
    any_good_no_c2 = any(
        row.get("status") == "ok"
        and row["nested_c2_test"]["c2_zero_reaches_target_good_fit"]
        and not row["nested_c2_test"]["c2_material_improvement"]
        for row in rows.values()
    )
    if any_good_no_c2:
        grid_classification = "some_fixed_L_reaches_good_fit_without_material_c2"
    elif any_no_c2:
        grid_classification = "some_fixed_L_removes_material_c2_need_but_not_residual"
    else:
        grid_classification = "high_L_grid_does_not_remove_material_c2_need"

    clean_screen = {
        key: value for key, value in screen.items() if not key.startswith("_")
    }
    return {
        "diagnostic_only": True,
        "production_default_topology_unchanged": True,
        "tested_question": (
            "Does allowing TES circuit inductance through the Elmer-configured "
            "12.3 nH scale and up to 100 nH remove the two-day normalized-noise "
            "need for phenomenological c2 under the same relaxed thermal backbone?"
        ),
        "L_scan": {
            "grid_nH": grid,
            "numerical_fixed_L_relative_half_width": rel_width,
            "implementation": (
                "Each row narrows the shared L profile to the nominal value "
                "+/- the stated relative half-width; all other relaxed thermal "
                "freedoms are unchanged."
            ),
        },
        "elmer_reference": {
            "project": str(elmer_path),
            "parameter": "parameter_expressions.L_tes",
            "expression": elmer_expression,
            "value_nH": elmer_L_nH,
            "is_independent_elmer_inductance_estimate": False,
            "semantics": str(elmer_cfg["semantics"]),
        },
        "fit_semantics": {
            "branches_per_L": ["c2_zero", "c2_reference"],
            "shared_thermal_parameters_except_L": [
                "C_tes",
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
            "magnicon": {
                "order": 2,
                "cutoff_Hz": float(filter_cfg["cutoff_Hz"]),
                "normalization": str(filter_cfg["normalization"]),
                "pole_Q": float(canonical["pole_Q"]),
            },
        },
        "rows": rows,
        "summary": {
            "grid_classification": grid_classification,
            "elmer_row_key": elmer_key,
            "elmer_row_classification": elmer_classification,
            "best_c2_zero_row_key": best_key,
            "best_c2_zero_joint_rms_dB": (
                float(best["c2_zero_metrics"]["joint_continuum_rms_dB"])
                if best is not None
                else None
            ),
            "best_c2_zero_improvement_vs_first_grid_row_dB": best_improvement,
            "materiality_screen": clean_screen,
        },
        "interpretation": {
            "lab_relevance": (
                "If the 12.3 nH row or a nearby physically plausible row removes "
                "the material c2 need, prioritize direct circuit-inductance and "
                "electrical-transfer measurements. If c2 remains material through "
                "the high-L grid, inductance alone is disfavored as the missing "
                "zero-like degree of freedom."
            ),
            "guardrail": str(config["guardrail"]),
        },
        "inputs": {
            "config": str(config_path),
            "base_thermal_backbone_config": str(thermal_path),
            "base_separated_substrate_config": str(base_path),
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
                "grid_classification": result["summary"][
                    "grid_classification"
                ],
                "elmer_row_classification": result["summary"][
                    "elmer_row_classification"
                ],
                "best_c2_zero_row_key": result["summary"][
                    "best_c2_zero_row_key"
                ],
                "best_c2_zero_joint_rms_dB": result["summary"][
                    "best_c2_zero_joint_rms_dB"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

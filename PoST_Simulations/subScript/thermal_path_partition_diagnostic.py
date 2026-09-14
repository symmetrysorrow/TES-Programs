"""Diagnostic-only TES/Stycast thermal-path sensitivity scan.

Consumes an Opt_noise.py summary.json, freezes every best-fit parameter except
the two existing Stycast-link conductances, and scans G_tes-stycast x
G_stycast-abs across the same physical bounds used by the production optimizer.

This adds no noise source, thermal state, or production fit parameter.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SIMULATION_ROOT = Path(__file__).resolve().parents[1]
if str(SIMULATION_ROOT) not in sys.path:
    sys.path.insert(0, str(SIMULATION_ROOT))

import Opt_noise as optimizer  # noqa: E402

MID_BAND = "5000-15000_Hz"
HIGH_BAND = "40000-100000_Hz"
TAIL_BAND = "100000-200000_Hz"


def fit_args_from_summary(summary: dict) -> SimpleNamespace:
    fit = summary["fit"]
    return SimpleNamespace(
        fit_min_hz=float(fit["min_hz"]),
        fit_max_hz=float(fit["max_hz"]),
        fit_weight_start_hz=float(fit["weight_start_hz"]),
        high_frequency_weight=float(fit["high_frequency_weight"]),
        robust_delta_dex=float(fit["robust_delta_dex"]),
        fit_points=int(fit["points"]),
        absolute_asd_weight=float(fit.get("absolute_asd_weight", 0.0)),
    )


def grid_with_baseline(
    lower: float,
    upper: float,
    baseline: float,
    points: int,
) -> tuple[float, ...]:
    if points < 2:
        raise ValueError("grid points must be at least 2")
    if lower <= 0.0 or upper <= lower:
        raise ValueError("conductance bounds must be positive and ordered")
    if not lower <= baseline <= upper:
        raise ValueError(
            f"baseline conductance {baseline:g} is outside physical bounds "
            f"{lower:g}--{upper:g}"
        )
    values = list(np.geomspace(lower, upper, points))
    values.append(float(baseline))
    return tuple(sorted(set(float(value) for value in values)))


def band_rows(
    model: np.ndarray,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> dict:
    diagnostics = optimizer.band_fit_diagnostics(model, target, fit_freq, args)
    rows = {}
    for name, values in diagnostics.items():
        mean_log = float(values["mean_log10_ratio"])
        rows[name] = {
            "rms_log10_ratio": float(values["rms_log10_ratio"]),
            "mean_log10_ratio": mean_log,
            "geometric_mean_model_over_measurement": float(10.0**mean_log),
            "max_abs_log10_ratio": float(values["max_abs_log10_ratio"]),
        }
    return rows


def run_diagnostic(
    candidate: dict,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
    *,
    g_tes_stycast_grid: tuple[float, ...] | None = None,
    g_stycast_abs_grid: tuple[float, ...] | None = None,
    grid_points: int = 9,
) -> dict:
    """Scan only the existing TES--Stycast--absorber conductance path."""

    if str(candidate.get("thermal_link_model", "")).lower() != "stycast_node":
        raise ValueError(
            "thermal-path diagnostic requires thermal_link_model='stycast_node'"
        )

    baseline = dict(candidate)
    baseline.pop("thermal_extension", None)
    baseline_point = optimizer.tes_operating_point(baseline)
    if not baseline_point.get("valid") or not baseline_point.get("stable"):
        raise ValueError(
            f"baseline operating point is not stable: "
            f"{baseline_point.get('reason')}"
        )

    baseline_model, _ = optimizer.deterministic_simulated_spectrum(
        baseline.copy(), fit_freq
    )
    baseline_score = float(
        optimizer.fit_score(baseline_model, target, fit_freq, args)
    )
    baseline_bands = band_rows(baseline_model, target, fit_freq, args)
    for required in (MID_BAND, HIGH_BAND, TAIL_BAND):
        if required not in baseline_bands:
            raise ValueError(
                f"fit range does not contain required diagnostic band {required}"
            )

    if g_tes_stycast_grid is None:
        g_tes_stycast_grid = grid_with_baseline(
            optimizer.G_TES_STYCAST_FIT_MIN_W_PER_K,
            optimizer.G_TES_STYCAST_FIT_MAX_W_PER_K,
            float(baseline["G_tes-stycast"]),
            grid_points,
        )
    if g_stycast_abs_grid is None:
        g_stycast_abs_grid = grid_with_baseline(
            optimizer.G_STYCAST_ABS_FIT_MIN_W_PER_K,
            optimizer.G_STYCAST_ABS_FIT_MAX_W_PER_K,
            float(baseline["G_stycast-abs"]),
            grid_points,
        )

    base_mid = float(baseline_bands[MID_BAND]["mean_log10_ratio"])
    base_high = float(baseline_bands[HIGH_BAND]["mean_log10_ratio"])
    base_tail = float(baseline_bands[TAIL_BAND]["mean_log10_ratio"])
    base_current = float(baseline_point["current_A"])
    base_joule = float(baseline_point["joule_power_W"])

    rows = []
    for g_tes_stycast in g_tes_stycast_grid:
        for g_stycast_abs in g_stycast_abs_grid:
            trial = baseline.copy()
            trial["G_tes-stycast"] = float(g_tes_stycast)
            trial["G_stycast-abs"] = float(g_stycast_abs)
            point = optimizer.tes_operating_point(trial)
            row = {
                "G_tes_stycast_W_per_K": float(g_tes_stycast),
                "G_stycast_abs_W_per_K": float(g_stycast_abs),
                "G_tes_stycast_over_baseline": float(
                    g_tes_stycast / baseline["G_tes-stycast"]
                ),
                "G_stycast_abs_over_baseline": float(
                    g_stycast_abs / baseline["G_stycast-abs"]
                ),
                "valid": bool(point.get("valid", False)),
                "stable": bool(point.get("stable", False)),
                "reason": point.get("reason"),
            }
            if not row["valid"] or not row["stable"]:
                rows.append(row)
                continue

            model, _ = optimizer.deterministic_simulated_spectrum(
                trial.copy(), fit_freq
            )
            bands = band_rows(model, target, fit_freq, args)
            mid = float(bands[MID_BAND]["mean_log10_ratio"])
            high = float(bands[HIGH_BAND]["mean_log10_ratio"])
            tail = float(bands[TAIL_BAND]["mean_log10_ratio"])
            row.update(
                {
                    "shape_score": float(
                        optimizer.fit_score(model, target, fit_freq, args)
                    ),
                    "current_ratio_to_baseline": float(
                        point["current_A"] / base_current
                    ),
                    "joule_power_ratio_to_baseline": float(
                        point["joule_power_W"] / base_joule
                    ),
                    "bands": bands,
                    "moves_mid_and_high_in_correct_direction": bool(
                        mid > base_mid and high < base_high
                    ),
                    "improves_mid_and_high_absolute_error": bool(
                        abs(mid) < abs(base_mid) and abs(high) < abs(base_high)
                    ),
                    "improves_mid_and_high_without_worsening_tail": bool(
                        abs(mid) < abs(base_mid)
                        and abs(high) < abs(base_high)
                        and abs(tail) <= abs(base_tail) + 1.0e-12
                    ),
                }
            )
            rows.append(row)

    stable_rows = [
        row for row in rows if row.get("stable") and "shape_score" in row
    ]
    best_global = (
        min(stable_rows, key=lambda row: row["shape_score"])
        if stable_rows
        else None
    )
    correct_direction = [
        row
        for row in stable_rows
        if row["moves_mid_and_high_in_correct_direction"]
    ]
    best_correct_direction = (
        min(
            correct_direction,
            key=lambda row: (
                abs(row["bands"][MID_BAND]["mean_log10_ratio"])
                + abs(row["bands"][HIGH_BAND]["mean_log10_ratio"]),
                row["shape_score"],
            ),
        )
        if correct_direction
        else None
    )
    simultaneous = [
        row
        for row in stable_rows
        if row["improves_mid_and_high_absolute_error"]
    ]
    strict = [
        row
        for row in stable_rows
        if row["improves_mid_and_high_without_worsening_tail"]
    ]

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "added_noise_source": False,
        "added_thermal_state": False,
        "held_fixed_except": ["G_tes-stycast", "G_stycast-abs"],
        "purpose": (
            "Test whether the existing TES-Stycast-absorber path can move the "
            "5-15 kHz deficit upward and the 40-100 kHz excess downward at "
            "fixed best-fit operating/electrical parameters."
        ),
        "baseline": {
            "shape_score": baseline_score,
            "G_tes_stycast_W_per_K": float(baseline["G_tes-stycast"]),
            "G_stycast_abs_W_per_K": float(baseline["G_stycast-abs"]),
            "G_tes_bath_W_per_K": float(baseline["G_tes-bath"]),
            "internal_to_bath_ratio_G_tes_stycast_over_G_tes_bath": float(
                baseline["G_tes-stycast"] / baseline["G_tes-bath"]
            ),
            "current_A": base_current,
            "joule_power_W": base_joule,
            "bands": baseline_bands,
        },
        "grid": {
            "G_tes_stycast_W_per_K": list(map(float, g_tes_stycast_grid)),
            "G_stycast_abs_W_per_K": list(map(float, g_stycast_abs_grid)),
            "total_points": (
                len(g_tes_stycast_grid) * len(g_stycast_abs_grid)
            ),
        },
        "stable_points": len(stable_rows),
        "best_global_shape_score_row": best_global,
        "best_correct_direction_row": best_correct_direction,
        "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction": bool(
            correct_direction
        ),
        "can_improve_5_15k_and_40_100k_absolute_error": bool(simultaneous),
        "can_improve_5_15k_and_40_100k_without_worsening_100_200k": bool(
            strict
        ),
        "best_strict_row": (
            min(strict, key=lambda row: row["shape_score"])
            if strict
            else None
        ),
        "rows": rows,
        "interpretation_guardrail": (
            "A positive screen only shows leverage within the existing reduced "
            "thermal topology. It does not identify a physical conductance "
            "without independent pulse/DC/impedance constraints."
        ),
    }


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary",
        type=Path,
        required=True,
        help="summary.json written by Opt_noise.py",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON path; defaults beside --summary.",
    )
    parser.add_argument(
        "--grid-points",
        type=int,
        default=9,
        help="Log-spaced points per conductance axis before adding the baseline.",
    )
    return parser.parse_args()


def main() -> None:
    args = arguments()
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    candidate = dict(summary["best_case_parameters"])
    fit_args = fit_args_from_summary(summary)
    fit_freq, target, _ = optimizer.target_spectrum(fit_args)
    result = run_diagnostic(
        candidate,
        target,
        fit_freq,
        fit_args,
        grid_points=args.grid_points,
    )
    result["source_summary"] = str(args.summary)
    output = args.output or args.summary.with_name(
        "thermal_path_partition_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "baseline_shape_score": result["baseline"]["shape_score"],
                "best_global_shape_score": (
                    result["best_global_shape_score_row"]["shape_score"]
                    if result["best_global_shape_score_row"]
                    else None
                ),
                "correct_direction": result[
                    "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction"
                ],
                "strict_improvement": result[
                    "can_improve_5_15k_and_40_100k_without_worsening_100_200k"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

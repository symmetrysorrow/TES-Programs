"""Boundary-stress scan of ETF, TES heat capacity, and inductance."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Opt_noise as opt  # noqa: E402

MID = "5000-15000_Hz"
HIGH = "40000-100000_Hz"
TAIL = "100000-200000_Hz"


def fit_args(summary):
    f = summary["fit"]
    return SimpleNamespace(
        fit_min_hz=float(f["min_hz"]),
        fit_max_hz=float(f["max_hz"]),
        fit_weight_start_hz=float(f["weight_start_hz"]),
        high_frequency_weight=float(f["high_frequency_weight"]),
        robust_delta_dex=float(f["robust_delta_dex"]),
        fit_points=int(f["points"]),
        absolute_asd_weight=float(f.get("absolute_asd_weight", 0.0)),
    )


def grid(lower, upper, baseline, points, *, log=False):
    if points < 2 or upper <= lower or not lower <= baseline <= upper:
        raise ValueError("invalid grid or baseline")
    values = (
        np.geomspace(lower, upper, points)
        if log
        else np.linspace(lower, upper, points)
    )
    return tuple(sorted(set(map(float, [*values, baseline]))))


def alpha_inward_grid(baseline, points, minimum_fraction=0.05):
    if baseline <= 0 or not 0 < minimum_fraction < 1:
        raise ValueError("invalid alpha baseline or minimum fraction")
    lower = baseline * minimum_fraction
    return grid(lower, baseline, baseline, points, log=False)


def bands(model, target, freq, args):
    raw = opt.band_fit_diagnostics(model, target, freq, args)
    out = {}
    for name, row in raw.items():
        mean = float(row["mean_log10_ratio"])
        out[name] = {
            "mean_log10_ratio": mean,
            "geometric_mean_model_over_measurement": float(10**mean),
            "rms_log10_ratio": float(row["rms_log10_ratio"]),
        }
    return out


def run(candidate, target, freq, args, alpha_grid, ctes_grid, l_grid):
    base = dict(candidate)
    base.pop("thermal_extension", None)
    point0 = opt.tes_operating_point(base)
    if not point0.get("valid") or not point0.get("stable"):
        raise ValueError("unstable baseline")

    model0, _ = opt.deterministic_simulated_spectrum(base.copy(), freq)
    bands0 = bands(model0, target, freq, args)
    score0 = float(opt.fit_score(model0, target, freq, args))
    bm = float(bands0[MID]["mean_log10_ratio"])
    bh = float(bands0[HIGH]["mean_log10_ratio"])
    bt = float(bands0[TAIL]["mean_log10_ratio"])
    current0 = float(point0["current_A"])
    power0 = float(point0["joule_power_W"])

    rows = []
    for alpha in alpha_grid:
        for ctes in ctes_grid:
            for inductance in l_grid:
                trial = dict(base)
                trial.update(
                    {
                        "alpha": float(alpha),
                        "C_tes": float(ctes),
                        "L": float(inductance),
                    }
                )
                point = opt.tes_operating_point(trial)
                row = {
                    "alpha": float(alpha),
                    "C_tes_J_per_K": float(ctes),
                    "C_tes_over_material": float(
                        ctes / opt.C_TES_MATERIAL_J_PER_K
                    ),
                    "L_H": float(inductance),
                    "valid": bool(point.get("valid")),
                    "stable": bool(point.get("stable")),
                    "reason": point.get("reason"),
                }
                if not row["valid"] or not row["stable"]:
                    rows.append(row)
                    continue

                model, _ = opt.deterministic_simulated_spectrum(
                    trial.copy(), freq
                )
                b = bands(model, target, freq, args)
                mid = float(b[MID]["mean_log10_ratio"])
                high = float(b[HIGH]["mean_log10_ratio"])
                tail = float(b[TAIL]["mean_log10_ratio"])
                row.update(
                    {
                        "shape_score": float(
                            opt.fit_score(model, target, freq, args)
                        ),
                        "current_ratio_to_baseline": float(
                            point["current_A"] / current0
                        ),
                        "joule_power_ratio_to_baseline": float(
                            point["joule_power_W"] / power0
                        ),
                        "bands": b,
                        "correct_direction": bool(mid > bm and high < bh),
                        "improves_mid_high": bool(
                            abs(mid) < abs(bm) and abs(high) < abs(bh)
                        ),
                        "strict_improvement": bool(
                            abs(mid) < abs(bm)
                            and abs(high) < abs(bh)
                            and abs(tail) <= abs(bt) + 1e-12
                        ),
                    }
                )
                rows.append(row)

    stable = [r for r in rows if r.get("stable") and "shape_score" in r]
    correct = [r for r in stable if r["correct_direction"]]
    improve = [r for r in stable if r["improves_mid_high"]]
    strict = [r for r in stable if r["strict_improvement"]]
    max_i = max(
        (abs(r["current_ratio_to_baseline"] - 1.0) for r in stable),
        default=None,
    )
    max_p = max(
        (abs(r["joule_power_ratio_to_baseline"] - 1.0) for r in stable),
        default=None,
    )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "added_noise_source": False,
        "added_thermal_state": False,
        "varied": ["alpha", "C_tes", "L"],
        "held_fixed": [
            "T_c",
            "T_bath",
            "n",
            "R",
            "G_tes-bath",
            "G_tes-stycast",
            "G_stycast-abs",
            "C_stycast",
            "C_abs",
            "beta",
            "R_series_eff",
        ],
        "baseline": {
            "shape_score": score0,
            "alpha": float(base["alpha"]),
            "C_tes_J_per_K": float(base["C_tes"]),
            "C_tes_over_material": float(
                base["C_tes"] / opt.C_TES_MATERIAL_J_PER_K
            ),
            "L_H": float(base["L"]),
            "current_A": current0,
            "joule_power_W": power0,
            "bands": bands0,
        },
        "grid": {
            "alpha": list(map(float, alpha_grid)),
            "C_tes_J_per_K": list(map(float, ctes_grid)),
            "L_H": list(map(float, l_grid)),
            "total_points": len(alpha_grid) * len(ctes_grid) * len(l_grid),
        },
        "stable_points": len(stable),
        "max_abs_current_ratio_minus_one": max_i,
        "max_abs_joule_power_ratio_minus_one": max_p,
        "best_global_shape_score_row": min(
            stable, key=lambda r: r["shape_score"], default=None
        ),
        "best_correct_direction_row": min(
            correct,
            key=lambda r: (
                abs(r["bands"][MID]["mean_log10_ratio"])
                + abs(r["bands"][HIGH]["mean_log10_ratio"]),
                r["shape_score"],
            ),
            default=None,
        ),
        "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction": bool(
            correct
        ),
        "can_improve_5_15k_and_40_100k_absolute_error": bool(improve),
        "can_improve_5_15k_and_40_100k_without_worsening_100_200k": bool(
            strict
        ),
        "best_strict_row": min(
            strict, key=lambda r: r["shape_score"], default=None
        ),
        "rows": rows,
        "guardrail": (
            "A positive screen shows dynamical leverage only. alpha, C_tes, "
            "and L require independent pulse or complex-impedance constraints "
            "before physical interpretation."
        ),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--output", type=Path)
    p.add_argument("--alpha-points", type=int, default=7)
    p.add_argument("--ctes-points", type=int, default=7)
    p.add_argument("--l-points", type=int, default=7)
    p.add_argument("--alpha-min-fraction", type=float, default=0.05)
    a = p.parse_args()

    summary = json.loads(a.summary.read_text(encoding="utf-8"))
    candidate = dict(summary["best_case_parameters"])
    args = fit_args(summary)
    freq, target, _ = opt.target_spectrum(args)
    result = run(
        candidate,
        target,
        freq,
        args,
        alpha_inward_grid(
            float(candidate["alpha"]),
            a.alpha_points,
            a.alpha_min_fraction,
        ),
        grid(
            opt.C_TES_FIT_MIN_J_PER_K,
            opt.C_TES_FIT_MAX_J_PER_K,
            float(candidate["C_tes"]),
            a.ctes_points,
            log=True,
        ),
        grid(
            opt.L_FIT_MIN_H,
            opt.L_FIT_MAX_H,
            float(candidate["L"]),
            a.l_points,
            log=True,
        ),
    )
    result["source_summary"] = str(a.summary)
    output = a.output or a.summary.with_name(
        "etf_ctes_inductance_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(output),
        "stable_points": result["stable_points"],
        "correct_direction": result[
            "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction"
        ],
        "strict_improvement": result[
            "can_improve_5_15k_and_40_100k_without_worsening_100_200k"
        ],
        "best_global_shape_score": (
            result["best_global_shape_score_row"]["shape_score"]
            if result["best_global_shape_score_row"] else None
        ),
    }, indent=2))


if __name__ == "__main__":
    main()

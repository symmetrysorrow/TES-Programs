"""Scan IV-compatible TES-bath ETF and existing TES-Stycast coupling."""

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


def target_iv_power(summary):
    best_rsh = summary.get("best_case_R_SH_ohm")
    cases = summary.get("cases", [])
    for case in cases:
        if best_rsh is None or np.isclose(
            float(case["R_SH_ohm"]), float(best_rsh), rtol=0.0, atol=1e-12
        ):
            power = float(case["iv_operating_point"]["P_J_W"])
            if np.isfinite(power) and power > 0:
                return power
    raise ValueError("selected IV operating-point P_J_W is missing")


def grid(lower, upper, baseline, points, *, log=False):
    if points < 2 or upper <= lower:
        raise ValueError("invalid grid")
    if not lower <= baseline <= upper:
        raise ValueError("baseline outside grid")
    values = (
        np.geomspace(lower, upper, points)
        if log
        else np.linspace(lower, upper, points)
    )
    return tuple(sorted(set(map(float, [*values, baseline]))))


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


def run(candidate, iv_power, target, freq, args, tc_grid, n_grid, gst_grid):
    base = dict(candidate)
    base.pop("thermal_extension", None)
    if base.get("thermal_link_model") != "stycast_node":
        raise ValueError("stycast_node model required")

    base_point = opt.tes_operating_point(base)
    if not base_point.get("stable") or not base_point.get("valid"):
        raise ValueError("unstable baseline")
    base_model, _ = opt.deterministic_simulated_spectrum(base.copy(), freq)
    base_bands = bands(base_model, target, freq, args)
    base_score = float(opt.fit_score(base_model, target, freq, args))
    bm = float(base_bands[MID]["mean_log10_ratio"])
    bh = float(base_bands[HIGH]["mean_log10_ratio"])
    bt = float(base_bands[TAIL]["mean_log10_ratio"])

    rows = []
    for tc in tc_grid:
        for exponent in n_grid:
            if tc <= base["T_bath"]:
                continue
            gb = opt.g_tes_bath_from_joule_power(
                iv_power, tc, base["T_bath"], exponent
            )
            for gst in gst_grid:
                trial = dict(base)
                trial.update(
                    {
                        "T_c": float(tc),
                        "n": float(exponent),
                        "G_tes-bath": float(gb),
                        "G_tes-stycast": float(gst),
                    }
                )
                point = opt.tes_operating_point(trial)
                row = {
                    "T_c_K": float(tc),
                    "n": float(exponent),
                    "G_tes_bath_W_per_K": float(gb),
                    "G_tes_stycast_W_per_K": float(gst),
                    "G_tes_stycast_over_G_tes_bath": float(gst / gb),
                    "valid": bool(point.get("valid")),
                    "stable": bool(point.get("stable")),
                    "reason": point.get("reason"),
                }
                if not row["valid"] or not row["stable"]:
                    rows.append(row)
                    continue

                model, _ = opt.deterministic_simulated_spectrum(trial, freq)
                b = bands(model, target, freq, args)
                mid = float(b[MID]["mean_log10_ratio"])
                high = float(b[HIGH]["mean_log10_ratio"])
                tail = float(b[TAIL]["mean_log10_ratio"])
                pj = float(point["joule_power_W"])
                row.update(
                    {
                        "shape_score": float(opt.fit_score(model, target, freq, args)),
                        "iv_power_relative_error": float((pj - iv_power) / iv_power),
                        "current_ratio_to_baseline": float(
                            point["current_A"] / base_point["current_A"]
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
    max_power_error = max(
        (abs(r["iv_power_relative_error"]) for r in stable), default=None
    )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "added_noise_source": False,
        "added_thermal_state": False,
        "varied": ["T_c", "n", "G_tes-stycast"],
        "derived_each_trial": "G_tes-bath from fixed IV P_J",
        "fixed_IV_joule_power_W": float(iv_power),
        "baseline": {
            "shape_score": base_score,
            "T_c_K": float(base["T_c"]),
            "n": float(base["n"]),
            "G_tes_bath_W_per_K": float(base["G_tes-bath"]),
            "G_tes_stycast_W_per_K": float(base["G_tes-stycast"]),
            "bands": base_bands,
        },
        "grid": {
            "T_c_K": list(map(float, tc_grid)),
            "n": list(map(float, n_grid)),
            "G_tes_stycast_W_per_K": list(map(float, gst_grid)),
            "total_points": len(tc_grid) * len(n_grid) * len(gst_grid),
        },
        "stable_points": len(stable),
        "max_abs_iv_power_relative_error": max_power_error,
        "all_stable_points_preserve_iv_power": bool(
            stable and max_power_error is not None and max_power_error <= 1e-10
        ),
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
        "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction": bool(correct),
        "can_improve_5_15k_and_40_100k_absolute_error": bool(improve),
        "can_improve_5_15k_and_40_100k_without_worsening_100_200k": bool(strict),
        "best_strict_row": min(
            strict, key=lambda r: r["shape_score"], default=None
        ),
        "rows": rows,
        "guardrail": (
            "Positive leverage is not an independent parameter measurement; "
            "pulse/DC/complex-impedance constraints are still required."
        ),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--output", type=Path)
    p.add_argument("--tc-points", type=int, default=7)
    p.add_argument("--n-points", type=int, default=7)
    p.add_argument("--g-points", type=int, default=7)
    a = p.parse_args()

    summary = json.loads(a.summary.read_text(encoding="utf-8"))
    candidate = dict(summary["best_case_parameters"])
    if "T_c_fit_range_K" not in summary:
        raise ValueError("summary.json is missing T_c_fit_range_K")

    args = fit_args(summary)
    freq, target, _ = opt.target_spectrum(args)
    tc_lo, tc_hi = map(float, summary["T_c_fit_range_K"])
    result = run(
        candidate,
        target_iv_power(summary),
        target,
        freq,
        args,
        grid(tc_lo, tc_hi, candidate["T_c"], a.tc_points),
        grid(opt.N_FIT_MIN, opt.N_FIT_MAX, candidate["n"], a.n_points),
        grid(
            opt.G_TES_STYCAST_FIT_MIN_W_PER_K,
            opt.G_TES_STYCAST_FIT_MAX_W_PER_K,
            candidate["G_tes-stycast"],
            a.g_points,
            log=True,
        ),
    )
    result["source_summary"] = str(a.summary)
    output = a.output or a.summary.with_name("iv_power_thermal_path_diagnostic.json")
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "stable_points": result["stable_points"],
        "iv_power_preserved": result["all_stable_points_preserve_iv_power"],
        "correct_direction": result[
            "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction"
        ],
        "strict_improvement": result[
            "can_improve_5_15k_and_40_100k_without_worsening_100_200k"
        ],
    }, indent=2))


if __name__ == "__main__":
    main()

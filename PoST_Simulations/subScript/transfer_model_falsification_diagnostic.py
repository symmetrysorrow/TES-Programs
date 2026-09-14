"""Fit low-order stable pole/zero magnitudes to the required residual transfer."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import differential_evolution, least_squares

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Opt_noise as opt  # noqa: E402

FAMILIES = (
    ("unity", 0, 0),
    ("1p", 1, 0),
    ("1p1z", 1, 1),
    ("2p1z", 2, 1),
    ("2p2z", 2, 2),
    ("3p2z", 3, 2),
    ("3p3z", 3, 3),
)


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


def pole_zero_transfer(frequency, poles_hz, zeros_hz, reference_hz=1_000.0):
    """Magnitude of a stable minimum-phase real pole/zero transfer, normalized."""

    f = np.asarray(frequency, dtype=float)
    if np.any(f <= 0.0) or reference_hz <= 0.0:
        raise ValueError("frequencies and reference must be positive")
    poles = np.asarray(poles_hz, dtype=float)
    zeros = np.asarray(zeros_hz, dtype=float)
    if np.any(poles <= 0.0) or np.any(zeros <= 0.0):
        raise ValueError("pole and zero corners must be positive")

    log_mag = np.zeros_like(f)
    ref_log_mag = 0.0
    for corner in zeros:
        log_mag += 0.5 * np.log1p((f / corner) ** 2)
        ref_log_mag += 0.5 * np.log1p((reference_hz / corner) ** 2)
    for corner in poles:
        log_mag -= 0.5 * np.log1p((f / corner) ** 2)
        ref_log_mag -= 0.5 * np.log1p((reference_hz / corner) ** 2)
    return np.exp(log_mag - ref_log_mag)


def decode(log_corners, n_poles, n_zeros):
    x = np.asarray(log_corners, dtype=float)
    poles = np.sort(10.0 ** x[:n_poles])
    zeros = np.sort(10.0 ** x[n_poles:n_poles + n_zeros])
    return poles, zeros


def residual_db(log_corners, n_poles, n_zeros, frequency, required):
    poles, zeros = decode(log_corners, n_poles, n_zeros)
    fitted = pole_zero_transfer(frequency, poles, zeros)
    return 20.0 * np.log10(fitted / required)


def band_metrics(residual, frequency, args):
    rows = {}
    for low, high, _weight in opt.FIT_BANDS_HZ:
        lo = max(float(low), float(args.fit_min_hz))
        hi = min(float(high), float(args.fit_max_hz))
        mask = (frequency >= lo) & (frequency <= hi)
        if not np.any(mask):
            continue
        values = np.asarray(residual[mask], dtype=float)
        rows[f"{lo:g}-{hi:g}_Hz"] = {
            "mean_residual_dB": float(np.mean(values)),
            "rms_residual_dB": float(np.sqrt(np.mean(values**2))),
            "max_abs_residual_dB": float(np.max(np.abs(values))),
        }
    return rows


def fit_family(
    name,
    n_poles,
    n_zeros,
    frequency,
    required,
    model,
    target,
    args,
    *,
    corner_min_hz,
    corner_max_hz,
    seed,
    de_maxiter,
    rms_threshold_db,
    max_threshold_db,
):
    n_parameters = n_poles + n_zeros
    if n_parameters == 0:
        fitted = np.ones_like(frequency)
        poles = np.asarray([], dtype=float)
        zeros = np.asarray([], dtype=float)
        success = True
        evaluations = 0
    else:
        lower = np.log10(float(corner_min_hz))
        upper = np.log10(float(corner_max_hz))
        bounds = [(lower, upper)] * n_parameters

        def objective(x):
            r = residual_db(
                x, n_poles, n_zeros, frequency, required
            )
            return float(np.mean(r**2))

        de = differential_evolution(
            objective,
            bounds,
            seed=int(seed),
            maxiter=int(de_maxiter),
            popsize=8,
            polish=False,
            tol=1e-8,
            workers=1,
            updating="immediate",
        )
        ls = least_squares(
            lambda x: residual_db(
                x, n_poles, n_zeros, frequency, required
            ),
            de.x,
            bounds=(lower, upper),
            max_nfev=2000,
        )
        poles, zeros = decode(ls.x, n_poles, n_zeros)
        fitted = pole_zero_transfer(frequency, poles, zeros)
        success = bool(de.success or ls.success)
        evaluations = int(de.nfev + ls.nfev)

    r_db = 20.0 * np.log10(fitted / required)
    rms_db = float(np.sqrt(np.mean(r_db**2)))
    max_db = float(np.max(np.abs(r_db)))
    corrected_model = model * fitted
    score = float(opt.fit_score(corrected_model, target, frequency, args))
    passes = bool(
        rms_db <= float(rms_threshold_db)
        and max_db <= float(max_threshold_db)
    )
    return {
        "family": name,
        "n_poles": int(n_poles),
        "n_zeros": int(n_zeros),
        "n_parameters": int(n_parameters),
        "proper": bool(n_poles >= n_zeros),
        "stable_minimum_phase_real_sections": True,
        "success": success,
        "evaluations": evaluations,
        "poles_Hz": list(map(float, poles)),
        "zeros_Hz": list(map(float, zeros)),
        "rms_residual_dB": rms_db,
        "max_abs_residual_dB": max_db,
        "mean_residual_dB": float(np.mean(r_db)),
        "corrected_shape_score": score,
        "passes_screen_tolerance": passes,
        "band_residuals": band_metrics(r_db, frequency, args),
        "_fitted_transfer": fitted,
        "_residual_dB": r_db,
    }


def run(
    candidate,
    target,
    frequency,
    args,
    *,
    corner_min_hz=300.0,
    corner_max_hz=2_000_000.0,
    seed=20260914,
    de_maxiter=80,
    rms_threshold_db=1.0,
    max_threshold_db=3.0,
):
    model, _ = opt.deterministic_simulated_spectrum(
        dict(candidate), frequency
    )
    req_summary, req_curves = opt.required_transfer_diagnostics(
        candidate, model, target, frequency, args
    )
    required = np.asarray(
        req_curves["required_ASD_transfer"], dtype=float
    )
    baseline_score = float(opt.fit_score(model, target, frequency, args))

    rows = []
    for index, (name, n_poles, n_zeros) in enumerate(FAMILIES):
        rows.append(
            fit_family(
                name,
                n_poles,
                n_zeros,
                frequency,
                required,
                model,
                target,
                args,
                corner_min_hz=corner_min_hz,
                corner_max_hz=corner_max_hz,
                seed=seed + index,
                de_maxiter=de_maxiter,
                rms_threshold_db=rms_threshold_db,
                max_threshold_db=max_threshold_db,
            )
        )

    best_rms = min(rows, key=lambda row: row["rms_residual_dB"])
    best_score = min(rows, key=lambda row: row["corrected_shape_score"])
    passing = [
        row for row in rows if row["passes_screen_tolerance"]
    ]
    smallest_passing = (
        min(
            passing,
            key=lambda row: (
                row["n_parameters"],
                row["rms_residual_dB"],
            ),
        )
        if passing
        else None
    )

    sample_indices = np.unique(
        np.linspace(
            0,
            len(frequency) - 1,
            min(81, len(frequency)),
            dtype=int,
        )
    )
    sample = []
    best_fit = np.asarray(best_rms["_fitted_transfer"], dtype=float)
    for index in sample_indices:
        sample.append(
            {
                "frequency_Hz": float(frequency[index]),
                "required_ASD_transfer": float(required[index]),
                "required_correction_dB": float(
                    20.0 * np.log10(required[index])
                ),
                "best_low_order_transfer": float(best_fit[index]),
                "best_low_order_correction_dB": float(
                    20.0 * np.log10(best_fit[index])
                ),
                "best_low_order_residual_dB": float(
                    best_rms["_residual_dB"][index]
                ),
            }
        )

    clean_rows = []
    for row in rows:
        clean_rows.append(
            {
                key: value
                for key, value in row.items()
                if not key.startswith("_")
            }
        )
    clean_best_rms = next(
        row for row in clean_rows
        if row["family"] == best_rms["family"]
    )
    clean_best_score = next(
        row for row in clean_rows
        if row["family"] == best_score["family"]
    )
    clean_smallest = (
        next(
            row for row in clean_rows
            if row["family"] == smallest_passing["family"]
        )
        if smallest_passing
        else None
    )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "detector_parameters_held_fixed": True,
        "definition": req_summary["definition"],
        "model_class": (
            "post-model multiplicative ASD magnitude built from stable "
            "minimum-phase real first-order poles/zeros and normalized at 1 kHz"
        ),
        "important_limitation": (
            "Magnitude-only success is not identification of a physical "
            "readout location or proof of the true phase response. A failure "
            "does falsify these tested low-order single-path magnitude forms "
            "within the stated corner-frequency box and numerical tolerance."
        ),
        "corner_search_Hz": [
            float(corner_min_hz),
            float(corner_max_hz),
        ],
        "screen_tolerance": {
            "rms_residual_dB_max": float(rms_threshold_db),
            "max_abs_residual_dB_max": float(max_threshold_db),
            "configurable_not_physical_prior": True,
        },
        "baseline_shape_score": baseline_score,
        "required_transfer_context": req_summary,
        "families": clean_rows,
        "any_low_order_family_passes_screen_tolerance": bool(passing),
        "smallest_passing_family": clean_smallest,
        "best_by_rms_dB": clean_best_rms,
        "best_by_corrected_shape_score": clean_best_score,
        "curve_sample_best_rms_family": sample,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--output", type=Path)
    p.add_argument("--corner-min-hz", type=float, default=300.0)
    p.add_argument("--corner-max-hz", type=float, default=2_000_000.0)
    p.add_argument("--de-maxiter", type=int, default=80)
    p.add_argument("--seed", type=int, default=20260914)
    p.add_argument("--rms-threshold-db", type=float, default=1.0)
    p.add_argument("--max-threshold-db", type=float, default=3.0)
    a = p.parse_args()

    summary = json.loads(a.summary.read_text(encoding="utf-8"))
    candidate = dict(summary["best_case_parameters"])
    args = fit_args(summary)
    frequency, target, _ = opt.target_spectrum(args)
    result = run(
        candidate,
        target,
        frequency,
        args,
        corner_min_hz=a.corner_min_hz,
        corner_max_hz=a.corner_max_hz,
        seed=a.seed,
        de_maxiter=a.de_maxiter,
        rms_threshold_db=a.rms_threshold_db,
        max_threshold_db=a.max_threshold_db,
    )
    result["source_summary"] = str(a.summary)
    output = a.output or a.summary.with_name(
        "transfer_model_falsification_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(output),
        "baseline_shape_score": result["baseline_shape_score"],
        "any_family_passes": result[
            "any_low_order_family_passes_screen_tolerance"
        ],
        "smallest_passing_family": (
            result["smallest_passing_family"]["family"]
            if result["smallest_passing_family"] else None
        ),
        "best_rms_family": result["best_by_rms_dB"]["family"],
        "best_rms_dB": result["best_by_rms_dB"]["rms_residual_dB"],
        "best_max_abs_dB": result[
            "best_by_rms_dB"
        ]["max_abs_residual_dB"],
    }, indent=2))


if __name__ == "__main__":
    main()

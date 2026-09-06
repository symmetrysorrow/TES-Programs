"""Scan intermediate high-frequency targets between raw trend and lower baseline.

Exploratory/noise-guided only. Above 10 kHz, optimize a geometric interpolation
between the asymmetric lower baseline and a smoothed raw experimental trend:

    log T_lambda = (1-lambda) log T_lower + lambda log T_raw_smooth

Each lambda is optimized independently. Final candidates are selected against a
fixed lambda=0.5 midpoint target, so the scan explores parameter basins without
letting the convex feature itself define the final target.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import optimize, sparse
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations"))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))

from explore_high_frequency_noise_parameters import HARDWARE_BESSEL_ORDER, pulse_gate, rms_log_ratio  # noqa: E402
from explore_high_frequency_noise_parameters_adaptive import SEARCH_PARAMETERS, parameters_from_unit, specs_from_envelope, unit_from_parameters  # noqa: E402
from high_frequency_noise_comparison import experiment_asd, log_interp, normalized  # noqa: E402
from noise_measurement_model import HARDWARE_BESSEL_CUTOFF_HZ  # noqa: E402
from refine_hardware_bessel_lower_envelope import (  # noqa: E402
    HARDWARE_MODES,
    deterministic_normalized_curve,
    dump,
    lower_envelope_targets,
    realize_mode,
    seed_parameters,
)

STRICT = "C — exact target physical case remains unidentified"
DEFAULT_LAMBDAS = (0.20, 0.35, 0.50, 0.65, 0.80)
BAND_WEIGHTS = {"mid": 1.0, "high": 1.0, "edge": 1.25}


def scan_grids(edge_max_hz: float = 120_000.0) -> dict[str, np.ndarray]:
    edge_max_hz = float(edge_max_hz)
    if edge_max_hz <= 80_000.0:
        raise ValueError("edge_max_hz must exceed 80 kHz")
    return {
        "mid": np.logspace(np.log10(1_000.0), np.log10(10_000.0), 180),
        "high": np.logspace(np.log10(10_000.0), np.log10(80_000.0), 240),
        "edge": np.logspace(np.log10(80_000.0), np.log10(edge_max_hz), 140),
        "all": np.logspace(np.log10(1_000.0), np.log10(edge_max_hz), 520),
    }


def symmetric_log_trend(log_values: np.ndarray, smoothness: float = 6_000.0) -> np.ndarray:
    y = np.asarray(log_values, dtype=float)
    if y.ndim != 1 or len(y) < 5 or np.any(~np.isfinite(y)):
        raise ValueError("log_values must be a finite one-dimensional array")
    if smoothness <= 0.0:
        raise ValueError("smoothness must be positive")
    n = len(y)
    difference = sparse.diags(
        [np.ones(n - 2), -2.0 * np.ones(n - 2), np.ones(n - 2)],
        [0, 1, 2],
        shape=(n - 2, n),
        format="csc",
    )
    penalty = float(smoothness) * (difference.T @ difference)
    return np.asarray(spsolve(sparse.eye(n, format="csc") + penalty, y), dtype=float)


def geometric_intermediate(lower: np.ndarray, upper: np.ndarray, lam: float) -> np.ndarray:
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    if lower.shape != upper.shape:
        raise ValueError("lower and upper must have identical shapes")
    if np.any(lower <= 0.0) or np.any(upper <= 0.0):
        raise ValueError("ASD targets must be positive")
    lam = float(lam)
    if not 0.0 <= lam <= 1.0:
        raise ValueError("lambda must lie in [0, 1]")
    return np.exp((1.0 - lam) * np.log(lower) + lam * np.log(upper))


def target_from_support(raw_values, grid, support_frequency, support_lower, support_upper, lam, blend_end_hz):
    grid = np.asarray(grid, dtype=float)
    raw_values = np.asarray(raw_values, dtype=float)
    values = raw_values.copy()
    high = grid >= 10_000.0
    if not np.any(high):
        return values
    lower = np.exp(np.interp(np.log(grid[high]), np.log(support_frequency), np.log(support_lower)))
    upper = np.exp(np.interp(np.log(grid[high]), np.log(support_frequency), np.log(support_upper)))
    target = geometric_intermediate(lower, upper, lam)
    if blend_end_hz <= 10_000.0:
        blend = np.ones(np.count_nonzero(high))
    else:
        blend = np.clip(
            (np.log(grid[high]) - np.log(10_000.0))
            / (np.log(float(blend_end_hz)) - np.log(10_000.0)),
            0.0,
            1.0,
        )
    values[high] = np.exp((1.0 - blend) * np.log(raw_values[high]) + blend * np.log(target))
    return values


def build_target_family(exp_freq, exp_norm, grids, lambdas, lower_smoothness, lower_asymmetry, trend_smoothness, blend_end_hz, support_max_hz):
    raw_targets, lower_targets, metadata = lower_envelope_targets(
        exp_freq,
        exp_norm,
        grids,
        smoothness=lower_smoothness,
        asymmetry=lower_asymmetry,
        blend_end_hz=blend_end_hz,
        support_max_hz=support_max_hz,
    )
    support_frequency = np.asarray(metadata["support_frequency_Hz"], dtype=float)
    support_raw = np.asarray(metadata["support_raw_normalized_asd"], dtype=float)
    support_lower = np.asarray(metadata["support_lower_baseline_normalized_asd"], dtype=float)
    support_upper = np.exp(symmetric_log_trend(np.log(support_raw), smoothness=trend_smoothness))
    support_upper = np.maximum(support_upper, support_lower)

    family = {}
    all_lambdas = sorted(set(float(v) for v in lambdas) | {0.5})
    for lam in all_lambdas:
        family[f"{lam:.6f}"] = {
            name: target_from_support(
                raw_targets[name], grid, support_frequency, support_lower,
                support_upper, lam, blend_end_hz
            )
            for name, grid in grids.items()
        }
    upper_targets = {
        name: target_from_support(
            raw_targets[name], grid, support_frequency, support_lower,
            support_upper, 1.0, blend_end_hz
        )
        for name, grid in grids.items()
    }
    metadata = {
        **metadata,
        "raw_trend_smoothness": float(trend_smoothness),
        "interpolation": "geometric_log_ASD_between_lower_and_smoothed_raw_trend",
        "selection_lambda": 0.5,
        "support_smoothed_raw_normalized_asd": support_upper.tolist(),
        "scanned_lambdas": sorted(set(float(v) for v in lambdas)),
    }
    return {
        "raw": raw_targets,
        "lower": lower_targets,
        "upper": upper_targets,
        "midpoint": family["0.500000"],
    }, family, metadata


def combined_frequency(grids):
    return np.unique(np.concatenate([np.asarray(grids[name], dtype=float) for name in BAND_WEIGHTS]))


def band_scores(frequency, simulated, targets, grids):
    out = {}
    weighted = 0.0
    total_weight = 0.0
    for name, weight in BAND_WEIGHTS.items():
        curve = log_interp(frequency, simulated, grids[name])
        rms, maximum = rms_log_ratio(curve, targets[name])
        out[name] = {"rms_log_ratio": rms, "max_abs_log_ratio": maximum}
        weighted += float(weight) * rms**2
        total_weight += float(weight)
    out["weighted_rms_log_ratio"] = float(np.sqrt(weighted / total_weight))
    return out


def bracket_position(frequency, simulated, lower_targets, upper_targets, grids):
    out = {}
    for name in ("high", "edge"):
        curve = log_interp(frequency, simulated, grids[name])
        lower = np.asarray(lower_targets[name], dtype=float)
        upper = np.asarray(upper_targets[name], dtype=float)
        width = np.log(upper) - np.log(lower)
        valid = width > 1e-3
        if not np.any(valid):
            out[name] = {"mean": None, "median": None, "count": 0}
            continue
        z = (np.log(curve[valid]) - np.log(lower[valid])) / width[valid]
        out[name] = {"mean": float(np.mean(z)), "median": float(np.median(z)), "count": int(np.count_nonzero(valid))}
    return out


def seed_pool_for_mode(search_doc, lower_doc, mode, pareto_seed_count):
    seeds = seed_parameters(search_doc, pareto_seed_count)
    if lower_doc:
        row = lower_doc.get("best_by_hardware_mode", {}).get(mode)
        if row and "parameters" in row:
            seeds.append({"seed_id": f"previous_lower_{mode}", "parameters": dict(row["parameters"])})
    return list({seed["seed_id"]: seed for seed in seeds}.values())


def fit_one(seed, mode, lam, specs, template, envelope, constraints, grids, target, reference_targets, rate, cutoff, max_nfev):
    x0 = unit_from_parameters(seed["parameters"], specs)
    fit_frequency = combined_frequency(grids)

    def residual(unit):
        parameters = parameters_from_unit(unit, specs, template, envelope)
        ok, _gate = pulse_gate(parameters, constraints)
        if not ok:
            return np.full(sum(len(grids[name]) for name in BAND_WEIGHTS), 6.0)
        try:
            simulated = deterministic_normalized_curve(parameters, fit_frequency, rate, cutoff, mode)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError, OverflowError):
            return np.full(sum(len(grids[name]) for name in BAND_WEIGHTS), 6.0)
        pieces = []
        for name, weight in BAND_WEIGHTS.items():
            curve = log_interp(fit_frequency, simulated, grids[name])
            r = np.log(curve / target[name])
            pieces.append(r * np.sqrt(float(weight) / len(r)))
        return np.concatenate(pieces)

    result = optimize.least_squares(
        residual,
        x0,
        bounds=(np.zeros(len(specs)), np.ones(len(specs))),
        method="trf",
        jac="2-point",
        loss="soft_l1",
        f_scale=0.08,
        x_scale="jac",
        max_nfev=int(max_nfev),
    )
    parameters = parameters_from_unit(result.x, specs, template, envelope)
    gate_ok, gate = pulse_gate(parameters, constraints)
    simulated = deterministic_normalized_curve(parameters, fit_frequency, rate, cutoff, mode)
    return {
        "seed_id": seed["seed_id"],
        "hardware_mode": mode,
        "lambda": float(lam),
        "parameters": parameters,
        "pulse_gate_ok": bool(gate_ok),
        "pulse_consistency": gate,
        "optimizer": {
            "success": bool(result.success),
            "status": int(result.status),
            "cost": float(result.cost),
            "optimality": float(result.optimality),
            "nfev": int(result.nfev),
        },
        "fit_target_scores": band_scores(fit_frequency, simulated, target, grids),
        "midpoint_validation_scores": band_scores(fit_frequency, simulated, reference_targets["midpoint"], grids),
        "raw_experiment_scores": band_scores(fit_frequency, simulated, reference_targets["raw"], grids),
        "lower_target_scores": band_scores(fit_frequency, simulated, reference_targets["lower"], grids),
        "smoothed_raw_scores": band_scores(fit_frequency, simulated, reference_targets["upper"], grids),
        "bracket_position": bracket_position(
            fit_frequency, simulated, reference_targets["lower"], reference_targets["upper"], grids
        ),
    }


def best_lambda_candidates(rows, modes):
    by_mode_lambda = {}
    best_mode = {}
    for mode in modes:
        mode_rows = [r for r in rows if r["hardware_mode"] == mode and r["pulse_gate_ok"]]
        lambda_rows = {}
        for lam in sorted(set(r["lambda"] for r in mode_rows)):
            candidates = [r for r in mode_rows if np.isclose(r["lambda"], lam)]
            if candidates:
                lambda_rows[f"{lam:.6f}"] = min(
                    candidates,
                    key=lambda r: r["fit_target_scores"]["weighted_rms_log_ratio"],
                )
        by_mode_lambda[mode] = lambda_rows
        if lambda_rows:
            best_mode[mode] = min(
                lambda_rows.values(),
                key=lambda r: r["midpoint_validation_scores"]["weighted_rms_log_ratio"],
            )
    return by_mode_lambda, best_mode


def finite_metrics(frequency, values, reference_targets, grids):
    return {
        "midpoint_validation_scores": band_scores(frequency, values, reference_targets["midpoint"], grids),
        "raw_experiment_scores": band_scores(frequency, values, reference_targets["raw"], grids),
        "lower_target_scores": band_scores(frequency, values, reference_targets["lower"], grids),
        "smoothed_raw_scores": band_scores(frequency, values, reference_targets["upper"], grids),
        "bracket_position": bracket_position(
            frequency, values, reference_targets["lower"], reference_targets["upper"], grids
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--search-json", type=Path)
    ap.add_argument("--lower-json", type=Path)
    ap.add_argument("--hardware-modes", nargs="+", choices=HARDWARE_MODES, default=list(HARDWARE_MODES))
    ap.add_argument("--lambdas", nargs="+", type=float, default=list(DEFAULT_LAMBDAS))
    ap.add_argument("--pareto-seeds", type=int, default=3)
    ap.add_argument("--max-nfev", type=int, default=90)
    ap.add_argument("--lower-smoothness", type=float, default=80_000.0)
    ap.add_argument("--lower-asymmetry", type=float, default=0.08)
    ap.add_argument("--raw-trend-smoothness", type=float, default=6_000.0)
    ap.add_argument("--blend-end-hz", type=float, default=20_000.0)
    ap.add_argument("--support-max-hz", type=float, default=150_000.0)
    ap.add_argument("--edge-max-hz", type=float, default=120_000.0)
    ap.add_argument("--finite-records", type=int, default=0)
    ap.add_argument("--finite-seed", type=int, default=20260906)
    args = ap.parse_args()

    lambdas = sorted(set(float(v) for v in args.lambdas))
    if not lambdas or any(v < 0.0 or v > 1.0 for v in lambdas):
        raise ValueError("all lambdas must lie in [0, 1]")

    case = args.case_dir
    search_json = args.search_json or case / "correlated_high_frequency_parameter_search.json"
    search_doc = json.loads(search_json.read_text(encoding="utf-8"))
    lower_json = args.lower_json or case / "lower_envelope_hardware_refinement.json"
    lower_doc = json.loads(lower_json.read_text(encoding="utf-8")) if lower_json.exists() else None
    scenarios = json.loads((case / "proxy_scenarios.json").read_text(encoding="utf-8"))
    envelope = json.loads((case / "proxy_parameter_envelope.json").read_text(encoding="utf-8"))
    constraints = json.loads((case / "pulse_combination_constraints.json").read_text(encoding="utf-8"))
    frozen = scenarios["pulse_consistent_scenarios"]
    if not frozen:
        raise RuntimeError("no frozen pulse-consistent scenarios")

    exp_freq, exp_asd, exp_paths, acquisition = experiment_asd(args.target_root)
    rate = float(acquisition["rate_Hz"])
    sample = int(acquisition["samples"])
    cutoff = float(acquisition["analysis_bessel_cutoff_Hz"])
    if args.edge_max_hz >= rate / 2.0:
        raise ValueError("edge_max_hz must lie below Nyquist")
    exp_norm = normalized(exp_asd, exp_freq)
    grids = scan_grids(args.edge_max_hz)
    reference_targets, target_family, target_metadata = build_target_family(
        exp_freq,
        exp_norm,
        grids,
        lambdas,
        lower_smoothness=args.lower_smoothness,
        lower_asymmetry=args.lower_asymmetry,
        trend_smoothness=args.raw_trend_smoothness,
        blend_end_hz=args.blend_end_hz,
        support_max_hz=max(args.support_max_hz, args.edge_max_hz),
    )

    template = dict(frozen[0]["parameters"])
    specs = specs_from_envelope(envelope, template, excess_max=5.0)
    fit_rows = []
    for mode in args.hardware_modes:
        seeds = seed_pool_for_mode(search_doc, lower_doc, mode, args.pareto_seeds)
        for lam in lambdas:
            target = target_family[f"{lam:.6f}"]
            for seed in seeds:
                fit_rows.append(
                    fit_one(
                        seed, mode, lam, specs, template, envelope, constraints,
                        grids, target, reference_targets, rate, cutoff, args.max_nfev
                    )
                )
    if not fit_rows:
        raise RuntimeError("no intermediate-target fits were produced")

    by_mode_lambda, best_by_mode = best_lambda_candidates(fit_rows, list(args.hardware_modes))
    if not best_by_mode:
        raise RuntimeError("no pulse-consistent intermediate-target solution")

    records = args.finite_records if args.finite_records > 0 else len(exp_paths)
    finite = {}
    for mode, row in best_by_mode.items():
        frequency, values = realize_mode(row, mode, sample, rate, cutoff, records, args.finite_seed)
        row["finite_validation"] = finite_metrics(frequency, values, reference_targets, grids)
        finite[mode] = {"frequency": frequency, "normalized_asd": values}

    physical = {m: r for m, r in best_by_mode.items() if m != "bypass"}
    best_physical = min(
        physical.values(),
        key=lambda r: r["finite_validation"]["midpoint_validation_scores"]["weighted_rms_log_ratio"],
    ) if physical else None
    best_including_bypass = min(
        best_by_mode.values(),
        key=lambda r: r["finite_validation"]["midpoint_validation_scores"]["weighted_rms_log_ratio"],
    )

    output = {
        "stage": "intermediate_target_hardware_convention_parameter_refinement",
        "strict_target_conclusion": STRICT,
        "strict_target_parameter_estimate_allowed": False,
        "searched_parameters": list(SEARCH_PARAMETERS),
        "hardware_modes": list(args.hardware_modes),
        "hardware_bypass_is_diagnostic_only": True,
        "fixed": {
            "hardware_bessel_order": HARDWARE_BESSEL_ORDER,
            "hardware_bessel_cutoff_Hz": HARDWARE_BESSEL_CUTOFF_HZ,
            "analysis_bessel_cutoff_Hz": cutoff,
            "rate_Hz": rate,
            "samples": sample,
            "T_bath_K": envelope["parameters"]["T_bath"]["nominal"],
        },
        "band_definitions_Hz": {
            "mid": [1_000.0, 10_000.0],
            "high": [10_000.0, 80_000.0],
            "edge": [80_000.0, float(args.edge_max_hz)],
        },
        "band_weights": BAND_WEIGHTS,
        "target_family": target_metadata,
        "selection_rule": (
            "best fit per lambda, then minimum deterministic weighted error to fixed "
            "lambda=0.5 midpoint; finite-record validation is reported afterward"
        ),
        "raw_experiment_preserved_for_reporting": True,
        "noise_residual_fit": False,
        "additive_noise_parameter_fit": False,
        "simulation_amplitude_rescale": False,
        "fit_count": len(fit_rows),
        "best_by_mode_and_lambda": by_mode_lambda,
        "best_by_hardware_mode": best_by_mode,
        "best_physical_hardware_mode": best_physical,
        "best_including_bypass_diagnostic": best_including_bypass,
    }
    dump(case / "intermediate_target_hardware_refinement.json", output)

    import matplotlib.pyplot as plt

    plot_frequency = np.logspace(np.log10(rate / sample), np.log10(rate / 2.0), 1200)
    exp_plot = log_interp(exp_freq[1:], exp_norm[1:], plot_frequency)
    support_f = np.asarray(target_metadata["support_frequency_Hz"], dtype=float)
    support_lower = np.asarray(target_metadata["support_lower_baseline_normalized_asd"], dtype=float)
    support_upper = np.asarray(target_metadata["support_smoothed_raw_normalized_asd"], dtype=float)
    support_mid = geometric_intermediate(support_lower, support_upper, 0.5)
    colors = {"phase": "tab:blue", "mag": "tab:green", "bypass": "tab:orange"}
    labels = {
        "phase": "100 kHz Bessel phase norm",
        "mag": "100 kHz Bessel mag norm (-3 dB @ 100 kHz)",
        "bypass": "100 kHz hardware Bessel bypass (diagnostic)",
    }
    curves = {
        mode: log_interp(data["frequency"][1:], data["normalized_asd"][1:], plot_frequency)
        for mode, data in finite.items()
    }

    plt.figure(figsize=(9.8, 5.9))
    plt.plot(plot_frequency, exp_plot, color="black", lw=1.5, label=f"Experiment ({len(exp_paths)} accepted records)")
    plt.plot(support_f, support_lower, color="gray", ls=":", lw=1.1, label="Lower baseline")
    plt.plot(support_f, support_upper, color="gray", ls="-.", lw=1.1, label="Smoothed raw trend")
    plt.plot(support_f, support_mid, color="gray", ls="--", lw=1.4, label="Geometric midpoint target")
    for mode, values in curves.items():
        lam = best_by_mode[mode]["lambda"]
        plt.plot(plot_frequency, values, lw=1.2, color=colors[mode], label=f"{labels[mode]} (lambda={lam:.2f})")
    plt.xscale("log"); plt.yscale("log"); plt.xlim(rate / sample, rate / 2.0)
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("Intermediate-target parameter refinement")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=7); plt.tight_layout()
    plt.savefig(case / "intermediate_target_hardware_refinement.png", dpi=180); plt.close()

    comparison_mask = (plot_frequency >= 1_000.0) & (plot_frequency <= args.edge_max_hz)
    plt.figure(figsize=(9.6, 5.2))
    for mode, values in curves.items():
        plt.plot(plot_frequency[comparison_mask], values[comparison_mask] / exp_plot[comparison_mask], lw=1.2, color=colors[mode], label=labels[mode])
    plt.axhline(1.0, color="black", lw=1.0)
    plt.xscale("log"); plt.yscale("log"); plt.xlim(1_000.0, args.edge_max_hz)
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Simulation / raw experiment")
    plt.title("Intermediate-target refinement — residual to raw experiment")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "intermediate_target_hardware_refinement_raw_ratio.png", dpi=180); plt.close()

    midpoint_plot = exp_plot.copy()
    high_plot = (plot_frequency >= 10_000.0) & (plot_frequency <= support_f[-1])
    midpoint_plot[high_plot] = np.exp(np.interp(np.log(plot_frequency[high_plot]), np.log(support_f), np.log(support_mid)))
    plt.figure(figsize=(9.6, 5.2))
    for mode, values in curves.items():
        plt.plot(plot_frequency[comparison_mask], values[comparison_mask] / midpoint_plot[comparison_mask], lw=1.2, color=colors[mode], label=labels[mode])
    plt.axhline(1.0, color="black", lw=1.0)
    plt.xscale("log"); plt.yscale("log"); plt.xlim(1_000.0, args.edge_max_hz)
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Simulation / geometric-midpoint target")
    plt.title("Intermediate-target refinement — residual to midpoint")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "intermediate_target_hardware_refinement_midpoint_ratio.png", dpi=180); plt.close()

    plt.figure(figsize=(7.4, 4.8))
    for mode, lambda_rows in by_mode_lambda.items():
        x = np.array([float(key) for key in lambda_rows], dtype=float)
        order = np.argsort(x)
        y = np.array([
            lambda_rows[f"{value:.6f}"]["midpoint_validation_scores"]["weighted_rms_log_ratio"]
            for value in x
        ], dtype=float)
        plt.plot(x[order], y[order], marker="o", lw=1.2, color=colors[mode], label=labels[mode])
    plt.axvline(0.5, color="black", ls="--", lw=1.0, label="midpoint lambda=0.5")
    plt.xlabel("Optimization-target lambda")
    plt.ylabel("Deterministic weighted RMS log-ratio to fixed midpoint")
    plt.title("Lambda scan: target level vs midpoint validation")
    plt.grid(True, alpha=0.25); plt.legend(fontsize=7); plt.tight_layout()
    plt.savefig(case / "intermediate_target_hardware_refinement_lambda_scan.png", dpi=180); plt.close()

    plt.figure(figsize=(9.0, 4.8))
    for lam in sorted(set(lambdas) | {0.5}):
        target = geometric_intermediate(support_lower, support_upper, lam)
        plt.plot(support_f, target, lw=1.0, label=f"lambda={lam:.2f}")
    plt.plot(support_f, support_lower, color="black", ls=":", lw=1.1, label="lower")
    plt.plot(support_f, support_upper, color="black", ls="-.", lw=1.1, label="smoothed raw")
    plt.xscale("log"); plt.yscale("log")
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Normalized ASD")
    plt.title("High-frequency intermediate target family")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=7, ncol=2); plt.tight_layout()
    plt.savefig(case / "intermediate_target_hardware_refinement_target_family.png", dpi=180); plt.close()

    print(json.dumps({
        "lambdas": lambdas,
        "fit_count": len(fit_rows),
        "records": records,
        "best_by_mode": {
            mode: {
                "lambda": row["lambda"],
                "seed_id": row["seed_id"],
                "deterministic_midpoint_score": row["midpoint_validation_scores"]["weighted_rms_log_ratio"],
                "finite_midpoint_score": row["finite_validation"]["midpoint_validation_scores"]["weighted_rms_log_ratio"],
                "finite_edge_bracket_position": row["finite_validation"]["bracket_position"]["edge"],
            }
            for mode, row in best_by_mode.items()
        },
    }, indent=2))


if __name__ == "__main__":
    main()

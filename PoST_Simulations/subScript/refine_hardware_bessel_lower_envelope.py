"""Refine TES noise parameters against the lower high-frequency baseline.

This is an exploratory continuation of the correlated multivariate search.
It is meant for the situation where a broad/narrow convex feature in the
measured high-frequency ASD may pull an ordinary least-squares fit upward.
The raw experimental ASD is never altered or hidden: it remains the reporting
reference and is plotted explicitly.  Only the optimization target above
10 kHz is replaced by a smooth asymmetric lower-baseline estimate.

Three hardware conventions can be refined with identical TES seeds:
- phase: historical 4th-order 100 kHz analog Bessel normalization,
- mag:   4th-order 100 kHz Bessel with -3 dB at 100 kHz,
- bypass: diagnostic-only unity hardware response.

The confirmed 10 kHz software Bessel, finite-record estimator, sample rate,
record length, and T_bath stay fixed.  No additive residual/noise floor or
post-hoc amplitude scale is fitted.
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

from explore_high_frequency_noise_parameters import (  # noqa: E402
    HARDWARE_BESSEL_ORDER,
    evaluation_grids,
    experimental_targets,
    pulse_gate,
    rms_log_ratio,
)
from explore_high_frequency_noise_parameters_adaptive import (  # noqa: E402
    SEARCH_PARAMETERS,
    parameters_from_unit,
    specs_from_envelope,
    unit_from_parameters,
)
from high_frequency_noise_comparison import experiment_asd, log_interp, normalized  # noqa: E402
from noise_measurement_model import (  # noqa: E402
    HARDWARE_BESSEL_CUTOFF_HZ,
    expected_post_analysis_asd,
    finite_record_post_analysis_asd,
    hardware_sampled_asd,
)

STRICT = "C — exact target physical case remains unidentified"
HARDWARE_MODES = ("phase", "mag", "bypass")


def dump(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def asymmetric_lower_baseline(
    log_values: np.ndarray,
    smoothness: float = 80_000.0,
    asymmetry: float = 0.08,
    iterations: int = 18,
) -> np.ndarray:
    """Return a smooth lower-trend baseline using asymmetric least squares.

    The input grid is assumed to be uniform in log frequency.  ``asymmetry``
    below 0.5 gives points above the baseline less leverage than points below
    it, so broad convex excesses do not drag the target upward.
    """
    y = np.asarray(log_values, dtype=float)
    if y.ndim != 1 or len(y) < 5 or np.any(~np.isfinite(y)):
        raise ValueError("log_values must be a finite one-dimensional array")
    if smoothness <= 0.0:
        raise ValueError("smoothness must be positive")
    if not 0.0 < asymmetry < 0.5:
        raise ValueError("asymmetry must be between 0 and 0.5")
    if iterations <= 0:
        raise ValueError("iterations must be positive")

    n = len(y)
    difference = sparse.diags(
        [np.ones(n - 2), -2.0 * np.ones(n - 2), np.ones(n - 2)],
        [0, 1, 2],
        shape=(n - 2, n),
        format="csc",
    )
    penalty = float(smoothness) * (difference.T @ difference)
    baseline = y.copy()
    weights = np.ones(n)
    for _ in range(int(iterations)):
        weight_matrix = sparse.diags(weights, 0, shape=(n, n), format="csc")
        baseline = spsolve(weight_matrix + penalty, weights * y)
        weights = np.where(y > baseline, float(asymmetry), 1.0 - float(asymmetry))
    return np.asarray(baseline, dtype=float)


def lower_envelope_targets(
    exp_freq: np.ndarray,
    exp_norm: np.ndarray,
    grids: dict[str, np.ndarray],
    smoothness: float,
    asymmetry: float,
    blend_end_hz: float,
    support_max_hz: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict]:
    """Build raw targets and a lower-baseline optimization target above 10 kHz."""
    raw = experimental_targets(exp_freq, exp_norm, grids)
    positive = np.asarray(exp_freq, dtype=float) > 0.0
    max_available = float(np.max(np.asarray(exp_freq, dtype=float)[positive]))
    support_hi = min(float(support_max_hz), max_available)
    if support_hi <= 100_000.0:
        support_hi = min(max_available, 125_000.0)
    if support_hi <= 10_000.0:
        raise ValueError("experimental spectrum does not extend above 10 kHz")

    support_frequency = np.logspace(np.log10(10_000.0), np.log10(support_hi), 500)
    support_raw = log_interp(exp_freq[positive], exp_norm[positive], support_frequency)
    support_log_baseline = asymmetric_lower_baseline(
        np.log(support_raw),
        smoothness=smoothness,
        asymmetry=asymmetry,
    )
    support_baseline = np.exp(support_log_baseline)

    target: dict[str, np.ndarray] = {}
    for name, grid in grids.items():
        grid = np.asarray(grid, dtype=float)
        raw_values = np.asarray(raw[name], dtype=float)
        values = raw_values.copy()
        high = grid >= 10_000.0
        if np.any(high):
            baseline = np.exp(
                np.interp(
                    np.log(grid[high]),
                    np.log(support_frequency),
                    support_log_baseline,
                )
            )
            if blend_end_hz <= 10_000.0:
                blend = np.ones(np.count_nonzero(high))
            else:
                blend = np.clip(
                    (np.log(grid[high]) - np.log(10_000.0))
                    / (np.log(float(blend_end_hz)) - np.log(10_000.0)),
                    0.0,
                    1.0,
                )
            values[high] = np.exp(
                (1.0 - blend) * np.log(raw_values[high])
                + blend * np.log(baseline)
            )
        target[name] = values

    metadata = {
        "method": "asymmetric_least_squares_lower_baseline_in_log_frequency_log_ASD",
        "high_frequency_start_Hz": 10_000.0,
        "blend_end_Hz": float(blend_end_hz),
        "support_max_Hz": float(support_hi),
        "smoothness": float(smoothness),
        "asymmetry": float(asymmetry),
        "support_frequency_Hz": support_frequency.tolist(),
        "support_raw_normalized_asd": support_raw.tolist(),
        "support_lower_baseline_normalized_asd": support_baseline.tolist(),
    }
    return raw, target, metadata


def mode_kwargs(mode: str) -> dict:
    if mode == "phase":
        return {"hardware_norm": "phase", "bypass_hardware": False}
    if mode == "mag":
        return {"hardware_norm": "mag", "bypass_hardware": False}
    if mode == "bypass":
        return {"hardware_norm": "phase", "bypass_hardware": True}
    raise ValueError(f"unknown hardware mode: {mode}")


def deterministic_normalized_curve(
    parameters: dict,
    frequency: np.ndarray,
    rate: float,
    cutoff: float,
    mode: str,
) -> np.ndarray:
    query = np.unique(np.concatenate((np.asarray([1_000.0]), np.asarray(frequency, dtype=float))))
    kwargs = mode_kwargs(mode)
    expected = expected_post_analysis_asd(
        parameters,
        query,
        rate_hz=rate,
        hardware_cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
        analysis_cutoff_hz=cutoff,
        hardware_norm=kwargs["hardware_norm"],
        bypass_hardware=kwargs["bypass_hardware"],
    )
    norm_1k = log_interp(query, expected, np.asarray([1_000.0]))[0]
    return log_interp(query, expected / norm_1k, np.asarray(frequency, dtype=float))


def score_curve(simulated: np.ndarray, target: np.ndarray) -> dict:
    rms, maximum = rms_log_ratio(simulated, target)
    return {"rms_log_ratio": rms, "max_abs_log_ratio": maximum}


def seed_parameters(search_doc: dict, pareto_seed_count: int) -> list[dict]:
    selected: dict[str, dict] = {}
    for name, row in search_doc.get("best_deterministic", {}).items():
        if "parameters" in row:
            selected[f"best_{name}"] = row["parameters"]
    front = [row for row in search_doc.get("pareto_front", []) if "parameters" in row]
    if front and pareto_seed_count > 0:
        indices = np.unique(
            np.rint(
                np.linspace(0, len(front) - 1, min(int(pareto_seed_count), len(front)))
            ).astype(int)
        )
        for index in indices:
            selected[f"pareto_{int(index):03d}"] = front[int(index)]["parameters"]
    return [{"seed_id": key, "parameters": dict(value)} for key, value in selected.items()]


def fit_one_seed(
    seed: dict,
    mode: str,
    specs: list[dict],
    template: dict,
    envelope: dict,
    constraints: dict,
    grid: np.ndarray,
    target: np.ndarray,
    rate: float,
    cutoff: float,
    max_nfev: int,
) -> dict:
    x0 = unit_from_parameters(seed["parameters"], specs)

    def residual(unit: np.ndarray) -> np.ndarray:
        parameters = parameters_from_unit(unit, specs, template, envelope)
        ok, _gate = pulse_gate(parameters, constraints)
        if not ok:
            return np.full(len(grid), 6.0)
        try:
            simulated = deterministic_normalized_curve(parameters, grid, rate, cutoff, mode)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError, OverflowError):
            return np.full(len(grid), 6.0)
        return np.log(simulated / target)

    result = optimize.least_squares(
        residual,
        x0,
        bounds=(np.zeros(len(specs)), np.ones(len(specs))),
        method="trf",
        jac="2-point",
        loss="soft_l1",
        f_scale=0.10,
        x_scale="jac",
        max_nfev=int(max_nfev),
    )
    parameters = parameters_from_unit(result.x, specs, template, envelope)
    gate_ok, gate = pulse_gate(parameters, constraints)
    simulated = deterministic_normalized_curve(parameters, grid, rate, cutoff, mode)
    return {
        "seed_id": seed["seed_id"],
        "hardware_mode": mode,
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
        "target_score": score_curve(simulated, target),
    }


def realize_mode(
    row: dict,
    mode: str,
    sample: int,
    rate: float,
    cutoff: float,
    records: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    kwargs = mode_kwargs(mode)
    pre_analysis = hardware_sampled_asd(
        row["parameters"],
        frequency,
        rate_hz=rate,
        cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
        order=HARDWARE_BESSEL_ORDER,
        norm=kwargs["hardware_norm"],
        bypass=kwargs["bypass_hardware"],
    )
    finite = finite_record_post_analysis_asd(
        pre_analysis,
        sample,
        rate,
        analysis_cutoff_hz=cutoff,
        records=records,
        seed=seed,
    )
    return frequency, normalized(finite, frequency)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument(
        "--search-json",
        type=Path,
        help="defaults to correlated_high_frequency_parameter_search.json in case-dir",
    )
    ap.add_argument(
        "--hardware-modes",
        nargs="+",
        choices=HARDWARE_MODES,
        default=list(HARDWARE_MODES),
    )
    ap.add_argument("--pareto-seeds", type=int, default=6)
    ap.add_argument("--max-nfev", type=int, default=100)
    ap.add_argument("--lower-smoothness", type=float, default=80_000.0)
    ap.add_argument("--lower-asymmetry", type=float, default=0.08)
    ap.add_argument("--lower-blend-end-hz", type=float, default=20_000.0)
    ap.add_argument("--lower-support-max-hz", type=float, default=150_000.0)
    ap.add_argument("--finite-records", type=int, default=0)
    ap.add_argument("--finite-seed", type=int, default=20260906)
    args = ap.parse_args()

    case = args.case_dir
    search_json = args.search_json or case / "correlated_high_frequency_parameter_search.json"
    search_doc = json.loads(search_json.read_text(encoding="utf-8"))
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
    exp_norm = normalized(exp_asd, exp_freq)
    grids = evaluation_grids()
    raw_targets, lower_targets, lower_metadata = lower_envelope_targets(
        exp_freq,
        exp_norm,
        grids,
        smoothness=args.lower_smoothness,
        asymmetry=args.lower_asymmetry,
        blend_end_hz=args.lower_blend_end_hz,
        support_max_hz=args.lower_support_max_hz,
    )

    template = dict(frozen[0]["parameters"])
    specs = specs_from_envelope(envelope, template, excess_max=5.0)
    seeds = seed_parameters(search_doc, args.pareto_seeds)
    if not seeds:
        raise RuntimeError("no parameter seeds found in correlated search JSON")

    fit_grid = grids["all"]
    fit_target = lower_targets["all"]
    fit_rows = []
    for mode in args.hardware_modes:
        for seed in seeds:
            fit_rows.append(
                fit_one_seed(
                    seed,
                    mode,
                    specs,
                    template,
                    envelope,
                    constraints,
                    fit_grid,
                    fit_target,
                    rate,
                    cutoff,
                    args.max_nfev,
                )
            )

    best_by_mode = {}
    for mode in args.hardware_modes:
        candidates = [row for row in fit_rows if row["hardware_mode"] == mode and row["pulse_gate_ok"]]
        if not candidates:
            continue
        best_by_mode[mode] = min(candidates, key=lambda row: row["target_score"]["rms_log_ratio"])
    if not best_by_mode:
        raise RuntimeError("no pulse-consistent refined solution")

    records = args.finite_records if args.finite_records > 0 else len(exp_paths)
    finite = {}
    for mode, row in best_by_mode.items():
        frequency, values = realize_mode(
            row,
            mode,
            sample,
            rate,
            cutoff,
            records,
            args.finite_seed,
        )
        raw_all = log_interp(exp_freq[1:], exp_norm[1:], grids["all"])
        finite_all = log_interp(frequency[1:], values[1:], grids["all"])
        row["finite_raw_experiment_score"] = score_curve(finite_all, raw_all)
        row["finite_lower_target_score"] = score_curve(finite_all, lower_targets["all"])
        finite[mode] = {"frequency": frequency, "normalized_asd": values}

    output = {
        "stage": "lower_envelope_hardware_convention_parameter_refinement",
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
        "optimization_target": lower_metadata,
        "raw_experiment_preserved_for_reporting": True,
        "noise_residual_fit": False,
        "additive_noise_parameter_fit": False,
        "simulation_amplitude_rescale": False,
        "seed_count": len(seeds),
        "fit_count": len(fit_rows),
        "best_by_hardware_mode": best_by_mode,
    }
    dump(case / "lower_envelope_hardware_refinement.json", output)

    import matplotlib.pyplot as plt

    plot_frequency = np.logspace(np.log10(rate / sample), np.log10(rate / 2.0), 1200)
    exp_plot = log_interp(exp_freq[1:], exp_norm[1:], plot_frequency)
    support_f = np.asarray(lower_metadata["support_frequency_Hz"], dtype=float)
    support_lower = np.asarray(lower_metadata["support_lower_baseline_normalized_asd"], dtype=float)
    lower_plot = np.full_like(plot_frequency, np.nan)
    high_mask = (plot_frequency >= 10_000.0) & (plot_frequency <= support_f[-1])
    lower_plot[high_mask] = np.exp(
        np.interp(np.log(plot_frequency[high_mask]), np.log(support_f), np.log(support_lower))
    )

    colors = {"phase": "tab:blue", "mag": "tab:green", "bypass": "tab:orange"}
    labels = {
        "phase": "100 kHz Bessel phase norm",
        "mag": "100 kHz Bessel mag norm (-3 dB @ 100 kHz)",
        "bypass": "100 kHz hardware Bessel bypass (diagnostic)",
    }

    plt.figure(figsize=(9.5, 5.8))
    plt.plot(plot_frequency, exp_plot, color="black", lw=1.5, label=f"Experiment ({len(exp_paths)} accepted records)")
    plt.plot(plot_frequency, lower_plot, color="gray", ls="--", lw=1.2, label="Lower-baseline optimization target")
    plot_curves = {}
    for mode, data in finite.items():
        values = log_interp(data["frequency"][1:], data["normalized_asd"][1:], plot_frequency)
        plot_curves[mode] = values
        plt.plot(plot_frequency, values, lw=1.2, color=colors[mode], label=labels[mode])
    plt.xscale("log"); plt.yscale("log"); plt.xlim(rate / sample, rate / 2.0)
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("Lower-baseline parameter refinement by hardware convention")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "lower_envelope_hardware_refinement.png", dpi=180); plt.close()

    comparison_mask = (plot_frequency >= 1_000.0) & (plot_frequency <= 100_000.0)
    plt.figure(figsize=(9.5, 5.2))
    for mode, values in plot_curves.items():
        plt.plot(plot_frequency[comparison_mask], values[comparison_mask] / exp_plot[comparison_mask], lw=1.2, color=colors[mode], label=labels[mode])
    plt.axhline(1.0, color="black", lw=1.0)
    plt.xscale("log"); plt.yscale("log"); plt.xlim(1_000.0, 100_000.0)
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Simulation / raw experiment")
    plt.title("Refined residual ratio to raw experiment")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "lower_envelope_hardware_refinement_raw_ratio.png", dpi=180); plt.close()

    target_plot = exp_plot.copy()
    target_high = comparison_mask & (plot_frequency >= 10_000.0)
    target_plot[target_high] = np.exp(
        np.interp(np.log(plot_frequency[target_high]), np.log(support_f), np.log(support_lower))
    )
    plt.figure(figsize=(9.5, 5.2))
    for mode, values in plot_curves.items():
        plt.plot(plot_frequency[comparison_mask], values[comparison_mask] / target_plot[comparison_mask], lw=1.2, color=colors[mode], label=labels[mode])
    plt.axhline(1.0, color="black", lw=1.0)
    plt.xscale("log"); plt.yscale("log"); plt.xlim(1_000.0, 100_000.0)
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Simulation / optimization target")
    plt.title("Refined residual ratio to lower-baseline target")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "lower_envelope_hardware_refinement_target_ratio.png", dpi=180); plt.close()

    plt.figure(figsize=(9.0, 4.8))
    support_raw = np.asarray(lower_metadata["support_raw_normalized_asd"], dtype=float)
    plt.plot(support_f, support_raw, color="black", lw=1.1, label="Raw experiment trend")
    plt.plot(support_f, support_lower, color="gray", ls="--", lw=1.4, label="Asymmetric lower baseline")
    plt.xscale("log"); plt.yscale("log")
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Normalized ASD")
    plt.title("High-frequency target construction")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "lower_envelope_hardware_refinement_target.png", dpi=180); plt.close()

    print(json.dumps({
        "hardware_modes": list(best_by_mode),
        "seed_count": len(seeds),
        "records": records,
        "best": {
            mode: {
                "seed_id": row["seed_id"],
                "target_score": row["target_score"],
                "finite_raw_experiment_score": row["finite_raw_experiment_score"],
                "finite_lower_target_score": row["finite_lower_target_score"],
            }
            for mode, row in best_by_mode.items()
        },
    }, indent=2))


if __name__ == "__main__":
    main()

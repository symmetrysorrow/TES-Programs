"""Correlated multivariate TES noise search.

Exploratory/noise-guided only. The confirmed 4th-order 100 kHz hardware
Bessel, 10 kHz software Bessel, sample rate, record length, and T_bath remain
fixed. The same 13 TES model parameters as the adaptive search are varied.
This continuation adds covariance-adaptive proposals, robust residual-vector
least-squares, and pair-interaction maps.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import optimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations"))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))

from explore_high_frequency_noise_parameters import (  # noqa: E402
    HARDWARE_BESSEL_ORDER, SCORE_FIELDS, best_by_objective, evaluation_grids,
    experimental_targets, finite_best_by_objective, pareto_front, public_row,
    pulse_gate, realize_finite, score_parameters,
)
from explore_high_frequency_noise_parameters_adaptive import (  # noqa: E402
    SEARCH_PARAMETERS, boundary_hits, correlated_r_sh, dump, elite_pool,
    expand_bounds, latin_hypercube, parameters_from_unit, polish,
    specs_from_envelope, unit_from_parameters,
)
from high_frequency_noise_comparison import experiment_asd, log_interp, normalized  # noqa: E402
from noise_measurement_model import (  # noqa: E402
    HARDWARE_BESSEL_CUTOFF_HZ, expected_post_analysis_asd,
)

STRICT = "C — exact target physical case remains unidentified"
INTERACTION_PAIRS = (
    ("L", "R_l"), ("L", "beta"), ("L", "alpha"),
    ("alpha", "C_tes"), ("C_tes", "G_tes-bath"),
)


def row_id(row):
    return str(row.get("trial_id", row.get("scenario_id", "")))


def evaluate(unit, trial_id, stage, specs, template, envelope, constraints,
             grids, targets, rate, cutoff, index):
    p = parameters_from_unit(unit, specs, template, envelope)
    ok, gate = pulse_gate(p, constraints)
    if not ok:
        return None, "unstable" if not gate.get("stable", False) else "pulse_gate"
    try:
        scores = score_parameters(p, grids, targets, rate, cutoff)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError, OverflowError):
        return None, "model_failure"
    return {
        "trial_id": trial_id, "search_stage": stage, "evaluation_index": index,
        **scores, "parameters": p, "pulse_consistency": gate,
        "strict_target_allowed": False,
    }, None


def diverse_anchors(rows, count):
    selected = {row_id(r): r for r in best_by_objective(rows).values()}
    front = pareto_front(rows)
    if front and count > 0:
        indices = np.unique(np.rint(np.linspace(
            0, len(front) - 1, min(count, len(front))
        )).astype(int))
        for i in indices:
            selected[row_id(front[i])] = front[i]
    return list(selected.values())


def rank_weights(rows):
    n = len(rows)
    rank_sum = np.zeros(n)
    for field in SCORE_FIELDS.values():
        order = np.argsort([float(r[field]) for r in rows], kind="stable")
        ranks = np.empty(n)
        ranks[order] = np.arange(n)
        rank_sum += ranks / max(n - 1, 1)
    weights = np.exp(-2.0 * rank_sum)
    return weights / weights.sum()


def learned_covariance(rows, specs, floor_sigma=0.02):
    d = len(specs)
    if len(rows) < 2:
        return np.eye(d) * floor_sigma**2
    x = np.asarray([unit_from_parameters(r["parameters"], specs) for r in rows])
    w = rank_weights(rows)
    mean = np.sum(x * w[:, None], axis=0)
    centered = x - mean
    cov = (centered * w[:, None]).T @ centered
    cov = 0.5 * (cov + cov.T) + np.eye(d) * floor_sigma**2
    values, vectors = np.linalg.eigh(cov)
    values = np.maximum(values, floor_sigma**2)
    return (vectors * values) @ vectors.T


def covariance_correlations(cov, specs, limit=12):
    sigma = np.sqrt(np.maximum(np.diag(cov), 1e-30))
    corr = cov / np.outer(sigma, sigma)
    out = []
    for i in range(len(specs)):
        for j in range(i + 1, len(specs)):
            out.append({
                "parameter_1": specs[i]["name"], "parameter_2": specs[j]["name"],
                "correlation": float(corr[i, j]),
            })
    out.sort(key=lambda x: abs(x["correlation"]), reverse=True)
    return out[:limit]


def covariance_refine(rows, specs, template, envelope, constraints, grids, targets,
                      rate, cutoff, rng, trials, round_index, start_index,
                      elite_count, scale):
    geometry = elite_pool(rows, max(20, 2 * elite_count))
    anchors = diverse_anchors(rows, max(8, elite_count))
    cov = learned_covariance(geometry, specs)
    round_scale = scale * 0.62 ** max(round_index - 1, 0)
    iso = max(0.008, 0.22 * round_scale)
    proposal_cov = cov * round_scale**2 + np.eye(len(specs)) * iso**2
    accepted, rejected, index = [], {"unstable": 0, "pulse_gate": 0, "model_failure": 0}, start_index
    for i in range(trials):
        center = unit_from_parameters(anchors[int(rng.integers(len(anchors)))]["parameters"], specs)
        delta = rng.multivariate_normal(np.zeros(len(specs)), proposal_cov, check_valid="ignore")
        index += 1
        row, reason = evaluate(
            np.clip(center + delta, 0, 1), f"cov_{round_index}_{i:06d}",
            f"covariance_adaptive_refinement_{round_index}", specs, template,
            envelope, constraints, grids, targets, rate, cutoff, index,
        )
        if row is None:
            rejected[reason] += 1
        else:
            accepted.append(row)
    diagnostics = {
        "round": round_index, "anchor_count": len(anchors),
        "geometry_elite_count": len(geometry), "proposal_scale": float(round_scale),
        "proposal_trace": float(np.trace(proposal_cov)),
        "top_correlations": covariance_correlations(proposal_cov, specs),
    }
    return accepted, index, rejected, diagnostics


def deterministic_residual(parameters, grid, target, rate, cutoff):
    f = np.unique(np.concatenate(([1000.0], np.asarray(grid, dtype=float))))
    p = dict(parameters)
    p["hardware_bessel_order"] = HARDWARE_BESSEL_ORDER
    expected = expected_post_analysis_asd(
        p, f, rate_hz=rate, hardware_cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
        analysis_cutoff_hz=cutoff,
    )
    norm = log_interp(f, expected, np.asarray([1000.0]))[0]
    simulated = log_interp(f, expected / norm, grid)
    residual = np.log(simulated / target)
    if np.any(~np.isfinite(residual)):
        raise FloatingPointError("non-finite residual")
    return residual


def residual_from_unit(unit, objective, specs, template, envelope, constraints,
                       grids, targets, rate, cutoff):
    p = parameters_from_unit(unit, specs, template, envelope)
    grid = grids[objective]
    ok, _ = pulse_gate(p, constraints)
    if not ok:
        return np.full(len(grid), 6.0)
    try:
        return deterministic_residual(p, grid, targets[objective], rate, cutoff)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError, OverflowError):
        return np.full(len(grid), 6.0)


def least_squares_plan(rows, pareto_count):
    plan = [(obj, row) for obj, row in best_by_objective(rows).items()]
    front = pareto_front(rows)
    if front and pareto_count > 0:
        indices = np.unique(np.rint(np.linspace(
            0, len(front) - 1, min(pareto_count, len(front))
        )).astype(int))
        plan += [("all", front[i]) for i in indices]
    unique = {(obj, row_id(row)): (obj, row) for obj, row in plan}
    return list(unique.values())


def residual_least_squares(rows, specs, template, envelope, constraints, grids,
                           targets, rate, cutoff, pareto_seeds, max_nfev, index):
    added, diagnostics = [], []
    for seed_index, (objective, seed) in enumerate(least_squares_plan(rows, pareto_seeds)):
        x0 = unit_from_parameters(seed["parameters"], specs)

        def fun(x):
            return residual_from_unit(
                x, objective, specs, template, envelope, constraints,
                grids, targets, rate, cutoff,
            )

        result = optimize.least_squares(
            fun, x0, bounds=(np.zeros(len(specs)), np.ones(len(specs))),
            method="trf", jac="2-point", loss="soft_l1", f_scale=0.12,
            x_scale="jac", max_nfev=max_nfev,
        )
        index += 1
        row, reason = evaluate(
            result.x, f"least_squares_{objective}_{seed_index:03d}",
            f"residual_vector_least_squares_{objective}", specs, template,
            envelope, constraints, grids, targets, rate, cutoff, index,
        )
        record = {
            "objective": objective, "seed_trial_id": row_id(seed),
            "success": bool(result.success), "status": int(result.status),
            "cost": float(result.cost), "optimality": float(result.optimality),
            "nfev": int(result.nfev), "accepted_final_point": row is not None,
            "rejection_reason": reason,
        }
        diagnostics.append(record)
        if row is not None:
            row["least_squares"] = record
            added.append(row)
    return added, index, diagnostics


def finite_pool(rows, n):
    selected = {}
    for field in SCORE_FIELDS.values():
        ordered = sorted(rows, key=lambda r, f=field: (r[f], row_id(r)))
        for row in ordered[:n]:
            selected[row_id(row)] = row
    for row in diverse_anchors(rows, max(5, n)):
        selected[row_id(row)] = row
    return list(selected.values())


def derived_coordinates(p):
    r_eff = float(p["R_l"]) + float(p["R"]) * (1.0 + float(p["beta"]))
    tau_el = float(p["L"]) / r_eff
    g_eff = 1.0 / (1.0 / float(p["G_abs-tes"]) + 1.0 / (2.0 * float(p["G_abs-abs"])))
    tc, tb, n, r, g = map(float, (p["T_c"], p["T_bath"], p["n"], p["R"], p["G_tes-bath"]))
    current = math.sqrt(g * tc * (1.0 - (tb / tc) ** n) / (n * r))
    loop_gain = float(p["alpha"]) * current**2 * r / (g * tc)
    denom = (1.0 - loop_gain) * g
    tau_i = None if np.isclose(denom, 0.0) else float(p["C_tes"]) / denom
    return {
        "R_eff_ohm": r_eff, "tau_el_s": tau_el,
        "electrical_pole_Hz": 1.0 / (2.0 * np.pi * tau_el),
        "loop_gain": loop_gain, "tau_i_s": tau_i,
        "G_eff_W_per_K": g_eff, "current_A": current,
    }


def safe_matrix(a):
    return [[float(v) if np.isfinite(v) else None for v in row] for row in np.asarray(a)]


def interaction_scan(base, pair, specs, envelope, constraints, grids, targets,
                     rate, cutoff, points, half_width):
    dims = {name: i for i, name in enumerate(SEARCH_PARAMETERS)}
    ix, iy = dims[pair[0]], dims[pair[1]]
    base_unit = unit_from_parameters(base, specs)
    ux = np.linspace(max(0, base_unit[ix] - half_width), min(1, base_unit[ix] + half_width), points)
    uy = np.linspace(max(0, base_unit[iy] - half_width), min(1, base_unit[iy] + half_width), points)
    x_values, y_values = [], []
    for u in ux:
        v = base_unit.copy(); v[ix] = u
        x_values.append(parameters_from_unit(v, specs, base, envelope)[pair[0]])
    for u in uy:
        v = base_unit.copy(); v[iy] = u
        y_values.append(parameters_from_unit(v, specs, base, envelope)[pair[1]])
    scores = np.full((points, points), np.nan)
    high = np.full_like(scores, np.nan)
    for jy, vy in enumerate(uy):
        for jx, vx in enumerate(ux):
            v = base_unit.copy(); v[ix] = vx; v[iy] = vy
            p = parameters_from_unit(v, specs, base, envelope)
            if "R" in pair:
                p["R_SH"] = correlated_r_sh(envelope, p["R"])
            ok, _ = pulse_gate(p, constraints)
            if not ok:
                continue
            try:
                s = score_parameters(p, grids, targets, rate, cutoff)
            except Exception:
                continue
            scores[jy, jx] = s[SCORE_FIELDS["all"]]
            high[jy, jx] = s[SCORE_FIELDS["high"]]
    return {
        "pair": list(pair), "x_values": [float(x) for x in x_values],
        "y_values": [float(y) for y in y_values],
        "all_score": safe_matrix(scores), "high_score": safe_matrix(high),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--trials", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=20260906)
    ap.add_argument("--refine-rounds", type=int, default=4)
    ap.add_argument("--refine-trials", type=int, default=1000)
    ap.add_argument("--elite-count", type=int, default=20)
    ap.add_argument("--covariance-scale", type=float, default=1.0)
    ap.add_argument("--excess-johnson-max", type=float, default=5.0)
    ap.add_argument("--least-squares-pareto-seeds", type=int, default=8)
    ap.add_argument("--least-squares-max-nfev", type=int, default=120)
    ap.add_argument("--interaction-points", type=int, default=11)
    ap.add_argument("--interaction-half-width", type=float, default=0.28)
    ap.add_argument("--finite-top", type=int, default=4)
    ap.add_argument("--finite-records", type=int, default=0)
    ap.add_argument("--finite-seed", type=int, default=20260906)
    args = ap.parse_args()

    case = args.case_dir
    scenarios = json.loads((case / "proxy_scenarios.json").read_text(encoding="utf-8"))
    envelope = json.loads((case / "proxy_parameter_envelope.json").read_text(encoding="utf-8"))
    constraints = json.loads((case / "pulse_combination_constraints.json").read_text(encoding="utf-8"))
    if scenarios.get("freeze_status") != "frozen" or envelope.get("freeze_status") != "frozen":
        raise RuntimeError("correlated search requires frozen Stage-A inputs")
    frozen = scenarios["pulse_consistent_scenarios"]

    exp_f, exp_asd, exp_paths, acq = experiment_asd(args.target_root)
    rate, sample = float(acq["rate_Hz"]), int(acq["samples"])
    cutoff = float(acq["analysis_bessel_cutoff_Hz"])
    if rate / 2.0 < 100_000.0:
        raise RuntimeError("Nyquist below 100 kHz")
    exp_norm = normalized(exp_asd, exp_f)
    grids = evaluation_grids()
    targets = experimental_targets(exp_f, exp_norm, grids)

    template = dict(frozen[0]["parameters"])
    specs = specs_from_envelope(envelope, template, args.excess_johnson_max)
    frozen_rows = []
    for scenario in frozen:
        p = dict(scenario["parameters"])
        p["hardware_bessel_order"] = HARDWARE_BESSEL_ORDER
        p["excess_johnson_M"] = float(p.get("excess_johnson_M", 0.0))
        frozen_rows.append({
            "scenario_id": scenario["scenario_id"],
            **score_parameters(p, grids, targets, rate, cutoff),
            "parameters": p,
        })
    frozen_best = best_by_objective(frozen_rows)

    rng = np.random.default_rng(args.seed)
    rows, rejected, index = [], {"unstable": 0, "pulse_gate": 0, "model_failure": 0}, 0

    # Preserve the best frozen geometries as explicit starting islands.
    seen = set()
    for frozen_row in frozen_best.values():
        sid = frozen_row["scenario_id"]
        if sid in seen:
            continue
        seen.add(sid)
        index += 1
        row, reason = evaluate(
            unit_from_parameters(frozen_row["parameters"], specs),
            f"seed_{sid}", "frozen_geometry_seed", specs, template, envelope,
            constraints, grids, targets, rate, cutoff, index,
        )
        if row is None:
            rejected[reason] += 1
        else:
            rows.append(row)

    for trial, vector in enumerate(latin_hypercube(args.trials, len(specs), rng)):
        index += 1
        row, reason = evaluate(
            vector, f"global_{trial:06d}", "global_latin_hypercube", specs,
            template, envelope, constraints, grids, targets, rate, cutoff, index,
        )
        if row is None:
            rejected[reason] += 1
        else:
            rows.append(row)
    if not rows:
        raise RuntimeError("no accepted exploratory candidates")

    expansion_events, covariance_history = [], []
    for round_index in range(1, args.refine_rounds + 1):
        events = expand_bounds(specs, elite_pool(rows, args.elite_count))
        for event in events:
            event["round"] = round_index
        expansion_events.extend(events)
        added, index, local_rejected, diag = covariance_refine(
            rows, specs, template, envelope, constraints, grids, targets, rate,
            cutoff, rng, args.refine_trials, round_index, index,
            args.elite_count, args.covariance_scale,
        )
        rows.extend(added); covariance_history.append(diag)
        for key, value in local_rejected.items():
            rejected[key] += value

    ls_rows, index, ls_diagnostics = residual_least_squares(
        rows, specs, template, envelope, constraints, grids, targets, rate,
        cutoff, args.least_squares_pareto_seeds, args.least_squares_max_nfev, index,
    )
    rows.extend(ls_rows)
    for objective, seed in list(best_by_objective(rows).items()):
        added, _, index = polish(
            seed, objective, specs, template, envelope, constraints,
            grids, targets, rate, cutoff, index,
        )
        rows.extend(added)

    best, pareto = best_by_objective(rows), pareto_front(rows)
    interactions = {
        "__".join(pair): interaction_scan(
            best["all"]["parameters"], pair, specs, envelope, constraints,
            grids, targets, rate, cutoff, args.interaction_points,
            args.interaction_half_width,
        )
        for pair in INTERACTION_PAIRS
    }

    record_count = args.finite_records if args.finite_records > 0 else len(exp_paths)
    finite_rows = [
        realize_finite(r, sample, rate, cutoff, record_count, args.finite_seed, grids, targets)
        for r in finite_pool(rows, args.finite_top)
    ]
    finite_best = finite_best_by_objective(finite_rows)
    frozen_unique = {row_id(r): r for r in frozen_best.values()}
    frozen_finite = [
        realize_finite(r, sample, rate, cutoff, record_count, args.finite_seed, grids, targets)
        for r in frozen_unique.values()
    ]
    frozen_finite_best = finite_best_by_objective(frozen_finite)

    dump(case / "correlated_high_frequency_parameter_search.json", {
        "stage": "correlated_multivariate_all_parameter_noise_guided_search",
        "strict_target_conclusion": STRICT,
        "strict_target_parameter_estimate_allowed": False,
        "searched_parameters": list(SEARCH_PARAMETERS),
        "searched_parameter_count": len(SEARCH_PARAMETERS),
        "fixed": {
            "T_bath_K": envelope["parameters"]["T_bath"]["nominal"],
            "hardware_bessel_order": HARDWARE_BESSEL_ORDER,
            "hardware_bessel_cutoff_Hz": 100_000.0,
            "analysis_bessel_cutoff_Hz": cutoff, "rate_Hz": rate,
            "samples": sample, "F_LINK": 0.9,
        },
        "method": {
            "global": "Latin hypercube", "global_trials": args.trials,
            "multivariate_local": "Pareto-seeded covariance-adaptive Gaussian",
            "refine_rounds": args.refine_rounds,
            "refine_trials_per_round": args.refine_trials,
            "residual_local_optimizer": "scipy.optimize.least_squares soft_l1",
            "least_squares_pareto_seeds": args.least_squares_pareto_seeds,
            "interaction_pairs": [list(p) for p in INTERACTION_PAIRS],
        },
        "noise_residual_fit": False, "additive_noise_parameter_fit": False,
        "simulation_amplitude_rescale": False,
        "evaluation_count": index, "accepted_count": len(rows),
        "rejected": rejected, "bound_expansions": expansion_events,
        "covariance_history": covariance_history,
        "least_squares_diagnostics": ls_diagnostics,
        "interaction_scans": interactions,
        "best_deterministic": {
            key: {**public_row(row), "boundary_hits": boundary_hits(row, specs),
                  "derived_coordinates": derived_coordinates(row["parameters"])}
            for key, row in best.items()
        },
        "best_finite": {key: public_row(row) for key, row in finite_best.items()},
        "pareto_front": [public_row(row) for row in pareto[:300]],
    })

    import matplotlib.pyplot as plt
    pf = np.logspace(np.log10(rate / sample), np.log10(rate / 2.0), 1200)
    ep = log_interp(exp_f[1:], exp_norm[1:], pf)

    def curve(row):
        return log_interp(row["_finite_frequency_Hz"][1:], row["_finite_normalized_asd"][1:], pf)

    curves = {
        "Frozen balanced": curve(frozen_finite_best["all"]),
        "Mid best": curve(finite_best["mid"]),
        "High best": curve(finite_best["high"]),
        "Balanced best": curve(finite_best["all"]),
    }
    plt.figure(figsize=(9, 5.5))
    plt.plot(pf, ep, color="black", lw=1.5, label=f"Experiment ({len(exp_paths)} accepted records)")
    for label, values in curves.items():
        plt.plot(pf, values, lw=1.1, label=label)
    plt.xscale("log"); plt.yscale("log"); plt.xlim(rate / sample, rate / 2.0)
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("Correlated multivariate TES noise search")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "correlated_high_frequency_parameter_search.png", dpi=180); plt.close()

    mask = (pf >= 1_000) & (pf <= 100_000)
    plt.figure(figsize=(9, 5))
    for label, values in curves.items():
        plt.plot(pf[mask], values[mask] / ep[mask], lw=1.1, label=label)
    plt.axhline(1, color="black", lw=1); plt.xscale("log"); plt.yscale("log")
    plt.xlim(1_000, 100_000); plt.xlabel("Frequency [Hz]"); plt.ylabel("Simulation / experiment")
    plt.title("Correlated-search residual ratio")
    plt.grid(True, which="both", alpha=0.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "correlated_high_frequency_parameter_search_ratio.png", dpi=180); plt.close()

    plt.figure(figsize=(6.2, 5.2))
    plt.scatter([r[SCORE_FIELDS["mid"]] for r in rows], [r[SCORE_FIELDS["high"]] for r in rows], s=8, alpha=0.2)
    plt.plot([r[SCORE_FIELDS["mid"]] for r in pareto], [r[SCORE_FIELDS["high"]] for r in pareto], marker="o", ms=3, lw=1)
    plt.xlabel("1–10 kHz RMS log-ratio"); plt.ylabel("10–100 kHz RMS log-ratio")
    plt.title("Correlated-search Pareto front"); plt.grid(True, alpha=0.25); plt.tight_layout()
    plt.savefig(case / "correlated_high_frequency_parameter_pareto.png", dpi=180); plt.close()

    columns = 3
    fig, axes = plt.subplots(int(np.ceil(len(INTERACTION_PAIRS) / columns)), columns, figsize=(12.2, 7.6), squeeze=False)
    spec_by_name = {s["name"]: s for s in specs}
    for ax, pair in zip(axes.ravel(), INTERACTION_PAIRS):
        scan = interactions["__".join(pair)]
        x = np.asarray(scan["x_values"]); y = np.asarray(scan["y_values"])
        z = np.asarray([[np.nan if v is None else v for v in row] for row in scan["all_score"]])
        mesh = ax.pcolormesh(x, y, z, shading="auto")
        ax.set_xlabel(pair[0]); ax.set_ylabel(pair[1]); ax.set_title(f"{pair[0]} × {pair[1]}: 1–100 kHz RMS")
        if spec_by_name[pair[0]]["kind"] == "log": ax.set_xscale("log")
        if spec_by_name[pair[1]]["kind"] == "log": ax.set_yscale("log")
        fig.colorbar(mesh, ax=ax, label="RMS log-ratio")
    for ax in axes.ravel()[len(INTERACTION_PAIRS):]:
        ax.axis("off")
    fig.suptitle("Two-parameter interactions around final balanced solution")
    fig.tight_layout()
    fig.savefig(case / "correlated_high_frequency_parameter_interactions.png", dpi=180)
    plt.close(fig)

    print(json.dumps({
        "evaluations": index, "accepted": len(rows), "rejected": rejected,
        "least_squares_runs": len(ls_diagnostics),
        "best": {
            key: {"id": row["trial_id"], "score": row[SCORE_FIELDS[key]],
                  "derived_coordinates": derived_coordinates(row["parameters"])}
            for key, row in best.items()
        },
    }, indent=2))


if __name__ == "__main__":
    main()

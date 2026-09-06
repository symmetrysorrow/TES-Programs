"""Adaptive all-parameter TES noise-shape search.

Noise-guided and exploratory only. The confirmed 4th-order 100 kHz hardware
Bessel, target 10 kHz software Bessel, sample rate, record length, and T_bath
are fixed. All other parameters consumed by tes_noise_model are searched:
T_c, R, R_l, alpha, beta, L, n, C_tes, C_abs, G_tes-bath, G_abs-tes,
G_abs-abs, and excess_johnson_M.

The algorithm uses one-at-a-time sensitivity, Latin-hypercube global coverage,
elite local refinement, coordinate polish, multi-band/Pareto ranking, then
finite-record re-evaluation of only the best candidates.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations"))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))

from explore_high_frequency_noise_parameters import (  # noqa: E402
    BANDS_HZ,
    HARDWARE_BESSEL_ORDER,
    SCORE_FIELDS,
    best_by_objective,
    evaluation_grids,
    experimental_targets,
    finite_best_by_objective,
    pareto_front,
    public_row,
    pulse_gate,
    realize_finite,
    score_parameters,
)
from high_frequency_noise_comparison import experiment_asd, log_interp, normalized  # noqa: E402

STRICT = "C — exact target physical case remains unidentified"
SEARCH_PARAMETERS = (
    "T_c", "R", "R_l", "alpha", "beta", "L", "n", "C_tes", "C_abs",
    "G_tes-bath", "G_abs-tes", "G_abs-abs", "excess_johnson_M",
)


def dump(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def latin_hypercube(n: int, d: int, rng: np.random.Generator) -> np.ndarray:
    if n <= 0 or d <= 0:
        raise ValueError("n and d must be positive")
    out = np.empty((n, d))
    for j in range(d):
        column = (np.arange(n) + rng.random(n)) / n
        rng.shuffle(column)
        out[:, j] = column
    return out


def specs_from_envelope(envelope: dict, template: dict, excess_max: float) -> list[dict]:
    p = envelope["parameters"]
    ref = envelope["sensitivity_reference"]
    specs = [
        {"name": "T_c", "kind": "linear", "low": p["T_c"]["range"][0],
         "high": p["T_c"]["range"][1], "expand": False,
         "basis": "frozen RT proxy envelope"},
        {"name": "R", "kind": "linear", "low": p["R"]["range"][0],
         "high": p["R"]["range"][1], "expand": False,
         "basis": "frozen IV-derived proxy envelope"},
    ]
    factors = {
        "R_l": (0.25, 4.0), "alpha": (0.1, 10.0), "beta": (0.1, 10.0),
        "n": (0.5, 2.0), "C_tes": (0.1, 10.0), "C_abs": (0.1, 10.0),
        "G_tes-bath": (0.1, 10.0), "G_abs-tes": (0.1, 10.0),
        "G_abs-abs": (0.1, 10.0),
    }
    for name, (lo, hi) in factors.items():
        nominal = float(ref[name])
        low = nominal * lo
        if name == "n":
            low = max(1.01, low)
        specs.append({"name": name, "kind": "log", "low": low,
                      "high": nominal * hi, "expand": True,
                      "reference": nominal,
                      "basis": f"exploratory {lo:g}x--{hi:g}x simulation reference"})
    r_eff = float(p["R"]["nominal"]) * (1.0 + float(ref["beta"])) + float(ref["R_l"])
    l0 = float(ref["L"])
    l_lo = min(l0 / 10.0, r_eff / (2 * np.pi * 300_000.0) / 3.0)
    l_hi = max(l0 * 10.0, r_eff / (2 * np.pi * 1_000.0) * 3.0)
    specs.append({"name": "L", "kind": "log", "low": l_lo, "high": l_hi,
                  "expand": True, "reference": l0,
                  "basis": "electrical-pole-informed exploratory range"})
    specs.append({"name": "excess_johnson_M", "kind": "linear", "low": 0.0,
                  "high": float(excess_max), "expand": True,
                  "reference": float(template.get("excess_johnson_M", 0.0)),
                  "basis": "existing optional tes_noise_model source parameter"})
    by_name = {s["name"]: s for s in specs}
    return [by_name[name] for name in SEARCH_PARAMETERS]


def decode(spec: dict, u: float) -> float:
    u = float(np.clip(u, 0, 1))
    lo, hi = float(spec["low"]), float(spec["high"])
    if spec["kind"] == "log":
        return float(np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))))
    return float(lo + u * (hi - lo))


def encode(spec: dict, value: float) -> float:
    lo, hi, value = float(spec["low"]), float(spec["high"]), float(value)
    if spec["kind"] == "log":
        u = (np.log(max(value, 1e-300)) - np.log(lo)) / (np.log(hi) - np.log(lo))
    else:
        u = (value - lo) / (hi - lo)
    return float(np.clip(u, 0, 1))


def correlated_r_sh(envelope: dict, r: float) -> float:
    r0, r1 = map(float, envelope["parameters"]["R"]["range"])
    s0, s1 = map(float, envelope["parameters"]["R_SH"]["range"])
    u = 0.5 if np.isclose(r0, r1) else np.clip((r - r0) / (r1 - r0), 0, 1)
    return float(s0 + u * (s1 - s0))


def parameters_from_unit(unit: np.ndarray, specs: list[dict], template: dict, envelope: dict) -> dict:
    p = dict(template)
    for spec, u in zip(specs, unit):
        p[spec["name"]] = decode(spec, u)
    p["T_bath"] = float(envelope["parameters"]["T_bath"]["nominal"])
    p["R_SH"] = correlated_r_sh(envelope, p["R"])
    p["hardware_bessel_order"] = HARDWARE_BESSEL_ORDER
    return p


def unit_from_parameters(p: dict, specs: list[dict]) -> np.ndarray:
    return np.array([encode(s, p[s["name"]]) for s in specs], dtype=float)


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


def elite_pool(rows: list[dict], n: int) -> list[dict]:
    selected = {}
    for field in SCORE_FIELDS.values():
        for row in sorted(rows, key=lambda x, f=field: (x[f], x["trial_id"]))[:n]:
            selected[row["trial_id"]] = row
    for row in pareto_front(rows)[:max(n, 8)]:
        selected[row["trial_id"]] = row
    return list(selected.values())


def expand_bounds(specs: list[dict], elites: list[dict], factor=3.0) -> list[dict]:
    events = []
    if not elites:
        return events
    threshold = max(2, math.ceil(0.15 * len(elites)))
    for s in specs:
        if not s.get("expand"):
            continue
        u = np.array([encode(s, r["parameters"][s["name"]]) for r in elites])
        low_press, high_press = int((u < .05).sum()), int((u > .95).sum())
        before = [float(s["low"]), float(s["high"])]
        if low_press >= threshold and s["kind"] == "log":
            s["low"] = max(float(s["low"]) / factor, 1e-30)
        if high_press >= threshold:
            s["high"] = float(s["high"]) * (factor if s["kind"] == "log" else 2.0)
        if s["name"] == "n":
            s["low"] = max(1.01, float(s["low"]))
        if s["name"] == "excess_johnson_M":
            s["low"], s["high"] = 0.0, min(float(s["high"]), 20.0)
        after = [float(s["low"]), float(s["high"])]
        if before != after:
            events.append({"parameter": s["name"], "before": before, "after": after,
                           "low_pressure": low_press, "high_pressure": high_press})
    return events


def sensitivity_scan(base: dict, specs: list[dict], envelope: dict, constraints: dict,
                     grids, targets, rate, cutoff, points: int) -> dict:
    base_scores = score_parameters(base, grids, targets, rate, cutoff)
    result = {}
    for s in specs:
        rows = []
        for u in np.linspace(0, 1, points):
            p = dict(base)
            p[s["name"]] = decode(s, u)
            if s["name"] == "R":
                p["R_SH"] = correlated_r_sh(envelope, p["R"])
            ok, _ = pulse_gate(p, constraints)
            if not ok:
                continue
            try:
                scores = score_parameters(p, grids, targets, rate, cutoff)
            except Exception:
                continue
            rows.append({"value": p[s["name"]], **scores})
        summary = {}
        for band, field in SCORE_FIELDS.items():
            if rows:
                best = min(rows, key=lambda r, f=field: r[f])
                summary[band] = {"baseline": base_scores[field], "best": best[field],
                                 "best_value": best["value"],
                                 "improvement": base_scores[field] - best[field]}
            else:
                summary[band] = {"baseline": base_scores[field], "best": None,
                                 "best_value": None, "improvement": None}
        result[s["name"]] = {"rows": rows, "summary": summary}
    return {"baseline": base_scores, "parameters": result}


def polish(seed: dict, objective: str, specs, template, envelope, constraints,
           grids, targets, rate, cutoff, start_index: int, passes=2, points=7):
    current, out, idx = seed, [], start_index
    field = SCORE_FIELDS[objective]
    for pass_i in range(passes):
        width = .12 * (.55 ** pass_i)
        for dim, s in enumerate(specs):
            center = encode(s, current["parameters"][s["name"]])
            for u in np.linspace(max(0, center-width), min(1, center+width), points):
                vec = unit_from_parameters(current["parameters"], specs)
                vec[dim] = u
                idx += 1
                row, _ = evaluate(vec, f"polish_{objective}_{idx}", f"polish_{objective}",
                                  specs, template, envelope, constraints,
                                  grids, targets, rate, cutoff, idx)
                if row is not None:
                    out.append(row)
                    if row[field] < current[field]:
                        current = row
    return out, current, idx


def finite_pool(rows: list[dict], n: int) -> list[dict]:
    selected = {}
    for field in SCORE_FIELDS.values():
        for r in sorted(rows, key=lambda x, f=field: (x[f], x["trial_id"]))[:n]:
            selected[r["trial_id"]] = r
    return list(selected.values())


def boundary_hits(row: dict, specs: list[dict]) -> list[dict]:
    hits = []
    for s in specs:
        u = encode(s, row["parameters"][s["name"]])
        if u < .03 or u > .97:
            hits.append({"parameter": s["name"], "side": "low" if u < .03 else "high",
                         "value": row["parameters"][s["name"]],
                         "range": [s["low"], s["high"]]})
    return hits


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--trials", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=20260906)
    ap.add_argument("--refine-rounds", type=int, default=3)
    ap.add_argument("--refine-trials", type=int, default=600)
    ap.add_argument("--elite-count", type=int, default=12)
    ap.add_argument("--local-sigma", type=float, default=.15)
    ap.add_argument("--sensitivity-points", type=int, default=9)
    ap.add_argument("--excess-johnson-max", type=float, default=5.0)
    ap.add_argument("--finite-top", type=int, default=3)
    ap.add_argument("--finite-records", type=int, default=0)
    ap.add_argument("--finite-seed", type=int, default=20260906)
    args = ap.parse_args()

    case = args.case_dir
    scenarios = json.loads((case / "proxy_scenarios.json").read_text(encoding="utf-8"))
    envelope = json.loads((case / "proxy_parameter_envelope.json").read_text(encoding="utf-8"))
    constraints = json.loads((case / "pulse_combination_constraints.json").read_text(encoding="utf-8"))
    if scenarios.get("freeze_status") != "frozen" or envelope.get("freeze_status") != "frozen":
        raise RuntimeError("adaptive search requires frozen Stage-A inputs")
    frozen = scenarios["pulse_consistent_scenarios"]
    if not frozen:
        raise RuntimeError("no frozen pulse-consistent scenarios")

    exp_f, exp_asd, exp_paths, acq = experiment_asd(args.target_root)
    rate, sample = float(acq["rate_Hz"]), int(acq["samples"])
    cutoff = float(acq["analysis_bessel_cutoff_Hz"])
    if rate / 2 < 100_000:
        raise RuntimeError("Nyquist below 100 kHz")
    exp_norm = normalized(exp_asd, exp_f)
    grids = evaluation_grids()
    targets = experimental_targets(exp_f, exp_norm, grids)

    template = dict(frozen[0]["parameters"])
    specs = specs_from_envelope(envelope, template, args.excess_johnson_max)

    frozen_rows = []
    for s in frozen:
        p = dict(s["parameters"])
        p["hardware_bessel_order"] = HARDWARE_BESSEL_ORDER
        p["excess_johnson_M"] = float(p.get("excess_johnson_M", 0.0))
        frozen_rows.append({"scenario_id": s["scenario_id"],
                            **score_parameters(p, grids, targets, rate, cutoff),
                            "parameters": p})
    frozen_best = best_by_objective(frozen_rows)
    sensitivity = sensitivity_scan(dict(frozen_best["all"]["parameters"]), specs,
                                   envelope, constraints, grids, targets, rate, cutoff,
                                   args.sensitivity_points)

    rng = np.random.default_rng(args.seed)
    rows, rejected, idx = [], {"unstable": 0, "pulse_gate": 0, "model_failure": 0}, 0

    seen = set()
    for fr in frozen_best.values():
        sid = fr["scenario_id"]
        if sid in seen:
            continue
        seen.add(sid)
        idx += 1
        row, reason = evaluate(unit_from_parameters(fr["parameters"], specs),
                               f"seed_{sid}", "frozen_geometry_seed", specs, template,
                               envelope, constraints, grids, targets, rate, cutoff, idx)
        if row is None:
            rejected[reason] += 1
        else:
            rows.append(row)

    for i, vec in enumerate(latin_hypercube(args.trials, len(specs), rng)):
        idx += 1
        row, reason = evaluate(vec, f"global_{i:06d}", "global_latin_hypercube",
                               specs, template, envelope, constraints, grids, targets,
                               rate, cutoff, idx)
        if row is None:
            rejected[reason] += 1
        else:
            rows.append(row)
    if not rows:
        raise RuntimeError("no accepted exploratory candidates")

    expansion_events = []
    for round_i in range(args.refine_rounds):
        elites = elite_pool(rows, args.elite_count)
        events = expand_bounds(specs, elites)
        for e in events:
            e["round"] = round_i + 1
        expansion_events.extend(events)
        sigma = args.local_sigma * (.55 ** round_i)
        for i in range(args.refine_trials):
            center = unit_from_parameters(elites[int(rng.integers(len(elites)))]["parameters"], specs)
            vec = np.clip(center + rng.normal(0, sigma, len(specs)), 0, 1)
            idx += 1
            row, reason = evaluate(vec, f"refine_{round_i+1}_{i:06d}",
                                   f"elite_refinement_{round_i+1}", specs, template,
                                   envelope, constraints, grids, targets, rate, cutoff, idx)
            if row is None:
                rejected[reason] += 1
            else:
                rows.append(row)

    for objective, seed in list(best_by_objective(rows).items()):
        added, _, idx = polish(seed, objective, specs, template, envelope, constraints,
                               grids, targets, rate, cutoff, idx)
        rows.extend(added)

    best = best_by_objective(rows)
    pareto = pareto_front(rows)
    nrecords = args.finite_records if args.finite_records > 0 else len(exp_paths)
    finite_rows = [realize_finite(r, sample, rate, cutoff, nrecords, args.finite_seed,
                                  grids, targets) for r in finite_pool(rows, args.finite_top)]
    finite_best = finite_best_by_objective(finite_rows)

    frozen_unique = {r["scenario_id"]: r for r in frozen_best.values()}
    frozen_finite = [realize_finite(r, sample, rate, cutoff, nrecords, args.finite_seed,
                                    grids, targets) for r in frozen_unique.values()]
    frozen_finite_best = finite_best_by_objective(frozen_finite)

    result = {
        "stage": "adaptive_all_parameter_noise_guided_search",
        "strict_target_conclusion": STRICT,
        "strict_target_parameter_estimate_allowed": False,
        "searched_parameters": list(SEARCH_PARAMETERS),
        "searched_parameter_count": len(SEARCH_PARAMETERS),
        "fixed": {"T_bath_K": envelope["parameters"]["T_bath"]["nominal"],
                  "hardware_bessel_order": HARDWARE_BESSEL_ORDER,
                  "hardware_bessel_cutoff_Hz": 100_000.0,
                  "analysis_bessel_cutoff_Hz": cutoff, "rate_Hz": rate,
                  "samples": sample, "F_LINK": 0.9},
        "R_SH_note": "provenance-only; tes_noise_model does not consume R_SH",
        "parameter_specs_final": specs,
        "bound_expansions": expansion_events,
        "method": {"global": "Latin hypercube", "global_trials": args.trials,
                   "refine_rounds": args.refine_rounds,
                   "refine_trials_per_round": args.refine_trials,
                   "coordinate_polish": True},
        "noise_residual_fit": False, "simulation_amplitude_rescale": False,
        "evaluation_count": idx, "accepted_count": len(rows), "rejected": rejected,
        "sensitivity": sensitivity,
        "frozen_best": {k: public_row(v) for k, v in frozen_best.items()},
        "best_deterministic": {k: {**public_row(v), "boundary_hits": boundary_hits(v, specs)}
                               for k, v in best.items()},
        "best_finite": {k: public_row(v) for k, v in finite_best.items()},
        "pareto_front": [public_row(v) for v in pareto[:200]],
    }
    dump(case / "adaptive_high_frequency_parameter_search.json", result)

    import matplotlib.pyplot as plt

    pf = np.logspace(np.log10(rate / sample), np.log10(rate / 2), 1200)
    ep = log_interp(exp_f[1:], exp_norm[1:], pf)
    def curve(row):
        return log_interp(row["_finite_frequency_Hz"][1:],
                          row["_finite_normalized_asd"][1:], pf)

    curves = {
        "Frozen balanced": curve(frozen_finite_best["all"]),
        "Mid best": curve(finite_best["mid"]),
        "High best": curve(finite_best["high"]),
        "Balanced best": curve(finite_best["all"]),
    }
    plt.figure(figsize=(9, 5.5))
    plt.plot(pf, ep, color="black", lw=1.5,
             label=f"Experiment ({len(exp_paths)} accepted records)")
    for label, values in curves.items():
        plt.plot(pf, values, lw=1.1, label=label)
    plt.xscale("log"); plt.yscale("log"); plt.xlim(rate/sample, rate/2)
    plt.xlabel("Frequency [Hz]"); plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("Adaptive all-parameter TES noise search")
    plt.grid(True, which="both", alpha=.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "adaptive_high_frequency_parameter_search.png", dpi=180); plt.close()

    mask = (pf >= 1_000) & (pf <= 100_000)
    plt.figure(figsize=(9, 5))
    for label, values in curves.items():
        plt.plot(pf[mask], values[mask] / ep[mask], lw=1.1, label=label)
    plt.axhline(1, color="black", lw=1); plt.xscale("log"); plt.yscale("log")
    plt.xlim(1_000, 100_000); plt.xlabel("Frequency [Hz]")
    plt.ylabel("Simulation / experiment"); plt.title("Adaptive search residual ratio")
    plt.grid(True, which="both", alpha=.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "adaptive_high_frequency_parameter_search_ratio.png", dpi=180); plt.close()

    names = list(SEARCH_PARAMETERS)
    high_imp = [sensitivity["parameters"][n]["summary"]["high"]["improvement"] or 0 for n in names]
    all_imp = [sensitivity["parameters"][n]["summary"]["all"]["improvement"] or 0 for n in names]
    y = np.arange(len(names)); h = .38
    plt.figure(figsize=(9.5, 6.2))
    plt.barh(y-h/2, high_imp, h, label="10–100 kHz")
    plt.barh(y+h/2, all_imp, h, label="1–100 kHz")
    plt.yticks(y, names); plt.axvline(0, color="black", lw=.8)
    plt.xlabel("Best one-at-a-time RMS-log improvement")
    plt.title("Parameter sensitivity around frozen balanced baseline")
    plt.grid(True, axis="x", alpha=.25); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(case / "adaptive_high_frequency_parameter_sensitivity.png", dpi=180); plt.close()

    plt.figure(figsize=(6.2, 5.2))
    plt.scatter([r[SCORE_FIELDS["mid"]] for r in rows],
                [r[SCORE_FIELDS["high"]] for r in rows], s=8, alpha=.2)
    plt.plot([r[SCORE_FIELDS["mid"]] for r in pareto],
             [r[SCORE_FIELDS["high"]] for r in pareto], marker="o", ms=3, lw=1)
    plt.xlabel("1–10 kHz RMS log-ratio"); plt.ylabel("10–100 kHz RMS log-ratio")
    plt.title("Adaptive-search Pareto front"); plt.grid(True, alpha=.25); plt.tight_layout()
    plt.savefig(case / "adaptive_high_frequency_parameter_pareto.png", dpi=180); plt.close()

    print(json.dumps({
        "searched_parameters": list(SEARCH_PARAMETERS),
        "evaluations": idx, "accepted": len(rows), "rejected": rejected,
        "bound_expansions": len(expansion_events),
        "best": {k: {"id": v["trial_id"], "score": v[SCORE_FIELDS[k]],
                     "boundary_hits": boundary_hits(v, specs)} for k, v in best.items()},
    }, indent=2))


if __name__ == "__main__":
    main()

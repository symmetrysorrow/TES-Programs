"""Exploratory parameter search for the TES 1--10 kHz noise shape.

This script is deliberately separate from the frozen proxy-ensemble workflow.
It generates *new* trial points inside the same Stage-A proxy coordinate ranges,
applies the same pulse-time-constant gate, and uses the measured noise spectrum
only to rank those trials.  Therefore its output is exploratory and must never
be promoted to a strict target parameter estimate.

Search strategy:
1. random draws from the original Stage-A parameter ranges,
2. stability + pulse slow-pole gate,
3. deterministic 1--10 kHz post-analysis shape score using
   100 kHz analog hardware Bessel + 10 kHz software Bessel expectation,
4. finite-record time-domain re-evaluation of only the top candidates using
   the same accepted-record count and a fixed common random seed.

No additive residual/noise-floor parameter is fitted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations"))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))

from high_frequency_noise_comparison import (  # noqa: E402
    experiment_asd,
    log_interp,
    normalized,
)
from noise_measurement_model import (  # noqa: E402
    DEFAULT_FINITE_RECORD_SEED,
    HARDWARE_BESSEL_CUTOFF_HZ,
    expected_post_analysis_asd,
    finite_record_post_analysis_asd,
    hardware_sampled_asd,
)
from proxy_physics import linear_modes, operating_point  # noqa: E402


STRICT = "C — exact target physical case remains unidentified"
SEARCH_VARY = (
    "R_l",
    "alpha",
    "beta",
    "L",
    "n",
    "C_tes",
    "C_abs",
    "G_tes-bath",
    "G_abs-tes",
    "G_abs-abs",
)


def dump(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def rms_log_ratio(simulated: np.ndarray, experimental: np.ndarray) -> tuple[float, float]:
    log_ratio = np.log(
        np.asarray(simulated, dtype=float) / np.asarray(experimental, dtype=float)
    )
    return (
        float(np.sqrt(np.mean(log_ratio**2))),
        float(np.max(np.abs(log_ratio))),
    )


def sample_parameters(
    rng: np.random.Generator,
    template: dict,
    envelope: dict,
) -> tuple[dict, dict]:
    """Draw one new point from exactly the original Stage-A proxy ranges."""
    generic = envelope["sensitivity_reference"]
    params = dict(template)

    tc_low, tc_high = map(float, envelope["parameters"]["T_c"]["range"])
    r_low, r_high = map(float, envelope["parameters"]["R"]["range"])
    params["T_c"] = float(rng.uniform(tc_low, tc_high))
    params["T_bath"] = float(envelope["parameters"]["T_bath"]["nominal"])

    branch = int(rng.integers(0, 2))
    params["R"] = (r_low, r_high)[branch]
    params["R_SH"] = (0.0038, 0.0039)[branch]

    factors = {}
    for name in SEARCH_VARY:
        factor = float(np.exp(rng.uniform(np.log(0.5), np.log(2.0))))
        factors[name] = factor
        params[name] = float(generic[name]) * factor

    # Acquisition metadata are not searched.  In particular, ``cutoff`` remains
    # provenance for the 100 kHz hardware configuration and is not reused as
    # the 10 kHz analysis filter.
    params["rate"] = float(template["rate"])
    params["samples"] = int(template["samples"])
    return params, factors


def pulse_gate(parameters: dict, constraints: dict) -> tuple[bool, dict]:
    point = operating_point(parameters)
    if not point.get("stable", False):
        return False, {
            "stable": False,
            "reason": point.get("reason", "unstable"),
            "slowest_model_time_constant_s": None,
        }
    modes = linear_modes(parameters)
    slow_tau = modes[0]["time_constant_s"] if modes else None
    low, high = map(
        float,
        constraints["scenario_slow_pole_time_constant_range_s"],
    )
    accepted = bool(slow_tau is not None and low <= slow_tau <= high)
    return accepted, {
        "stable": True,
        "slowest_model_time_constant_s": (
            float(slow_tau) if slow_tau is not None else None
        ),
        "allowed_range_s": [low, high],
        "accepted": accepted,
    }


def score_parameters(
    parameters: dict,
    eval_freq: np.ndarray,
    exp_eval: np.ndarray,
    rate_hz: float,
    analysis_cutoff_hz: float,
) -> tuple[float, float]:
    expected = expected_post_analysis_asd(
        parameters,
        eval_freq,
        rate_hz=rate_hz,
        hardware_cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
        analysis_cutoff_hz=analysis_cutoff_hz,
    )
    expected_norm = expected / expected[0]
    return rms_log_ratio(expected_norm, exp_eval)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--scenario-file", type=Path)
    ap.add_argument("--envelope-file", type=Path)
    ap.add_argument("--constraints-file", type=Path)
    ap.add_argument("--trials", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260906)
    ap.add_argument(
        "--finite-top",
        type=int,
        default=3,
        help="number of deterministic top candidates to re-evaluate in time domain",
    )
    ap.add_argument(
        "--finite-records",
        type=int,
        default=0,
        help="0 uses the experimental accepted-record count",
    )
    ap.add_argument(
        "--finite-seed",
        type=int,
        default=DEFAULT_FINITE_RECORD_SEED,
        help="common seed used for every finite candidate for fair re-ranking",
    )
    args = ap.parse_args()

    if args.trials <= 0:
        raise ValueError("--trials must be positive")
    if args.finite_top <= 0:
        raise ValueError("--finite-top must be positive")

    args.scenario_file = args.scenario_file or args.case_dir / "proxy_scenarios.json"
    args.envelope_file = (
        args.envelope_file or args.case_dir / "proxy_parameter_envelope.json"
    )
    args.constraints_file = (
        args.constraints_file or args.case_dir / "pulse_combination_constraints.json"
    )

    scenarios_doc = json.loads(args.scenario_file.read_text(encoding="utf-8"))
    envelope = json.loads(args.envelope_file.read_text(encoding="utf-8"))
    constraints = json.loads(args.constraints_file.read_text(encoding="utf-8"))
    frozen = scenarios_doc["pulse_consistent_scenarios"]
    if scenarios_doc.get("freeze_status") != "frozen":
        raise RuntimeError("search requires the already-frozen Stage-A scenario document")
    if envelope.get("freeze_status") != "frozen":
        raise RuntimeError("search requires the frozen Stage-A parameter envelope")
    if not frozen:
        raise RuntimeError("no frozen pulse-consistent scenario is available")

    exp_freq, exp_asd, exp_paths, acquisition = experiment_asd(args.target_root)
    rate = float(acquisition["rate_Hz"])
    sample = int(acquisition["samples"])
    analysis_cutoff = float(acquisition["analysis_bessel_cutoff_Hz"])
    exp_norm = normalized(exp_asd, exp_freq)
    eval_freq = np.logspace(np.log10(1000.0), np.log10(10000.0), 200)
    exp_eval = log_interp(exp_freq[1:], exp_norm[1:], eval_freq)

    # Baseline: re-score the frozen candidates through the corrected deterministic
    # measurement chain before creating any new trial.
    frozen_rows = []
    for scenario in frozen:
        score, max_abs = score_parameters(
            scenario["parameters"],
            eval_freq,
            exp_eval,
            rate,
            analysis_cutoff,
        )
        frozen_rows.append(
            {
                "scenario_id": scenario["scenario_id"],
                "rms_log_ratio": score,
                "max_abs_log_ratio": max_abs,
                "parameters": scenario["parameters"],
            }
        )
    frozen_best = min(
        frozen_rows,
        key=lambda row: (row["rms_log_ratio"], row["scenario_id"]),
    )

    rng = np.random.default_rng(int(args.seed))
    template = dict(frozen[0]["parameters"])
    accepted_rows = []
    rejected_unstable = 0
    rejected_pulse_gate = 0

    for index in range(int(args.trials)):
        params, factors = sample_parameters(rng, template, envelope)
        accepted, gate = pulse_gate(params, constraints)
        if not gate.get("stable", False):
            rejected_unstable += 1
            continue
        if not accepted:
            rejected_pulse_gate += 1
            continue
        try:
            score, max_abs = score_parameters(
                params,
                eval_freq,
                exp_eval,
                rate,
                analysis_cutoff,
            )
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            rejected_unstable += 1
            continue
        accepted_rows.append(
            {
                "trial_id": f"search_{index:06d}",
                "rms_log_ratio": score,
                "max_abs_log_ratio": max_abs,
                "parameters": params,
                "factor_from_generic_reference": factors,
                "pulse_consistency": gate,
                "origin": "new_noise_ranked_draw_inside_original_stage_A_ranges",
                "strict_target_allowed": False,
            }
        )

    if not accepted_rows:
        raise RuntimeError("no exploratory trial passed stability and pulse gates")
    accepted_rows.sort(key=lambda row: (row["rms_log_ratio"], row["trial_id"]))
    deterministic_best = accepted_rows[0]

    finite_count = (
        int(args.finite_records) if args.finite_records > 0 else len(exp_paths)
    )
    finite_top = min(int(args.finite_top), len(accepted_rows))
    full_freq = np.fft.rfftfreq(sample, d=1.0 / rate)

    # Realize the frozen baseline with the same finite-record estimator and the
    # same random seed used for exploratory candidates.  This keeps the main
    # comparison figure apples-to-apples and avoids plotting the exact Nyquist
    # zero of the analytic digital-IIR response as if it were measured ASD.
    frozen_pre_analysis = hardware_sampled_asd(
        frozen_best["parameters"],
        full_freq,
        rate_hz=rate,
        cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
    )
    frozen_finite_asd = finite_record_post_analysis_asd(
        frozen_pre_analysis,
        sample,
        rate,
        analysis_cutoff_hz=analysis_cutoff,
        records=finite_count,
        seed=int(args.finite_seed),
    )
    frozen_finite_norm = normalized(frozen_finite_asd, full_freq)
    frozen_finite_eval = log_interp(full_freq[1:], frozen_finite_norm[1:], eval_freq)
    frozen_finite_score, frozen_finite_max = rms_log_ratio(
        frozen_finite_eval,
        exp_eval,
    )

    finite_rows = []
    for row in accepted_rows[:finite_top]:
        pre_analysis = hardware_sampled_asd(
            row["parameters"],
            full_freq,
            rate_hz=rate,
            cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
        )
        finite_asd = finite_record_post_analysis_asd(
            pre_analysis,
            sample,
            rate,
            analysis_cutoff_hz=analysis_cutoff,
            records=finite_count,
            seed=int(args.finite_seed),
        )
        finite_norm = normalized(finite_asd, full_freq)
        finite_eval = log_interp(full_freq[1:], finite_norm[1:], eval_freq)
        finite_score, finite_max = rms_log_ratio(finite_eval, exp_eval)
        finite_rows.append(
            {
                **row,
                "finite_record_rms_log_ratio": finite_score,
                "finite_record_max_abs_log_ratio": finite_max,
                "finite_record_count": finite_count,
                "finite_record_seed": int(args.finite_seed),
                "_finite_frequency_Hz": full_freq,
                "_finite_normalized_asd": finite_norm,
            }
        )

    finite_rows.sort(
        key=lambda row: (
            row["finite_record_rms_log_ratio"],
            row["trial_id"],
        )
    )
    finite_best = finite_rows[0]
    score_improvement = float(
        frozen_best["rms_log_ratio"] - deterministic_best["rms_log_ratio"]
    )
    finite_improvement_vs_frozen_deterministic = float(
        frozen_best["rms_log_ratio"] - finite_best["finite_record_rms_log_ratio"]
    )
    finite_improvement_vs_frozen_finite = float(
        frozen_finite_score - finite_best["finite_record_rms_log_ratio"]
    )

    result = {
        "stage": "exploratory_noise_guided_parameter_search",
        "strict_target_conclusion": STRICT,
        "strict_target_parameter_estimate_allowed": False,
        "separate_from_frozen_ensemble": True,
        "warning": (
            "These new trial parameters are selected using the experimental noise "
            "shape and therefore are exploratory only. They must not be written back "
            "into proxy_scenarios.json or the strict target input."
        ),
        "selection_band_Hz": [1000.0, 10000.0],
        "low_frequency_excluded_from_score": True,
        "normalization_frequency_Hz": 1000.0,
        "search_method": (
            "fixed-seed random draws inside original Stage-A ranges; stability and "
            "the original pulse slow-pole gate are applied before noise ranking"
        ),
        "noise_used_for_ranking": True,
        "noise_residual_fit": False,
        "additive_noise_parameter_fit": False,
        "simulation_amplitude_rescale": False,
        "hardware_bessel_cutoff_Hz": float(HARDWARE_BESSEL_CUTOFF_HZ),
        "analysis_bessel_cutoff_Hz": analysis_cutoff,
        "trial_count_requested": int(args.trials),
        "trial_count_gate_accepted": len(accepted_rows),
        "rejected_unstable_or_model_failure": rejected_unstable,
        "rejected_pulse_gate": rejected_pulse_gate,
        "random_seed": int(args.seed),
        "finite_record_common_seed": int(args.finite_seed),
        "finite_record_count": finite_count,
        "finite_rerank_count": finite_top,
        "frozen_baseline_best": frozen_best,
        "frozen_baseline_finite_record": {
            "finite_record_rms_log_ratio": frozen_finite_score,
            "finite_record_max_abs_log_ratio": frozen_finite_max,
            "finite_record_count": finite_count,
            "finite_record_seed": int(args.finite_seed),
        },
        "best_exploratory_deterministic": deterministic_best,
        "best_exploratory_finite_record": {
            key: value
            for key, value in finite_best.items()
            if not key.startswith("_")
        },
        "deterministic_score_improvement_vs_frozen": score_improvement,
        "finite_score_improvement_vs_frozen_deterministic_baseline": (
            finite_improvement_vs_frozen_deterministic
        ),
        "finite_score_improvement_vs_frozen_finite_baseline": (
            finite_improvement_vs_frozen_finite
        ),
        "plot_curve_semantics": (
            "main PNG compares experiment, frozen baseline, and exploratory best "
            "through the same finite-record estimator; the deterministic analytic "
            "curve is intentionally not drawn through the exact Nyquist endpoint"
        ),
        "top_exploratory_deterministic": accepted_rows[: min(100, len(accepted_rows))],
        "top_exploratory_finite_record": [
            {key: value for key, value in row.items() if not key.startswith("_")}
            for row in finite_rows
        ],
    }
    dump(args.case_dir / "high_frequency_parameter_search.json", result)

    import matplotlib.pyplot as plt

    plot_freq = np.logspace(
        np.log10(rate / sample),
        np.log10(rate / 2.0),
        1200,
    )
    exp_plot = log_interp(exp_freq[1:], exp_norm[1:], plot_freq)
    frozen_finite_plot = log_interp(
        full_freq[1:],
        frozen_finite_norm[1:],
        plot_freq,
    )
    finite_plot = log_interp(
        finite_best["_finite_frequency_Hz"][1:],
        finite_best["_finite_normalized_asd"][1:],
        plot_freq,
    )

    plt.figure(figsize=(9, 5.5))
    plt.plot(
        plot_freq,
        exp_plot,
        color="black",
        lw=1.5,
        label=f"Experiment post-analysis ({len(exp_paths)} accepted records)",
    )
    plt.plot(
        plot_freq,
        frozen_finite_plot,
        lw=1.3,
        ls="--",
        label=f"Frozen best finite ({frozen_best['scenario_id']})",
    )
    plt.plot(
        plot_freq,
        finite_plot,
        lw=1.2,
        label=f"Exploratory finite best ({finite_best['trial_id']})",
    )
    plt.scatter([1000.0], [1.0], color="black", s=28, zorder=5)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(rate / sample, rate / 2.0)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("Exploratory TES parameter search — corrected measurement chain")
    plt.suptitle(
        "Ranking uses 1–10 kHz only; plotted simulation curves use the same "
        "finite-record estimator as the experiment.",
        fontsize=9,
        y=0.94,
    )
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout(rect=(0, 0, 0.99, 0.90))
    plt.savefig(args.case_dir / "high_frequency_parameter_search.png", dpi=180)
    plt.close()

    print(
        json.dumps(
            {
                "frozen_best": {
                    "id": frozen_best["scenario_id"],
                    "deterministic_rms_log_ratio": frozen_best["rms_log_ratio"],
                    "finite_record_rms_log_ratio": frozen_finite_score,
                },
                "exploratory_expected_best": {
                    "id": deterministic_best["trial_id"],
                    "rms_log_ratio": deterministic_best["rms_log_ratio"],
                },
                "exploratory_finite_best": {
                    "id": finite_best["trial_id"],
                    "rms_log_ratio": finite_best["finite_record_rms_log_ratio"],
                },
                "gate_accepted": len(accepted_rows),
                "trials": int(args.trials),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

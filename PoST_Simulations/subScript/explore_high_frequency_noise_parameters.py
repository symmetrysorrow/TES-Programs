"""Exploratory multi-band parameter search for TES noise shape.

This script is deliberately separate from the frozen proxy-ensemble workflow.
It generates *new* trial points inside the same Stage-A proxy coordinate ranges,
applies the same pulse-time-constant gate, and uses the measured noise spectrum
only to rank those trials. Therefore its output is exploratory and must never
be promoted to a strict target parameter estimate.

The confirmed readout chain is held fixed during the search:
- 4th-order, 100 kHz analog hardware Bessel before the ADC,
- measured sample rate and record length,
- target 10 kHz second-order software Bessel ``filtfilt``,
- Hann -> rFFT -> record power average -> one-sided ASD.

Search strategy:
1. random draws from the original Stage-A parameter ranges,
2. stability + pulse slow-pole gate,
3. deterministic post-analysis scoring in three bands:
   1--10 kHz, 10--100 kHz, and 1--100 kHz,
4. Pareto analysis of the 1--10 kHz / 10--100 kHz trade-off,
5. finite-record time-domain re-evaluation of the top candidates from each
   deterministic objective using a common seed and the experimental record count.

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
HARDWARE_BESSEL_ORDER = 4
BANDS_HZ = {
    "mid": (1_000.0, 10_000.0),
    "high": (10_000.0, 100_000.0),
    "all": (1_000.0, 100_000.0),
}
SCORE_FIELDS = {
    "mid": "rms_log_ratio_1_10_kHz",
    "high": "rms_log_ratio_10_100_kHz",
    "all": "rms_log_ratio_1_100_kHz",
}
MAX_FIELDS = {
    "mid": "max_abs_log_ratio_1_10_kHz",
    "high": "max_abs_log_ratio_10_100_kHz",
    "all": "max_abs_log_ratio_1_100_kHz",
}
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


def rms_log_ratio(
    simulated: np.ndarray,
    experimental: np.ndarray,
) -> tuple[float, float]:
    log_ratio = np.log(
        np.asarray(simulated, dtype=float) / np.asarray(experimental, dtype=float)
    )
    return (
        float(np.sqrt(np.mean(log_ratio**2))),
        float(np.max(np.abs(log_ratio))),
    )


def evaluation_grids() -> dict[str, np.ndarray]:
    """Use equal log-frequency density: 200 points per decade."""
    return {
        "mid": np.logspace(np.log10(1_000.0), np.log10(10_000.0), 200),
        "high": np.logspace(np.log10(10_000.0), np.log10(100_000.0), 200),
        "all": np.logspace(np.log10(1_000.0), np.log10(100_000.0), 400),
    }


def experimental_targets(
    exp_freq: np.ndarray,
    exp_norm: np.ndarray,
    grids: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    return {
        name: log_interp(exp_freq[1:], exp_norm[1:], grid)
        for name, grid in grids.items()
    }


def _row_id(row: dict) -> str:
    return str(row.get("trial_id", row.get("scenario_id", "")))


def pareto_front(
    rows: list[dict],
    mid_field: str = SCORE_FIELDS["mid"],
    high_field: str = SCORE_FIELDS["high"],
) -> list[dict]:
    """Return the non-dominated mid/high-error front, sorted by mid error."""
    ordered = sorted(
        rows,
        key=lambda row: (
            float(row[mid_field]),
            float(row[high_field]),
            _row_id(row),
        ),
    )
    front = []
    best_high = np.inf
    for row in ordered:
        high = float(row[high_field])
        if high < best_high:
            front.append(row)
            best_high = high
    return front


def best_by_objective(rows: list[dict]) -> dict[str, dict]:
    return {
        name: min(
            rows,
            key=lambda row, field=field: (float(row[field]), _row_id(row)),
        )
        for name, field in SCORE_FIELDS.items()
    }


def finite_best_by_objective(rows: list[dict]) -> dict[str, dict]:
    return {
        name: min(
            rows,
            key=lambda row, field=field: (
                float(row[f"finite_{field}"]),
                _row_id(row),
            ),
        )
        for name, field in SCORE_FIELDS.items()
    }


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

    # Acquisition/hardware metadata are fixed, not searched.
    params["rate"] = float(template["rate"])
    params["samples"] = int(template["samples"])
    params["hardware_bessel_order"] = HARDWARE_BESSEL_ORDER
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
    grids: dict[str, np.ndarray],
    exp_targets: dict[str, np.ndarray],
    rate_hz: float,
    analysis_cutoff_hz: float,
) -> dict[str, float]:
    model_frequency = np.unique(
        np.concatenate(
            (
                np.asarray([1_000.0]),
                *[np.asarray(grid, dtype=float) for grid in grids.values()],
            )
        )
    )
    fixed = dict(parameters)
    fixed["hardware_bessel_order"] = HARDWARE_BESSEL_ORDER
    expected = expected_post_analysis_asd(
        fixed,
        model_frequency,
        rate_hz=rate_hz,
        hardware_cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
        analysis_cutoff_hz=analysis_cutoff_hz,
    )
    norm_1k = log_interp(
        model_frequency,
        expected,
        np.asarray([1_000.0]),
    )[0]
    expected_norm = expected / norm_1k
    scores: dict[str, float] = {}
    for name, grid in grids.items():
        simulated = log_interp(model_frequency, expected_norm, grid)
        rms, maximum = rms_log_ratio(simulated, exp_targets[name])
        scores[SCORE_FIELDS[name]] = rms
        scores[MAX_FIELDS[name]] = maximum
    return scores


def finite_scores(
    finite_norm: np.ndarray,
    full_freq: np.ndarray,
    grids: dict[str, np.ndarray],
    exp_targets: dict[str, np.ndarray],
) -> dict[str, float]:
    scores: dict[str, float] = {}
    for name, grid in grids.items():
        simulated = log_interp(full_freq[1:], finite_norm[1:], grid)
        rms, maximum = rms_log_ratio(simulated, exp_targets[name])
        scores[f"finite_{SCORE_FIELDS[name]}"] = rms
        scores[f"finite_{MAX_FIELDS[name]}"] = maximum
    return scores


def deterministic_candidate_pool(
    rows: list[dict],
    per_objective: int,
) -> list[dict]:
    """Union of the top N deterministic candidates from each objective."""
    selected: dict[str, dict] = {}
    for field in SCORE_FIELDS.values():
        ordered = sorted(
            rows,
            key=lambda row, score_field=field: (
                float(row[score_field]),
                _row_id(row),
            ),
        )
        for row in ordered[:per_objective]:
            selected[_row_id(row)] = row
    return [selected[key] for key in sorted(selected)]


def realize_finite(
    row: dict,
    sample: int,
    rate: float,
    analysis_cutoff: float,
    records: int,
    seed: int,
    grids: dict[str, np.ndarray],
    exp_targets: dict[str, np.ndarray],
) -> dict:
    parameters = dict(row["parameters"])
    parameters["hardware_bessel_order"] = HARDWARE_BESSEL_ORDER
    full_freq = np.fft.rfftfreq(sample, d=1.0 / rate)
    pre_analysis = hardware_sampled_asd(
        parameters,
        full_freq,
        rate_hz=rate,
        cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
        order=HARDWARE_BESSEL_ORDER,
    )
    finite_asd = finite_record_post_analysis_asd(
        pre_analysis,
        sample,
        rate,
        analysis_cutoff_hz=analysis_cutoff,
        records=records,
        seed=seed,
    )
    finite_norm = normalized(finite_asd, full_freq)
    scores = finite_scores(finite_norm, full_freq, grids, exp_targets)
    return {
        **row,
        **scores,
        "finite_record_count": records,
        "finite_record_seed": seed,
        "_finite_frequency_Hz": full_freq,
        "_finite_normalized_asd": finite_norm,
    }


def public_row(row: dict) -> dict:
    return {key: value for key, value in row.items() if not key.startswith("_")}


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
        help=(
            "top deterministic candidates per objective to re-evaluate in time "
            "domain; the union of mid/high/all lists is used"
        ),
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
    if rate / 2.0 < BANDS_HZ["high"][1]:
        raise RuntimeError(
            "Nyquist frequency is below 100 kHz; 10--100 kHz score is unavailable"
        )
    exp_norm = normalized(exp_asd, exp_freq)
    grids = evaluation_grids()
    exp_targets = experimental_targets(exp_freq, exp_norm, grids)

    # Re-score the frozen candidates through the fixed measurement chain.
    frozen_rows = []
    for scenario in frozen:
        scores = score_parameters(
            scenario["parameters"],
            grids,
            exp_targets,
            rate,
            analysis_cutoff,
        )
        frozen_rows.append(
            {
                "scenario_id": scenario["scenario_id"],
                **scores,
                "parameters": scenario["parameters"],
            }
        )
    frozen_best = best_by_objective(frozen_rows)
    frozen_pareto = pareto_front(frozen_rows)

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
            scores = score_parameters(
                params,
                grids,
                exp_targets,
                rate,
                analysis_cutoff,
            )
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            rejected_unstable += 1
            continue
        accepted_rows.append(
            {
                "trial_id": f"search_{index:06d}",
                **scores,
                "parameters": params,
                "factor_from_generic_reference": factors,
                "pulse_consistency": gate,
                "origin": "new_noise_ranked_draw_inside_original_stage_A_ranges",
                "strict_target_allowed": False,
            }
        )

    if not accepted_rows:
        raise RuntimeError("no exploratory trial passed stability and pulse gates")

    deterministic_best = best_by_objective(accepted_rows)
    deterministic_pareto = pareto_front(accepted_rows)

    finite_count = (
        int(args.finite_records) if args.finite_records > 0 else len(exp_paths)
    )

    # Finite-record evaluation is performed for the top N of each objective.
    candidate_pool = deterministic_candidate_pool(accepted_rows, int(args.finite_top))
    finite_rows = [
        realize_finite(
            row,
            sample,
            rate,
            analysis_cutoff,
            finite_count,
            int(args.finite_seed),
            grids,
            exp_targets,
        )
        for row in candidate_pool
    ]
    finite_best = finite_best_by_objective(finite_rows)
    finite_pareto = pareto_front(
        finite_rows,
        mid_field=f"finite_{SCORE_FIELDS['mid']}",
        high_field=f"finite_{SCORE_FIELDS['high']}",
    )

    # Evaluate each distinct frozen deterministic objective winner with the same
    # finite-record estimator for objective-by-objective baselines.
    frozen_unique: dict[str, dict] = {}
    for row in frozen_best.values():
        frozen_unique[_row_id(row)] = row
    frozen_finite_rows = [
        realize_finite(
            row,
            sample,
            rate,
            analysis_cutoff,
            finite_count,
            int(args.finite_seed),
            grids,
            exp_targets,
        )
        for row in frozen_unique.values()
    ]
    frozen_finite_best = finite_best_by_objective(frozen_finite_rows)

    deterministic_improvement = {
        name: float(frozen_best[name][field] - deterministic_best[name][field])
        for name, field in SCORE_FIELDS.items()
    }
    finite_improvement = {
        name: float(
            frozen_finite_best[name][f"finite_{field}"]
            - finite_best[name][f"finite_{field}"]
        )
        for name, field in SCORE_FIELDS.items()
    }

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
        "score_bands_Hz": {name: list(bounds) for name, bounds in BANDS_HZ.items()},
        "low_frequency_below_1kHz_excluded_from_score": True,
        "normalization_frequency_Hz": 1_000.0,
        "score_definition": (
            "RMS natural-log ASD ratio on log-spaced frequencies after independent "
            "1 kHz normalization; 200 points/decade"
        ),
        "ranking_policy": (
            "mid (1-10 kHz), high (10-100 kHz), and balanced/all (1-100 kHz) "
            "are ranked independently; Pareto front minimizes mid and high scores"
        ),
        "search_method": (
            "fixed-seed random draws inside original Stage-A ranges; stability and "
            "the original pulse slow-pole gate are applied before noise ranking"
        ),
        "noise_used_for_ranking": True,
        "noise_residual_fit": False,
        "additive_noise_parameter_fit": False,
        "simulation_amplitude_rescale": False,
        "hardware_bessel_cutoff_Hz": float(HARDWARE_BESSEL_CUTOFF_HZ),
        "hardware_bessel_order": HARDWARE_BESSEL_ORDER,
        "hardware_parameters_searched": False,
        "analysis_bessel_cutoff_Hz": analysis_cutoff,
        "trial_count_requested": int(args.trials),
        "trial_count_gate_accepted": len(accepted_rows),
        "rejected_unstable_or_model_failure": rejected_unstable,
        "rejected_pulse_gate": rejected_pulse_gate,
        "random_seed": int(args.seed),
        "finite_record_common_seed": int(args.finite_seed),
        "finite_record_count": finite_count,
        "finite_top_per_objective": int(args.finite_top),
        "finite_candidate_pool_size": len(candidate_pool),
        "finite_candidate_pool_ids": [_row_id(row) for row in candidate_pool],
        "frozen_best_deterministic_by_objective": {
            name: public_row(row) for name, row in frozen_best.items()
        },
        "frozen_pareto_count": len(frozen_pareto),
        "frozen_pareto_front": [public_row(row) for row in frozen_pareto],
        "frozen_best_finite_by_objective": {
            name: public_row(row) for name, row in frozen_finite_best.items()
        },
        "best_exploratory_deterministic_by_objective": {
            name: public_row(row) for name, row in deterministic_best.items()
        },
        "deterministic_pareto_count": len(deterministic_pareto),
        "deterministic_pareto_front": [
            public_row(row) for row in deterministic_pareto[:100]
        ],
        "best_exploratory_finite_by_objective": {
            name: public_row(row) for name, row in finite_best.items()
        },
        "finite_pareto_count": len(finite_pareto),
        "finite_pareto_front": [public_row(row) for row in finite_pareto],
        "deterministic_score_improvement_vs_frozen": deterministic_improvement,
        "finite_score_improvement_vs_frozen": finite_improvement,
        "top_exploratory_deterministic": {
            name: [
                public_row(row)
                for row in sorted(
                    accepted_rows,
                    key=lambda row, score_field=field: (
                        float(row[score_field]),
                        _row_id(row),
                    ),
                )[:100]
            ]
            for name, field in SCORE_FIELDS.items()
        },
        "top_exploratory_finite_record": [
            public_row(row)
            for row in sorted(
                finite_rows,
                key=lambda row: (
                    float(row[f"finite_{SCORE_FIELDS['all']}"]),
                    _row_id(row),
                ),
            )
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

    # Main figure: one frozen balanced baseline plus unique finite winners from
    # each exploratory objective.
    frozen_plot_row = frozen_finite_best["all"]
    frozen_plot = log_interp(
        frozen_plot_row["_finite_frequency_Hz"][1:],
        frozen_plot_row["_finite_normalized_asd"][1:],
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
        frozen_plot,
        lw=1.3,
        ls="--",
        label=f"Frozen balanced finite ({_row_id(frozen_plot_row)})",
    )

    labels = {
        "mid": "Exploratory mid-band best",
        "high": "Exploratory high-band best",
        "all": "Exploratory balanced best",
    }
    plotted_ids = set()
    for name in ("mid", "high", "all"):
        row = finite_best[name]
        identity = _row_id(row)
        if identity in plotted_ids:
            continue
        plotted_ids.add(identity)
        curve = log_interp(
            row["_finite_frequency_Hz"][1:],
            row["_finite_normalized_asd"][1:],
            plot_freq,
        )
        plt.plot(
            plot_freq,
            curve,
            lw=1.2,
            label=f"{labels[name]} ({identity})",
        )

    plt.scatter([1_000.0], [1.0], color="black", s=28, zorder=5)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(rate / sample, rate / 2.0)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("Exploratory TES parameter search — multi-band comparison")
    plt.suptitle(
        "Hardware fixed at 4th-order 100 kHz Bessel; rankings: 1–10 kHz, "
        "10–100 kHz, and 1–100 kHz.",
        fontsize=9,
        y=0.94,
    )
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout(rect=(0, 0, 0.99, 0.90))
    plt.savefig(args.case_dir / "high_frequency_parameter_search.png", dpi=180)
    plt.close()

    # Ratio figure directly exposes where each selected simulation is high/low.
    ratio_freq = grids["all"]
    exp_ratio_reference = exp_targets["all"]
    plt.figure(figsize=(9, 4.8))
    plt.axhline(1.0, color="black", lw=1.0)
    plotted_ids.clear()
    for name in ("mid", "high", "all"):
        row = finite_best[name]
        identity = _row_id(row)
        if identity in plotted_ids:
            continue
        plotted_ids.add(identity)
        simulated = log_interp(
            row["_finite_frequency_Hz"][1:],
            row["_finite_normalized_asd"][1:],
            ratio_freq,
        )
        plt.plot(
            ratio_freq,
            simulated / exp_ratio_reference,
            lw=1.2,
            label=f"{labels[name]} ({identity})",
        )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(BANDS_HZ["all"])
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Simulation / experiment normalized ASD")
    plt.title("Selected exploratory candidates — 1–100 kHz shape ratio")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(
        args.case_dir / "high_frequency_parameter_search_ratio.png",
        dpi=180,
    )
    plt.close()

    # Pareto figure: every accepted deterministic trial plus the non-dominated front.
    pareto_mid = [row[SCORE_FIELDS["mid"]] for row in deterministic_pareto]
    pareto_high = [row[SCORE_FIELDS["high"]] for row in deterministic_pareto]
    plt.figure(figsize=(6.5, 5.5))
    plt.scatter(
        [row[SCORE_FIELDS["mid"]] for row in accepted_rows],
        [row[SCORE_FIELDS["high"]] for row in accepted_rows],
        s=8,
        alpha=0.25,
        label="accepted exploratory trials",
    )
    plt.plot(
        pareto_mid,
        pareto_high,
        marker="o",
        ms=3,
        lw=1.2,
        label="mid/high Pareto front",
    )
    balanced = deterministic_best["all"]
    plt.scatter(
        [balanced[SCORE_FIELDS["mid"]]],
        [balanced[SCORE_FIELDS["high"]]],
        s=45,
        marker="*",
        label=f"balanced best ({_row_id(balanced)})",
    )
    plt.xlabel("RMS log ratio: 1–10 kHz")
    plt.ylabel("RMS log ratio: 10–100 kHz")
    plt.title("TES parameter-search Pareto trade-off")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(
        args.case_dir / "high_frequency_parameter_search_pareto.png",
        dpi=180,
    )
    plt.close()

    print(
        json.dumps(
            {
                "hardware": {
                    "bessel_order": HARDWARE_BESSEL_ORDER,
                    "bessel_cutoff_Hz": float(HARDWARE_BESSEL_CUTOFF_HZ),
                    "searched": False,
                },
                "frozen_best": {
                    name: {
                        "id": _row_id(row),
                        "rms_log_ratio": row[SCORE_FIELDS[name]],
                    }
                    for name, row in frozen_best.items()
                },
                "exploratory_expected_best": {
                    name: {
                        "id": _row_id(row),
                        "rms_log_ratio": row[SCORE_FIELDS[name]],
                    }
                    for name, row in deterministic_best.items()
                },
                "exploratory_finite_best": {
                    name: {
                        "id": _row_id(row),
                        "rms_log_ratio": row[f"finite_{SCORE_FIELDS[name]}"],
                    }
                    for name, row in finite_best.items()
                },
                "deterministic_pareto_count": len(deterministic_pareto),
                "gate_accepted": len(accepted_rows),
                "trials": int(args.trials),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

"""High-frequency experiment vs frozen proxy-ensemble comparison.

The frozen scenario list is still selected without changing any parameter.
The comparison now applies the actual measurement chain consistently:

simulation:
    intrinsic TES ASD
    -> 100 kHz analog hardware Bessel + ADC alias fold
    -> finite 500 kS/s, 100000-sample time records
    -> target 10 kHz software Bessel ``filtfilt``
    -> Hann -> rFFT power average -> one-sided ASD

experiment:
    accepted raw CH0 records (hardware filtering is already in the acquisition)
    -> target 10 kHz software Bessel ``filtfilt``
    -> Hann -> rFFT power average -> one-sided ASD

Frozen-scenario selection remains deterministic and uses only 1--10 kHz shape
error.  The selected member is then rendered through finite time records so the
displayed simulation has the same finite-record estimator as the experiment.
No noise residual, amplitude rescale, or parameter optimization is performed.
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

from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    estimate_one_sided_asd,
)
from lib.tes_noise_model import SOURCE_CLASS_INDICES  # noqa: E402
from noise_measurement_model import (  # noqa: E402
    ANALYSIS_BESSEL_CUTOFF_HZ,
    DEFAULT_FINITE_RECORD_SEED,
    HARDWARE_BESSEL_CUTOFF_HZ,
    expected_post_analysis_asd,
    finite_record_post_analysis_asd,
    hardware_sampled_asd,
)
from proxy_physics import noise_components  # noqa: E402
from pulse_contamination_common import read_record  # noqa: E402
from pulse_contamination_v5 import RATE_HZ, SAMPLES, accepted_noise  # noqa: E402


STRICT = "C — exact target physical case remains unidentified"
ANCHORS = (1000, 2000, 3000, 5000, 7000, 10000)
CLASSES = ("TES_Johnson", "load_Johnson", "TES_bath_TFN", "TES_absorber_TFN")
DISPLAY_MIN_HZ = RATE_HZ / SAMPLES
DISPLAY_MAX_HZ = RATE_HZ / 2.0


def dump(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def log_interp(x: np.ndarray, y: np.ndarray, query: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    query = np.asarray(query, dtype=float)
    good = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < 2:
        raise ValueError("log interpolation requires at least two positive finite points")
    return np.exp(np.interp(np.log(query), np.log(x[good]), np.log(y[good])))


def target_acquisition(target: Path) -> dict:
    config = json.loads((target / "PulseConfig.json").read_text(encoding="utf-8"))
    rate = float(config["Readout"]["Rate"])
    samples = int(config["Readout"]["Sample"])
    cutoff = float(config["Analysis"]["CutoffFrequency"])
    if not np.isclose(rate, RATE_HZ):
        raise RuntimeError(f"target rate {rate} Hz does not match analysis rate {RATE_HZ} Hz")
    if samples != SAMPLES:
        raise RuntimeError(
            f"target samples {samples} does not match analysis samples {SAMPLES}"
        )
    return {
        "rate_Hz": rate,
        "samples": samples,
        "analysis_bessel_cutoff_Hz": cutoff,
    }


def experiment_asd(
    target: Path,
) -> tuple[np.ndarray, np.ndarray, list[Path], dict]:
    acquisition = target_acquisition(target)
    paths = accepted_noise(target, "CH0")
    if not paths:
        raise RuntimeError("no accepted CH0 noise records")

    def records():
        for path in paths:
            yield read_record(path)

    asd, count = estimate_one_sided_asd(
        records(),
        acquisition["samples"],
        acquisition["rate_Hz"],
        cutoff=acquisition["analysis_bessel_cutoff_Hz"],
        remove_mean=True,
    )
    if count != len(paths):
        raise RuntimeError(
            f"canonical ASD estimator accepted {count} records, expected {len(paths)}"
        )
    freq = np.fft.rfftfreq(
        acquisition["samples"],
        d=1.0 / acquisition["rate_Hz"],
    )
    return freq, asd, paths, acquisition


def component_labels() -> dict[str, tuple[int, ...]]:
    return {name: tuple(indices) for name, indices in SOURCE_CLASS_INDICES.items()}


def intervals(freq: np.ndarray, mask: np.ndarray) -> list[list[float]]:
    mask = np.asarray(mask, dtype=bool)
    values = np.asarray(freq)[mask]
    if not len(values):
        return []
    groups = np.split(values, np.where(np.diff(np.where(mask)[0]) > 1)[0] + 1)
    return [[float(group[0]), float(group[-1])] for group in groups if len(group)]


def normalized(values: np.ndarray, frequency: np.ndarray, reference_hz: float = 1000.0):
    values = np.asarray(values, dtype=float)
    index = int(np.argmin(np.abs(np.asarray(frequency) - reference_hz)))
    reference = float(values[index])
    if not np.isfinite(reference) or reference <= 0.0:
        raise ValueError("normalization reference must be positive and finite")
    return values / reference


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--scenario-file", type=Path)
    ap.add_argument(
        "--finite-records",
        type=int,
        default=0,
        help="0 uses the same number of accepted CH0 records as the experiment",
    )
    ap.add_argument("--finite-seed", type=int, default=DEFAULT_FINITE_RECORD_SEED)
    args = ap.parse_args()

    args.scenario_file = args.scenario_file or args.case_dir / "proxy_scenarios.json"
    scenarios_doc = json.loads(args.scenario_file.read_text(encoding="utf-8"))
    scenarios = sorted(
        scenarios_doc["pulse_consistent_scenarios"],
        key=lambda row: row["scenario_id"],
    )
    if (
        scenarios_doc.get("freeze_status") != "frozen"
        or scenarios_doc.get("pulse_consistent_count") != len(scenarios)
    ):
        raise RuntimeError("comparison requires the frozen pulse-consistent scenario list")

    exp_freq, exp_asd, exp_paths, acquisition = experiment_asd(args.target_root)
    rate = acquisition["rate_Hz"]
    sample = acquisition["samples"]
    analysis_cutoff = acquisition["analysis_bessel_cutoff_Hz"]
    exp_norm = normalized(exp_asd, exp_freq)

    eval_freq = np.logspace(np.log10(1000.0), np.log10(10000.0), 200)
    plot_freq = np.logspace(
        np.log10(DISPLAY_MIN_HZ),
        np.log10(DISPLAY_MAX_HZ),
        1200,
    )
    model_freq = np.unique(
        np.r_[1000.0, eval_freq, plot_freq, np.asarray(ANCHORS, dtype=float)]
    )
    exp_eval = log_interp(exp_freq[1:], exp_norm[1:], eval_freq)
    exp_plot = log_interp(exp_freq[1:], exp_norm[1:], plot_freq)

    rows = []
    curves = {}
    for scenario in scenarios:
        expected = expected_post_analysis_asd(
            scenario["parameters"],
            model_freq,
            rate_hz=rate,
            hardware_cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
            analysis_cutoff_hz=analysis_cutoff,
        )
        expected_norm = normalized(expected, model_freq)
        sim_eval = log_interp(model_freq, expected_norm, eval_freq)
        log_ratio = np.log(sim_eval / exp_eval)
        score = float(np.sqrt(np.mean(log_ratio**2)))
        rows.append(
            {
                "scenario_id": scenario["scenario_id"],
                "rms_log_ratio": score,
                "max_abs_log_ratio": float(np.max(np.abs(log_ratio))),
                "parameter_source_classification": scenario.get(
                    "source_class_by_parameter", {}
                ),
                "parameters": scenario["parameters"],
                "pulse_consistency": scenario.get("pulse_consistency"),
            }
        )
        curves[scenario["scenario_id"]] = expected_norm

    best = min(rows, key=lambda row: (row["rms_log_ratio"], row["scenario_id"]))
    best_expected = curves[best["scenario_id"]]
    best_expected_plot = log_interp(model_freq, best_expected, plot_freq)

    finite_records = int(args.finite_records) if args.finite_records > 0 else len(exp_paths)
    full_freq = np.fft.rfftfreq(sample, d=1.0 / rate)
    best_pre_analysis = hardware_sampled_asd(
        best["parameters"],
        full_freq,
        rate_hz=rate,
        cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
    )
    best_finite_asd = finite_record_post_analysis_asd(
        best_pre_analysis,
        sample,
        rate,
        analysis_cutoff_hz=analysis_cutoff,
        records=finite_records,
        seed=args.finite_seed,
    )
    best_finite_norm = normalized(best_finite_asd, full_freq)
    best_plot = log_interp(full_freq[1:], best_finite_norm[1:], plot_freq)
    best_eval = log_interp(full_freq[1:], best_finite_norm[1:], eval_freq)
    finite_log_ratio = np.log(best_eval / exp_eval)
    finite_rms = float(np.sqrt(np.mean(finite_log_ratio**2)))
    finite_max = float(np.max(np.abs(finite_log_ratio)))

    anchor_values = {}
    for hz in ANCHORS:
        query = np.asarray([float(hz)])
        anchor_values[str(hz)] = {
            "simulation_over_experiment": float(
                log_interp(full_freq[1:], best_finite_norm[1:], query)[0]
                / log_interp(exp_freq[1:], exp_norm[1:], query)[0]
            ),
            "simulation_finite_record_normalized": float(
                log_interp(full_freq[1:], best_finite_norm[1:], query)[0]
            ),
            "simulation_expected_normalized": float(
                log_interp(model_freq, best_expected, query)[0]
            ),
            "experiment_normalized": float(
                log_interp(exp_freq[1:], exp_norm[1:], query)[0]
            ),
        }

    ensemble = np.asarray(
        [
            log_interp(model_freq, curves[row["scenario_id"]], plot_freq)
            for row in rows
        ]
    )
    q05, median, q95 = (
        np.quantile(ensemble, q, axis=0) for q in (0.05, 0.5, 0.95)
    )
    inside = (exp_plot >= q05) & (exp_plot <= q95)
    scoring_mask = (plot_freq >= 1000.0) & (plot_freq <= 10000.0)
    inside_scoring = inside & scoring_mask

    comparison = {
        "stage": "high_frequency_best_simulation_comparison",
        "best_scenario_id": best["scenario_id"],
        "selection_band_Hz": [1000.0, 10000.0],
        "low_frequency_excluded_from_selection": True,
        "display_band_Hz": [float(DISPLAY_MIN_HZ), float(DISPLAY_MAX_HZ)],
        "display_excludes_DC": True,
        "normalization_frequency_Hz": 1000.0,
        "selection_metric_definition": (
            "deterministic expected post-analysis RMS log ratio on 200 log-spaced "
            "frequencies; sqrt(mean(log(ASD_sim/ASD_exp)^2))"
        ),
        "selection_rms_log_ratio": best["rms_log_ratio"],
        "selection_max_abs_log_ratio": best["max_abs_log_ratio"],
        "finite_record_rms_log_ratio": finite_rms,
        "finite_record_max_abs_log_ratio": finite_max,
        "evaluation_frequency_count": len(eval_freq),
        "evaluation_frequencies_Hz": eval_freq.tolist(),
        "anchor_values": anchor_values,
        "best_scenario_parameters": best["parameters"],
        "parameter_source_classification": best["parameter_source_classification"],
        "best_scenario_label": "best_pre_generated_frozen_scenario",
        "interpretation_label": (
            "best high-frequency matching member of the frozen proxy ensemble"
        ),
        "ensemble_size": int(
            scenarios_doc.get("sample_count", len(scenarios_doc.get("scenarios", [])))
        ),
        "pulse_consistent_ensemble_size": len(scenarios),
        "accepted_CH0_record_count": len(exp_paths),
        "simulation_finite_record_count": finite_records,
        "finite_record_seed": int(args.finite_seed),
        "hardware_bessel_cutoff_Hz": float(HARDWARE_BESSEL_CUTOFF_HZ),
        "analysis_bessel_cutoff_Hz": float(analysis_cutoff),
        "hardware_bessel_semantics": (
            "physical analog Bessel before ADC; applied to simulation only because "
            "the experimental raw records already passed through the hardware"
        ),
        "analysis_bessel_semantics": (
            "second-order digital Bessel filtfilt applied in time domain to both "
            "experimental and simulated records before Hann/FFT ASD estimation"
        ),
        "experimental_processing": (
            "accepted raw CH0 -> mean removal -> 10 kHz digital Bessel filtfilt "
            "-> Hann -> rFFT -> |FFT|^2 -> record power average -> one-sided ASD"
        ),
        "simulation_processing": (
            "intrinsic TES ASD -> 100 kHz analog hardware Bessel + first ADC alias "
            "fold -> Gaussian finite time records -> mean removal -> 10 kHz digital "
            "Bessel filtfilt -> Hann -> rFFT -> |FFT|^2 -> record power average "
            "-> one-sided ASD"
        ),
        "simulation_time_domain_realization": True,
        "experimental_subset": (
            "all accepted CH0 records; pulse_free_candidate mask not used"
        ),
        "simulation_source": (
            "proxy_scenarios.json:pulse_consistent_scenarios; frozen before comparison"
        ),
        "noise_residual_fit": False,
        "parameter_optimization": False,
        "simulation_amplitude_rescale": False,
        "strict_target_conclusion": STRICT,
        "ensemble_curve_semantics": (
            "deterministic expected post-analysis spectra; finite-record realization "
            "is generated only for the selected frozen scenario"
        ),
        "ensemble_curve_frequency_Hz": plot_freq.tolist(),
        "ensemble_normalized_asd": {
            "min": np.min(ensemble, axis=0).tolist(),
            "q05": q05.tolist(),
            "median": median.tolist(),
            "q95": q95.tolist(),
            "max": np.max(ensemble, axis=0).tolist(),
        },
        "experiment_normalized_asd": {
            "frequency_Hz": plot_freq.tolist(),
            "values": exp_plot.tolist(),
        },
        "best_finite_record_normalized_asd": {
            "frequency_Hz": plot_freq.tolist(),
            "values": best_plot.tolist(),
        },
        "best_expected_normalized_asd": {
            "frequency_Hz": plot_freq.tolist(),
            "values": best_expected_plot.tolist(),
        },
        "experiment_inside_ensemble_q05_q95_intervals_Hz": intervals(
            plot_freq, inside
        ),
        "experiment_inside_ensemble_q05_q95_intervals_in_selection_band_Hz": intervals(
            plot_freq, inside_scoring
        ),
        "experiment_inside_ensemble_q05_q95_fraction_in_selection_band": float(
            np.mean(inside_scoring[scoring_mask])
        ),
        "all_scenario_scores": sorted(
            rows, key=lambda row: (row["rms_log_ratio"], row["scenario_id"])
        ),
    }
    dump(args.case_dir / "high_frequency_best_simulation_comparison.json", comparison)

    import matplotlib.pyplot as plt

    args.case_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 5.5))
    plt.fill_between(
        plot_freq,
        np.min(ensemble, axis=0),
        np.max(ensemble, axis=0),
        color="C1",
        alpha=0.10,
        label="frozen ensemble expected min–max",
    )
    plt.fill_between(
        plot_freq,
        q05,
        q95,
        color="C1",
        alpha=0.25,
        label="frozen pulse-consistent expected q05–q95",
    )
    plt.plot(
        plot_freq,
        median,
        color="C1",
        lw=1.5,
        label="frozen ensemble expected median",
    )
    plt.plot(
        plot_freq,
        exp_plot,
        color="black",
        lw=1.5,
        label=f"Experiment CH0 post-analysis ({len(exp_paths)} accepted records)",
    )
    plt.plot(
        plot_freq,
        best_expected_plot,
        color="C0",
        lw=1.0,
        ls="--",
        alpha=0.75,
        label=f"{best['scenario_id']} expected post-analysis",
    )
    plt.plot(
        plot_freq,
        best_plot,
        color="C0",
        lw=1.2,
        label=f"{best['scenario_id']} finite-record simulation",
    )
    plt.scatter(
        [1000],
        [1],
        color="black",
        s=28,
        zorder=5,
        label="1 kHz normalization",
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(DISPLAY_MIN_HZ, DISPLAY_MAX_HZ)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("TES CH0 Noise: Experiment vs Frozen Intrinsic Simulation")
    plt.suptitle(
        "100 kHz hardware Bessel + finite records + 10 kHz software Bessel; "
        "selection uses only 1–10 kHz shape error.\n"
        "No noise-residual parameter fitting.",
        fontsize=9,
        y=0.95,
    )
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8, loc="best")
    plt.tight_layout(rect=(0, 0, 0.99, 0.90))
    plt.savefig(args.case_dir / "high_frequency_best_simulation_loglog.png", dpi=180)
    plt.close()

    ratio = best_plot / exp_plot
    plt.figure(figsize=(9, 4.5))
    plt.axhspan(0.8, 1.2, color="gray", alpha=0.10, label="±20%")
    plt.axhspan(0.9, 1.1, color="gray", alpha=0.16, label="±10%")
    plt.axhline(1, color="black", lw=1)
    plt.plot(
        plot_freq,
        ratio,
        color="C0",
        lw=1.2,
        label=f"{best['scenario_id']} finite / experiment",
    )
    plt.xscale("log")
    plt.xlim(DISPLAY_MIN_HZ, DISPLAY_MAX_HZ)
    plt.ylim(
        max(0.1, float(np.nanmin(ratio) * 0.8)),
        min(10, float(np.nanmax(ratio) * 1.2)),
    )
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Simulation / Experiment ASD")
    plt.title("Full-band finite-record shape ratio (selection used 1–10 kHz only)")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.case_dir / "high_frequency_best_simulation_ratio.png", dpi=180)
    plt.close()

    intrinsic_components, intrinsic_meta = noise_components(
        best["parameters"], model_freq
    )
    intrinsic_total = np.asarray(intrinsic_meta["total_asd"], dtype=float)
    norm_1k = intrinsic_total[np.argmin(np.abs(model_freq - 1000.0))].item()
    plt.figure(figsize=(9, 5.5))
    for name, indices in component_labels().items():
        component = np.sqrt(
            np.sum(np.asarray(intrinsic_components[:, indices]) ** 2, axis=1)
        )
        plt.plot(model_freq, component / norm_1k, lw=1.5, label=name)
    plt.plot(
        model_freq,
        intrinsic_total / norm_1k,
        color="black",
        lw=2,
        label="total intrinsic",
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(DISPLAY_MIN_HZ, DISPLAY_MAX_HZ)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Intrinsic ASD / total intrinsic ASD at 1 kHz")
    plt.title(f"Frozen intrinsic noise components before measurement chain — {best['scenario_id']}")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.case_dir / "high_frequency_best_simulation_components.png", dpi=180)
    plt.close()


if __name__ == "__main__":
    main()

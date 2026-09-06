"""High-frequency experiment vs frozen proxy-ensemble comparison.

This is a descriptive shape comparison only.  It consumes the already frozen
``pulse_consistent_scenarios`` and never changes parameters or fits a noise
residual.  The experiment is the complete accepted CH0 noise set; no
pulse-free mask is used.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "PoST_Simulations"))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))
from lib.tes_noise_model import SOURCE_CLASS_INDICES  # noqa: E402
from proxy_physics import noise_components  # noqa: E402
from pulse_contamination_common import read_record  # noqa: E402
from pulse_contamination_v5 import RATE_HZ, SAMPLES, accepted_noise  # noqa: E402
from Analyze_Experimental_Data.tes_analysis.noise_utils import one_sided_asd_from_power  # noqa: E402

STRICT = "C — exact target physical case remains unidentified"
ANCHORS = (1000, 2000, 3000, 5000, 7000, 10000)
CLASSES = ("TES_Johnson", "load_Johnson", "TES_bath_TFN", "TES_absorber_TFN")
# The comparison score is intentionally restricted to 1--10 kHz.  Figures,
# however, show the complete positive-frequency range available in the
# acquisition (the DC bin cannot be shown on a logarithmic axis).
DISPLAY_MIN_HZ = RATE_HZ / SAMPLES
DISPLAY_MAX_HZ = RATE_HZ / 2.0


def dump(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def log_interp(x: np.ndarray, y: np.ndarray, query: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float); query = np.asarray(query, dtype=float)
    good = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    return np.exp(np.interp(np.log(query), np.log(x[good]), np.log(y[good])))


def experiment_asd(target: Path) -> tuple[np.ndarray, np.ndarray, list[Path]]:
    paths = accepted_noise(target, "CH0")
    if not paths:
        raise RuntimeError("no accepted CH0 noise records")
    window = np.hanning(SAMPLES)
    total_power = np.zeros(SAMPLES // 2 + 1, dtype=float)
    for path in paths:
        raw = read_record(path)
        total_power += np.abs(np.fft.rfft((raw - np.mean(raw)) * window)) ** 2
    total_power /= len(paths)
    asd = one_sided_asd_from_power(total_power, SAMPLES, RATE_HZ, np.sqrt(np.mean(window ** 2)))
    freq = np.fft.rfftfreq(SAMPLES, 1 / RATE_HZ)
    return freq, asd, paths


def component_labels() -> dict[str, tuple[int, ...]]:
    return {name: tuple(indices) for name, indices in SOURCE_CLASS_INDICES.items()}


def intervals(freq: np.ndarray, mask: np.ndarray) -> list[list[float]]:
    values = np.asarray(freq)[np.asarray(mask, dtype=bool)]
    if not len(values): return []
    groups = np.split(values, np.where(np.diff(np.where(mask)[0]) > 1)[0] + 1)
    return [[float(group[0]), float(group[-1])] for group in groups if len(group)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--scenario-file", type=Path)
    args = ap.parse_args()
    args.scenario_file = args.scenario_file or args.case_dir / "proxy_scenarios.json"
    scenarios_doc = json.loads(args.scenario_file.read_text(encoding="utf-8"))
    scenarios = sorted(scenarios_doc["pulse_consistent_scenarios"], key=lambda row: row["scenario_id"])
    if scenarios_doc.get("freeze_status") != "frozen" or scenarios_doc.get("pulse_consistent_count") != len(scenarios):
        raise RuntimeError("comparison requires the frozen pulse-consistent scenario list")
    exp_freq, exp_asd, exp_paths = experiment_asd(args.target_root)
    norm_i = int(np.argmin(np.abs(exp_freq - 1000.0)))
    exp_norm = exp_asd / exp_asd[norm_i]
    eval_freq = np.logspace(np.log10(1000.0), np.log10(10000.0), 200)
    # Display the complete positive-frequency target range.  The 1--10 kHz
    # band remains the only band used for scenario selection.  DC is omitted
    # because both axes are logarithmic.
    plot_freq = np.logspace(np.log10(DISPLAY_MIN_HZ), np.log10(DISPLAY_MAX_HZ), 1200)
    model_freq = np.unique(np.r_[1000.0, eval_freq, plot_freq, np.asarray(ANCHORS, dtype=float)])
    exp_eval = log_interp(exp_freq[1:], exp_norm[1:], eval_freq)
    exp_plot = log_interp(exp_freq[1:], exp_norm[1:], plot_freq)
    rows = []
    curves = {}
    for scenario in scenarios:
        components, meta = noise_components(scenario["parameters"], model_freq)
        total = np.asarray(meta["total_asd"], dtype=float)
        total_norm = total / total[np.argmin(np.abs(model_freq - 1000.0))]
        sim_eval = log_interp(model_freq, total_norm, eval_freq)
        log_ratio = np.log(sim_eval / exp_eval)
        score = float(np.sqrt(np.mean(log_ratio ** 2)))
        rows.append({"scenario_id": scenario["scenario_id"], "rms_log_ratio": score, "max_abs_log_ratio": float(np.max(np.abs(log_ratio))), "parameter_source_classification": scenario.get("source_class_by_parameter", {}), "parameters": scenario["parameters"], "pulse_consistency": scenario.get("pulse_consistency")})
        norm_1k = total[np.argmin(np.abs(model_freq - 1000.0))].item()
        curves[scenario["scenario_id"]] = {"total": total_norm, "components": {name: np.sqrt(np.sum(np.asarray(components[:, indices]) ** 2, axis=1)) / norm_1k for name, indices in component_labels().items()}}
    best = min(rows, key=lambda row: (row["rms_log_ratio"], row["scenario_id"]))
    best_curve = curves[best["scenario_id"]]["total"]
    best_plot = log_interp(model_freq, best_curve, plot_freq)
    best_eval = log_interp(model_freq, best_curve, eval_freq)
    best_anchor = log_interp(model_freq, best_curve, np.asarray(ANCHORS, dtype=float))
    anchor_values = {str(hz): {"simulation_over_experiment": float(log_interp(model_freq, best_curve, np.asarray([hz]))[0] / log_interp(exp_freq[1:], exp_norm[1:], np.asarray([hz]))[0]), "simulation_normalized": float(best_anchor[i]), "experiment_normalized": float(log_interp(exp_freq[1:], exp_norm[1:], np.asarray([hz]))[0])} for i, hz in enumerate(ANCHORS)}
    ensemble = np.asarray([log_interp(model_freq, curves[row["scenario_id"]]["total"], plot_freq) for row in rows])
    q05, median, q95 = (np.quantile(ensemble, q, axis=0) for q in (.05, .5, .95))
    inside = (exp_plot >= q05) & (exp_plot <= q95)
    scoring_mask = (plot_freq >= 1000.0) & (plot_freq <= 10000.0)
    inside_scoring = inside & scoring_mask
    comparison = {
        "stage": "high_frequency_best_simulation_comparison",
        "best_scenario_id": best["scenario_id"],
        "selection_band_Hz": [1000.0, 10000.0],
        "display_band_Hz": [float(DISPLAY_MIN_HZ), float(DISPLAY_MAX_HZ)],
        "display_excludes_DC": True,
        "normalization_frequency_Hz": 1000.0,
        "metric_definition": "RMS log ratio on 200 log-spaced frequencies; sqrt(mean(log(ASD_sim/ASD_exp)^2))",
        "evaluation_frequency_count": len(eval_freq),
        "evaluation_frequencies_Hz": eval_freq.tolist(),
        "rms_log_ratio": best["rms_log_ratio"],
        "max_abs_log_ratio": best["max_abs_log_ratio"],
        "anchor_values": anchor_values,
        "best_scenario_parameters": best["parameters"],
        "parameter_source_classification": best["parameter_source_classification"],
        "best_scenario_label": "best_pre_generated_frozen_scenario",
        "interpretation_label": "best high-frequency matching member of the frozen proxy ensemble",
        "ensemble_size": int(scenarios_doc.get("sample_count", len(scenarios_doc.get("scenarios", [])))),
        "pulse_consistent_ensemble_size": len(scenarios),
        "accepted_CH0_record_count": len(exp_paths),
        "experimental_processing": "raw -> mean removal -> Hann -> rFFT -> |FFT|^2 -> record power average -> one-sided ASD",
        "pre_analysis_has_bessel": False,
        "experimental_subset": "all accepted CH0 records; pulse_free_candidate mask not used",
        "simulation_source": "proxy_scenarios.json:pulse_consistent_scenarios; frozen before comparison",
        "noise_residual_fit": False,
        "parameter_optimization": False,
        "simulation_amplitude_rescale": False,
        "strict_target_conclusion": STRICT,
        "ensemble_curve_frequency_Hz": plot_freq.tolist(),
        "ensemble_normalized_asd": {"min": np.min(ensemble, axis=0).tolist(), "q05": q05.tolist(), "median": median.tolist(), "q95": q95.tolist(), "max": np.max(ensemble, axis=0).tolist()},
        "experiment_normalized_asd": {"frequency_Hz": plot_freq.tolist(), "values": exp_plot.tolist()},
        "experiment_inside_ensemble_q05_q95_intervals_Hz": intervals(plot_freq, inside),
        "experiment_inside_ensemble_q05_q95_intervals_in_selection_band_Hz": intervals(plot_freq, inside_scoring),
        "experiment_inside_ensemble_q05_q95_fraction_in_selection_band": float(np.mean(inside_scoring[scoring_mask])),
        "all_scenario_scores": sorted(rows, key=lambda row: (row["rms_log_ratio"], row["scenario_id"])),
    }
    dump(args.case_dir / "high_frequency_best_simulation_comparison.json", comparison)
    import matplotlib.pyplot as plt
    args.case_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 5.5))
    plt.fill_between(plot_freq, np.min(ensemble, axis=0), np.max(ensemble, axis=0), color="C1", alpha=.10, label="frozen ensemble min–max")
    plt.fill_between(plot_freq, q05, q95, color="C1", alpha=.25, label="frozen pulse-consistent ensemble q05–q95")
    plt.plot(plot_freq, median, color="C1", lw=1.5, label="frozen ensemble median")
    plt.plot(plot_freq, exp_plot, color="black", lw=1.7, label=f"Experiment CH0 ({len(exp_paths)} accepted records)")
    plt.plot(plot_freq, best_plot, color="C0", lw=2.0, label=f"{best['scenario_id']} — best_pre-generated_frozen_scenario")
    plt.scatter([1000], [1], color="black", s=28, zorder=5, label="1 kHz normalization")
    plt.xscale("log"); plt.yscale("log"); plt.xlim(DISPLAY_MIN_HZ, DISPLAY_MAX_HZ); plt.xlabel("Frequency [Hz]"); plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("TES CH0 Noise: Experiment vs Best Frozen Intrinsic Simulation")
    plt.suptitle("Best scenario selected only from the pre-generated noise-blind frozen proxy ensemble using 1–10 kHz shape error.\nNo noise-residual parameter fitting.", fontsize=9, y=.94)
    plt.grid(True, which="both", alpha=.25); plt.legend(fontsize=8, loc="best"); plt.tight_layout(rect=(0,0,.99,.90)); plt.savefig(args.case_dir / "high_frequency_best_simulation_loglog.png", dpi=180); plt.close()
    ratio = best_plot / exp_plot
    plt.figure(figsize=(9, 4.5)); plt.axhspan(.8, 1.2, color="gray", alpha=.10, label="±20%"); plt.axhspan(.9, 1.1, color="gray", alpha=.16, label="±10%"); plt.axhline(1, color="black", lw=1); plt.plot(plot_freq, ratio, color="C0", lw=2, label=f"{best['scenario_id']} / experiment"); plt.xscale("log"); plt.xlim(DISPLAY_MIN_HZ, DISPLAY_MAX_HZ); plt.ylim(max(.1,float(np.nanmin(ratio)*.8)), min(10,float(np.nanmax(ratio)*1.2))); plt.xlabel("Frequency [Hz]"); plt.ylabel("Simulation / Experiment ASD"); plt.title("Full-band shape ratio (selection used 1–10 kHz only)"); plt.grid(True, which="both", alpha=.25); plt.legend(); plt.tight_layout(); plt.savefig(args.case_dir / "high_frequency_best_simulation_ratio.png", dpi=180); plt.close()
    components = curves[best["scenario_id"]]["components"]
    plt.figure(figsize=(9, 5.5));
    for name in CLASSES: plt.plot(model_freq, components[name], lw=1.5, label=name)
    plt.plot(model_freq, best_curve, color="black", lw=2, label="total intrinsic")
    plt.xscale("log"); plt.yscale("log"); plt.xlim(DISPLAY_MIN_HZ, DISPLAY_MAX_HZ); plt.xlabel("Frequency [Hz]"); plt.ylabel("ASD / best scenario total ASD at 1 kHz"); plt.title(f"Frozen intrinsic noise components — {best['scenario_id']}"); plt.grid(True, which="both", alpha=.25); plt.legend(fontsize=8); plt.tight_layout(); plt.savefig(args.case_dir / "high_frequency_best_simulation_components.png", dpi=180); plt.close()


if __name__ == "__main__": main()

"""Stage B model runner: consume only the frozen Stage-A scenarios.

This script still does not read experimental data.  It produces physical
source ensembles, performs stability-first exclusion, and records normalized
pre-analysis model spectra.  Experimental comparison is a separate script.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from proxy_physics import SOURCE_NAMES, noise_components, operating_point


FREQUENCIES = np.array([10.0, 100.0, 1000.0, 3000.0, 5000.0, 7000.0, 10000.0])
SENSITIVITY_PARAMETERS = ("T_c", "T_bath", "R", "R_l", "alpha", "beta", "L", "n", "C_tes", "C_abs", "G_tes-bath", "G_abs-tes", "G_abs-abs")
BANDS = {"10-100_Hz": (10.0, 100.0), "100-1000_Hz": (100.0, 1000.0), "1-3_kHz": (1000.0, 3000.0), "3-10_kHz": (3000.0, 10000.0)}


def _dump(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _normalized(values: np.ndarray) -> np.ndarray:
    reference = values[2]
    return values / reference if reference > 0.0 else np.full_like(values, np.nan)


def _sensitivity(baseline: dict, baseline_source: str) -> dict:
    rows = {}
    epsilon = 0.01
    for name in SENSITIVITY_PARAMETERS:
        base = float(baseline[name])
        low = dict(baseline)
        high = dict(baseline)
        low[name] = base * (1.0 - epsilon)
        high[name] = base * (1.0 + epsilon)
        try:
            low_components, low_meta = noise_components(low, FREQUENCIES)
            high_components, high_meta = noise_components(high, FREQUENCIES)
            low_total = _normalized(low_meta["total_asd"])
            high_total = _normalized(high_meta["total_asd"])
            derivative = (np.log(high_total) - np.log(low_total)) / (np.log(1.0 + epsilon) - np.log(1.0 - epsilon))
            band_values = {}
            for band, (left, right) in BANDS.items():
                mask = (FREQUENCIES >= left) & (FREQUENCIES <= right)
                band_values[band] = float(np.mean(np.abs(derivative[mask])))
            rows[name] = {"source_class": "simulation_reference_only", "baseline": base, "fractional_step": epsilon, "d_ln_ASD_d_ln_parameter": [float(x) for x in derivative], "frequencies_Hz": FREQUENCIES.tolist(), "band_mean_absolute_sensitivity": band_values, "normalization": "each spectrum normalized at 1 kHz", "stable_plus_minus": True}
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            rows[name] = {"source_class": "simulation_reference_only", "baseline": base, "fractional_step": epsilon, "d_ln_ASD_d_ln_parameter": None, "frequencies_Hz": FREQUENCIES.tolist(), "band_mean_absolute_sensitivity": None, "stable_plus_minus": False, "reason": str(exc)}
    return {"stage": "B_model_sensitivity", "experimental_spectrum_read": False, "baseline_source": baseline_source, "method": "central finite difference around pulse-consistent ensemble median; no optimizer", "frequencies_Hz": FREQUENCIES.tolist(), "parameters": rows}


def _plot(output_dir: Path, stable_rows: list[dict], sensitivity: dict) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    if stable_rows:
        values = np.asarray([row["normalized_total_asd"] for row in stable_rows])
        plt.figure(figsize=(7, 4))
        plt.fill_between(FREQUENCIES, np.min(values, axis=0), np.max(values, axis=0), alpha=0.25, label="stable scenario min/max")
        plt.plot(FREQUENCIES, np.median(values, axis=0), label="sampled median")
        plt.xscale("log")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Normalized model ASD")
        plt.grid(True, which="both")
        plt.legend()
        plt.tight_layout()
        plt.savefig(plot_dir / "proxy_noise_envelope_pre_analysis.png", dpi=160)
        plt.close()
        source_values = np.asarray([row["normalized_components"] for row in stable_rows])
        plt.figure(figsize=(7, 4))
        for index, name in enumerate(SOURCE_NAMES):
            plt.plot(FREQUENCIES, np.median(source_values[:, :, index], axis=0), label=name)
        plt.xscale("log")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Normalized source ASD")
        plt.grid(True, which="both")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(plot_dir / "source_decomposition.png", dpi=160)
        plt.close()
    names = list(sensitivity["parameters"])
    matrix = np.array([[sensitivity["parameters"][name]["band_mean_absolute_sensitivity"][band] if sensitivity["parameters"][name]["band_mean_absolute_sensitivity"] else np.nan for band in BANDS] for name in names])
    plt.figure(figsize=(9, 5))
    plt.imshow(matrix, aspect="auto", interpolation="nearest")
    plt.colorbar(label="mean |d ln ASD / d ln p|")
    plt.xticks(range(len(BANDS)), list(BANDS), rotation=35, ha="right")
    plt.yticks(range(len(names)), names)
    plt.tight_layout()
    plt.savefig(plot_dir / "parameter_sensitivity.png", dpi=160)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-a-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stage_a_manifest = json.loads((args.stage_a_dir / "noise_blind_sweep_manifest.json").read_text(encoding="utf-8"))
    envelope = json.loads((args.stage_a_dir / "proxy_parameter_envelope.json").read_text(encoding="utf-8"))
    scenarios = json.loads((args.stage_a_dir / "proxy_scenarios.json").read_text(encoding="utf-8"))
    generic = envelope["sensitivity_reference"]
    if stage_a_manifest.get("freeze_status") != "frozen" or not scenarios.get("freeze_status") == "frozen":
        raise RuntimeError("Stage-A parameter ranges are not frozen")
    candidate_scenarios = scenarios.get("pulse_consistent_scenarios", scenarios["scenarios"])
    stable_rows = []
    excluded_rows = []
    for scenario in candidate_scenarios:
        params = scenario["parameters"]
        try:
            point = operating_point(params)
        except Exception as exc:
            point = {"stable": False, "reason": f"operating_point_exception:{exc}"}
        if not point.get("stable", False):
            excluded_rows.append({"scenario_id": scenario["scenario_id"], "reason": point.get("reason", "unstable"), "excluded_before_noise": True})
            continue
        try:
            components, meta = noise_components(params, FREQUENCIES)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            excluded_rows.append({"scenario_id": scenario["scenario_id"], "reason": f"noise_matrix_failed:{exc}", "excluded_before_noise": False})
            continue
        source_scale = meta["total_asd"][2]
        stable_rows.append({"scenario_id": scenario["scenario_id"], "scenario_class": scenario["scenario_class"], "source_class_by_parameter": scenario["source_class_by_parameter"], "pulse_consistency": scenario.get("pulse_consistency"), "operating_point": meta["operating_point"], "frequencies_Hz": FREQUENCIES.tolist(), "total_asd": [float(x) for x in meta["total_asd"]], "normalized_total_asd": [float(x) for x in _normalized(meta["total_asd"])], "normalized_components": [[float(x) for x in row] for row in (components / source_scale)], "normalized_source_class_components": {name: [float(x) for x in values / source_scale] for name, values in meta["source_class_components"].items()}, "source_names": list(SOURCE_NAMES)})
    values = np.asarray([row["normalized_total_asd"] for row in stable_rows], dtype=float) if stable_rows else np.empty((0, len(FREQUENCIES)))
    component_values = np.asarray([row["normalized_components"] for row in stable_rows], dtype=float) if stable_rows else np.empty((0, len(FREQUENCIES), len(SOURCE_NAMES)))
    envelope = {
        "stage": "B_model_ensemble",
        "input_stage": "A_frozen_proxy_scenarios_only",
        "experimental_spectrum_read": False,
        "comparison_performed": False,
        "frequencies_Hz": FREQUENCIES.tolist(),
        "sampled_quantiles_are_not_probabilities": True,
        "scenario_selection": "pulse_consistent_scenarios_only; reference_sensitivity_scenarios are not included in this physical envelope",
        "scenario_count": len(candidate_scenarios),
        "stable_count": len(stable_rows),
        "excluded_count": len(excluded_rows),
        "stability_exclusion_policy": "Exclude only nonphysical or linearly unstable scenarios before noise calculation; never exclude by experimental mismatch.",
        "excluded_scenario_schema": ["scenario_id", "reason", "excluded_before_noise"],
        "normalized_total_asd": {"min": np.min(values, axis=0).tolist() if len(values) else None, "q05": np.quantile(values, 0.05, axis=0).tolist() if len(values) else None, "q50": np.quantile(values, 0.50, axis=0).tolist() if len(values) else None, "q95": np.quantile(values, 0.95, axis=0).tolist() if len(values) else None, "max": np.max(values, axis=0).tolist() if len(values) else None},
        "normalized_source_components": {name: {"min": np.min(component_values[:, :, i], axis=0).tolist() if len(values) else None, "q50": np.quantile(component_values[:, :, i], 0.50, axis=0).tolist() if len(values) else None, "max": np.max(component_values[:, :, i], axis=0).tolist() if len(values) else None} for i, name in enumerate(SOURCE_NAMES)},
        "source_class_aggregation": "PSD aggregation after eight independent source transfer ASDs; classes are TES_Johnson, load_Johnson, TES_bath_TFN, TES_absorber_TFN",
        "stable_scenarios": stable_rows,
        "excluded_scenarios": excluded_rows,
        "source_decomposition": list(SOURCE_NAMES),
        "empirical_white_noise_A_rtHz": 0.0,
        "empirical_readout_floor_A_rtHz": 0.0,
        "resistance_fluctuation_model": "none",
        "target_hanging_model": "disabled",
        "strict_target_conclusion": "C — exact target physical case remains unidentified",
    }
    _dump(args.output_dir / "proxy_noise_envelope.json", envelope)
    if candidate_scenarios:
        baseline = {name: float(np.median([row["parameters"][name] for row in candidate_scenarios])) for name in SENSITIVITY_PARAMETERS}
        baseline.update({"rate": 500000.0, "samples": 100000})
        baseline_source = "median_of_pulse_consistent_scenarios"
    else:
        baseline = dict(generic)
        baseline_source = "generic_reference_fallback_because_pulse_gate_empty"
    sensitivity = _sensitivity(baseline, baseline_source)
    _dump(args.output_dir / "sensitivity_summary.json", sensitivity)
    reference_sensitivity = _sensitivity(generic, "generic_simulation_reference_only")
    _dump(args.output_dir / "reference_sensitivity_summary.json", reference_sensitivity)
    markdown = ["# Noise-blind target-like sensitivity", "", f"This is a deterministic sensitivity study around `{baseline_source}`. It is not a target parameter estimate and uses no experimental spectrum.", "A separate generic 0.5x–2x reference sensitivity is stored in `reference_sensitivity_summary.json` and is not used for the physical envelope.", "", "| parameter | 10-100 Hz | 100-1000 Hz | 1-3 kHz | 3-10 kHz |", "|---|---:|---:|---:|---:|"]
    for name, row in sensitivity["parameters"].items():
        vals = row["band_mean_absolute_sensitivity"]
        markdown.append(f"| {name} | " + " | ".join(f"{vals[band]:.4g}" if vals else "n/a" for band in BANDS) + " |")
    markdown.extend(["", "Interpretation: values are finite-difference shape sensitivities around the generic nominal reference. They do not establish target provenance."])
    (args.output_dir / "sensitivity_summary.md").write_text("\n".join(markdown) + "\n", encoding="utf-8")
    _plot(args.output_dir, stable_rows, sensitivity)


if __name__ == "__main__":
    main()

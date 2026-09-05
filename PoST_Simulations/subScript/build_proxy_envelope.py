"""Stage A: construct and freeze non-noise proxy physics inputs.

This script never opens the target experimental spectrum.  It reads RT, IV,
pulse, configuration, and generic simulation-reference files only.  Its
outputs are immutable inputs for Stage B; no comparison code is imported.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np


CASE_NAME = "tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"
TARGET_PATH = Path(r"G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2")
RT_FILES = [
    Path(r"G:/tagawa/20241202/room1-ch1-rt/output/RT_10uA.csv"),
    Path(r"G:/tagawa/20241202/room1-ch1-rt/output/RT_20uA.csv"),
]
IV_FILES = sorted(Path(r"G:/tagawa/20241203/room1-ch1-iv/calibration").glob("IV_*mK.txt"))
TARGET_IV = Path(r"G:/tagawa/20241206/room1-ch1-iv3/calibration/IV_215mK.txt")
GENERIC_INPUT = Path(__file__).resolve().parents[1] / "input.json"


def _dump(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rt_summary() -> dict:
    rows = []
    for path in RT_FILES:
        data = np.array([(float(row["T"]), float(row["R"])) for row in csv.DictReader(path.open(newline=""))])
        rn = float(np.median(data[data[:, 0] >= np.percentile(data[:, 0], 80), 1]))
        crossings = {}
        for fraction in (0.1, 0.5, 0.9):
            index = int(np.where(data[:, 1] >= fraction * rn)[0][0])
            if index == 0:
                temperature = float(data[0, 0])
            else:
                x0, y0 = data[index - 1]
                x1, y1 = data[index]
                temperature = float(x0 + (fraction * rn - y0) * (x1 - x0) / (y1 - y0))
            crossings[str(fraction)] = temperature / 1000.0
        rows.append({
            "source_file": path.as_posix(),
            "bias_label_uA": int(path.stem.split("_")[1].removesuffix("uA")),
            "Rn_proxy_ohm": rn / 1000.0,
            "crossings_K": crossings,
            "transition_width_10_to_90_K": crossings["0.9"] - crossings["0.1"],
        })
    midpoint = [row["crossings_K"]["0.5"] for row in rows]
    onset = [row["crossings_K"]["0.1"] for row in rows]
    high = [row["crossings_K"]["0.9"] for row in rows]
    return {
        "source_class": "nearby_run_proxy",
        "definition": "Rn is median plateau above the 80th temperature percentile; crossings are linear interpolation of R/Rn.",
        "rows": rows,
        "T_c_range_K": [min(onset), max(high)],
        "T_c_nominal_K": float(np.median(midpoint)),
        "method_dependent_definitions": {
            "onset_R_over_Rn_0.1_K": [min(onset), max(onset)],
            "midpoint_R_over_Rn_0.5_K": [min(midpoint), max(midpoint)],
            "high_R_over_Rn_0.9_K": [min(high), max(high)],
        },
        "not_selected_using_spectrum": True,
    }


def _iv_summary() -> dict:
    rows = []
    for path in IV_FILES:
        data = np.loadtxt(path)
        ibias, vout = data[:2]
        eta = float(1.0 / np.polyfit(ibias[:10], vout[:10], 1)[0])
        per_shunt = []
        for r_sh in (0.0038, 0.0039):
            ites = eta * vout
            valid = np.isfinite(ites) & (ites > 0.0) & (ibias > ites)
            resistance = np.full_like(ites, np.nan, dtype=float)
            resistance[valid] = r_sh * (ibias[valid] - ites[valid]) / ites[valid]
            usable = np.isfinite(resistance) & (resistance > 0.0) & (resistance < 0.2)
            representative = None
            if np.count_nonzero(usable) >= 10:
                rn_proxy = float(np.median(np.sort(resistance[usable])[-10:]))
                midpoint_indices = np.where(usable & (resistance >= 0.5 * rn_proxy))[0]
                if len(midpoint_indices):
                    midpoint = int(midpoint_indices[0])
                    vtes = float((ibias[midpoint] - ites[midpoint]) * 1e-6 * r_sh)
                    representative = {
                        "I_bias_uA": float(ibias[midpoint]),
                        "I_TES_candidate_uA": float(ites[midpoint]),
                        "V_TES_V": vtes,
                        "R_TES_ohm": float(resistance[midpoint]),
                        "P_J_W": float(vtes * ites[midpoint] * 1e-6),
                        "selection": "first usable point at or above 0.5*R_n_proxy",
                    }
            per_shunt.append({
                "R_SH_ohm": r_sh,
                "usable_points": int(np.count_nonzero(usable)),
                "R_n_proxy_ohm": rn_proxy if np.count_nonzero(usable) >= 10 else None,
                "R_min_ohm": float(np.nanmin(resistance[usable])) if np.any(usable) else None,
                "R_max_ohm": float(np.nanmax(resistance[usable])) if np.any(usable) else None,
                "representative_midpoint": representative,
                "data_quality": "usable_for_context_only" if np.count_nonzero(usable) >= 10 else "insufficient_or_nonmonotonic",
            })
        rows.append({"source_file": path.as_posix(), "temperature_label_K": int(path.stem.split("_")[1].removesuffix("mK")) / 1000.0, "eta_uA_per_V": eta, "shunt_branches": per_shunt})

    target = np.loadtxt(TARGET_IV)
    eta = float(1.0 / np.polyfit(target[0, :10], target[1, :10], 1)[0])
    index = int(np.argmin(np.abs(target[0] - 1400.0)))
    ites = float(eta * target[1, index])
    target_rows = []
    for r_sh in (0.0038, 0.0039):
        ish = float(target[0, index] - ites)
        resistance = r_sh * ish / ites
        vtes = ish * 1e-6 * r_sh
        power = vtes * ites * 1e-6
        target_rows.append({"R_SH_ohm": r_sh, "I_bias_uA": float(target[0, index]), "I_TES_candidate_uA": ites, "V_TES_V": vtes, "R_candidate_ohm": resistance, "P_J_candidate_W": power, "source_class": "same_campaign_unlinked_proxy"})
    return {
        "source_class": "same_campaign_unlinked_proxy",
        "files": [path.as_posix() for path in IV_FILES],
        "temperature_series": rows,
        "R_SH_branches_are_reference_only": True,
        "target_same_day_candidate": {"source_file": TARGET_IV.as_posix(), "eta_uA_per_V": eta, "rows_at_1400uA": target_rows, "admissible_for_strict_target": False},
        "thermal_law_status": "not_identifiable_from_unlinked_and_nonmonotonic_series",
        "thermal_law_note": "DC points are retained for provenance and quality audit; no noise likelihood or residual is used.",
    }


def _pulse_constraints() -> dict:
    root = TARGET_PATH / "CH0_pulse/rawdata"
    files = sorted(root.glob("CH0_*.dat"), key=lambda p: int(p.stem.split("_")[1]))[:30]
    waveforms = []
    for path in files:
        values = np.fromfile(path, dtype=np.float64, offset=4)
        if values.size == 100000:
            waveforms.append(values - np.mean(values[:1000]))
    if not waveforms:
        return {"status": "unavailable", "source_class": "target_confirmed", "reason": "no pulse records"}
    waveform = np.median(np.asarray(waveforms), axis=0)
    peak = int(np.argmax(waveform))
    amplitude = float(waveform[peak])
    rise = {}
    for fraction in (0.2, 0.9):
        indices = np.where(waveform[:peak] >= fraction * amplitude)[0]
        rise[f"{fraction}"] = int(indices[0]) if len(indices) else None
    decay = {}
    for fraction in (0.9, 0.1):
        indices = np.where(waveform[peak:] <= fraction * amplitude)[0]
        decay[f"{fraction}"] = peak + int(indices[0]) if len(indices) else None
    return {
        "status": "combination_constraints_only",
        "source_class": "target_confirmed",
        "source_files": [path.as_posix() for path in files],
        "sample_rate_Hz": 500000.0,
        "records_summarized": len(waveforms),
        "median_waveform_peak_index": peak,
        "median_waveform_peak_time_s": peak / 500000.0,
        "median_waveform_peak_amplitude_raw": amplitude,
        "rise_crossing_indices": rise,
        "decay_crossing_indices": decay,
        "observed_decay_90_to_10_s": ((decay["0.1"] - decay["0.9"]) / 500000.0) if decay["0.1"] is not None and decay["0.9"] is not None else None,
        "interpretation": "Observed time scales constrain combinations such as tau_el and effective thermal poles; alpha, beta, C, G, and L are not individually fitted.",
        "identifiable_combinations": ["tau_el=L/[R_l+R(1+beta)]", "effective thermal time constant"],
        "individual_parameter_fit": False,
    }


def _generic_reference() -> dict:
    return json.loads(GENERIC_INPUT.read_text(encoding="utf-8"))


def _parameter_envelope(rt: dict, iv: dict, generic: dict) -> dict:
    tc_min, tc_max = rt["T_c_range_K"]
    tc_nominal = rt["T_c_nominal_K"]
    candidate_r = [row["R_candidate_ohm"] for row in iv["target_same_day_candidate"]["rows_at_1400uA"]]
    params = {
        "T_c": {"range": [tc_min, tc_max], "nominal": tc_nominal, "unit": "K", "source": [path.as_posix() for path in RT_FILES], "source_class": "nearby_run_proxy", "derivation": "RT transition envelope from onset/midpoint/high R/Rn definitions", "confidence": "conditional", "correlations": ["T_c definition with Rn and bias", "T_c-G-n thermal-law degeneracy"], "allowed_for_proxy_sweep": True, "allowed_for_strict_target": False},
        "T_bath": {"range": [0.215, 0.215], "nominal": 0.215, "unit": "K", "source": ["target directory label 215mK"], "source_class": "target_setpoint", "derivation": "setpoint label only", "confidence": "conditional", "correlations": ["operating-point power balance"], "allowed_for_proxy_sweep": True, "allowed_for_strict_target": False},
        "R": {"range": [min(candidate_r), max(candidate_r)], "nominal": float(np.mean(candidate_r)), "unit": "ohm", "source": [TARGET_IV.as_posix()], "source_class": "same_campaign_unlinked_proxy", "derivation": "R_SH-dependent deterministic IV conversion at 1400 uA candidate", "confidence": "conditional", "correlations": ["R_SH", "eta", "I_TES", "P_J"], "allowed_for_proxy_sweep": True, "allowed_for_strict_target": False},
        "R_SH": {"range": [0.0038, 0.0039], "nominal": 0.00385, "unit": "ohm", "source": ["Analyze_Experimental_Data/tes_analysis/iv.py", "getpara/thermalconductivity.py"], "source_class": "generic_code_reference", "derivation": "reference branches only; no board evidence", "confidence": "none_for_target", "correlations": ["R", "V_TES", "P_J"], "allowed_for_proxy_sweep": False, "allowed_for_strict_target": False},
    }
    generic_only = {"R_l", "alpha", "beta", "L", "n", "C_tes", "C_abs", "G_tes-bath", "G_abs-tes", "G_abs-abs"}
    for name in generic_only:
        value = generic[name]
        params[name] = {"range": None, "nominal": value, "unit": "SI_or_dimensionless", "source": [GENERIC_INPUT.as_posix()], "source_class": "simulation_reference_only", "derivation": "generic simulation nominal used only for sensitivity/reference scenarios", "confidence": "none_for_target", "correlations": ["multi-parameter electrothermal degeneracy"], "allowed_for_proxy_sweep": False, "allowed_for_strict_target": False, "allowed_for_sensitivity_reference": True}
    params["K"] = {"range": None, "nominal": None, "unit": "W/K^n", "source": iv["files"], "source_class": "unresolved", "derivation": "P=K(T_c^n-T_bath^n) cannot be target-linked from the unlinked/nonmonotonic series", "confidence": "none", "correlations": ["G_tes-bath=n*K*T_c^(n-1)", "T_c", "n"], "allowed_for_proxy_sweep": False, "allowed_for_strict_target": False}
    return {
        "stage": "A",
        "stage_name": "noise_blind_physics_construction",
        "freeze_status": "frozen",
        "noise_spectrum_read": False,
        "noise_likelihood_used": False,
        "parameters": params,
        "sensitivity_reference": {name: generic[name] for name in ("T_c", "T_bath", "R", "R_l", "alpha", "beta", "L", "n", "C_tes", "C_abs", "G_tes-bath", "G_abs-tes", "G_abs-abs")},
        "thermal_convention": {"formula": "P=K*(T_c^n-T_bath^n)=G_tes-bath*T_c/n*(1-(T_bath/T_c)^n)", "G_definition": "G_tes-bath=n*K*T_c^(n-1) at T_c", "repository_sources": ["PoST_Simulations/PoST_Simulation.py", "getpara/thermalconductivity.py"]},
        "strict_target_conclusion": "C — exact target physical case remains unidentified",
    }


def _scenarios(envelope: dict, generic: dict, count: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    tc_low, tc_high = envelope["parameters"]["T_c"]["range"]
    r_low, r_high = envelope["parameters"]["R"]["range"]
    records = []
    vary = ["R_l", "alpha", "beta", "L", "n", "C_tes", "C_abs", "G_tes-bath", "G_abs-tes", "G_abs-abs"]
    for index in range(count):
        params = dict(generic)
        params.update({"T_c": float(rng.uniform(tc_low, tc_high)), "T_bath": 0.215, "R": float(rng.choice([r_low, r_high])), "R_SH": float(rng.choice([0.0038, 0.0039])), "rate": 500000.0, "samples": 100000})
        factors = {}
        for name in vary:
            factor = float(np.exp(rng.uniform(np.log(0.5), np.log(2.0))))
            factors[name] = factor
            params[name] = float(generic[name]) * factor
        records.append({"scenario_id": f"proxy_{index:04d}", "scenario_class": "conditional_proxy_plus_simulation_reference_only", "seed": seed, "parameters": params, "factor_from_generic_reference": factors, "source_class_by_parameter": {"T_c": "nearby_run_proxy", "T_bath": "target_setpoint", "R": "same_campaign_unlinked_proxy", "R_SH": "generic_code_reference", **{name: "simulation_reference_only" for name in vary}}, "strict_target_allowed": False, "parameter_range_frozen_before_comparison": True})
    return {"stage": "A", "freeze_status": "frozen", "generation_method": "fixed_seed_log_uniform_reference_factors", "seed": seed, "sample_count": count, "noise_spectrum_read": False, "scenarios": records}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--count", type=int, default=96)
    parser.add_argument("--seed", type=int, default=20260905)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rt = _rt_summary()
    iv = _iv_summary()
    generic = _generic_reference()
    envelope = _parameter_envelope(rt, iv, generic)
    scenarios = _scenarios(envelope, generic, args.count, args.seed)
    constraints = _pulse_constraints()
    _dump(args.output_dir / "proxy_parameter_envelope.json", envelope)
    _dump(args.output_dir / "proxy_scenarios.json", scenarios)
    _dump(args.output_dir / "pulse_combination_constraints.json", constraints)
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[2], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    manifest = {
        "stage": "A",
        "stage_name": "noise_blind_physics_construction",
        "freeze_status": "frozen",
        "generation_method": scenarios["generation_method"],
        "seed": args.seed,
        "sample_count": args.count,
        "source_commit": commit,
        "parameter_ranges_frozen_before_comparison": True,
        "experimental_spectrum_read": False,
        "experimental_noise_path_read": False,
        "inputs": {"rt": [path.as_posix() for path in RT_FILES], "multi_temperature_iv": [path.as_posix() for path in IV_FILES], "same_day_iv_candidate": TARGET_IV.as_posix(), "pulse": (TARGET_PATH / "CH0_pulse/rawdata").as_posix(), "generic_reference": GENERIC_INPUT.as_posix()},
        "strict_target_input": {"path": (args.output_dir / "input.json").as_posix(), "exists_at_generation": False, "modified": False, "note": "Current main checkout has no case input.json; no strict target input was created or changed."},
        "outputs": ["proxy_parameter_envelope.json", "proxy_scenarios.json", "pulse_combination_constraints.json"],
        "prohibited_in_stage_A": ["experimental spectrum", "noise likelihood", "residual fit", "parameter optimizer using noise"],
        "strict_target_conclusion": "C — exact target physical case remains unidentified",
    }
    _dump(args.output_dir / "noise_blind_sweep_manifest.json", manifest)


if __name__ == "__main__":
    main()

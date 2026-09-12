"""Standalone PoST optimizer: differential evolution -> Powell.

The optimizer is specialized for the 2024-12-06 215 mK / 1400 uA target
case.  During optimization it evaluates the deterministic five-state TES noise
model through the confirmed 100 kHz fourth-order analog Bessel, sampling/alias
fold, and 10 kHz digital analysis Bessel.  This avoids fitting 4096-sample
finite-record artifacts.  After each R_SH branch is optimized, the best point
is re-run with the full experimental record length so the generated comparison
plot uses the same finite-record processing as the measurement.
The TES resistance is no longer fitted freely.  It is calculated from the
same-campaign 1400 uA IV point for each R_SH value in the configured sweep.
The .dat file must include the 100 kHz hardware Bessel and a 10 kHz analysis
Bessel, so each candidate fixes input.json["cutoff"] to 10000.

By default input.json is restored after optimization.  Add --apply-final to
keep the best parameters (including the best R_SH branch) and run PoST one
final time.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import differential_evolution, minimize

from lib.tes_noise_model import operating_point as tes_operating_point
from subScript.noise_measurement_model import (
    HARDWARE_BESSEL_CUTOFF_HZ,
    analysis_filter_magnitude,
    hardware_sampled_asd,
)


# ---------- Paths ----------
# Keep the defaults relocatable.  The old H: drive paths were specific to the
# original workstation and make the script fail before optimization starts.
SCRIPT_DIR = Path(__file__).resolve().parent
POST_SCRIPT = SCRIPT_DIR / "PoST_Simulation.py"
INPUT_PATH = SCRIPT_DIR / "input.json"
TARGET_CASE_DIR = (
    SCRIPT_DIR
    / "cases"
    / "tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"
)
TARGET_ENVELOPE_PATH = TARGET_CASE_DIR / "proxy_parameter_envelope.json"
TARGET_SCENARIOS_PATH = TARGET_CASE_DIR / "proxy_scenarios.json"
NOISE_DAT_PATH = SCRIPT_DIR / "noise_total-bessel100k.dat"
EXPERIMENT_ROOT = Path(
    r"G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_"
    r"difftrig5e-5_rate500k_samples100k_gain5_day2"
)
MODELNOISE_PATH = EXPERIMENT_ROOT / "CH0_noise" / "modelnoise.txt"
PULSE_CONFIG_PATH = EXPERIMENT_ROOT / "PulseConfig.json"
TARGET_IV_PATH = Path(r"G:\tagawa\20241206\room1-ch1-iv3\calibration\IV_215mK.txt")
TARGET_BIAS_UA = 1400.0

# ---------- Comparison settings ----------
# Fit the full requested band.  Log-spaced samples give comparable weight per
# frequency decade; the default extra high-frequency multiplier is disabled.
FIT_MIN_HZ = 1_000.0
FIT_MAX_HZ = 200_000.0
FIT_WEIGHT_START_HZ = 30_000.0
FIT_HIGH_FREQUENCY_WEIGHT = 1.0
FIT_ROBUST_DELTA_DEX = 0.12
FIT_POINTS = 401
SIM_ANALYSIS_CUTOFF_HZ = 10_000.0  # experimental NoiseAnalysis cutoff.
TARGET_HARDWARE_BESSEL_ORDER = 4
EXCESS_JOHNSON_MAX = 5.0


@dataclass(frozen=True)
class Bound:
    lower: float
    upper: float
    logarithmic: bool = True


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-maxiter", type=int, default=8)
    parser.add_argument("--de-popsize", type=int, default=3)
    parser.add_argument("--powell-maxfev", type=int, default=250)
    parser.add_argument("--skip-de", action="store_true")
    parser.add_argument("--apply-final", action="store_true")
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument(
        "--fit-min-hz",
        type=float,
        default=FIT_MIN_HZ,
        help="Lowest frequency included in the fit. Frequencies below this are ignored.",
    )
    parser.add_argument(
        "--fit-max-hz",
        type=float,
        default=FIT_MAX_HZ,
        help="Highest frequency included in the fit.",
    )
    parser.add_argument(
        "--fit-weight-start-hz",
        type=float,
        default=FIT_WEIGHT_START_HZ,
        help=(
            "Frequency above which the log-frequency fit weight ramps from 1 "
            "toward --high-frequency-weight."
        ),
    )
    parser.add_argument(
        "--high-frequency-weight",
        type=float,
        default=FIT_HIGH_FREQUENCY_WEIGHT,
        help="Fit-weight multiplier reached at --fit-max-hz.",
    )
    parser.add_argument(
        "--robust-delta-dex",
        type=float,
        default=FIT_ROBUST_DELTA_DEX,
        help=(
            "Huber transition for log10(model/measurement) residuals. "
            "Set 0 for ordinary weighted least squares."
        ),
    )
    parser.add_argument(
        "--fit-points",
        type=int,
        default=FIT_POINTS,
        help="Number of log-spaced frequencies used by the objective.",
    )
    parser.add_argument(
        "--case-dir",
        type=Path,
        default=TARGET_CASE_DIR,
        help=(
            "Target-case directory containing proxy_parameter_envelope.json "
            "and proxy_scenarios.json."
        ),
    )
    parser.add_argument(
        "--reference-input",
        type=Path,
        default=None,
        help=(
            "Optional explicit starting input. By default the first frozen, "
            "pulse-consistent target-case proxy scenario is used."
        ),
    )
    parser.add_argument("--target-iv", type=Path, default=TARGET_IV_PATH,
                        help="Same-element IV file used to derive TES R at the target bias.")
    parser.add_argument("--target-bias-ua", type=float, default=TARGET_BIAS_UA,
                        help="Bias current in uA at which TES R is derived from the IV.")
    parser.add_argument("--r-shunt-start-mohm", type=float, default=3.8,
                        help="First shunt-resistance branch in mOhm.")
    parser.add_argument("--r-shunt-stop-mohm", type=float, default=3.9,
                        help="Last shunt-resistance branch in mOhm.")
    parser.add_argument("--r-shunt-count", type=int, default=3,
                        help="Number of equally spaced R_SH branches, including endpoints.")
    return parser.parse_args()


def load_json(path: Path):
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def target_case_reference(case_dir: Path, explicit_reference: Path | None = None):
    """Return a target-consistent starting point and the frozen envelope."""

    envelope_path = case_dir / "proxy_parameter_envelope.json"
    scenarios_path = case_dir / "proxy_scenarios.json"
    envelope = load_json(envelope_path)
    scenarios = load_json(scenarios_path)

    frozen = scenarios.get("pulse_consistent_scenarios", [])
    if not frozen:
        raise RuntimeError(f"No pulse-consistent proxy scenarios in {scenarios_path}")

    if explicit_reference is None:
        reference = dict(frozen[0]["parameters"])
        reference_source = f"{scenarios_path}:pulse_consistent_scenarios[0]"
    else:
        reference = load_json(explicit_reference)
        reference_source = str(explicit_reference)

    # The target acquisition is the 215 mK case.  Keep the measured/setpoint
    # bath condition fixed instead of inheriting the unrelated 136 mK generic
    # PoST input.
    reference["T_bath"] = float(envelope["parameters"]["T_bath"]["nominal"])
    reference["rate"] = 500_000.0
    reference["samples"] = 100_000
    reference["hardware_bessel_order"] = TARGET_HARDWARE_BESSEL_ORDER
    reference["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    reference.setdefault("excess_johnson_M", 0.0)
    return reference, envelope, reference_source


def tes_resistance_from_iv(iv_path: Path, shunt_resistance_ohm: float,
                           target_bias_uA: float) -> dict:
    """Derive the TES operating resistance from a calibrated IV point.

    The IV files store two rows: bias current in uA and measured output
    voltage.  The first ten points define the superconducting calibration
    slope, matching the existing IV analysis workflow.
    """

    data = np.asarray(np.loadtxt(iv_path), dtype=float)
    if data.ndim != 2 or data.shape[0] < 2:
        raise ValueError(f"Expected two-row IV data in {iv_path}")

    bias_uA = data[0]
    output_voltage = data[1]
    finite = np.isfinite(bias_uA) & np.isfinite(output_voltage)
    if np.count_nonzero(finite) < 10:
        raise ValueError(f"Not enough finite IV points in {iv_path}")

    bias_uA = bias_uA[finite]
    output_voltage = output_voltage[finite]
    order = np.argsort(bias_uA)
    bias_uA = bias_uA[order]
    output_voltage = output_voltage[order]

    calibration_count = min(10, len(bias_uA))
    slope, intercept = np.polyfit(
        bias_uA[:calibration_count],
        output_voltage[:calibration_count],
        1,
    )
    if not np.isfinite(slope) or slope == 0.0:
        raise ValueError(f"Invalid superconducting IV slope in {iv_path}")

    eta_uA_per_V = 1.0 / slope
    index = int(np.argmin(np.abs(bias_uA - float(target_bias_uA))))
    i_bias_uA = float(bias_uA[index])
    i_tes_uA = float(eta_uA_per_V * output_voltage[index])
    i_shunt_uA = i_bias_uA - i_tes_uA
    if i_tes_uA <= 0.0 or i_shunt_uA <= 0.0:
        raise ValueError(
            f"Invalid operating point in {iv_path}: "
            f"I_bias={i_bias_uA:g} uA, I_TES={i_tes_uA:g} uA"
        )

    v_tes = i_shunt_uA * 1.0e-6 * float(shunt_resistance_ohm)
    r_tes = v_tes / (i_tes_uA * 1.0e-6)
    return {
        "R_SH_ohm": float(shunt_resistance_ohm),
        "I_bias_uA": i_bias_uA,
        "I_TES_uA": i_tes_uA,
        "V_TES_V": float(v_tes),
        "R_TES_ohm": float(r_tes),
        "P_J_W": float(v_tes * i_tes_uA * 1.0e-6),
        "eta_uA_per_V": float(eta_uA_per_V),
        "iv_point_index": index,
    }


def write_json_atomically(path: Path, data: dict):
    temporary = path.with_suffix(path.suffix + ".optimize-tmp")
    with temporary.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)
        f.write("\n")
    temporary.replace(path)


def normalize_at(freq, asd, reference_hz=1000.0):
    value = asd[np.argmin(np.abs(freq - reference_hz))]
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"Invalid ASD at {reference_hz} Hz")
    return asd / value


def fit_frequency_weights(frequency: np.ndarray, args) -> np.ndarray:
    """Return a smooth high-frequency emphasis on a log-frequency axis."""

    frequency = np.asarray(frequency, dtype=float)
    weights = np.ones_like(frequency)
    high_weight = float(args.high_frequency_weight)
    weight_start = max(float(args.fit_weight_start_hz), float(args.fit_min_hz))
    fit_max = float(args.fit_max_hz)

    if high_weight <= 1.0 or fit_max <= weight_start:
        return weights

    denominator = np.log10(fit_max) - np.log10(weight_start)
    if denominator <= 0.0:
        return weights

    ramp = (
        np.log10(np.maximum(frequency, weight_start)) - np.log10(weight_start)
    ) / denominator
    ramp = np.clip(ramp, 0.0, 1.0)
    return 1.0 + (high_weight - 1.0) * ramp


def fit_score(
    model: np.ndarray,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> float:
    """Weighted robust loss in log-ASD ratio across the configured fit band."""

    model = np.asarray(model, dtype=float)
    target = np.asarray(target, dtype=float)
    if (
        np.any(~np.isfinite(model))
        or np.any(~np.isfinite(target))
        or np.any(model <= 0.0)
        or np.any(target <= 0.0)
    ):
        raise ValueError("Fit spectra must be finite and strictly positive")

    residual = np.log10(model) - np.log10(target)
    delta = float(args.robust_delta_dex)
    if delta > 0.0:
        absolute = np.abs(residual)
        # Huber-like loss with the same quadratic scale as residual**2 near zero.
        loss = np.where(
            absolute <= delta,
            residual**2,
            2.0 * delta * absolute - delta**2,
        )
    else:
        loss = residual**2

    weights = fit_frequency_weights(fit_freq, args)
    return float(np.sum(weights * loss) / np.sum(weights))


def target_spectrum(args):
    pulse_config = load_json(PULSE_CONFIG_PATH)
    exp_asd = np.loadtxt(MODELNOISE_PATH)
    exp_rate = float(pulse_config["Readout"]["Rate"])
    exp_freq = np.arange(len(exp_asd)) * (exp_rate / 2.0) / len(exp_asd)
    fit_freq = np.geomspace(
        float(args.fit_min_hz),
        float(args.fit_max_hz),
        int(args.fit_points),
    )
    target = np.interp(fit_freq, exp_freq, normalize_at(exp_freq, exp_asd))
    return fit_freq, target


def run_post(timeout_seconds: int, output_dir: Path, noise_path: Path):
    before = noise_path.stat().st_mtime_ns if noise_path.exists() else -1
    result = subprocess.run(
        [
            sys.executable,
            str(POST_SCRIPT),
            "--noise-only",
            "--output",
            str(output_dir),
        ],
        cwd=POST_SCRIPT.parent,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    if result.returncode != 0:
        raise RuntimeError("PoST failed:\n" + result.stderr[-4000:])
    if not noise_path.exists() or noise_path.stat().st_mtime_ns <= before:
        raise RuntimeError(
            "noise_total-bessel100k.dat was not updated. "
            "PoST_Simulation.py must end with MakeNoise() followed by _ShowNoiseSpectrum()."
        )


def simulated_spectrum(candidate: dict, fit_freq: np.ndarray, noise_path: Path):
    sim_asd = np.loadtxt(noise_path)
    sim_rate = float(candidate["rate"])
    sim_freq = np.arange(len(sim_asd)) * (sim_rate / 2.0) / len(sim_asd)
    return np.interp(fit_freq, sim_freq, normalize_at(sim_freq, sim_asd))


def plot_sweep_comparison(summary: dict, output_path: Path) -> Path:
    """Plot normalized measured ASD against every fixed-R_SH result."""

    pulse_config = load_json(PULSE_CONFIG_PATH)
    experimental_rate = float(pulse_config["Readout"]["Rate"])
    experimental_asd = np.asarray(np.loadtxt(MODELNOISE_PATH), dtype=float)
    experimental_frequency = (
        np.arange(len(experimental_asd))
        * (experimental_rate / 2.0)
        / len(experimental_asd)
    )
    experimental_normalized = normalize_at(
        experimental_frequency,
        experimental_asd,
    )
    experimental_mask = (
        (experimental_frequency > 0.0)
        & np.isfinite(experimental_normalized)
        & (experimental_normalized > 0.0)
    )
    plot_min_frequency = float(np.min(experimental_frequency[experimental_mask]))
    plot_max_frequency = float(np.max(experimental_frequency[experimental_mask]))

    fig, (spectrum_axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    spectrum_axis.loglog(
        experimental_frequency[experimental_mask],
        experimental_normalized[experimental_mask],
        color="black",
        linewidth=2.0,
        label="Measured CH0",
    )

    colors = ("tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple")
    cases = summary["cases"]
    fit_min_hz = float(summary["fit"]["min_hz"])
    fit_max_hz = float(summary["fit"]["max_hz"])
    best_r_sh = float(summary["best_case_R_SH_ohm"])
    for index, case in enumerate(cases):
        noise_path = Path(case["work_dir"]) / NOISE_DAT_PATH.name
        model_asd = np.asarray(np.loadtxt(noise_path), dtype=float)
        model_frequency = (
            np.arange(len(model_asd))
            * (experimental_rate / 2.0)
            / len(model_asd)
        )
        model_normalized = normalize_at(model_frequency, model_asd)
        model_mask = (
            (model_frequency > 0.0)
            & np.isfinite(model_normalized)
            & (model_normalized > 0.0)
        )
        label = f"R_SH={float(case['R_SH_ohm']) * 1e3:.2f} mOhm"
        is_best = np.isclose(float(case["R_SH_ohm"]), best_r_sh)
        spectrum_axis.loglog(
            model_frequency[model_mask],
            model_normalized[model_mask],
            color=colors[index % len(colors)],
            linewidth=2.4 if is_best else 1.3,
            linestyle="-" if is_best else "--",
            label=label + (" (best)" if is_best else ""),
        )

        measured_on_model_grid = np.interp(
            model_frequency[model_mask],
            experimental_frequency[experimental_mask],
            experimental_normalized[experimental_mask],
        )
        ratio = model_normalized[model_mask] / measured_on_model_grid
        ratio_axis.semilogx(
            model_frequency[model_mask],
            ratio,
            color=colors[index % len(colors)],
            linewidth=2.0 if is_best else 1.1,
            linestyle="-" if is_best else "--",
            label=label,
        )

    spectrum_axis.set_ylabel("Normalized ASD")
    spectrum_axis.set_title("PoST noise model vs measured CH0 noise (full frequency range)")
    spectrum_axis.grid(True, which="both", alpha=0.25)
    spectrum_axis.axvspan(
        fit_min_hz,
        fit_max_hz,
        color="gray",
        alpha=0.08,
        label="fit band",
    )
    spectrum_axis.legend(fontsize=9)
    ratio_axis.axhline(1.0, color="black", linewidth=1.0)
    ratio_axis.fill_between(
        [fit_min_hz, fit_max_hz],
        [0.9, 0.9],
        [1.1, 1.1],
        color="gray",
        alpha=0.15,
        label="±10%",
    )
    ratio_axis.set_xlabel("Frequency [Hz]")
    ratio_axis.set_ylabel("Model / measured")
    ratio_axis.set_xlim(plot_min_frequency, plot_max_frequency)
    ratio_axis.set_yscale("log")
    ratio_axis.grid(True, which="both", alpha=0.25)
    ratio_axis.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def parameter_bounds(reference: dict, envelope: dict, fixed_r_ohm: float):
    """Target-case search box derived from the frozen proxy envelope.

    T_bath and R are not fitted here: T_bath is fixed to the 215 mK target
    setpoint and R is fixed by the same-campaign IV conversion for each R_SH
    branch.  The remaining parameters cover the stable target-case region
    already exercised by the repository's adaptive/correlated searches.
    """

    sensitivity = envelope["sensitivity_reference"]
    tc_low, tc_high = map(float, envelope["parameters"]["T_c"]["range"])

    def ref(name):
        return float(sensitivity[name])

    r_eff = (
        float(fixed_r_ohm) * (1.0 + ref("beta"))
        + ref("R_l")
    )
    l_low = min(ref("L") / 10.0, r_eff / (2.0 * np.pi * 300_000.0) / 3.0)
    l_high = max(ref("L") * 10.0, r_eff / (2.0 * np.pi * 1_000.0) * 3.0)

    return {
        "T_c": Bound(tc_low, tc_high, logarithmic=False),
        "R_l": Bound(ref("R_l") * 0.25, ref("R_l") * 6.0),
        "alpha": Bound(ref("alpha") * 0.05, ref("alpha") * 15.0),
        "beta": Bound(0.0, 12.0, logarithmic=False),
        "L": Bound(max(l_low, 1e-12), l_high),
        "n": Bound(max(1.01, ref("n") * 0.5), ref("n") * 2.0),
        "C_tes": Bound(ref("C_tes") * 0.1, ref("C_tes") * 10.0),
        "C_abs": Bound(ref("C_abs") * 0.1, ref("C_abs") * 10.0),
        "G_tes-bath": Bound(ref("G_tes-bath") * 0.1, ref("G_tes-bath") * 10.0),
        "G_abs-tes": Bound(ref("G_abs-tes") * 0.1, ref("G_abs-tes") * 20.0),
        "G_abs-abs": Bound(ref("G_abs-abs") * 0.1, ref("G_abs-abs") * 10.0),
        "excess_johnson_M": Bound(0.0, EXCESS_JOHNSON_MAX, logarithmic=False),
    }


def vector_bounds(bounds: dict):
    keys = list(bounds)
    transformed = []
    for key in keys:
        bound = bounds[key]
        if bound.logarithmic:
            transformed.append((np.log10(bound.lower), np.log10(bound.upper)))
        else:
            transformed.append((bound.lower, bound.upper))
    return keys, transformed


def encode(candidate: dict, keys, bounds):
    return np.asarray([
        np.log10(candidate[key]) if bounds[key].logarithmic else candidate[key]
        for key in keys
    ], dtype=float)


def decode(
    vector,
    original: dict,
    keys,
    bounds,
    sample_rate: float,
    post_filter_white_asd: float,
):
    candidate = original.copy()
    for key, value in zip(keys, vector):
        candidate[key] = float(10.0 ** value if bounds[key].logarithmic else value)
    # Match the experimental modelnoise.txt analysis filter; this is not fitted.
    candidate["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    # Filter coefficients and the frequency grid must use the acquisition rate
    # of modelnoise.txt, even when an old optimization reference used 300 kHz.
    candidate["rate"] = sample_rate
    candidate["post_filter_white_asd_A_rtHz"] = post_filter_white_asd
    candidate.pop("readout_white_asd_A_rtHz", None)
    candidate["hardware_bessel_order"] = TARGET_HARDWARE_BESSEL_ORDER
    return candidate


def deterministic_simulated_spectrum(
    candidate: dict,
    fit_freq: np.ndarray,
) -> np.ndarray:
    """Expected normalized post-analysis ASD without finite-record Monte Carlo."""

    frequency = np.unique(
        np.concatenate((np.asarray([1_000.0]), np.asarray(fit_freq, dtype=float)))
    )
    rate = float(candidate["rate"])
    pre_analysis = hardware_sampled_asd(
        candidate,
        frequency,
        rate_hz=rate,
        cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
        order=TARGET_HARDWARE_BESSEL_ORDER,
    )
    white = float(candidate.get("post_filter_white_asd_A_rtHz", 0.0))
    if white < 0.0:
        raise ValueError("post_filter_white_asd_A_rtHz must be non-negative")
    pre_analysis = np.sqrt(pre_analysis**2 + white**2)
    expected = pre_analysis * analysis_filter_magnitude(
        frequency,
        rate,
        cutoff_hz=SIM_ANALYSIS_CUTOFF_HZ,
    )
    normalized = normalize_at(frequency, expected)
    return np.interp(fit_freq, frequency, normalized)


def optimize_case(
    args,
    original: dict,
    reference: dict,
    envelope: dict,
    fit_freq: np.ndarray,
    target: np.ndarray,
    experimental_rate: float,
    experimental_samples: int,
    post_filter_white_asd: float,
    operating_point: dict,
    work_dir: Path,
):
    """Optimize one fixed-R case derived from one R_SH branch."""

    fixed_r_ohm = float(operating_point["R_TES_ohm"])
    shunt_resistance_ohm = float(operating_point["R_SH_ohm"])

    # Start from the frozen 215 mK target-case proxy, not the generic
    # PoST_Simulations/input.json (which belongs to a different thermal point).
    case_original = reference.copy()
    case_original["R"] = fixed_r_ohm
    case_original["R_SH"] = shunt_resistance_ohm
    case_original["T_bath"] = float(envelope["parameters"]["T_bath"]["nominal"])
    case_original["rate"] = float(experimental_rate)
    case_original["samples"] = int(experimental_samples)
    case_original["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    case_original["hardware_bessel_order"] = TARGET_HARDWARE_BESSEL_ORDER

    work_dir.mkdir(parents=True, exist_ok=True)
    work_input_path = work_dir / "input.json"
    work_noise_path = work_dir / NOISE_DAT_PATH.name
    bounds = parameter_bounds(reference, envelope, fixed_r_ohm=fixed_r_ohm)
    keys, scipy_bounds = vector_bounds(bounds)

    initial = case_original.copy()
    initial["post_filter_white_asd_A_rtHz"] = post_filter_white_asd
    initial.pop("readout_white_asd_A_rtHz", None)
    initial_x = encode(initial, keys, bounds)

    cache = {}
    evaluation_count = 0
    stability_rejection_count = 0
    simulation_failure_count = 0
    best_score = np.inf
    best_candidate = initial.copy()

    def objective(vector):
        nonlocal evaluation_count, stability_rejection_count
        nonlocal simulation_failure_count, best_score, best_candidate

        cache_key = tuple(np.round(vector, 12))
        if cache_key in cache:
            return cache[cache_key]

        candidate = decode(
            vector,
            case_original,
            keys,
            bounds,
            experimental_rate,
            post_filter_white_asd,
        )
        candidate["R"] = fixed_r_ohm
        candidate["R_SH"] = shunt_resistance_ohm
        candidate["T_bath"] = float(envelope["parameters"]["T_bath"]["nominal"])
        candidate["samples"] = int(experimental_samples)
        evaluation_count += 1

        point = tes_operating_point(candidate)
        if not point.get("stable", False):
            stability_rejection_count += 1
            max_real = point.get("max_eigenvalue_real_s_inv")
            detail = (
                f", max Re(lambda)={float(max_real):.6g} 1/s"
                if max_real is not None and np.isfinite(max_real)
                else ""
            )
            print(
                f"R_SH={shunt_resistance_ohm * 1e3:.4f} mOhm, "
                f"evaluation {evaluation_count:4d}: rejected "
                f"({point.get('reason', 'unstable')}{detail})"
            )
            score = 1e12
            cache[cache_key] = score
            return score

        try:
            model = deterministic_simulated_spectrum(candidate, fit_freq)
            score = fit_score(model, target, fit_freq, args)
            print(
                f"R_SH={shunt_resistance_ohm * 1e3:.4f} mOhm, "
                f"evaluation {evaluation_count:4d}: {score:.6g}"
            )
        except Exception as error:
            simulation_failure_count += 1
            print(
                f"R_SH={shunt_resistance_ohm * 1e3:.4f} mOhm, "
                f"evaluation {evaluation_count:4d}: rejected ({error})"
            )
            score = 1e12

        cache[cache_key] = score
        if score < best_score and score < 1e11:
            best_score = score
            best_candidate = candidate.copy()
            print("  new best score:", best_score)
        return score

    print(
        f"\n=== R_SH={shunt_resistance_ohm * 1e3:.4f} mOhm: "
        f"R_TES={fixed_r_ohm * 1e3:.6f} mOhm ==="
    )
    print(
        f"Target assumptions: T_bath={case_original['T_bath']:.6g} K, "
        f"T_c search={bounds['T_c'].lower:.6g}--{bounds['T_c'].upper:.6g} K"
    )
    print("Initial deterministic evaluation")
    objective(initial_x)

    if args.skip_de:
        starting_x = initial_x
    else:
        print("\nStage 1/2: differential_evolution")
        global_result = differential_evolution(
            objective,
            scipy_bounds,
            maxiter=args.de_maxiter,
            popsize=args.de_popsize,
            seed=args.seed,
            workers=1,
            updating="immediate",
            polish=False,
        )
        starting_x = global_result.x
        print("DE result:", global_result.fun)

    print("\nStage 2/2: Powell")
    local_result = minimize(
        objective,
        starting_x,
        method="Powell",
        bounds=scipy_bounds,
        options={"maxfev": args.powell_maxfev, "disp": True},
    )
    objective(local_result.x)

    if not np.isfinite(best_score) or best_score >= 1e11:
        raise RuntimeError(
            f"No stable valid candidate found for R_SH={shunt_resistance_ohm * 1e3:.4f} mOhm"
        )

    # Validate only the winner through the full production finite-record path.
    # This keeps the optimizer fast while ensuring the comparison plot is made
    # with the same 100k-sample record length as the experiment.
    validation_candidate = best_candidate.copy()
    validation_candidate["samples"] = int(experimental_samples)
    validation_candidate["rate"] = float(experimental_rate)
    write_json_atomically(work_input_path, validation_candidate)
    run_post(args.timeout, work_dir, work_noise_path)
    finite_model = simulated_spectrum(
        validation_candidate,
        fit_freq,
        work_noise_path,
    )
    finite_validation_score = fit_score(
        finite_model,
        target,
        fit_freq,
        args,
    )

    print("\nObjective evaluations:", evaluation_count)
    print("Stability rejections:", stability_rejection_count)
    print("Model-evaluation failures:", simulation_failure_count)
    print("Best deterministic score:", best_score)
    print("Full-record validation score:", finite_validation_score)
    print("Best fitted parameters:")
    print(json.dumps({key: best_candidate[key] for key in keys}, indent=2))

    return {
        "R_SH_ohm": shunt_resistance_ohm,
        "R_TES_ohm": fixed_r_ohm,
        "iv_operating_point": operating_point,
        "best_score": float(best_score),
        "finite_validation_score": float(finite_validation_score),
        "evaluations": evaluation_count,
        "stability_rejections": stability_rejection_count,
        "simulation_failures": simulation_failure_count,
        "searched_parameters": list(keys),
        "search_bounds": {
            key: {
                "lower": float(bounds[key].lower),
                "upper": float(bounds[key].upper),
                "logarithmic": bool(bounds[key].logarithmic),
            }
            for key in keys
        },
        "best_candidate": best_candidate,
        "work_dir": str(work_dir),
    }

def main():
    args = arguments()

    required_paths = [
        POST_SCRIPT,
        INPUT_PATH,
        MODELNOISE_PATH,
        PULSE_CONFIG_PATH,
        args.target_iv,
        args.case_dir / "proxy_parameter_envelope.json",
        args.case_dir / "proxy_scenarios.json",
    ]
    if args.reference_input is not None:
        required_paths.append(args.reference_input)
    for required_path in required_paths:
        if not required_path.is_file():
            raise FileNotFoundError(f"Missing required file: {required_path}")

    if args.r_shunt_count < 1:
        raise ValueError("--r-shunt-count must be at least 1")
    if args.r_shunt_start_mohm <= 0.0 or args.r_shunt_stop_mohm <= 0.0:
        raise ValueError("R_SH values must be positive")
    if args.fit_min_hz <= 0.0 or args.fit_max_hz <= args.fit_min_hz:
        raise ValueError("Fit band must satisfy 0 < --fit-min-hz < --fit-max-hz")
    if args.fit_weight_start_hz <= 0.0:
        raise ValueError("--fit-weight-start-hz must be positive")
    if args.high_frequency_weight < 1.0:
        raise ValueError("--high-frequency-weight must be at least 1")
    if args.robust_delta_dex < 0.0:
        raise ValueError("--robust-delta-dex must be non-negative")
    if args.fit_points < 16:
        raise ValueError("--fit-points must be at least 16")

    original = load_json(INPUT_PATH)
    reference, envelope, reference_source = target_case_reference(
        args.case_dir,
        args.reference_input,
    )
    pulse_config = load_json(PULSE_CONFIG_PATH)
    experimental_rate = float(pulse_config["Readout"]["Rate"])
    experimental_samples = int(pulse_config["Readout"]["Sample"])
    if args.fit_max_hz >= experimental_rate / 2.0:
        raise ValueError(
            f"--fit-max-hz={args.fit_max_hz:g} must be below Nyquist "
            f"({experimental_rate / 2.0:g} Hz)"
        )

    reference["rate"] = experimental_rate
    reference["samples"] = experimental_samples
    reference["T_bath"] = float(envelope["parameters"]["T_bath"]["nominal"])
    reference["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    reference["hardware_bessel_order"] = TARGET_HARDWARE_BESSEL_ORDER

    post_filter_white_asd = float(
        reference.get(
            "post_filter_white_asd_A_rtHz",
            reference.get("readout_white_asd_A_rtHz", 0.0),
        )
    )

    backup = INPUT_PATH.with_name(
        f"input.before_scipy_optimization.{time.strftime('%Y%m%d-%H%M%S')}.json"
    )
    shutil.copy2(INPUT_PATH, backup)
    print("backup:", backup)
    print("target-case reference:", reference_source)
    print(
        "fixed acquisition:",
        {
            "T_bath_K": reference["T_bath"],
            "rate_Hz": experimental_rate,
            "samples": experimental_samples,
            "hardware_bessel_order": TARGET_HARDWARE_BESSEL_ORDER,
            "hardware_bessel_cutoff_Hz": HARDWARE_BESSEL_CUTOFF_HZ,
            "analysis_bessel_cutoff_Hz": SIM_ANALYSIS_CUTOFF_HZ,
        },
    )
    fit_freq, target = target_spectrum(args)

    try:
        shunt_values_ohm = np.linspace(
            args.r_shunt_start_mohm * 1e-3,
            args.r_shunt_stop_mohm * 1e-3,
            args.r_shunt_count,
        )
        sweep_root = INPUT_PATH.parent / ".noise_optimization_work_rsh_sweep"
        sweep_root.mkdir(exist_ok=True)
        cases = []

        for shunt_resistance_ohm in shunt_values_ohm:
            operating_point = tes_resistance_from_iv(
                args.target_iv,
                float(shunt_resistance_ohm),
                args.target_bias_ua,
            )
            case_name = (
                f"rsh_{shunt_resistance_ohm * 1e3:.4f}mohm".replace(".", "p")
            )
            result = optimize_case(
                args,
                original,
                reference,
                envelope,
                fit_freq,
                target,
                experimental_rate,
                experimental_samples,
                post_filter_white_asd,
                operating_point,
                sweep_root / case_name,
            )
            cases.append(result)

        # All branches use the same finite-record seed in production, so this
        # is a fair branch-to-branch re-ranking after deterministic fitting.
        best_case = min(cases, key=lambda row: row["finite_validation_score"])
        summary = {
            "target_case_dir": str(args.case_dir),
            "target_case_reference": reference_source,
            "target_iv": str(args.target_iv),
            "target_bias_uA": args.target_bias_ua,
            "r_shunt_sweep_mOhm": [
                float(value * 1e3) for value in shunt_values_ohm
            ],
            "tes_resistance_is_fitted": False,
            "fixed_target_assumptions": {
                "T_bath_K": float(reference["T_bath"]),
                "rate_Hz": experimental_rate,
                "samples": experimental_samples,
                "hardware_bessel_order": TARGET_HARDWARE_BESSEL_ORDER,
                "hardware_bessel_cutoff_Hz": float(HARDWARE_BESSEL_CUTOFF_HZ),
                "analysis_bessel_cutoff_Hz": float(SIM_ANALYSIS_CUTOFF_HZ),
            },
            "T_c_proxy_range_K": [
                float(envelope["parameters"]["T_c"]["range"][0]),
                float(envelope["parameters"]["T_c"]["range"][1]),
            ],
            "optimization_evaluator": (
                "deterministic expected post-analysis ASD; full-record "
                "production validation only for each branch winner"
            ),
            "fit": {
                "min_hz": float(args.fit_min_hz),
                "max_hz": float(args.fit_max_hz),
                "weight_start_hz": float(args.fit_weight_start_hz),
                "high_frequency_weight": float(args.high_frequency_weight),
                "robust_delta_dex": float(args.robust_delta_dex),
                "points": int(args.fit_points),
                "loss": "weighted robust log10(model/measurement) residual",
                "frequencies_outside_fit_band_ignored": True,
            },
            "cases": [
                {
                    key: value
                    for key, value in case.items()
                    if key != "best_candidate"
                }
                for case in cases
            ],
            "best_case_R_SH_ohm": best_case["R_SH_ohm"],
            "best_case_R_TES_ohm": best_case["R_TES_ohm"],
            "best_case_score": best_case["best_score"],
            "best_case_finite_validation_score": best_case[
                "finite_validation_score"
            ],
        }

        summary_path = sweep_root / "summary.json"
        write_json_atomically(summary_path, summary)
        comparison_path = plot_sweep_comparison(
            summary,
            sweep_root / "noise_comparison.png",
        )
        summary["comparison_plot"] = str(comparison_path)
        write_json_atomically(summary_path, summary)
        print("\nR_SH sweep summary:")
        print(json.dumps(summary, indent=2))

        if args.apply_final:
            final_candidate = best_case["best_candidate"].copy()
            final_candidate["samples"] = experimental_samples
            final_candidate["rate"] = experimental_rate
            final_candidate["T_bath"] = float(reference["T_bath"])
            final_candidate["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
            final_candidate["hardware_bessel_order"] = (
                TARGET_HARDWARE_BESSEL_ORDER
            )
            write_json_atomically(INPUT_PATH, final_candidate)
            run_post(args.timeout, INPUT_PATH.parent, NOISE_DAT_PATH)
            print(
                "Best input.json applied and final .dat regenerated: "
                f"R_SH={best_case['R_SH_ohm'] * 1e3:.4f} mOhm, "
                f"R_TES={best_case['R_TES_ohm'] * 1e3:.6f} mOhm"
            )
        else:
            print(
                "Original input.json was not changed; use --apply-final to "
                "keep the best branch."
            )
    except Exception:
        write_json_atomically(INPUT_PATH, original)
        print("Original input.json restored after error.")
        raise


if __name__ == "__main__":
    main()

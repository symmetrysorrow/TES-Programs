"""Standalone PoST optimizer: differential evolution -> Powell.

For each trial, this script writes an isolated input.json, runs
PoST_Simulation.py, then compares the newly generated
noise_total-bessel100k.dat with tagawa CH0_noise/modelnoise.txt.  The default
fit scores the full 1 kHz--200 kHz band on a log-frequency grid.  The loss is
computed from log-ASD ratios with a robust Huber-like penalty so isolated
spectral spikes do not dominate the fit.  The generated comparison plot still
shows the full positive frequency range.
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


# ---------- Paths ----------
# Keep the defaults relocatable.  The old H: drive paths were specific to the
# original workstation and make the script fail before optimization starts.
SCRIPT_DIR = Path(__file__).resolve().parent
POST_SCRIPT = SCRIPT_DIR / "PoST_Simulation.py"
INPUT_PATH = SCRIPT_DIR / "input.json"
# This is the original, constraint-defining input.  Do not use a previously
# optimized input.json as the reference, or R's 50--200% range will compound.
DEFAULT_REFERENCE_INPUT_PATH = SCRIPT_DIR / "input.json"
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
SIM_ANALYSIS_CUTOFF_HZ = 10000  # modelnoise.txt is filtered at this cutoff.
OPTIMIZATION_SAMPLES = 4096


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
    parser.add_argument("--reference-input", type=Path, default=DEFAULT_REFERENCE_INPUT_PATH,
                        help="Original JSON that defines all parameter bounds and fixed values.")
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
    """Weighted robust loss in log-ASD ratio, focused on high frequency."""

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


def parameter_bounds(original: dict, fixed_r_ohm: float | None = None):
    """Only parameters the user has allowed to vary are listed here."""
    bounds = {
        "C_abs": Bound(original["C_abs"] * 1e-3, original["C_abs"] * 1e3),
        "C_tes": Bound(original["C_tes"] * 1e-3, original["C_tes"] * 1e3),
        "G_abs-abs": Bound(original["G_abs-abs"] * 1e-3, original["G_abs-abs"] * 1e3),
        "G_abs-tes": Bound(original["G_abs-tes"] * 1e-3, original["G_abs-tes"] * 1e3),
        "G_tes-bath": Bound(original["G_tes-bath"] * 1e-3, original["G_tes-bath"] * 1e3),
        "T_bath": Bound(1e-6, original["T_c"]),
        # The 10/20 uA RT data is not the target operating current, so it is
        # used only as a broad condition-specific prior for the ~252 uA TES
        # operating point.  Do not let noise-only fitting push alpha to 1000.
        "alpha": Bound(50.0, 100.0),
        "beta": Bound(0.0, 5.0, logarithmic=False),
        # L is now intentionally variable.  Tighten this upper bound if known.
        "L": Bound(original["L"], 1e-5),
    }
    if fixed_r_ohm is None:
        bounds["R"] = Bound(original["R"] * 0.5, original["R"] * 2.0)
    return bounds


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
    return candidate


def optimize_case(
    args,
    original: dict,
    reference: dict,
    fit_freq: np.ndarray,
    target: np.ndarray,
    experimental_rate: float,
    post_filter_white_asd: float,
    operating_point: dict,
    work_dir: Path,
):
    """Optimize one fixed-R case derived from one R_SH branch."""

    fixed_r_ohm = float(operating_point["R_TES_ohm"])
    shunt_resistance_ohm = float(operating_point["R_SH_ohm"])
    case_original = original.copy()
    case_original["R"] = fixed_r_ohm
    # R_SH is metadata for the IV-derived R only; the TES noise model does not
    # consume it directly.
    case_original["R_SH"] = shunt_resistance_ohm

    work_dir.mkdir(parents=True, exist_ok=True)
    work_input_path = work_dir / "input.json"
    work_noise_path = work_dir / NOISE_DAT_PATH.name
    bounds = parameter_bounds(reference, fixed_r_ohm=fixed_r_ohm)
    keys, scipy_bounds = vector_bounds(bounds)

    initial = case_original.copy()
    initial["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    initial["rate"] = experimental_rate
    initial["post_filter_white_asd_A_rtHz"] = post_filter_white_asd
    initial.pop("readout_white_asd_A_rtHz", None)
    initial["samples"] = OPTIMIZATION_SAMPLES
    initial_x = encode(initial, keys, bounds)

    cache = {}
    evaluation_count = 0
    best_score = np.inf
    best_candidate = initial.copy()

    def objective(vector):
        nonlocal evaluation_count, best_score, best_candidate
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
        candidate["samples"] = OPTIMIZATION_SAMPLES
        evaluation_count += 1
        try:
            write_json_atomically(work_input_path, candidate)
            run_post(args.timeout, work_dir, work_noise_path)
            model = simulated_spectrum(candidate, fit_freq, work_noise_path)
            score = fit_score(model, target, fit_freq, args)
            print(
                f"R_SH={shunt_resistance_ohm * 1e3:.4f} mOhm, "
                f"evaluation {evaluation_count:4d}: {score:.6g}"
            )
        except Exception as error:
            print(
                f"R_SH={shunt_resistance_ohm * 1e3:.4f} mOhm, "
                f"evaluation {evaluation_count:4d}: rejected ({error})"
            )
            score = 1e12
        cache[cache_key] = score
        if score < best_score:
            best_score = score
            best_candidate = candidate.copy()
            print("  new best score:", best_score)
        return score

    print(
        f"\n=== R_SH={shunt_resistance_ohm * 1e3:.4f} mOhm: "
        f"R_TES={fixed_r_ohm * 1e3:.6f} mOhm ==="
    )
    print("Initial evaluation")
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

    print("\nPoST evaluations:", evaluation_count)
    print("Best score:", best_score)
    print("Best fitted parameters:")
    print(json.dumps({key: best_candidate[key] for key in keys}, indent=2))
    return {
        "R_SH_ohm": shunt_resistance_ohm,
        "R_TES_ohm": fixed_r_ohm,
        "iv_operating_point": operating_point,
        "best_score": float(best_score),
        "evaluations": evaluation_count,
        "best_candidate": best_candidate,
        "work_dir": str(work_dir),
    }


def main():
    args = arguments()
    for path in (
        POST_SCRIPT,
        INPUT_PATH,
        MODELNOISE_PATH,
        PULSE_CONFIG_PATH,
        args.reference_input,
        args.target_iv,
    ):
        if not path.is_file():
            raise FileNotFoundError(f"Missing required file: {path}")

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

    original = load_json(INPUT_PATH)  # state to restore when --apply-final is absent
    reference = load_json(args.reference_input)
    pulse_config = load_json(PULSE_CONFIG_PATH)
    experimental_rate = float(pulse_config["Readout"]["Rate"])
    if args.fit_max_hz >= experimental_rate / 2.0:
        raise ValueError(
            f"--fit-max-hz={args.fit_max_hz:g} must be below Nyquist "
            f"({experimental_rate / 2.0:g} Hz)"
        )
    post_filter_white_asd = float(
        original.get(
            "post_filter_white_asd_A_rtHz",
            original.get("readout_white_asd_A_rtHz", 0.0),
        )
    )
    backup = INPUT_PATH.with_name(f"input.before_scipy_optimization.{time.strftime('%Y%m%d-%H%M%S')}.json")
    shutil.copy2(INPUT_PATH, backup)
    print("backup:", backup)
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
            case_name = f"rsh_{shunt_resistance_ohm * 1e3:.4f}mohm".replace(".", "p")
            result = optimize_case(
                args,
                original,
                reference,
                fit_freq,
                target,
                experimental_rate,
                post_filter_white_asd,
                operating_point,
                sweep_root / case_name,
            )
            cases.append(result)

        best_case = min(cases, key=lambda row: row["best_score"])
        summary = {
            "target_iv": str(args.target_iv),
            "target_bias_uA": args.target_bias_ua,
            "r_shunt_sweep_mOhm": [float(value * 1e3) for value in shunt_values_ohm],
            "tes_resistance_is_fitted": False,
            "alpha_bounds": [50.0, 100.0],
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
            final_candidate["samples"] = original["samples"]
            write_json_atomically(INPUT_PATH, final_candidate)
            run_post(args.timeout, INPUT_PATH.parent, NOISE_DAT_PATH)
            print(
                "Best input.json applied and final .dat regenerated: "
                f"R_SH={best_case['R_SH_ohm'] * 1e3:.4f} mOhm, "
                f"R_TES={best_case['R_TES_ohm'] * 1e3:.6f} mOhm"
            )
        else:
            print("Original input.json was not changed; use --apply-final to keep the best branch.")
    except Exception:
        write_json_atomically(INPUT_PATH, original)
        print("Original input.json restored after error.")
        raise


if __name__ == "__main__":
    main()

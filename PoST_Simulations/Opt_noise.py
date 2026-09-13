"""Standalone PoST optimizer: differential evolution -> Powell.

The optimizer is specialized for the 2024-12-06 215 mK / 1400 uA target
case.  During optimization it evaluates the deterministic five-state TES noise
model through the fixed fourth-order 100 kHz analog Bessel, sampling/alias
fold, and the fixed 10 kHz digital analysis Bessel.
This avoids fitting 4096-sample
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
from scipy.optimize import differential_evolution, least_squares, minimize

from lib.tes_noise_model import (
    ELECTRICAL_LINK_MODEL_RC,
    linearized_matrix as tes_linearized_matrix,
    noise_components as tes_noise_components,
    operating_point as tes_operating_point,
    tes_johnson_voltage_asd,
)
from subScript.noise_measurement_model import (
    HARDWARE_BESSEL_CUTOFF_HZ,
    analysis_filter_magnitude,
    hardware_filter_magnitude,
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
# frequency decade, with an extra ramp above 40 kHz so both the 50--100 kHz
# shoulder and the 100--200 kHz tail contribute strongly to the solution.
FIT_MIN_HZ = 1_000.0
FIT_MAX_HZ = 200_000.0
FIT_WEIGHT_START_HZ = 40_000.0
FIT_HIGH_FREQUENCY_WEIGHT = 4.0
FIT_ROBUST_DELTA_DEX = 0.30
FIT_POINTS = 601
FIT_BAND_MEAN_PENALTY = 3.0
ABSOLUTE_ASD_REFERENCE_HZ = 1_000.0
# The exact target CH0 readout gain/current calibration is unresolved in the
# target-case provenance, so absolute ASD must not influence the production fit
# by default.  Keep the diagnostic path available for a future independently
# calibrated CH0 scale.
ABSOLUTE_ASD_WEIGHT_DEFAULT = 0.0
MEASURED_ASD_PA_TO_A = 1.0e-12
PRODUCTION_ASD_UA_TO_A = 1.0e-6
BETA_DIAGNOSTIC_GRID = (
    0.0,
    0.5,
    1.0,
    1.5,
    2.0,
    3.0,
    4.0,
    6.0,
    8.0,
    12.0,
)
JOHNSON_SOURCE_SCALE_DIAGNOSTIC_GRID = tuple(
    float(value) for value in np.linspace(0.0, 1.0, 11)
)
SOURCE_DIAGNOSTIC_FREQUENCIES_HZ = (
    1_000.0,
    5_000.0,
    10_000.0,
    20_000.0,
    40_000.0,
    70_000.0,
    100_000.0,
    150_000.0,
    200_000.0,
)
FIT_BANDS_HZ = (
    (1_000.0, 5_000.0, 1.0),
    (5_000.0, 15_000.0, 1.25),
    (15_000.0, 40_000.0, 1.25),
    (40_000.0, 100_000.0, 1.5),
    (100_000.0, 200_000.0, 2.0),
)
SIM_ANALYSIS_CUTOFF_HZ = 10_000.0  # experimental NoiseAnalysis cutoff.
TARGET_HARDWARE_BESSEL_ORDER = 4
# A quoted analog "100 kHz cutoff" is interpreted as the -3 dB frequency.
# The target-case hardware refinement in this repository also preferred
# scipy's magnitude-normalized Bessel convention over the legacy phase norm.
TARGET_HARDWARE_BESSEL_NORM = "mag"
TARGET_HARDWARE_BESSEL_CUTOFF_HZ = 100_000.0
T_BATH_FIT_HALF_WIDTH_K = 0.002
ALPHA_FIT_MAX = 200.0
L_FIT_MIN_H = 1.0e-10
L_FIT_MAX_H = 12.3e-9
R_RC_FIT_MIN_OHM = 1.0e-6
R_RC_FIT_MAX_OHM = 20.0e-3
F_RC_FIT_MIN_HZ = 5_000.0
F_RC_FIT_MAX_HZ = 500_000.0
R_RC_INITIAL_OHM = 1.0e-3
F_RC_INITIAL_HZ = 100_000.0

# Diagnostic-only profile used to distinguish a genuine finite-frequency RC
# relaxation from a near-static redistribution between R_l and R_rc.
RC_DEGENERACY_FIXED_RL_GRID_OHM = tuple(
    value_mohm * 1.0e-3
    for value_mohm in (2.5, 4.0, 6.0, 8.0, 10.0, 12.0)
)
RC_DEGENERACY_PARTITION_POINTS = 9
RC_DEGENERACY_PROFILE_MAXFEV = 120

# Broad effective-parameter priors from the Elmer single/dual-pixel geometry
# and material table.  They are intentionally permissive: the hand-built
# absorber/glue geometry is not precise, sub-kelvin Pb/Stycast transport can
# be boundary/ballistic limited, and the reduced five-state model has no
# explicit Stycast heat-capacity node.
R_L_FIT_MIN_OHM = 2.5e-3
R_L_FIT_MAX_OHM = 12.0e-3

# Static effective-series model used by default after the RC profile showed
# that f_rc runs to its 500 kHz ceiling and R_l/R_rc are strongly degenerate.
# This parameter is deliberately named separately from the physical load:
# it is an effective total series resistance, not a claimed decomposition of
# shunt/load/wiring/readout contributions.
R_SERIES_EFF_FIT_MIN_OHM = 2.5e-3
R_SERIES_EFF_FIT_MAX_OHM = 30.0e-3
N_FIT_MIN = 1.5
N_FIT_MAX = 6.0

# TES: 500 um x 500 um x (40+120) nm, rho=15695 kg/m3,
# cp=1.25e-3 J/(kg K) -> 7.85e-13 J/K.
C_TES_MATERIAL_J_PER_K = 7.8475e-13

# One Stycast pad in the Elmer geometry: d=498 um, t=20 um,
# rho=2400 kg/m3, cp=1.22e-3 J/(kg K) -> about 1.14e-11 J/K.
# Since the five-state model has no explicit glue node, C_tes is allowed to
# act as an effective local heat capacity ranging from mostly-bare TES to
# TES plus roughly one glue pad.
C_STYCAST_PAD_MATERIAL_J_PER_K = 1.1407e-11
C_STYCAST_FIT_MIN_J_PER_K = 0.03 * C_STYCAST_PAD_MATERIAL_J_PER_K
C_STYCAST_FIT_MAX_J_PER_K = 10.0 * C_STYCAST_PAD_MATERIAL_J_PER_K

# With Stycast represented explicitly, keep the TES heat capacity tied more
# closely to the bilayer rather than letting it absorb the glue heat capacity.
C_TES_FIT_MIN_J_PER_K = 0.25 * C_TES_MATERIAL_J_PER_K
C_TES_FIT_MAX_J_PER_K = 5.0 * C_TES_MATERIAL_J_PER_K

# Pb absorber: 20 mm x 1 mm x 0.7 mm, rho=9860 kg/m3,
# cp=3.26e-5 J/(kg K) -> 4.50e-9 J/K in the current Elmer table.
# The absorber machining/assembly uncertainty is large enough that roughly a
# factor-of-two geometric heat-capacity variation is credible.  Keep this
# quantity physical instead of allowing the normalized-noise fit to use a
# 0.02x--6x effective heat capacity as a nuisance degree of freedom.
C_ABS_MATERIAL_J_PER_K = 4.500104e-9
C_ABS_FIT_MIN_J_PER_K = 0.50 * C_ABS_MATERIAL_J_PER_K
C_ABS_FIT_MAX_J_PER_K = 2.00 * C_ABS_MATERIAL_J_PER_K

# End-to-end Pb conductance from k*A/L using k=1.68e-2 W/(m K),
# A=1 mm*0.7 mm and L=20 mm -> 5.88e-7 W/K.  At sub-kelvin temperature the
# effective conductance can be dominated by boundary-limited/ballistic paths,
# contact area, and nonuniform geometry, so use this only as a scale.
G_ABS_ABS_MATERIAL_W_PER_K = 5.88e-7
G_ABS_ABS_FIT_MIN_W_PER_K = 0.03 * G_ABS_ABS_MATERIAL_W_PER_K
G_ABS_ABS_FIT_MAX_W_PER_K = 30.0 * G_ABS_ABS_MATERIAL_W_PER_K

# One Stycast pad: diameter 498 um, thickness 20 um,
# k=2.69094e-6 W/(m K) -> about 2.62e-8 W/K.  The real cured thickness,
# contact area, cracks/voids, and low-T transport are poorly known, so permit
# a similarly broad effective TES--absorber conductance.
G_ABS_TES_MATERIAL_W_PER_K = 2.622e-8
G_TES_STYCAST_FIT_MIN_W_PER_K = 0.03 * G_ABS_TES_MATERIAL_W_PER_K
G_TES_STYCAST_FIT_MAX_W_PER_K = 30.0 * G_ABS_TES_MATERIAL_W_PER_K
G_STYCAST_ABS_FIT_MIN_W_PER_K = 0.03 * G_ABS_TES_MATERIAL_W_PER_K
G_STYCAST_ABS_FIT_MAX_W_PER_K = 30.0 * G_ABS_TES_MATERIAL_W_PER_K

POST_FILTER_WHITE_FRACTION_MIN = 1.0e-6
POST_FILTER_WHITE_FRACTION_MAX = 1.0
POST_FILTER_WHITE_FRACTION_INITIAL = 1.0e-4
EXCESS_JOHNSON_MAX = 5.0


@dataclass(frozen=True)
class Bound:
    lower: float
    upper: float
    logarithmic: bool = True


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-maxiter", type=int, default=45)
    parser.add_argument("--de-popsize", type=int, default=10)
    parser.add_argument("--powell-maxfev", type=int, default=1000)
    parser.add_argument("--least-squares-max-nfev", type=int, default=350)
    parser.add_argument("--skip-de", action="store_true")
    parser.add_argument("--apply-final", action="store_true")
    parser.add_argument(
        "--use-rc-relaxation",
        action="store_true",
        help=(
            "Use the previous two-parameter passive RC relaxation branch. "
            "Default is the simpler static effective-series-resistance model."
        ),
    )
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
        "--absolute-asd-weight",
        type=float,
        default=ABSOLUTE_ASD_WEIGHT_DEFAULT,
        help=(
            "Optional penalty weight for the reported CH0 ASD level at 1 kHz. "
            "The target-case absolute readout calibration is unresolved, so the "
            "default is 0 (diagnostic only). Use a positive value only after an "
            "independent CH0 current calibration has been established."
        ),
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
    parser.add_argument("--r-shunt-start-mohm", type=float, default=3.85,
                        help="First shunt-resistance branch in mOhm.")
    parser.add_argument("--r-shunt-stop-mohm", type=float, default=3.85,
                        help="Last shunt-resistance branch in mOhm.")
    parser.add_argument("--r-shunt-count", type=int, default=1,
                        help=(
                            "Number of equally spaced R_SH branches. The default "
                            "focuses the deeper search on the repeatedly selected "
                            "3.85 mOhm branch; pass 3.8/3.9/3 to restore the sweep."
                        ))
    return parser.parse_args()


def load_json(path: Path):
    with path.open(encoding="utf-8") as f:
        return json.load(f)


BEST_PARAMETER_KEYS = (
    "T_c",
    "T_bath",
    "R",
    "R_SH",
    "R_l",
    "R_series_eff",
    "alpha",
    "beta",
    "L",
    "R_rc",
    "f_rc_Hz",
    "n",
    "C_tes",
    "C_stycast",
    "C_abs",
    "G_tes-bath",
    "G_tes-stycast",
    "G_stycast-abs",
    "G_abs-abs",
    "excess_johnson_M",
    "post_filter_white_fraction",
    "post_filter_white_asd_A_rtHz",
    "hardware_bessel_order",
    "hardware_bessel_norm",
    "hardware_bessel_cutoff_Hz",
    "thermal_link_model",
    "electrical_link_model",
    "rate",
    "samples",
    "cutoff",
)


def best_parameters_for_summary(candidate: dict) -> dict:
    """Return the physically relevant best-fit candidate fields as JSON scalars."""

    result = {}
    for key in BEST_PARAMETER_KEYS:
        if key not in candidate:
            continue
        value = candidate[key]
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, (str, bool, int, float)) or value is None:
            result[key] = value
    return result


def target_case_reference(
    case_dir: Path,
    explicit_reference: Path | None = None,
    use_rc_relaxation: bool = False,
):
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

    # The target acquisition is the 215 mK case. Use the setpoint as the
    # starting value, while the optimizer may move it slightly within the
    # configured target-stage uncertainty.
    reference["T_bath"] = float(envelope["parameters"]["T_bath"]["nominal"])
    reference["rate"] = 500_000.0
    reference["samples"] = 100_000
    reference["hardware_bessel_order"] = TARGET_HARDWARE_BESSEL_ORDER
    reference["hardware_bessel_norm"] = TARGET_HARDWARE_BESSEL_NORM
    reference["hardware_bessel_cutoff_Hz"] = TARGET_HARDWARE_BESSEL_CUTOFF_HZ
    reference.setdefault(
        "post_filter_white_fraction",
        POST_FILTER_WHITE_FRACTION_INITIAL,
    )
    reference["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    reference.setdefault("excess_johnson_M", 0.0)
    if use_rc_relaxation:
        reference.setdefault("R_rc", R_RC_INITIAL_OHM)
        reference.setdefault("f_rc_Hz", F_RC_INITIAL_HZ)
        reference["electrical_link_model"] = ELECTRICAL_LINK_MODEL_RC
    else:
        reference["R_series_eff"] = float(reference["R_l"])
        reference.pop("R_rc", None)
        reference.pop("f_rc_Hz", None)
        reference["electrical_link_model"] = "rl"
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


def log_ratio_residual(
    model: np.ndarray,
    target: np.ndarray,
) -> np.ndarray:
    model = np.asarray(model, dtype=float)
    target = np.asarray(target, dtype=float)
    if (
        np.any(~np.isfinite(model))
        or np.any(~np.isfinite(target))
        or np.any(model <= 0.0)
        or np.any(target <= 0.0)
    ):
        raise ValueError("Fit spectra must be finite and strictly positive")
    return np.log10(model) - np.log10(target)


def robust_point_loss(residual: np.ndarray, delta: float) -> np.ndarray:
    residual = np.asarray(residual, dtype=float)
    if delta <= 0.0:
        return residual**2
    absolute = np.abs(residual)
    return np.where(
        absolute <= delta,
        residual**2,
        2.0 * delta * absolute - delta**2,
    )


def band_balanced_residual_weights(
    fit_freq: np.ndarray,
    args,
) -> np.ndarray:
    """Weights whose squared residuals reproduce an equalized band objective."""

    fit_freq = np.asarray(fit_freq, dtype=float)
    base = fit_frequency_weights(fit_freq, args)
    weights = np.zeros_like(fit_freq)

    active = []
    for low, high, band_weight in FIT_BANDS_HZ:
        low = max(float(low), float(args.fit_min_hz))
        high = min(float(high), float(args.fit_max_hz))
        mask = (fit_freq >= low) & (
            fit_freq <= high if np.isclose(high, args.fit_max_hz) else fit_freq < high
        )
        count = int(np.count_nonzero(mask))
        if count:
            active.append((mask, float(band_weight), count))

    if not active:
        return base / max(float(np.sum(base)), 1e-30)

    total_band_weight = sum(item[1] for item in active)
    for mask, band_weight, _count in active:
        local = base[mask]
        local_sum = float(np.sum(local))
        if local_sum <= 0.0:
            continue
        weights[mask] = (
            band_weight / total_band_weight
        ) * local / local_sum

    missing = weights <= 0.0
    if np.any(missing):
        fallback = base[missing]
        fallback_sum = float(np.sum(fallback))
        if fallback_sum > 0.0:
            weights[missing] = 1e-6 * fallback / fallback_sum
    weights /= np.sum(weights)
    return weights


def band_mean_residuals(
    residual: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> np.ndarray:
    residual = np.asarray(residual, dtype=float)
    fit_freq = np.asarray(fit_freq, dtype=float)
    means = []
    for low, high, _band_weight in FIT_BANDS_HZ:
        low_eff = max(float(low), float(args.fit_min_hz))
        high_eff = min(float(high), float(args.fit_max_hz))
        mask = (fit_freq >= low_eff) & (fit_freq <= high_eff)
        if np.any(mask):
            means.append(float(np.mean(residual[mask])))
    return np.asarray(means, dtype=float)


def weighted_residual_vector(
    model: np.ndarray,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> np.ndarray:
    residual = log_ratio_residual(model, target)
    weights = band_balanced_residual_weights(fit_freq, args)
    # Preserve the residual's natural dex scale for scipy least_squares while
    # retaining the same relative weighting as the scalar band-balanced loss.
    point_residual = residual * np.sqrt(weights * len(weights))
    # Add broad-band mean residuals as pseudo-observations. This specifically
    # suppresses the smooth, alternating over/under-shoot seen in the ratio
    # plot without trying to chase narrow experimental spikes.
    means = band_mean_residuals(residual, fit_freq, args)
    mean_residual = means * np.sqrt(FIT_BAND_MEAN_PENALTY)
    return np.concatenate((point_residual, mean_residual))


def fit_score(
    model: np.ndarray,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> float:
    """Band-balanced robust loss in log-ASD ratio."""

    residual = log_ratio_residual(model, target)
    loss = robust_point_loss(residual, float(args.robust_delta_dex))
    weights = band_balanced_residual_weights(fit_freq, args)
    point_score = float(np.sum(weights * loss))
    means = band_mean_residuals(residual, fit_freq, args)
    mean_score = (
        float(np.mean(means**2)) * FIT_BAND_MEAN_PENALTY
        if len(means)
        else 0.0
    )
    return point_score + mean_score


def absolute_level_residual(
    model_reference_asd_A: float,
    target_reference_asd_A: float,
) -> float:
    """Return log10(model/measurement) for the absolute ASD anchor."""

    model = float(model_reference_asd_A)
    target = float(target_reference_asd_A)
    if (
        not np.isfinite(model)
        or not np.isfinite(target)
        or model <= 0.0
        or target <= 0.0
    ):
        raise ValueError("Absolute ASD levels must be finite and positive")
    return float(np.log10(model) - np.log10(target))


def combined_fit_score(
    model: np.ndarray,
    target: np.ndarray,
    fit_freq: np.ndarray,
    model_reference_asd_A: float,
    target_reference_asd_A: float,
    args,
) -> float:
    """Shape loss plus one independent absolute-ASD level anchor."""

    shape_score = fit_score(model, target, fit_freq, args)
    level = absolute_level_residual(
        model_reference_asd_A,
        target_reference_asd_A,
    )
    return shape_score + float(args.absolute_asd_weight) * level**2


def combined_residual_vector(
    model: np.ndarray,
    target: np.ndarray,
    fit_freq: np.ndarray,
    model_reference_asd_A: float,
    target_reference_asd_A: float,
    args,
) -> np.ndarray:
    """Least-squares residuals with a single absolute-level pseudo-observation."""

    shape = weighted_residual_vector(model, target, fit_freq, args)
    level = absolute_level_residual(
        model_reference_asd_A,
        target_reference_asd_A,
    )
    level_residual = np.asarray(
        [level * np.sqrt(float(args.absolute_asd_weight))],
        dtype=float,
    )
    return np.concatenate((shape, level_residual))


def band_fit_diagnostics(
    model: np.ndarray,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> dict:
    residual = log_ratio_residual(model, target)
    diagnostics = {}
    for low, high, band_weight in FIT_BANDS_HZ:
        low_eff = max(float(low), float(args.fit_min_hz))
        high_eff = min(float(high), float(args.fit_max_hz))
        mask = (fit_freq >= low_eff) & (fit_freq <= high_eff)
        if not np.any(mask):
            continue
        key = f"{low_eff:g}-{high_eff:g}_Hz"
        diagnostics[key] = {
            "rms_log10_ratio": float(np.sqrt(np.mean(residual[mask] ** 2))),
            "mean_log10_ratio": float(np.mean(residual[mask])),
            "max_abs_log10_ratio": float(np.max(np.abs(residual[mask]))),
            "band_weight": float(band_weight),
        }
    return diagnostics


def target_spectrum(args):
    pulse_config = load_json(PULSE_CONFIG_PATH)
    # Analyze_Experimental_Data/noise_main.py writes modelnoise.txt in
    # pA/rtHz after eta calibration.  Preserve that scale for normalization,
    # but convert the absolute-level anchor explicitly to A/rtHz.
    exp_asd_pA = np.asarray(np.loadtxt(MODELNOISE_PATH), dtype=float)
    exp_rate = float(pulse_config["Readout"]["Rate"])
    exp_freq = np.arange(len(exp_asd_pA)) * (exp_rate / 2.0) / len(exp_asd_pA)
    fit_freq = np.geomspace(
        float(args.fit_min_hz),
        float(args.fit_max_hz),
        int(args.fit_points),
    )
    target = np.interp(
        fit_freq,
        exp_freq,
        normalize_at(exp_freq, exp_asd_pA),
    )
    target_reference_asd_A = float(
        np.interp(ABSOLUTE_ASD_REFERENCE_HZ, exp_freq, exp_asd_pA)
        * MEASURED_ASD_PA_TO_A
    )
    if not np.isfinite(target_reference_asd_A) or target_reference_asd_A <= 0.0:
        raise ValueError("Invalid measured absolute ASD near 1 kHz")
    return fit_freq, target, target_reference_asd_A


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
    # PoST_Simulation.py persists noise_total-bessel100k.dat in uA/rtHz.
    sim_asd_uA = np.asarray(np.loadtxt(noise_path), dtype=float)
    sim_rate = float(candidate["rate"])
    sim_freq = np.arange(len(sim_asd_uA)) * (sim_rate / 2.0) / len(sim_asd_uA)
    normalized = np.interp(
        fit_freq,
        sim_freq,
        normalize_at(sim_freq, sim_asd_uA),
    )
    reference_asd_A = float(
        np.interp(ABSOLUTE_ASD_REFERENCE_HZ, sim_freq, sim_asd_uA)
        * PRODUCTION_ASD_UA_TO_A
    )
    return normalized, reference_asd_A


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


def parameter_bounds(
    reference: dict,
    envelope: dict,
    fixed_r_ohm: float,
    use_rc_relaxation: bool = False,
):
    """Target-case search box derived from the frozen proxy envelope.

    R is fixed by the same-campaign IV conversion for each R_SH branch.
    T_bath is allowed a small offset around the 215 mK setpoint, while T_c is
    kept inside the independent nearby-run RT envelope.  The remaining
    parameters cover the stable target-case region already exercised by the
    repository's adaptive/correlated searches.
    """

    sensitivity = envelope["sensitivity_reference"]
    tc_low, tc_high = map(float, envelope["parameters"]["T_c"]["range"])
    tbath_nominal = float(envelope["parameters"]["T_bath"]["nominal"])

    def ref(name):
        return float(sensitivity[name])

    bounds = {
        # Keep T_c inside the independently motivated nearby-run RT envelope.
        "T_c": Bound(tc_low, tc_high, logarithmic=False),
        "T_bath": Bound(
            tbath_nominal - T_BATH_FIT_HALF_WIDTH_K,
            tbath_nominal + T_BATH_FIT_HALF_WIDTH_K,
            logarithmic=False,
        ),
        "alpha": Bound(ref("alpha") * 0.05, ALPHA_FIT_MAX, logarithmic=False),
        "beta": Bound(0.0, 12.0, logarithmic=False),
        "L": Bound(L_FIT_MIN_H, L_FIT_MAX_H),
        "n": Bound(N_FIT_MIN, N_FIT_MAX, logarithmic=False),
        "C_tes": Bound(C_TES_FIT_MIN_J_PER_K, C_TES_FIT_MAX_J_PER_K),
        "C_stycast": Bound(
            C_STYCAST_FIT_MIN_J_PER_K,
            C_STYCAST_FIT_MAX_J_PER_K,
        ),
        "C_abs": Bound(C_ABS_FIT_MIN_J_PER_K, C_ABS_FIT_MAX_J_PER_K),
        "G_tes-stycast": Bound(
            G_TES_STYCAST_FIT_MIN_W_PER_K,
            G_TES_STYCAST_FIT_MAX_W_PER_K,
        ),
        "G_stycast-abs": Bound(
            G_STYCAST_ABS_FIT_MIN_W_PER_K,
            G_STYCAST_ABS_FIT_MAX_W_PER_K,
        ),
        "G_abs-abs": Bound(
            G_ABS_ABS_FIT_MIN_W_PER_K,
            G_ABS_ABS_FIT_MAX_W_PER_K,
        ),
        "excess_johnson_M": Bound(0.0, EXCESS_JOHNSON_MAX, logarithmic=False),
        "post_filter_white_fraction": Bound(
            POST_FILTER_WHITE_FRACTION_MIN,
            POST_FILTER_WHITE_FRACTION_MAX,
        ),
    }
    if use_rc_relaxation:
        bounds["R_l"] = Bound(R_L_FIT_MIN_OHM, R_L_FIT_MAX_OHM)
        bounds["R_rc"] = Bound(R_RC_FIT_MIN_OHM, R_RC_FIT_MAX_OHM)
        bounds["f_rc_Hz"] = Bound(F_RC_FIT_MIN_HZ, F_RC_FIT_MAX_HZ)
    else:
        bounds["R_series_eff"] = Bound(
            R_SERIES_EFF_FIT_MIN_OHM,
            R_SERIES_EFF_FIT_MAX_OHM,
        )
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


def clip_candidate_to_bounds(candidate: dict, keys, bounds):
    """Move an initial/warm-start point inside the configured physical box."""

    for key in keys:
        bound = bounds[key]
        value = float(candidate[key])
        if value < bound.lower or value > bound.upper:
            clipped = float(np.clip(value, bound.lower, bound.upper))
            print(
                f"Clipping initial {key}: {value:.6g} -> {clipped:.6g} "
                f"(allowed {bound.lower:.6g}--{bound.upper:.6g})"
            )
            candidate[key] = clipped
    return candidate


def encode(candidate: dict, keys, bounds):
    return np.asarray([
        np.log10(candidate[key]) if bounds[key].logarithmic else candidate[key]
        for key in keys
    ], dtype=float)


def parameter_boundary_diagnostics(
    candidate: dict,
    keys,
    bounds: dict,
    near_fraction: float = 0.01,
) -> dict:
    """Report where each fitted parameter lands inside its configured bound."""

    diagnostics = {}
    for key in keys:
        bound = bounds[key]
        value = float(candidate[key])
        if bound.logarithmic:
            value_t = np.log10(value)
            lower_t = np.log10(bound.lower)
            upper_t = np.log10(bound.upper)
        else:
            value_t = value
            lower_t = float(bound.lower)
            upper_t = float(bound.upper)
        span = upper_t - lower_t
        position = (value_t - lower_t) / span if span > 0.0 else 0.5
        position = float(np.clip(position, 0.0, 1.0))
        lower_distance = position
        upper_distance = 1.0 - position
        nearest = "lower" if lower_distance <= upper_distance else "upper"
        nearest_distance = min(lower_distance, upper_distance)
        diagnostics[key] = {
            "value": value,
            "position_fraction": position,
            "distance_to_lower_fraction": lower_distance,
            "distance_to_upper_fraction": upper_distance,
            "nearest_bound": nearest,
            "distance_to_nearest_bound_fraction": nearest_distance,
            "within_1pct_of_bound": bool(nearest_distance <= near_fraction),
            "logarithmic_position": bool(bound.logarithmic),
        }
    return diagnostics


def post_analysis_source_class_asd(
    candidate: dict,
    frequencies_hz: np.ndarray,
) -> dict:
    """Decompose CH0 ASD by physical source class through the full filters.

    Each intrinsic source class is propagated through the same 100 kHz analog
    Bessel, first ADC alias fold, and 10 kHz digital analysis response used by
    the deterministic optimizer. The empirical post-filter white term is kept
    as its own class so the PSD fractions reconstruct the full model.
    """

    frequency = np.asarray(frequencies_hz, dtype=float)
    rate = float(candidate["rate"])
    if np.any(frequency < 0.0) or np.any(frequency > rate / 2.0):
        raise ValueError("source-diagnostic frequencies must lie below Nyquist")

    diagnostic_candidate = candidate.copy()
    apply_post_filter_white_fraction(diagnostic_candidate)
    alias_frequency = rate - frequency
    query = np.unique(np.concatenate((frequency, alias_frequency)))
    intrinsic = tes_noise_components(diagnostic_candidate, query)

    main_response = hardware_filter_magnitude(
        frequency,
        cutoff_hz=TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=TARGET_HARDWARE_BESSEL_ORDER,
        norm=TARGET_HARDWARE_BESSEL_NORM,
    )
    alias_response = hardware_filter_magnitude(
        alias_frequency,
        cutoff_hz=TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=TARGET_HARDWARE_BESSEL_ORDER,
        norm=TARGET_HARDWARE_BESSEL_NORM,
    )
    analysis_response = analysis_filter_magnitude(
        frequency,
        rate,
        cutoff_hz=SIM_ANALYSIS_CUTOFF_HZ,
    )
    same_bin = np.isclose(
        alias_frequency,
        frequency,
        rtol=0.0,
        atol=max(rate, 1.0) * 1e-12,
    )

    class_asd = {}
    for name, intrinsic_asd in intrinsic["aggregated_components_ch0"].items():
        intrinsic_asd = np.asarray(intrinsic_asd, dtype=float)
        main = np.interp(frequency, query, intrinsic_asd) * main_response
        alias = np.interp(alias_frequency, query, intrinsic_asd) * alias_response
        alias = np.where(same_bin, 0.0, alias)
        class_asd[name] = np.sqrt(main**2 + alias**2) * analysis_response

    white = float(diagnostic_candidate.get("post_filter_white_asd_A_rtHz", 0.0))
    class_asd["post_filter_white"] = np.full_like(
        frequency,
        white,
        dtype=float,
    ) * analysis_response

    reconstructed = np.sqrt(
        np.sum(
            np.asarray([curve**2 for curve in class_asd.values()], dtype=float),
            axis=0,
        )
    )
    detector = hardware_sampled_asd(
        diagnostic_candidate,
        frequency,
        rate_hz=rate,
        cutoff_hz=TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=TARGET_HARDWARE_BESSEL_ORDER,
        norm=TARGET_HARDWARE_BESSEL_NORM,
    )
    expected = np.sqrt(detector**2 + white**2) * analysis_response
    denominator = np.maximum(expected, np.finfo(float).tiny)
    consistency = float(np.max(np.abs(reconstructed - expected) / denominator))

    return {
        "frequencies_Hz": frequency,
        "class_asd_A_rtHz": class_asd,
        "total_asd_A_rtHz": reconstructed,
        "reconstruction_max_relative_error": consistency,
    }


def source_class_diagnostics(
    candidate: dict,
    fit_freq: np.ndarray,
    args,
) -> tuple[dict, dict]:
    """Summarize which physical noise classes dominate the fitted spectrum."""

    curves = post_analysis_source_class_asd(candidate, fit_freq)
    frequency = curves["frequencies_Hz"]
    class_asd = curves["class_asd_A_rtHz"]
    total_asd = curves["total_asd_A_rtHz"]
    total_psd = np.maximum(total_asd**2, np.finfo(float).tiny)
    fractions = {
        name: np.asarray(asd, dtype=float) ** 2 / total_psd
        for name, asd in class_asd.items()
    }

    points = {}
    for requested in SOURCE_DIAGNOSTIC_FREQUENCIES_HZ:
        if requested < frequency[0] or requested > frequency[-1]:
            continue
        source_rows = {}
        for name in class_asd:
            asd_value = float(np.interp(requested, frequency, class_asd[name]))
            fraction_value = float(np.interp(requested, frequency, fractions[name]))
            source_rows[name] = {
                "asd_A_rtHz": asd_value,
                "psd_fraction": fraction_value,
            }
        dominant = max(
            source_rows,
            key=lambda name: source_rows[name]["psd_fraction"],
        )
        points[f"{requested:g}_Hz"] = {
            "total_model_asd_A_rtHz": float(
                np.interp(requested, frequency, total_asd)
            ),
            "dominant_source_class": dominant,
            "dominant_psd_fraction": float(
                source_rows[dominant]["psd_fraction"]
            ),
            "sources": source_rows,
        }

    bands = {}
    for low, high, _weight in FIT_BANDS_HZ:
        low_eff = max(float(low), float(args.fit_min_hz))
        high_eff = min(float(high), float(args.fit_max_hz))
        mask = (frequency >= low_eff) & (frequency <= high_eff)
        if not np.any(mask):
            continue
        source_rows = {}
        for name in class_asd:
            source_rows[name] = {
                "mean_psd_fraction": float(np.mean(fractions[name][mask])),
                "median_psd_fraction": float(np.median(fractions[name][mask])),
                "rms_asd_A_rtHz": float(
                    np.sqrt(np.mean(class_asd[name][mask] ** 2))
                ),
            }
        dominant = max(
            source_rows,
            key=lambda name: source_rows[name]["mean_psd_fraction"],
        )
        bands[f"{low_eff:g}-{high_eff:g}_Hz"] = {
            "dominant_source_class": dominant,
            "dominant_mean_psd_fraction": float(
                source_rows[dominant]["mean_psd_fraction"]
            ),
            "sources": source_rows,
        }

    summary = {
        "measurement_chain": (
            "intrinsic source class -> fixed 100 kHz 4th-order mag Bessel "
            "-> first ADC alias fold -> fixed 10 kHz digital analysis response"
        ),
        "source_classes": list(class_asd),
        "reconstruction_max_relative_error": float(
            curves["reconstruction_max_relative_error"]
        ),
        "representative_frequencies": points,
        "fit_bands": bands,
    }
    return summary, curves


def source_ablation_diagnostics(
    curves: dict,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> dict:
    """Report the normalized fit if one source class is removed in PSD."""

    frequency = np.asarray(curves["frequencies_Hz"], dtype=float)
    total_asd = np.asarray(curves["total_asd_A_rtHz"], dtype=float)
    target = np.asarray(target, dtype=float)
    fit_freq = np.asarray(fit_freq, dtype=float)
    if frequency.shape != fit_freq.shape or not np.allclose(
        frequency, fit_freq, rtol=0.0, atol=1e-9
    ):
        raise ValueError("source-ablation curves must use the optimizer fit grid")
    if target.shape != fit_freq.shape:
        raise ValueError("source-ablation target must use the optimizer fit grid")

    full_normalized = normalize_at(
        frequency, total_asd, reference_hz=ABSOLUTE_ASD_REFERENCE_HZ
    )
    full_score = float(fit_score(full_normalized, target, fit_freq, args))
    result = {
        "diagnostic_only": True,
        "normalization_Hz": float(ABSOLUTE_ASD_REFERENCE_HZ),
        "full_model_shape_score": full_score,
        "remove_one_source_class": {},
    }

    total_psd = total_asd**2
    tiny = np.finfo(float).tiny
    for name, class_asd in curves["class_asd_A_rtHz"].items():
        class_asd = np.asarray(class_asd, dtype=float)
        ablated_asd = np.sqrt(np.maximum(total_psd - class_asd**2, tiny))
        ablated_normalized = normalize_at(
            frequency, ablated_asd, reference_hz=ABSOLUTE_ASD_REFERENCE_HZ
        )
        residual = log_ratio_residual(ablated_normalized, target)
        score = float(fit_score(ablated_normalized, target, fit_freq, args))
        bands = {}
        for low, high, _weight in FIT_BANDS_HZ:
            low_eff = max(float(low), float(args.fit_min_hz))
            high_eff = min(float(high), float(args.fit_max_hz))
            mask = (fit_freq >= low_eff) & (fit_freq <= high_eff)
            if not np.any(mask):
                continue
            mean_log = float(np.mean(residual[mask]))
            bands[f"{low_eff:g}-{high_eff:g}_Hz"] = {
                "rms_log10_ratio": float(np.sqrt(np.mean(residual[mask] ** 2))),
                "mean_log10_ratio": mean_log,
                "geometric_mean_model_over_measurement": float(10.0 ** mean_log),
            }
        result["remove_one_source_class"][name] = {
            "shape_score": score,
            "score_change_vs_full": float(score - full_score),
            "bands": bands,
        }
    return result


def required_transfer_diagnostics(
    candidate: dict,
    model: np.ndarray,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> tuple[dict, dict]:
    """Infer the multiplicative ASD correction required after the current model.

    The diagnostic compares the normalized measured target with the normalized
    deterministic best-fit model.  It does not fit or apply a new transfer
    function.  H_required = measured / model, so values below unity mean that
    the present model would need additional attenuation at that frequency.

    Existing hardware/analysis responses and the first alias-fold fraction are
    reported only as context; they are already included in the model.
    """

    frequency = np.asarray(fit_freq, dtype=float)
    model = np.asarray(model, dtype=float)
    target = np.asarray(target, dtype=float)
    if model.shape != frequency.shape or target.shape != frequency.shape:
        raise ValueError(
            "required-transfer model, target, and frequency must have identical shapes"
        )
    if (
        np.any(~np.isfinite(model))
        or np.any(~np.isfinite(target))
        or np.any(model <= 0.0)
        or np.any(target <= 0.0)
    ):
        raise ValueError(
            "required-transfer model and target must be finite and positive"
        )

    required = target / model
    correction_db = 20.0 * np.log10(required)
    model_over_measurement = model / target

    rate = float(candidate["rate"])
    hardware_main = hardware_filter_magnitude(
        frequency,
        cutoff_hz=TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=TARGET_HARDWARE_BESSEL_ORDER,
        norm=TARGET_HARDWARE_BESSEL_NORM,
    )
    analysis = analysis_filter_magnitude(
        frequency,
        rate,
        cutoff_hz=SIM_ANALYSIS_CUTOFF_HZ,
    )
    combined_direct = hardware_main * analysis

    def normalized_response(response: np.ndarray) -> np.ndarray:
        response = np.asarray(response, dtype=float)
        reference = float(
            np.interp(ABSOLUTE_ASD_REFERENCE_HZ, frequency, response)
        )
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValueError("invalid transfer response at normalization frequency")
        return response / reference

    hardware_main_normalized = normalized_response(hardware_main)
    analysis_normalized = normalized_response(analysis)
    combined_direct_normalized = normalized_response(combined_direct)

    alias_frequency = rate - frequency
    query = np.unique(np.concatenate((frequency, alias_frequency)))
    intrinsic = tes_noise_components(candidate, query)
    intrinsic_total = np.asarray(intrinsic["total_ch0"], dtype=float)
    main_intrinsic = np.interp(frequency, query, intrinsic_total)
    alias_intrinsic = np.interp(alias_frequency, query, intrinsic_total)
    alias_hardware = hardware_filter_magnitude(
        alias_frequency,
        cutoff_hz=TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=TARGET_HARDWARE_BESSEL_ORDER,
        norm=TARGET_HARDWARE_BESSEL_NORM,
    )
    main_after_hardware = main_intrinsic * hardware_main
    alias_after_hardware = alias_intrinsic * alias_hardware
    same_bin = np.isclose(
        alias_frequency,
        frequency,
        rtol=0.0,
        atol=max(rate, 1.0) * 1e-12,
    )
    alias_after_hardware = np.where(same_bin, 0.0, alias_after_hardware)
    folded_psd = main_after_hardware**2 + alias_after_hardware**2
    alias_psd_fraction = np.divide(
        alias_after_hardware**2,
        folded_psd,
        out=np.zeros_like(folded_psd),
        where=folded_psd > 0.0,
    )

    points = {}
    for requested in SOURCE_DIAGNOSTIC_FREQUENCIES_HZ:
        if requested < frequency[0] or requested > frequency[-1]:
            continue
        correction = float(np.interp(requested, frequency, required))
        points[f"{requested:g}_Hz"] = {
            "required_ASD_transfer_measurement_over_model": correction,
            "required_correction_dB": float(20.0 * np.log10(correction)),
            "model_over_measurement": float(
                np.interp(requested, frequency, model_over_measurement)
            ),
            "already_modeled_hardware_main_normalized": float(
                np.interp(
                    requested,
                    frequency,
                    hardware_main_normalized,
                )
            ),
            "already_modeled_analysis_normalized": float(
                np.interp(
                    requested,
                    frequency,
                    analysis_normalized,
                )
            ),
            "already_modeled_combined_direct_normalized": float(
                np.interp(
                    requested,
                    frequency,
                    combined_direct_normalized,
                )
            ),
            "first_alias_fold_psd_fraction": float(
                np.interp(
                    requested,
                    frequency,
                    alias_psd_fraction,
                )
            ),
        }

    bands = {}
    log_required = np.log10(required)
    for low, high, _weight in FIT_BANDS_HZ:
        low_eff = max(float(low), float(args.fit_min_hz))
        high_eff = min(float(high), float(args.fit_max_hz))
        mask = (frequency >= low_eff) & (frequency <= high_eff)
        if not np.any(mask):
            continue
        mean_log = float(np.mean(log_required[mask]))
        bands[f"{low_eff:g}-{high_eff:g}_Hz"] = {
            "geometric_mean_required_ASD_transfer": float(10.0 ** mean_log),
            "mean_required_correction_dB": float(20.0 * mean_log),
            "median_required_ASD_transfer": float(np.median(required[mask])),
            "min_required_ASD_transfer": float(np.min(required[mask])),
            "max_required_ASD_transfer": float(np.max(required[mask])),
            "mean_first_alias_fold_psd_fraction": float(
                np.mean(alias_psd_fraction[mask])
            ),
        }

    deviation_onset = {}
    for fractional_deviation in (0.10, 0.20, 0.50):
        mask = np.abs(required - 1.0) >= fractional_deviation
        indices = np.flatnonzero(mask)
        deviation_onset[
            f"{int(round(fractional_deviation * 100.0))}pct"
        ] = (
            float(frequency[indices[0]]) if len(indices) else None
        )

    sample_indices = np.unique(
        np.linspace(
            0,
            len(frequency) - 1,
            min(81, len(frequency)),
            dtype=int,
        )
    )
    curve_sample = [
        {
            "frequency_Hz": float(frequency[index]),
            "required_ASD_transfer": float(required[index]),
            "required_correction_dB": float(correction_db[index]),
            "first_alias_fold_psd_fraction": float(
                alias_psd_fraction[index]
            ),
        }
        for index in sample_indices
    ]

    summary = {
        "diagnostic_only": True,
        "definition": "H_required(f) = measured_normalized_ASD / model_normalized_ASD",
        "interpretation": (
            "H_required < 1 means additional attenuation would be required; "
            "H_required > 1 means additional gain would be required. "
            "No correction is applied to the optimizer."
        ),
        "normalization_Hz": float(ABSOLUTE_ASD_REFERENCE_HZ),
        "known_chain_context": (
            "100 kHz 4th-order mag Bessel, first ADC alias fold, and "
            "10 kHz digital analysis response are already included"
        ),
        "representative_frequencies": points,
        "fit_bands": bands,
        "first_frequency_exceeding_fractional_deviation_from_unity_Hz": (
            deviation_onset
        ),
        "curve_sample": curve_sample,
    }
    curves = {
        "frequencies_Hz": frequency,
        "required_ASD_transfer": required,
        "required_correction_dB": correction_db,
        "hardware_main_normalized": hardware_main_normalized,
        "analysis_normalized": analysis_normalized,
        "combined_direct_normalized": combined_direct_normalized,
        "first_alias_fold_psd_fraction": alias_psd_fraction,
    }
    return summary, curves


def plot_required_transfer_diagnostics(
    curves: dict,
    output_path: Path,
) -> Path:
    """Plot the inferred residual transfer beside already-modeled responses."""

    frequency = np.asarray(curves["frequencies_Hz"], dtype=float)
    correction_db = np.asarray(
        curves["required_correction_dB"],
        dtype=float,
    )
    hardware = np.asarray(
        curves["hardware_main_normalized"],
        dtype=float,
    )
    analysis = np.asarray(
        curves["analysis_normalized"],
        dtype=float,
    )
    combined = np.asarray(
        curves["combined_direct_normalized"],
        dtype=float,
    )
    alias_fraction = np.asarray(
        curves["first_alias_fold_psd_fraction"],
        dtype=float,
    )

    fig, (correction_axis, context_axis) = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )
    correction_axis.semilogx(
        frequency,
        correction_db,
        linewidth=2.0,
        label="required residual ASD transfer",
    )
    correction_axis.axhline(0.0, linewidth=1.0)
    correction_axis.set_ylabel("Required correction [dB]")
    correction_axis.set_title(
        "Required residual transfer: measured / best-fit model"
    )
    correction_axis.grid(True, which="both", alpha=0.25)
    correction_axis.legend(fontsize=9)

    context_axis.semilogx(
        frequency,
        20.0 * np.log10(np.maximum(hardware, np.finfo(float).tiny)),
        linewidth=1.4,
        label="modeled 100 kHz hardware Bessel",
    )
    context_axis.semilogx(
        frequency,
        20.0 * np.log10(np.maximum(analysis, np.finfo(float).tiny)),
        linewidth=1.4,
        label="modeled 10 kHz analysis response",
    )
    context_axis.semilogx(
        frequency,
        20.0 * np.log10(np.maximum(combined, np.finfo(float).tiny)),
        linewidth=1.8,
        label="modeled direct-path product",
    )
    context_axis.set_xlabel("Frequency [Hz]")
    context_axis.set_ylabel("Normalized modeled response [dB]")
    context_axis.grid(True, which="both", alpha=0.25)
    context_axis.legend(fontsize=8, loc="lower left")

    alias_axis = context_axis.twinx()
    alias_axis.semilogx(
        frequency,
        alias_fraction,
        linewidth=1.2,
        linestyle="--",
        label="first alias PSD fraction",
    )
    alias_axis.set_ylabel("Alias PSD fraction")
    alias_axis.set_ylim(0.0, 1.0)
    alias_axis.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def johnson_source_scale_diagnostics(
    curves: dict,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> dict:
    """Scale only TES Johnson ASD after transfer, for diagnosis only.

    A scale s multiplies the TES Johnson ASD, so its PSD contribution is
    multiplied by s^2.  Every other source class and every transfer function
    remains exactly at the best-fit value.  The resulting total is normalized
    at 1 kHz exactly like the production shape objective.
    """

    frequency = np.asarray(curves["frequencies_Hz"], dtype=float)
    fit_freq = np.asarray(fit_freq, dtype=float)
    target = np.asarray(target, dtype=float)
    total_asd = np.asarray(curves["total_asd_A_rtHz"], dtype=float)
    johnson_asd = np.asarray(
        curves["class_asd_A_rtHz"]["TES_Johnson"],
        dtype=float,
    )
    if frequency.shape != fit_freq.shape or not np.allclose(
        frequency,
        fit_freq,
        rtol=0.0,
        atol=1e-9,
    ):
        raise ValueError(
            "Johnson-scale diagnostic curves must use the optimizer fit grid"
        )
    if target.shape != fit_freq.shape:
        raise ValueError(
            "Johnson-scale diagnostic target must use the optimizer fit grid"
        )

    total_psd = total_asd**2
    johnson_psd = johnson_asd**2
    other_psd = np.maximum(
        total_psd - johnson_psd,
        0.0,
    )
    tiny = np.finfo(float).tiny
    rows = []

    for scale in JOHNSON_SOURCE_SCALE_DIAGNOSTIC_GRID:
        scaled_johnson_psd = (float(scale) ** 2) * johnson_psd
        scaled_total_psd = np.maximum(
            other_psd + scaled_johnson_psd,
            tiny,
        )
        scaled_asd = np.sqrt(scaled_total_psd)
        normalized = normalize_at(
            frequency,
            scaled_asd,
            reference_hz=ABSOLUTE_ASD_REFERENCE_HZ,
        )
        residual = log_ratio_residual(normalized, target)
        score = float(fit_score(normalized, target, fit_freq, args))

        bands = {}
        for low, high, _weight in FIT_BANDS_HZ:
            low_eff = max(float(low), float(args.fit_min_hz))
            high_eff = min(float(high), float(args.fit_max_hz))
            mask = (fit_freq >= low_eff) & (fit_freq <= high_eff)
            if not np.any(mask):
                continue
            mean_log = float(np.mean(residual[mask]))
            bands[f"{low_eff:g}-{high_eff:g}_Hz"] = {
                "rms_log10_ratio": float(
                    np.sqrt(np.mean(residual[mask] ** 2))
                ),
                "mean_log10_ratio": mean_log,
                "geometric_mean_model_over_measurement": float(
                    10.0 ** mean_log
                ),
                "TES_Johnson_mean_psd_fraction_after_scaling": float(
                    np.mean(
                        scaled_johnson_psd[mask]
                        / scaled_total_psd[mask]
                    )
                ),
            }

        rows.append(
            {
                "TES_Johnson_asd_scale": float(scale),
                "TES_Johnson_psd_scale": float(scale) ** 2,
                "shape_score": score,
                "bands": bands,
            }
        )

    best_global = min(rows, key=lambda row: row["shape_score"])
    active_band_names = list(rows[0]["bands"]) if rows else []
    best_by_band = {}
    for band_name in active_band_names:
        best_row = min(
            rows,
            key=lambda row: abs(
                row["bands"][band_name]["mean_log10_ratio"]
            ),
        )
        best_by_band[band_name] = {
            "TES_Johnson_asd_scale": float(
                best_row["TES_Johnson_asd_scale"]
            ),
            "TES_Johnson_psd_scale": float(
                best_row["TES_Johnson_psd_scale"]
            ),
            "geometric_mean_model_over_measurement": float(
                best_row["bands"][band_name][
                    "geometric_mean_model_over_measurement"
                ]
            ),
            "mean_log10_ratio": float(
                best_row["bands"][band_name]["mean_log10_ratio"]
            ),
        }

    return {
        "diagnostic_only": True,
        "held_fixed_except": "TES_Johnson source amplitude after transfer",
        "scale_semantics": (
            "ASD_J -> s_J * ASD_J; PSD_J -> s_J^2 * PSD_J"
        ),
        "normalization_Hz": float(ABSOLUTE_ASD_REFERENCE_HZ),
        "rows": rows,
        "best_global_shape_score_row": best_global,
        "best_scale_by_band_mean_ratio": best_by_band,
    }


def _electrical_rc_diagnostic_row(
    trial: dict,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> dict:
    """Evaluate one diagnostic RC/load decomposition at fixed other parameters."""

    point = tes_operating_point(trial)
    row = {
        "R_l_ohm": float(trial["R_l"]),
        "R_rc_ohm": float(trial["R_rc"]),
        "f_rc_Hz": float(trial["f_rc_Hz"]),
        "R_l_plus_R_rc_ohm": float(trial["R_l"] + trial["R_rc"]),
        "stable": bool(point.get("stable", False)),
        "valid": bool(point.get("valid", False)),
        "reason": point.get("reason"),
    }
    if not row["stable"] or not row["valid"]:
        return row

    try:
        model, _ = deterministic_simulated_spectrum(
            trial.copy(),
            fit_freq,
        )
    except Exception as error:
        row["stable"] = False
        row["valid"] = False
        row["reason"] = f"evaluation_failed:{error}"
        return row

    diagnostics = band_fit_diagnostics(
        model,
        target,
        fit_freq,
        args,
    )
    bands = {}
    for band_name, values in diagnostics.items():
        mean_log = float(values["mean_log10_ratio"])
        bands[band_name] = {
            "rms_log10_ratio": float(values["rms_log10_ratio"]),
            "mean_log10_ratio": mean_log,
            "geometric_mean_model_over_measurement": float(
                10.0 ** mean_log
            ),
        }

    row.update(
        {
            "shape_score": float(
                fit_score(model, target, fit_freq, args)
            ),
            "bands": bands,
        }
    )
    return row


def electrical_rc_degeneracy_diagnostics(
    candidate: dict,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
    fixed_rl_grid_ohm=None,
    profile_maxfev: int = RC_DEGENERACY_PROFILE_MAXFEV,
) -> dict:
    """Diagnose whether the RC branch is dynamic or merely extra resistance.

    Two complementary shape-only diagnostics are performed while every
    non-electrical best-fit parameter is held fixed:

    1. Keep R_l + R_rc fixed at the best-fit DC value, keep f_rc fixed, and
       repartition the resistance between R_l and R_rc.
    2. Fix R_l to several values and re-optimize only R_rc and f_rc.

    If the score is insensitive to the first repartition, or the second profile
    repeatedly drives f_rc to its upper bound, the fitted RC branch is acting
    predominantly like an additional static series resistance rather than a
    resolved in-band relaxation.
    """

    if candidate.get("electrical_link_model") != ELECTRICAL_LINK_MODEL_RC:
        raise ValueError(
            "RC degeneracy diagnostics require the rl_rc_relaxation model"
        )

    best_rl = float(candidate["R_l"])
    best_rrc = float(candidate["R_rc"])
    best_frc = float(candidate["f_rc_Hz"])
    best_total = best_rl + best_rrc

    lower_partition_rl = max(
        R_L_FIT_MIN_OHM,
        best_total - R_RC_FIT_MAX_OHM,
    )
    upper_partition_rl = min(
        R_L_FIT_MAX_OHM,
        best_total - R_RC_FIT_MIN_OHM,
    )
    partition_rows = []
    if lower_partition_rl <= upper_partition_rl:
        partition_values = list(
            np.linspace(
                lower_partition_rl,
                upper_partition_rl,
                RC_DEGENERACY_PARTITION_POINTS,
            )
        )
        if lower_partition_rl <= best_rl <= upper_partition_rl:
            partition_values.append(best_rl)
        partition_values = sorted(
            set(float(value) for value in partition_values)
        )
        for r_l in partition_values:
            trial = candidate.copy()
            trial["R_l"] = float(r_l)
            trial["R_rc"] = float(best_total - r_l)
            trial["f_rc_Hz"] = best_frc
            partition_rows.append(
                _electrical_rc_diagnostic_row(
                    trial,
                    target,
                    fit_freq,
                    args,
                )
            )

    if fixed_rl_grid_ohm is None:
        fixed_rl_grid_ohm = RC_DEGENERACY_FIXED_RL_GRID_OHM
    profile_values = [
        float(value)
        for value in fixed_rl_grid_ohm
        if R_L_FIT_MIN_OHM <= float(value) <= R_L_FIT_MAX_OHM
    ]
    profile_values.append(best_rl)
    profile_values = sorted(set(profile_values))

    transformed_bounds = [
        (np.log10(R_RC_FIT_MIN_OHM), np.log10(R_RC_FIT_MAX_OHM)),
        (np.log10(F_RC_FIT_MIN_HZ), np.log10(F_RC_FIT_MAX_HZ)),
    ]
    seed = np.asarray(
        [np.log10(best_rrc), np.log10(best_frc)],
        dtype=float,
    )
    profile_rows = []

    for r_l in profile_values:
        def profile_objective(vector):
            trial = candidate.copy()
            trial["R_l"] = float(r_l)
            trial["R_rc"] = float(10.0 ** vector[0])
            trial["f_rc_Hz"] = float(10.0 ** vector[1])
            row = _electrical_rc_diagnostic_row(
                trial,
                target,
                fit_freq,
                args,
            )
            return float(row.get("shape_score", 1.0e12))

        seed_score = float(profile_objective(seed))
        result = minimize(
            profile_objective,
            seed,
            method="Powell",
            bounds=transformed_bounds,
            options={
                "maxfev": int(profile_maxfev),
                "xtol": 1.0e-4,
                "ftol": 1.0e-6,
                "disp": False,
            },
        )
        result_vector = np.asarray(result.x, dtype=float)
        result_score = float(profile_objective(result_vector))
        if seed_score <= result_score:
            chosen = seed
            chosen_from = "best_fit_seed"
        else:
            chosen = result_vector
            chosen_from = "profile_optimization"

        trial = candidate.copy()
        trial["R_l"] = float(r_l)
        trial["R_rc"] = float(10.0 ** chosen[0])
        trial["f_rc_Hz"] = float(10.0 ** chosen[1])
        row = _electrical_rc_diagnostic_row(
            trial,
            target,
            fit_freq,
            args,
        )
        row.update(
            {
                "profile_optimizer_success": bool(result.success),
                "profile_optimizer_status": int(result.status),
                "profile_optimizer_nfev": int(result.nfev),
                "chosen_from": chosen_from,
                "f_rc_at_lower_bound": bool(
                    np.isclose(
                        row["f_rc_Hz"],
                        F_RC_FIT_MIN_HZ,
                        rtol=0.0,
                        atol=F_RC_FIT_MIN_HZ * 1.0e-5,
                    )
                ),
                "f_rc_at_upper_bound": bool(
                    np.isclose(
                        row["f_rc_Hz"],
                        F_RC_FIT_MAX_HZ,
                        rtol=0.0,
                        atol=F_RC_FIT_MAX_HZ * 1.0e-5,
                    )
                ),
                "R_rc_at_lower_bound": bool(
                    np.isclose(
                        row["R_rc_ohm"],
                        R_RC_FIT_MIN_OHM,
                        rtol=0.0,
                        atol=R_RC_FIT_MIN_OHM * 1.0e-5,
                    )
                ),
                "R_rc_at_upper_bound": bool(
                    np.isclose(
                        row["R_rc_ohm"],
                        R_RC_FIT_MAX_OHM,
                        rtol=0.0,
                        atol=R_RC_FIT_MAX_OHM * 1.0e-5,
                    )
                ),
            }
        )
        profile_rows.append(row)

    finite_partition = [
        row for row in partition_rows if "shape_score" in row
    ]
    finite_profile = [
        row for row in profile_rows if "shape_score" in row
    ]
    best_profile = (
        min(finite_profile, key=lambda row: row["shape_score"])
        if finite_profile
        else None
    )
    partition_score_span = (
        float(
            max(row["shape_score"] for row in finite_partition)
            - min(row["shape_score"] for row in finite_partition)
        )
        if finite_partition
        else None
    )
    profile_score_span = (
        float(
            max(row["shape_score"] for row in finite_profile)
            - min(row["shape_score"] for row in finite_profile)
        )
        if finite_profile
        else None
    )
    upper_bound_count = sum(
        bool(row.get("f_rc_at_upper_bound"))
        for row in finite_profile
    )

    return {
        "diagnostic_only": True,
        "objective": "normalized shape score only; absolute ASD excluded",
        "held_fixed": (
            "all best-fit parameters except the explicitly profiled "
            "R_l, R_rc, and f_rc"
        ),
        "tes_resistance_response": "instantaneous alpha/beta (unchanged)",
        "best_fit": {
            "R_l_ohm": best_rl,
            "R_rc_ohm": best_rrc,
            "f_rc_Hz": best_frc,
            "R_l_plus_R_rc_ohm": best_total,
        },
        "constant_dc_sum_partition_sweep": {
            "definition": (
                "hold R_l + R_rc and f_rc at best-fit values; repartition "
                "the same DC resistance between R_l and R_rc"
            ),
            "rows": partition_rows,
            "shape_score_span": partition_score_span,
        },
        "fixed_R_l_profile": {
            "definition": (
                "fix R_l; re-optimize only R_rc and f_rc while all other "
                "best-fit parameters remain fixed"
            ),
            "profile_maxfev_per_R_l": int(profile_maxfev),
            "rows": profile_rows,
            "best_row": best_profile,
            "shape_score_span": profile_score_span,
            "rows_with_f_rc_at_upper_bound": int(upper_bound_count),
            "finite_row_count": int(len(finite_profile)),
        },
        "interpretation_guardrail": (
            "A flat constant-sum sweep or repeated f_rc upper-bound solutions "
            "indicates R_l/R_rc degeneracy or an effectively static added "
            "resistance; it is not evidence for a resolved TES resistance "
            "relaxation."
        ),
    }


def electrical_rc_ablation_diagnostics(
    candidate: dict,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> dict:
    """Compare the fitted RC model with the same parameters using plain R_l+L."""

    rc_model, _ = deterministic_simulated_spectrum(
        candidate.copy(),
        fit_freq,
    )
    plain = candidate.copy()
    plain["electrical_link_model"] = "rl"
    plain_model, _ = deterministic_simulated_spectrum(
        plain,
        fit_freq,
    )
    rc_score = float(fit_score(rc_model, target, fit_freq, args))
    plain_score = float(fit_score(plain_model, target, fit_freq, args))
    return {
        "diagnostic_only": True,
        "tes_resistance_response": "instantaneous alpha/beta in both models",
        "fitted_model": ELECTRICAL_LINK_MODEL_RC,
        "comparison_model": "rl",
        "R_rc_ohm": float(candidate["R_rc"]),
        "f_rc_Hz": float(candidate["f_rc_Hz"]),
        "rc_shape_score": rc_score,
        "same_parameters_without_rc_shape_score": plain_score,
        "score_improvement_from_rc_at_same_parameters": float(
            plain_score - rc_score
        ),
        "rc_bands": band_fit_diagnostics(
            rc_model,
            target,
            fit_freq,
            args,
        ),
        "same_parameters_without_rc_bands": band_fit_diagnostics(
            plain_model,
            target,
            fit_freq,
            args,
        ),
    }


def johnson_beta_audit(candidate: dict) -> dict:
    """Record the exact beta dependence used by the shared/production model."""

    beta = float(candidate["beta"])
    excess_m = float(candidate.get("excess_johnson_M", 0.0))
    temperature = float(candidate["T_c"])
    resistance = float(candidate["R"])
    voltage_asd = tes_johnson_voltage_asd(
        temperature,
        resistance,
        beta,
        excess_m,
    )
    point = tes_operating_point(candidate)
    tau_el = float(point["tau_el_s"])
    return {
        "status": "shared_production_convention",
        "source_voltage_psd_expression": (
            "4*k_B*T_c*R*(1+2*beta)*(1+excess_johnson_M^2)"
        ),
        "source_voltage_asd_expression": (
            "sqrt(4*k_B*T_c*R*(1+2*beta)*(1+excess_johnson_M^2))"
        ),
        "electrical_tau_denominator": "R_l + R*(1+beta)",
        "joule_coupling_factor": "2+beta",
        "beta": beta,
        "excess_johnson_M": excess_m,
        "source_psd_beta_factor": float(1.0 + 2.0 * beta),
        "source_asd_factor_vs_beta0_same_M": float(
            np.sqrt(1.0 + 2.0 * beta)
        ),
        "source_voltage_asd_V_rtHz": float(voltage_asd),
        "tau_el_s": tau_el,
        "electrical_corner_Hz": float(
            1.0 / (2.0 * np.pi * tau_el)
        ),
        "electrical_link_model": candidate.get(
            "electrical_link_model",
            "rl",
        ),
        "R_rc_ohm": (
            float(candidate["R_rc"]) if "R_rc" in candidate else None
        ),
        "f_rc_Hz": (
            float(candidate["f_rc_Hz"]) if "f_rc_Hz" in candidate else None
        ),
        "tau_rc_s": point.get("tau_rc_s"),
        "C_rc_equivalent_F": point.get("C_rc_equivalent_F"),
        "tes_resistance_response": "instantaneous alpha/beta",
        "note": (
            "This audit records the repository convention and code-path parity; "
            "it does not by itself establish the experimental beta convention."
        ),
    }


def beta_sweep_diagnostics(
    candidate: dict,
    target: np.ndarray,
    fit_freq: np.ndarray,
    args,
) -> dict:
    """Sweep beta only, holding every other best-fit parameter fixed."""

    beta_values = sorted(
        set(BETA_DIAGNOSTIC_GRID + (float(candidate["beta"]),))
    )
    rows = []
    for beta in beta_values:
        trial = candidate.copy()
        trial["beta"] = float(beta)
        point = tes_operating_point(trial)
        row = {
            "beta": float(beta),
            "stable": bool(point.get("stable", False)),
            "valid": bool(point.get("valid", False)),
            "reason": point.get("reason"),
        }
        if not row["stable"] or not row["valid"]:
            rows.append(row)
            continue

        model, _reference_asd = deterministic_simulated_spectrum(
            trial.copy(),
            fit_freq,
        )
        score = float(fit_score(model, target, fit_freq, args))
        band_diag = band_fit_diagnostics(model, target, fit_freq, args)
        source_summary, _curves = source_class_diagnostics(
            trial.copy(),
            fit_freq,
            args,
        )

        band_rows = {}
        for band_name, diagnostics in band_diag.items():
            mean_log = float(diagnostics["mean_log10_ratio"])
            source_band = source_summary["fit_bands"][band_name]["sources"]
            band_rows[band_name] = {
                "rms_log10_ratio": float(
                    diagnostics["rms_log10_ratio"]
                ),
                "mean_log10_ratio": mean_log,
                "geometric_mean_model_over_measurement": float(
                    10.0 ** mean_log
                ),
                "TES_Johnson_mean_psd_fraction": float(
                    source_band["TES_Johnson"]["mean_psd_fraction"]
                ),
            }

        audit = johnson_beta_audit(trial)
        row.update(
            {
                "shape_score": score,
                "source_psd_beta_factor": audit[
                    "source_psd_beta_factor"
                ],
                "source_asd_factor_vs_beta0_same_M": audit[
                    "source_asd_factor_vs_beta0_same_M"
                ],
                "TES_Johnson_voltage_asd_V_rtHz": audit[
                    "source_voltage_asd_V_rtHz"
                ],
                "tau_el_s": audit["tau_el_s"],
                "electrical_corner_Hz": audit["electrical_corner_Hz"],
                "bands": band_rows,
            }
        )
        rows.append(row)

    stable_rows = [
        row for row in rows
        if row.get("stable") and "shape_score" in row
    ]
    best_row = min(
        stable_rows,
        key=lambda row: row["shape_score"],
    ) if stable_rows else None

    return {
        "diagnostic_only": True,
        "held_fixed_except": "beta",
        "grid_includes_best_fit_beta": float(candidate["beta"]),
        "rows": rows,
        "best_fixed_other_parameters_row": best_row,
    }


def eigenmode_diagnostics(candidate: dict) -> dict:
    """Return eigenfrequencies and scaled state participation for each mode."""

    matrix = tes_linearized_matrix(candidate, 0.0)
    eigenvalues, eigenvectors = np.linalg.eig(-matrix)
    point = tes_operating_point(candidate)
    current_scale = float(point["current_A"])
    temperature_scale = float(candidate["T_c"])

    if matrix.shape[0] == 9:
        voltage_scale = current_scale * max(
            float(candidate.get("R_rc", 0.0)),
            float(candidate.get("R_l", 0.0)),
            1.0e-12,
        )
        state_names = (
            "I1", "Vrc1", "TES1", "Stycast1", "Pb_center",
            "Stycast2", "TES2", "Vrc2", "I2",
        )
        state_scales = np.asarray(
            [
                current_scale, voltage_scale, temperature_scale,
                temperature_scale, temperature_scale, temperature_scale,
                temperature_scale, voltage_scale, current_scale,
            ],
            dtype=float,
        )
    elif matrix.shape[0] == 7 and candidate.get(
        "electrical_link_model"
    ) == ELECTRICAL_LINK_MODEL_RC:
        voltage_scale = current_scale * max(
            float(candidate.get("R_rc", 0.0)),
            float(candidate.get("R_l", 0.0)),
            1.0e-12,
        )
        state_names = (
            "I1", "Vrc1", "TES1", "Pb_center", "TES2", "Vrc2", "I2",
        )
        state_scales = np.asarray(
            [
                current_scale, voltage_scale, temperature_scale,
                temperature_scale, temperature_scale, voltage_scale,
                current_scale,
            ],
            dtype=float,
        )
    elif matrix.shape[0] == 7:
        state_names = ("I1", "TES1", "Stycast1", "Pb_center", "Stycast2", "TES2", "I2")
        state_scales = np.asarray(
            [current_scale, temperature_scale, temperature_scale, temperature_scale,
             temperature_scale, temperature_scale, current_scale],
            dtype=float,
        )
    elif matrix.shape[0] == 5:
        state_names = ("I1", "TES1", "Pb_center", "TES2", "I2")
        state_scales = np.asarray(
            [current_scale, temperature_scale, temperature_scale,
             temperature_scale, current_scale],
            dtype=float,
        )
    else:
        raise ValueError(f"Unsupported TES state count: {matrix.shape[0]}")

    rows = []
    for index, value in enumerate(eigenvalues):
        real = float(np.real(value))
        imag = float(np.imag(value))
        magnitude = float(np.abs(value))
        decay_rate = max(-real, 0.0)
        scaled = eigenvectors[:, index] / state_scales
        raw = np.abs(scaled) ** 2
        norm = float(np.sum(raw))
        participation = raw / norm if np.isfinite(norm) and norm > 0.0 else np.zeros_like(raw)
        state_participation = {
            name: float(part) for name, part in zip(state_names, participation)
        }
        dominant_state = max(state_participation, key=state_participation.get)
        rows.append(
            {
                "real_s_inv": real,
                "imag_s_inv": imag,
                "time_constant_s": float(1.0 / decay_rate) if decay_rate > 0.0 else None,
                "decay_corner_Hz": float(decay_rate / (2.0 * np.pi)),
                "oscillation_Hz": float(abs(imag) / (2.0 * np.pi)),
                "natural_frequency_Hz": float(magnitude / (2.0 * np.pi)),
                "damping_ratio": float(decay_rate / magnitude) if magnitude > 0.0 else None,
                "dominant_state": dominant_state,
                "dominant_state_participation": float(state_participation[dominant_state]),
                "state_participation": state_participation,
            }
        )
    rows.sort(key=lambda row: row["natural_frequency_Hz"])
    return {
        "state_count": int(matrix.shape[0]),
        "state_order": list(state_names),
        "participation_scaling": {
            "current_states": "delta_I / operating_current",
            "thermal_states": "delta_T / T_c",
            "RC_voltage_states": "delta_Vrc / (operating_current * max(R_rc, R_l))",
            "normalization": "sum(abs(scaled_eigenvector)**2) = 1 per mode",
        },
        "stable": bool(all(row["real_s_inv"] < 0.0 for row in rows)),
        "modes_sorted_by_natural_frequency": rows,
    }


def plot_source_class_diagnostics(curves: dict, output_path: Path) -> Path:
    """Plot best-fit source-class ASD and their fractional PSD contributions."""

    frequency = np.asarray(curves["frequencies_Hz"], dtype=float)
    class_asd = curves["class_asd_A_rtHz"]
    total_asd = np.asarray(curves["total_asd_A_rtHz"], dtype=float)
    total_psd = np.maximum(total_asd**2, np.finfo(float).tiny)

    fig, (asd_axis, fraction_axis) = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )
    asd_axis.loglog(
        frequency,
        total_asd,
        linewidth=2.4,
        label="total model",
    )
    for name, asd in class_asd.items():
        asd = np.asarray(asd, dtype=float)
        asd_axis.loglog(frequency, asd, linewidth=1.4, label=name)
        fraction_axis.semilogx(
            frequency,
            asd**2 / total_psd,
            linewidth=1.4,
            label=name,
        )

    asd_axis.set_ylabel("Model ASD [A/rtHz]")
    asd_axis.set_title("Best-fit CH0 source-class decomposition")
    asd_axis.grid(True, which="both", alpha=0.25)
    asd_axis.legend(fontsize=8, ncol=2)
    fraction_axis.set_xlabel("Frequency [Hz]")
    fraction_axis.set_ylabel("PSD fraction")
    fraction_axis.set_ylim(0.0, 1.05)
    fraction_axis.grid(True, which="both", alpha=0.25)
    fraction_axis.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def g_tes_bath_from_joule_power(
    joule_power_w: float,
    t_c: float,
    t_bath: float,
    exponent: float,
) -> float:
    """Convert measured IV Joule power into differential G at T_c."""

    if not (t_c > t_bath > 0.0):
        raise ValueError("Require T_c > T_bath > 0 for IV power balance")
    if exponent <= 0.0 or joule_power_w <= 0.0:
        raise ValueError("n and target Joule power must be positive")
    factor = 1.0 - (t_bath / t_c) ** exponent
    if factor <= 0.0:
        raise ValueError("Invalid thermal power-law factor")
    return float(joule_power_w * exponent / (t_c * factor))


def decode(
    vector,
    original: dict,
    keys,
    bounds,
    sample_rate: float,
    post_filter_white_asd: float,
    electrical_link_model: str = "rl",
):
    candidate = original.copy()
    for key, value in zip(keys, vector):
        candidate[key] = float(10.0 ** value if bounds[key].logarithmic else value)
    # Match the experimental modelnoise.txt analysis filter; this is not fitted.
    candidate["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    # Filter coefficients and the frequency grid must use the acquisition rate
    # of modelnoise.txt, even when an old optimization reference used 300 kHz.
    candidate["rate"] = sample_rate

    target_joule_power_w = candidate.pop("_target_joule_power_W", None)
    if target_joule_power_w is not None:
        candidate["G_tes-bath"] = g_tes_bath_from_joule_power(
            float(target_joule_power_w),
            float(candidate["T_c"]),
            float(candidate["T_bath"]),
            float(candidate["n"]),
        )

    candidate["post_filter_white_asd_A_rtHz"] = post_filter_white_asd
    candidate.pop("readout_white_asd_A_rtHz", None)
    candidate["hardware_bessel_order"] = TARGET_HARDWARE_BESSEL_ORDER
    candidate["hardware_bessel_norm"] = TARGET_HARDWARE_BESSEL_NORM
    candidate["hardware_bessel_cutoff_Hz"] = TARGET_HARDWARE_BESSEL_CUTOFF_HZ
    candidate["thermal_link_model"] = "stycast_node"
    candidate["electrical_link_model"] = electrical_link_model
    if electrical_link_model == "rl":
        candidate["R_l"] = float(candidate["R_series_eff"])
        candidate.pop("R_rc", None)
        candidate.pop("f_rc_Hz", None)
    return candidate


def apply_post_filter_white_fraction(candidate: dict) -> float:
    """Convert the dimensionless fitted readout floor into absolute ASD."""

    fraction = float(
        candidate.get(
            "post_filter_white_fraction",
            POST_FILTER_WHITE_FRACTION_INITIAL,
        )
    )
    if fraction < 0.0:
        raise ValueError("post_filter_white_fraction must be non-negative")

    rate = float(candidate["rate"])
    cutoff_hz = TARGET_HARDWARE_BESSEL_CUTOFF_HZ
    candidate["hardware_bessel_cutoff_Hz"] = cutoff_hz
    detector_at_reference = hardware_sampled_asd(
        candidate,
        np.asarray([1_000.0]),
        rate_hz=rate,
        cutoff_hz=cutoff_hz,
        order=TARGET_HARDWARE_BESSEL_ORDER,
        norm=TARGET_HARDWARE_BESSEL_NORM,
    )[0]
    white_asd = fraction * float(detector_at_reference)
    candidate["post_filter_white_asd_A_rtHz"] = white_asd
    candidate.pop("readout_white_asd_A_rtHz", None)
    return white_asd


def deterministic_simulated_spectrum(
    candidate: dict,
    fit_freq: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Return normalized shape and absolute 1-kHz ASD in A/rtHz."""

    frequency = np.unique(
        np.concatenate(
            (
                np.asarray([ABSOLUTE_ASD_REFERENCE_HZ]),
                np.asarray(fit_freq, dtype=float),
            )
        )
    )
    rate = float(candidate["rate"])
    cutoff_hz = TARGET_HARDWARE_BESSEL_CUTOFF_HZ
    candidate["hardware_bessel_cutoff_Hz"] = cutoff_hz
    apply_post_filter_white_fraction(candidate)
    pre_analysis = hardware_sampled_asd(
        candidate,
        frequency,
        rate_hz=rate,
        cutoff_hz=cutoff_hz,
        order=TARGET_HARDWARE_BESSEL_ORDER,
        norm=TARGET_HARDWARE_BESSEL_NORM,
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
    normalized = normalize_at(
        frequency,
        expected,
        reference_hz=ABSOLUTE_ASD_REFERENCE_HZ,
    )
    reference_asd_A = float(
        np.interp(ABSOLUTE_ASD_REFERENCE_HZ, frequency, expected)
    )
    return np.interp(fit_freq, frequency, normalized), reference_asd_A


def optimize_case(
    args,
    original: dict,
    reference: dict,
    envelope: dict,
    fit_freq: np.ndarray,
    target: np.ndarray,
    target_reference_asd_A: float,
    experimental_rate: float,
    experimental_samples: int,
    post_filter_white_asd: float,
    operating_point: dict,
    work_dir: Path,
):
    """Optimize one fixed-R case derived from one R_SH branch."""

    fixed_r_ohm = float(operating_point["R_TES_ohm"])
    shunt_resistance_ohm = float(operating_point["R_SH_ohm"])
    use_rc_relaxation = bool(args.use_rc_relaxation)
    electrical_link_model = (
        ELECTRICAL_LINK_MODEL_RC if use_rc_relaxation else "rl"
    )

    # Start from the frozen 215 mK target-case proxy, not the generic
    # PoST_Simulations/input.json (which belongs to a different thermal point).
    case_original = reference.copy()
    case_original["R"] = fixed_r_ohm
    case_original["R_SH"] = shunt_resistance_ohm
    case_original["T_bath"] = float(envelope["parameters"]["T_bath"]["nominal"])
    case_original["_target_joule_power_W"] = float(
        operating_point["P_J_W"]
    )
    case_original["rate"] = float(experimental_rate)
    case_original["samples"] = int(experimental_samples)
    case_original["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    case_original["hardware_bessel_order"] = TARGET_HARDWARE_BESSEL_ORDER
    case_original["hardware_bessel_norm"] = TARGET_HARDWARE_BESSEL_NORM
    case_original["hardware_bessel_cutoff_Hz"] = TARGET_HARDWARE_BESSEL_CUTOFF_HZ
    case_original["thermal_link_model"] = "stycast_node"
    case_original["electrical_link_model"] = electrical_link_model
    if use_rc_relaxation:
        case_original.setdefault("R_rc", R_RC_INITIAL_OHM)
        case_original.setdefault("f_rc_Hz", F_RC_INITIAL_HZ)
        case_original.pop("R_series_eff", None)
    else:
        case_original["R_series_eff"] = float(
            case_original.get("R_series_eff", case_original["R_l"])
        )
        case_original["R_l"] = float(case_original["R_series_eff"])
        case_original.pop("R_rc", None)
        case_original.pop("f_rc_Hz", None)
    case_original.setdefault("C_stycast", C_STYCAST_PAD_MATERIAL_J_PER_K)
    case_original.setdefault("G_tes-stycast", G_ABS_TES_MATERIAL_W_PER_K)
    case_original.setdefault("G_stycast-abs", G_ABS_TES_MATERIAL_W_PER_K)
    case_original.setdefault(
        "post_filter_white_fraction",
        POST_FILTER_WHITE_FRACTION_INITIAL,
    )

    work_dir.mkdir(parents=True, exist_ok=True)
    work_input_path = work_dir / "input.json"
    work_noise_path = work_dir / NOISE_DAT_PATH.name
    bounds = parameter_bounds(
        reference,
        envelope,
        fixed_r_ohm=fixed_r_ohm,
        use_rc_relaxation=use_rc_relaxation,
    )
    keys, scipy_bounds = vector_bounds(bounds)

    initial = case_original.copy()
    initial["post_filter_white_asd_A_rtHz"] = post_filter_white_asd
    initial["post_filter_white_fraction"] = max(
        float(initial.get("post_filter_white_fraction", POST_FILTER_WHITE_FRACTION_INITIAL)),
        POST_FILTER_WHITE_FRACTION_MIN,
    )
    initial["hardware_bessel_cutoff_Hz"] = TARGET_HARDWARE_BESSEL_CUTOFF_HZ
    initial.pop("readout_white_asd_A_rtHz", None)

    # Reuse the last good target-case solution as a warm start when available.
    # The hardware cutoff itself is deliberately not inherited: it is forced
    # back to the known physical 100 kHz value before re-optimization.
    if work_input_path.is_file():
        try:
            previous = load_json(work_input_path)
            previous_tbath = float(previous.get("T_bath", np.nan))
            if np.isfinite(previous_tbath) and abs(
                previous_tbath - case_original["T_bath"]
            ) < 5e-3:
                reused = []
                for key in keys:
                    if key == "R_series_eff" and key not in previous:
                        if "R_rc" in previous and "R_l" in previous:
                            value = float(previous["R_l"]) + float(previous["R_rc"])
                        elif "R_l" in previous:
                            value = float(previous["R_l"])
                        else:
                            continue
                    elif key in previous:
                        value = float(previous[key])
                    else:
                        continue
                    bound = bounds[key]
                    if (
                        np.isfinite(value)
                        and value >= bound.lower
                        and value <= bound.upper
                    ):
                        initial[key] = value
                        reused.append(key)
                if reused:
                    print(
                        "Warm-starting from previous target-case solution:",
                        ", ".join(reused),
                    )
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            pass

    initial["hardware_bessel_cutoff_Hz"] = TARGET_HARDWARE_BESSEL_CUTOFF_HZ
    clip_candidate_to_bounds(initial, keys, bounds)
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
            electrical_link_model=electrical_link_model,
        )
        candidate["R"] = fixed_r_ohm
        candidate["R_SH"] = shunt_resistance_ohm
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
            model, model_reference_asd_A = deterministic_simulated_spectrum(
                candidate,
                fit_freq,
            )
            score = combined_fit_score(
                model,
                target,
                fit_freq,
                model_reference_asd_A,
                target_reference_asd_A,
                args,
            )
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
        f"Target assumptions: T_bath search="
        f"{bounds['T_bath'].lower:.6g}--{bounds['T_bath'].upper:.6g} K, "
        f"T_c search={bounds['T_c'].lower:.6g}--{bounds['T_c'].upper:.6g} K, "
        f"P_J={case_original['_target_joule_power_W']:.6g} W"
    )
    print(
        "Physical fit bounds:",
        {
            "T_c_K": [bounds["T_c"].lower, bounds["T_c"].upper],
            "T_bath_K": [bounds["T_bath"].lower, bounds["T_bath"].upper],
            "alpha": [bounds["alpha"].lower, bounds["alpha"].upper],
            "L_H": [L_FIT_MIN_H, L_FIT_MAX_H],
            "electrical_model": (
                "passive_rc_relaxation"
                if use_rc_relaxation
                else "static_effective_series_resistance"
            ),
            "series_resistance_search_ohm": (
                [R_L_FIT_MIN_OHM, R_L_FIT_MAX_OHM]
                if use_rc_relaxation
                else [R_SERIES_EFF_FIT_MIN_OHM, R_SERIES_EFF_FIT_MAX_OHM]
            ),
            "R_rc_ohm": (
                [R_RC_FIT_MIN_OHM, R_RC_FIT_MAX_OHM]
                if use_rc_relaxation
                else None
            ),
            "f_rc_Hz": (
                [F_RC_FIT_MIN_HZ, F_RC_FIT_MAX_HZ]
                if use_rc_relaxation
                else None
            ),
            "tes_resistance_response": "instantaneous alpha/beta (unchanged)",
            "n": [N_FIT_MIN, N_FIT_MAX],
            "C_tes_J_per_K": [C_TES_FIT_MIN_J_PER_K, C_TES_FIT_MAX_J_PER_K],
            "C_stycast_J_per_K": [
                C_STYCAST_FIT_MIN_J_PER_K,
                C_STYCAST_FIT_MAX_J_PER_K,
            ],
            "C_stycast_pad_reference_J_per_K": C_STYCAST_PAD_MATERIAL_J_PER_K,
            "C_abs_J_per_K": [C_ABS_FIT_MIN_J_PER_K, C_ABS_FIT_MAX_J_PER_K],
            "G_tes_stycast_W_per_K": [
                G_TES_STYCAST_FIT_MIN_W_PER_K,
                G_TES_STYCAST_FIT_MAX_W_PER_K,
            ],
            "G_stycast_abs_W_per_K": [
                G_STYCAST_ABS_FIT_MIN_W_PER_K,
                G_STYCAST_ABS_FIT_MAX_W_PER_K,
            ],
            "G_abs_abs_W_per_K": [
                G_ABS_ABS_FIT_MIN_W_PER_K,
                G_ABS_ABS_FIT_MAX_W_PER_K,
            ],
            "G_tes_bath": "derived from IV Joule power",
        },
    )
    print("Initial deterministic evaluation")
    objective(initial_x)

    if args.skip_de:
        starting_x = initial_x
    else:
        print("\nStage 1/3: differential_evolution")
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

    print("\nStage 2/3: Powell")
    local_result = minimize(
        objective,
        starting_x,
        method="Powell",
        bounds=scipy_bounds,
        options={"maxfev": args.powell_maxfev, "disp": True},
    )
    objective(local_result.x)

    print("\nStage 3/3: robust residual least_squares")
    ls_cache = {}
    least_squares_residual_size = 1 + len(fit_freq) + len(
        band_mean_residuals(
            np.zeros(len(fit_freq), dtype=float),
            fit_freq,
            args,
        )
    )
    lower = np.asarray([bound[0] for bound in scipy_bounds], dtype=float)
    upper = np.asarray([bound[1] for bound in scipy_bounds], dtype=float)

    def least_squares_residual(vector):
        nonlocal evaluation_count, stability_rejection_count
        nonlocal simulation_failure_count, best_score, best_candidate

        cache_key = tuple(np.round(vector, 12))
        cached = ls_cache.get(cache_key)
        if cached is not None:
            return cached

        candidate = decode(
            vector,
            case_original,
            keys,
            bounds,
            experimental_rate,
            post_filter_white_asd,
            electrical_link_model=electrical_link_model,
        )
        candidate["R"] = fixed_r_ohm
        candidate["R_SH"] = shunt_resistance_ohm
        candidate["samples"] = int(experimental_samples)
        evaluation_count += 1

        point = tes_operating_point(candidate)
        if not point.get("stable", False):
            stability_rejection_count += 1
            residual_vector = np.full(
                least_squares_residual_size,
                3.0,
                dtype=float,
            )
            ls_cache[cache_key] = residual_vector
            return residual_vector

        try:
            model, model_reference_asd_A = deterministic_simulated_spectrum(
                candidate,
                fit_freq,
            )
            residual_vector = combined_residual_vector(
                model,
                target,
                fit_freq,
                model_reference_asd_A,
                target_reference_asd_A,
                args,
            )
            score = combined_fit_score(
                model,
                target,
                fit_freq,
                model_reference_asd_A,
                target_reference_asd_A,
                args,
            )
            if score < best_score:
                best_score = score
                best_candidate = candidate.copy()
                print(
                    f"  least_squares new best: {best_score:.6g} "
                    f"(evaluation {evaluation_count})"
                )
        except Exception:
            simulation_failure_count += 1
            residual_vector = np.full(
                least_squares_residual_size,
                3.0,
                dtype=float,
            )

        ls_cache[cache_key] = residual_vector
        return residual_vector

    least_squares_start = encode(best_candidate, keys, bounds)
    least_squares_result = least_squares(
        least_squares_residual,
        least_squares_start,
        bounds=(lower, upper),
        method="trf",
        loss="soft_l1",
        f_scale=max(float(args.robust_delta_dex), 0.05),
        x_scale="jac",
        max_nfev=int(args.least_squares_max_nfev),
        verbose=1,
    )
    objective(least_squares_result.x)

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
    finite_model, finite_reference_asd_A = simulated_spectrum(
        validation_candidate,
        fit_freq,
        work_noise_path,
    )
    finite_validation_score = combined_fit_score(
        finite_model,
        target,
        fit_freq,
        finite_reference_asd_A,
        target_reference_asd_A,
        args,
    )
    (
        deterministic_best_model,
        deterministic_reference_asd_A,
    ) = deterministic_simulated_spectrum(
        best_candidate.copy(),
        fit_freq,
    )
    deterministic_band_diagnostics = band_fit_diagnostics(
        deterministic_best_model,
        target,
        fit_freq,
        args,
    )
    finite_band_diagnostics = band_fit_diagnostics(
        finite_model,
        target,
        fit_freq,
        args,
    )
    absolute_asd_diagnostics = {
        "reference_Hz": float(ABSOLUTE_ASD_REFERENCE_HZ),
        "weight": float(args.absolute_asd_weight),
        "measurement_A_rtHz": float(target_reference_asd_A),
        "measurement_pA_rtHz": float(
            target_reference_asd_A / MEASURED_ASD_PA_TO_A
        ),
        "deterministic_model_A_rtHz": float(deterministic_reference_asd_A),
        "deterministic_model_pA_rtHz": float(
            deterministic_reference_asd_A / MEASURED_ASD_PA_TO_A
        ),
        "deterministic_model_over_measurement": float(
            deterministic_reference_asd_A / target_reference_asd_A
        ),
        "finite_model_A_rtHz": float(finite_reference_asd_A),
        "finite_model_pA_rtHz": float(
            finite_reference_asd_A / MEASURED_ASD_PA_TO_A
        ),
        "finite_model_over_measurement": float(
            finite_reference_asd_A / target_reference_asd_A
        ),
    }
    boundary_diagnostics = parameter_boundary_diagnostics(
        best_candidate,
        keys,
        bounds,
    )
    source_diagnostics, source_curves = source_class_diagnostics(
        best_candidate.copy(),
        fit_freq,
        args,
    )
    required_transfer_summary, required_transfer_curves = (
        required_transfer_diagnostics(
            best_candidate.copy(),
            deterministic_best_model,
            target,
            fit_freq,
            args,
        )
    )
    source_ablation_summary = source_ablation_diagnostics(
        source_curves,
        target,
        fit_freq,
        args,
    )
    johnson_scale_summary = johnson_source_scale_diagnostics(
        source_curves,
        target,
        fit_freq,
        args,
    )
    electrical_rc_ablation_summary = electrical_rc_ablation_diagnostics(
        best_candidate.copy(),
        target,
        fit_freq,
        args,
    )
    electrical_rc_degeneracy_summary = electrical_rc_degeneracy_diagnostics(
        best_candidate.copy(),
        target,
        fit_freq,
        args,
    )
    johnson_audit_summary = johnson_beta_audit(best_candidate.copy())
    beta_sweep_summary = beta_sweep_diagnostics(
        best_candidate.copy(),
        target,
        fit_freq,
        args,
    )
    eigenmode_summary = eigenmode_diagnostics(best_candidate.copy())
    source_plot_path = plot_source_class_diagnostics(
        source_curves,
        work_dir / "source_contributions.png",
    )
    required_transfer_plot_path = plot_required_transfer_diagnostics(
        required_transfer_curves,
        work_dir / "required_transfer.png",
    )

    print("\nObjective evaluations:", evaluation_count)
    print("Stability rejections:", stability_rejection_count)
    print("Model-evaluation failures:", simulation_failure_count)
    print("Best deterministic score:", best_score)
    print("Full-record validation score:", finite_validation_score)
    print("Deterministic band diagnostics:")
    print(json.dumps(deterministic_band_diagnostics, indent=2))
    print("Full-record band diagnostics:")
    print(json.dumps(finite_band_diagnostics, indent=2))
    print("Absolute ASD diagnostics:")
    print(json.dumps(absolute_asd_diagnostics, indent=2))
    print("Parameter boundary diagnostics:")
    print(json.dumps(boundary_diagnostics, indent=2))
    print("Required-transfer diagnostics:")
    print(json.dumps(required_transfer_summary, indent=2))
    print("Source-class diagnostics:")
    print(json.dumps(source_diagnostics, indent=2))
    print("Source-ablation diagnostics:")
    print(json.dumps(source_ablation_summary, indent=2))
    print("Johnson source-scale diagnostics:")
    print(json.dumps(johnson_scale_summary, indent=2))
    print("Electrical RC ablation diagnostics:")
    print(json.dumps(electrical_rc_ablation_summary, indent=2))
    print("Electrical RC degeneracy diagnostics:")
    print(json.dumps(electrical_rc_degeneracy_summary, indent=2))
    print("Johnson/beta audit:")
    print(json.dumps(johnson_audit_summary, indent=2))
    print("Beta sweep diagnostics:")
    print(json.dumps(beta_sweep_summary, indent=2))
    print("Eigenmode diagnostics:")
    print(json.dumps(eigenmode_summary, indent=2))
    print("Best fitted parameters:")
    print(json.dumps({key: best_candidate[key] for key in keys}, indent=2))
    print(
        "Best absolute post-filter white ASD:",
        best_candidate.get("post_filter_white_asd_A_rtHz", 0.0),
        "A/rtHz",
    )

    return {
        "R_SH_ohm": shunt_resistance_ohm,
        "R_TES_ohm": fixed_r_ohm,
        "iv_operating_point": operating_point,
        "best_score": float(best_score),
        "finite_validation_score": float(finite_validation_score),
        "deterministic_band_diagnostics": deterministic_band_diagnostics,
        "finite_band_diagnostics": finite_band_diagnostics,
        "absolute_asd_diagnostics": absolute_asd_diagnostics,
        "parameter_boundary_diagnostics": boundary_diagnostics,
        "required_transfer_diagnostics": required_transfer_summary,
        "required_transfer_plot": str(required_transfer_plot_path),
        "source_class_diagnostics": source_diagnostics,
        "source_ablation_diagnostics": source_ablation_summary,
        "johnson_source_scale_diagnostics": johnson_scale_summary,
        "electrical_rc_ablation_diagnostics": electrical_rc_ablation_summary,
        "electrical_rc_degeneracy_diagnostics": electrical_rc_degeneracy_summary,
        "johnson_beta_audit": johnson_audit_summary,
        "beta_sweep_diagnostics": beta_sweep_summary,
        "eigenmode_diagnostics": eigenmode_summary,
        "source_contribution_plot": str(source_plot_path),
        "least_squares": {
            "success": bool(least_squares_result.success),
            "status": int(least_squares_result.status),
            "cost": float(least_squares_result.cost),
            "optimality": float(least_squares_result.optimality),
            "nfev": int(least_squares_result.nfev),
        },
        "evaluations": evaluation_count,
        "stability_rejections": stability_rejection_count,
        "simulation_failures": simulation_failure_count,
        "best_hardware_bessel_cutoff_Hz": float(
            best_candidate["hardware_bessel_cutoff_Hz"]
        ),
        "best_post_filter_white_fraction": float(
            best_candidate["post_filter_white_fraction"]
        ),
        "best_post_filter_white_asd_A_rtHz": float(
            best_candidate.get("post_filter_white_asd_A_rtHz", 0.0)
        ),
        "target_joule_power_W": float(operating_point["P_J_W"]),
        "derived_G_tes_bath_W_per_K": float(best_candidate["G_tes-bath"]),
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
    if args.absolute_asd_weight < 0.0:
        raise ValueError("--absolute-asd-weight must be non-negative")

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
    reference["hardware_bessel_norm"] = TARGET_HARDWARE_BESSEL_NORM
    reference["hardware_bessel_cutoff_Hz"] = TARGET_HARDWARE_BESSEL_CUTOFF_HZ
    reference.setdefault(
        "post_filter_white_fraction",
        POST_FILTER_WHITE_FRACTION_INITIAL,
    )

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
        "fixed acquisition and calibration:",
        {
            "T_bath_nominal_K": reference["T_bath"],
            "rate_Hz": experimental_rate,
            "samples": experimental_samples,
            "hardware_bessel_order": TARGET_HARDWARE_BESSEL_ORDER,
            "hardware_bessel_norm": TARGET_HARDWARE_BESSEL_NORM,
            "hardware_bessel_cutoff_Hz": TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
            "post_filter_white_fraction_search": [
                POST_FILTER_WHITE_FRACTION_MIN,
                POST_FILTER_WHITE_FRACTION_MAX,
            ],
            "analysis_bessel_cutoff_Hz": SIM_ANALYSIS_CUTOFF_HZ,
        },
    )
    fit_freq, target, target_reference_asd_A = target_spectrum(args)

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
                target_reference_asd_A,
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
            "T_bath_is_fitted": True,
            "T_bath_nominal_K": float(reference["T_bath"]),
            "T_bath_search_K": [
                float(reference["T_bath"] - T_BATH_FIT_HALF_WIDTH_K),
                float(reference["T_bath"] + T_BATH_FIT_HALF_WIDTH_K),
            ],
            "fixed_target_assumptions": {
                "rate_Hz": experimental_rate,
                "samples": experimental_samples,
                "hardware_bessel_order": TARGET_HARDWARE_BESSEL_ORDER,
                "hardware_bessel_norm": TARGET_HARDWARE_BESSEL_NORM,
                "hardware_bessel_cutoff_Hz": float(
                    TARGET_HARDWARE_BESSEL_CUTOFF_HZ
                ),
                "post_filter_white_fraction_search": [
                    float(POST_FILTER_WHITE_FRACTION_MIN),
                    float(POST_FILTER_WHITE_FRACTION_MAX),
                ],
                "analysis_bessel_cutoff_Hz": float(SIM_ANALYSIS_CUTOFF_HZ),
                "electrical_link_model": ELECTRICAL_LINK_MODEL_RC,
                "tes_resistance_response": "instantaneous alpha/beta (unchanged)",
            },
            "T_c_fit_range_K": [
                float(envelope["parameters"]["T_c"]["range"][0]),
                float(envelope["parameters"]["T_c"]["range"][1]),
            ],
            "alpha_fit_max": float(ALPHA_FIT_MAX),
            "thermal_link_model": "stycast_node",
            "electrical_link_model": ELECTRICAL_LINK_MODEL_RC,
            "electrical_rc_relaxation": {
                "topology": "series R_l + L plus one passive R||C relaxation element per TES branch",
                "fitted_parameters": ["R_rc", "f_rc_Hz"],
                "R_rc_search_ohm": [
                    float(R_RC_FIT_MIN_OHM),
                    float(R_RC_FIT_MAX_OHM),
                ],
                "f_rc_search_Hz": [
                    float(F_RC_FIT_MIN_HZ),
                    float(F_RC_FIT_MAX_HZ),
                ],
                "tes_resistance_response": "instantaneous alpha/beta (unchanged)",
            },
            "stycast_node": {
                "topology": "TES <-> Stycast <-> Pb absorber center",
                "symmetric_nodes": 2,
                "fitted_parameters": [
                    "C_stycast",
                    "G_tes-stycast",
                    "G_stycast-abs",
                ],
            },
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
                "loss": (
                    "band-balanced robust normalized log10(model/measurement) "
                    "shape residual; equalized broad bands plus high-frequency "
                    "ramp and broad-band mean-ratio penalty; plus one absolute "
                    "ASD level anchor at 1 kHz"
                ),
                "band_mean_penalty": float(FIT_BAND_MEAN_PENALTY),
                "absolute_asd_reference_Hz": float(
                    ABSOLUTE_ASD_REFERENCE_HZ
                ),
                "absolute_asd_weight": float(args.absolute_asd_weight),
                "absolute_asd_used_in_objective": bool(
                    args.absolute_asd_weight > 0.0
                ),
                "absolute_asd_calibration_status": (
                    "unresolved_for_exact_target; diagnostic only unless "
                    "--absolute-asd-weight is explicitly set > 0"
                ),
                "measurement_modelnoise_reported_units": "pA/rtHz",
                "internal_absolute_asd_units": "A/rtHz",
                "bands_Hz": [
                    {
                        "min": float(low),
                        "max": float(high),
                        "weight": float(weight),
                    }
                    for low, high, weight in FIT_BANDS_HZ
                ],
                "frequencies_outside_fit_band_ignored": True,
            },
            "cases": [
                {
                    **{
                        key: value
                        for key, value in case.items()
                        if key != "best_candidate"
                    },
                    "best_parameters": best_parameters_for_summary(
                        case["best_candidate"]
                    ),
                }
                for case in cases
            ],
            "best_case_R_SH_ohm": best_case["R_SH_ohm"],
            "best_case_R_TES_ohm": best_case["R_TES_ohm"],
            "best_case_score": best_case["best_score"],
            "best_case_finite_validation_score": best_case[
                "finite_validation_score"
            ],
            "best_case_parameters": best_parameters_for_summary(
                best_case["best_candidate"]
            ),
            "best_case_required_transfer_diagnostics": best_case[
                "required_transfer_diagnostics"
            ],
            "best_case_required_transfer_plot": best_case[
                "required_transfer_plot"
            ],
            "best_case_source_class_diagnostics": best_case[
                "source_class_diagnostics"
            ],
            "best_case_source_ablation_diagnostics": best_case[
                "source_ablation_diagnostics"
            ],
            "best_case_johnson_source_scale_diagnostics": best_case[
                "johnson_source_scale_diagnostics"
            ],
            "best_case_electrical_rc_ablation_diagnostics": best_case[
                "electrical_rc_ablation_diagnostics"
            ],
            "best_case_electrical_rc_degeneracy_diagnostics": best_case[
                "electrical_rc_degeneracy_diagnostics"
            ],
            "best_case_johnson_beta_audit": best_case[
                "johnson_beta_audit"
            ],
            "best_case_beta_sweep_diagnostics": best_case[
                "beta_sweep_diagnostics"
            ],
            "best_case_eigenmode_diagnostics": best_case[
                "eigenmode_diagnostics"
            ],
            "best_case_source_contribution_plot": best_case[
                "source_contribution_plot"
            ],
            "hardware_bessel_cutoff_Hz": float(
                TARGET_HARDWARE_BESSEL_CUTOFF_HZ
            ),
            "best_case_post_filter_white_fraction": best_case[
                "best_post_filter_white_fraction"
            ],
            "best_case_post_filter_white_asd_A_rtHz": best_case[
                "best_post_filter_white_asd_A_rtHz"
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
            final_candidate["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
            final_candidate["hardware_bessel_order"] = (
                TARGET_HARDWARE_BESSEL_ORDER
            )
            final_candidate["hardware_bessel_norm"] = TARGET_HARDWARE_BESSEL_NORM
            final_candidate["hardware_bessel_cutoff_Hz"] = (
                TARGET_HARDWARE_BESSEL_CUTOFF_HZ
            )
            final_candidate["thermal_link_model"] = "stycast_node"
            final_candidate["electrical_link_model"] = ELECTRICAL_LINK_MODEL_RC
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

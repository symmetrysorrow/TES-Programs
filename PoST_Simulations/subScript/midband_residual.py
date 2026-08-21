"""Quantify the CH0 mid-band residual after the best five-state model.

The comparison is performed before the digital 10 kHz analysis Bessel.  The
experimental spectrum is rebuilt from the selected raw records, while the
five-state model is generated in an isolated directory.  Several simple
independent-noise shapes are fitted in PSD, with a free scale for the detector
model because the published comparison normalizes each spectrum at 1 kHz.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, signal


SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))
import PoST_Simulation as post  # noqa: E402


EXPERIMENT_PATH = Path(
    r"G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_"
    r"difftrig5e-5_rate500k_samples100k_gain5_day2"
)
BASE_INPUT_PATH = Path(r"H:\hata2025\new\input.json")

# Best point in post_noise_bath_lr_sweep.csv (normalized 1--30 kHz score).
BEST_CHANGES = {
    "T_bath": 0.200,
    "L": 7.5e-7,
    "R": 0.022849268976,
}

# Stable internal-source-only candidate from a broad simultaneous search.  It
# is a mathematical existence proof, not yet a physically accepted parameter
# set; several values sit on the deliberately broad search bounds.
STABLE_INTERNAL_CHANGES = {
    "T_bath": 0.210,
    "L": 1.4477885504337302e-7,
    "R": 0.010,
    "C_tes": 1.3929355055125838e-12,
    "G_abs-tes": 4.4279364878930646e-8,
    "G_tes-bath": 2.2722137421653298e-9,
    "alpha": 861.8553493459841,
    "beta": 0.39938582819592644,
    "post_filter_white_asd_A_rtHz": 0.0,
    "readout_white_asd_A_rtHz": 0.0,
}


def experimental_asd(rate: float, sample: int, cutoff: float):
    window = np.hanning(sample)
    power_gain = np.sqrt(np.mean(window**2))
    numerator, denominator = signal.bessel(2, cutoff / (rate / 2), "low")
    amplitude = np.zeros(sample // 2 + 1)
    count = 0
    for path in (EXPERIMENT_PATH / "CH0_noise" / "rawdata").glob("CH0_*.dat"):
        values = np.frombuffer(path.read_bytes()[4:], dtype=np.float64).copy()
        if len(values) != sample or values.max() - values.min() > 0.04:
            continue
        values -= values.mean()
        values = signal.filtfilt(numerator, denominator, values)
        if values.max() - values.min() > 0.04:
            continue
        amplitude += np.abs(np.fft.rfft(values * window) / power_gain)
        count += 1
    if count == 0:
        raise RuntimeError("No experimental records passed selection")
    return amplitude / count, count


def analysis_magnitude(frequency: np.ndarray, rate: float, cutoff: float):
    numerator, denominator = signal.bessel(2, cutoff / (rate / 2), "low")
    _, response = signal.freqz(
        numerator,
        denominator,
        worN=2 * np.pi * frequency / rate,
    )
    return np.abs(response) ** 2  # forward/backward filtfilt magnitude


def simulated_pre_analysis_asd(
    work_dir: Path,
    frequency: np.ndarray,
    rate: float,
    parameters: dict,
):
    post.output = str(work_dir)
    post.MakeNoise()
    plt.close("all")
    with h5py.File(work_dir / "noise.h5", "r") as file:
        source_frequency = file["frequency"][:]
        detector = file["total"][:]
    order = int(parameters.get("hardware_bessel_order", 4))
    main = post.general.AnalogBesselMagnitudeResponse(
        frequency, 100_000.0, order=order
    )
    alias_frequency = rate - frequency
    folded = post.general.AnalogBesselMagnitudeResponse(
        alias_frequency, 100_000.0, order=order
    )
    detector_main = np.interp(frequency, source_frequency, detector) * main
    detector_fold = np.interp(alias_frequency, source_frequency, detector) * folded
    white = float(
        parameters.get(
            "post_filter_white_asd_A_rtHz",
            parameters.get("readout_white_asd_A_rtHz", 0.0),
        )
    )
    return np.sqrt(detector_main**2 + detector_fold**2 + white**2)


def normalize_at_1khz(values: np.ndarray, frequency: np.ndarray):
    return values / values[np.abs(frequency - 1_000.0).argmin()]


def binned_median(frequency, arrays, low=1_000.0, high=40_000.0, width=250.0):
    """Reduce FFT-bin scatter without changing the plotted band shape."""
    edges = np.arange(low, high + width, width)
    centers = 0.5 * (edges[:-1] + edges[1:])
    reduced = [np.full(len(centers), np.nan) for _ in arrays]
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (frequency >= left) & (frequency < right)
        for output_values, input_values in zip(reduced, arrays):
            output_values[index] = np.median(input_values[mask])
    return centers, reduced


def source_shape(name: str, frequency: np.ndarray, values: np.ndarray):
    if name == "white":
        amplitude = np.exp(values[1])
        return np.full_like(frequency, amplitude)
    if name == "lowpass":
        amplitude, corner, order = np.exp(values[1]), np.exp(values[2]), np.exp(values[3])
        return amplitude / np.sqrt(1.0 + (frequency / corner) ** (2.0 * order))
    if name == "bandpass":
        amplitude = np.exp(values[1])
        low = np.exp(values[2])
        high = low + np.exp(values[3])
        highpass = (frequency / low) / np.sqrt(1.0 + (frequency / low) ** 2)
        lowpass = 1.0 / np.sqrt(1.0 + (frequency / high) ** 2)
        return amplitude * highpass * lowpass
    raise ValueError(name)


def fit_candidate(name, frequency, experimental, simulated, fit_mask):
    if name == "white":
        initial = np.log([0.75, 0.5])
        lower = np.log([0.05, 1e-4])
        upper = np.log([2.0, 5.0])
    elif name == "lowpass":
        initial = np.log([0.75, 0.6, 15_000.0, 2.0])
        lower = np.log([0.05, 1e-4, 500.0, 0.25])
        upper = np.log([2.0, 5.0, 100_000.0, 8.0])
    else:
        initial = np.log([0.75, 0.8, 2_000.0, 20_000.0])
        lower = np.log([0.05, 1e-4, 100.0, 500.0])
        upper = np.log([2.0, 5.0, 50_000.0, 200_000.0])

    def residual(values):
        detector_scale = np.exp(values[0])
        extra = source_shape(name, frequency, values)
        combined = np.sqrt((detector_scale * simulated) ** 2 + extra**2)
        return np.log(combined[fit_mask] / experimental[fit_mask])

    result = optimize.least_squares(residual, initial, bounds=(lower, upper))
    values = result.x
    detector_scale = np.exp(values[0])
    extra = source_shape(name, frequency, values)
    combined = np.sqrt((detector_scale * simulated) ** 2 + extra**2)
    rms_log = float(np.sqrt(np.mean(residual(values) ** 2)))
    result_values = {
        "name": name,
        "detector_scale": detector_scale,
        "source_amplitude_at_dc": float(extra[0]),
        "rms_log": rms_log,
        "rms_factor": float(np.exp(rms_log)),
        "combined": combined,
        "source": extra,
    }
    if name == "lowpass":
        result_values["corner_hz"] = float(np.exp(values[2]))
        result_values["order"] = float(np.exp(values[3]))
    elif name == "bandpass":
        low = float(np.exp(values[2]))
        result_values["low_corner_hz"] = low
        result_values["high_corner_hz"] = low + float(np.exp(values[3]))
    return result_values


def fit_anchored_lowpass(
    frequency,
    experimental,
    simulated,
    detector_scale,
    fit_mask,
):
    """Fit the minimum extra source after anchoring the model near 30 kHz."""

    # Both spectra equal one at 1 kHz.  Fix the source amplitude required by
    # PSD closure there; otherwise an unconstrained fit can suppress the
    # five-state contribution entirely because absolute calibration was lost
    # during independent normalization.
    amplitude = np.sqrt(max(1.0 - detector_scale**2, 0.0))
    fit_frequency = frequency[fit_mask]
    log_grid = np.geomspace(fit_frequency[0], fit_frequency[-1], 300)
    fit_indices = np.unique(
        [np.abs(frequency - value).argmin() for value in log_grid]
    )

    def evaluate(values):
        corner, order = np.exp(values)
        source = amplitude / np.sqrt(
            1.0 + (frequency / corner) ** (2.0 * order)
        )
        combined = np.sqrt((detector_scale * simulated) ** 2 + source**2)
        return source, combined

    def residual(values):
        _, combined = evaluate(values)
        return np.log(combined[fit_indices] / experimental[fit_indices])

    result = optimize.least_squares(
        residual,
        np.log([12_000.0, 2.5]),
        bounds=(
            np.log([500.0, 0.25]),
            np.log([100_000.0, 8.0]),
        ),
    )
    corner, order = np.exp(result.x)
    source, combined = evaluate(result.x)
    rms_log = float(np.sqrt(np.mean(residual(result.x) ** 2)))
    return {
        "name": "anchored_lowpass",
        "detector_scale": detector_scale,
        "source_amplitude_at_dc": float(amplitude),
        "corner_hz": float(corner),
        "order": float(order),
        "rms_factor": float(np.exp(rms_log)),
        "source": source,
        "combined": combined,
    }


def fit_multiplicative_biquad(frequency, experimental, simulated, fit_mask):
    """Fit the smallest pole/zero form that reproduces the ASD ratio hump.

    H(s) = (1 + s/wz)^2 / (1 + s/(Q*w0) + (s/w0)^2)

    Its two real zeros prevent the resonant pole pair from continuing as a
    1/f^2 roll-off.  The response is normalized at 1 kHz like the data.
    """

    ratio = experimental / simulated
    reference = np.abs(frequency - 1_000.0).argmin()

    def response(values):
        pole_frequency, quality, zero_frequency = np.exp(values)
        numerator = 1.0 + (frequency / zero_frequency) ** 2
        denominator = np.sqrt(
            (1.0 - (frequency / pole_frequency) ** 2) ** 2
            + (frequency / (quality * pole_frequency)) ** 2
        )
        magnitude = numerator / denominator
        return magnitude / magnitude[reference]

    def residual(values):
        return np.log(response(values)[fit_mask] / ratio[fit_mask])

    result = optimize.least_squares(
        residual,
        np.log([13_000.0, 0.8, 18_000.0]),
        bounds=(
            np.log([2_000.0, 0.1, 1_000.0]),
            np.log([100_000.0, 10.0, 1_000_000.0]),
        ),
    )
    pole_frequency, quality, zero_frequency = np.exp(result.x)
    rms_log = float(np.sqrt(np.mean(residual(result.x) ** 2)))
    return {
        "pole_frequency_hz": float(pole_frequency),
        "quality_factor": float(quality),
        "zero_frequency_hz": float(zero_frequency),
        "high_frequency_gain": float((pole_frequency / zero_frequency) ** 2),
        "rms_factor": float(np.exp(rms_log)),
        "response": response(result.x),
    }


def five_state_poles_and_zeros(parameters: dict):
    c_abs = float(parameters["C_abs"])
    c_tes = float(parameters["C_tes"])
    g_abs_abs = float(parameters["G_abs-abs"])
    g_abs_tes = float(parameters["G_abs-tes"])
    g_bath = float(parameters["G_tes-bath"])
    resistance = float(parameters["R"])
    r_load = float(parameters["R_l"])
    t_c = float(parameters["T_c"])
    t_bath = float(parameters["T_bath"])
    alpha = float(parameters["alpha"])
    beta = float(parameters["beta"])
    inductance = float(parameters["L"])
    exponent = float(parameters["n"])
    current = np.sqrt(
        g_bath * t_c * (1.0 - (t_bath / t_c) ** exponent)
        / (exponent * resistance)
    )
    loop_gain = alpha * current**2 * resistance / (g_bath * t_c)
    tau_el = inductance / (r_load + resistance * (1.0 + beta))
    tau_i = c_tes / ((1.0 - loop_gain) * g_bath)
    g_eff = 1.0 / (1.0 / g_abs_tes + 1.0 / (2.0 * g_abs_abs))

    matrix = np.zeros((5, 5))
    matrix[0, 0] = 1.0 / tau_el
    matrix[0, 1] = loop_gain * g_bath / (current * inductance)
    matrix[1, 0] = -current * resistance * (2.0 + beta) / c_tes
    matrix[1, 1] = 1.0 / tau_i + g_eff / c_tes
    matrix[1, 2] = -g_eff / c_tes
    matrix[2, 1] = -g_eff / c_abs
    matrix[2, 2] = 2.0 * g_eff / c_abs
    matrix[2, 3] = -g_eff / c_abs
    matrix[3, 2] = -g_eff / c_tes
    matrix[3, 3] = 1.0 / tau_i + g_eff / c_tes
    matrix[3, 4] = -current * resistance * (2.0 + beta) / c_tes
    matrix[4, 3] = loop_gain * g_bath / (current * inductance)
    matrix[4, 4] = 1.0 / tau_el

    # Unit heat-current source at TES1 and current readout at CH0.
    source = np.zeros((5, 1))
    source[1, 0] = 1.0 / c_tes
    readout = np.zeros((1, 5))
    readout[0, 0] = 1.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", signal.BadCoefficients)
        zeros, poles, _ = signal.ss2zpk(
            -matrix, source, readout, np.zeros((1, 1))
        )
    electrical_rate = 1.0 / tau_el
    open_loop_thermal_rate = (
        (1.0 - loop_gain) * g_bath + g_eff
    ) / c_tes
    critical_inductance = (
        r_load + resistance * (1.0 + beta)
    ) / max(-open_loop_thermal_rate, np.finfo(float).tiny)
    critical_c_tes = max(
        ((loop_gain - 1.0) * g_bath - g_eff) / electrical_rate,
        0.0,
    )
    stability = {
        "loop_gain": float(loop_gain),
        "electrical_rate_per_s": float(electrical_rate),
        "open_loop_thermal_rate_per_s": float(open_loop_thermal_rate),
        "critical_inductance_h": float(critical_inductance),
        "critical_c_tes_j_per_k": float(critical_c_tes),
        "inductance_over_critical": float(inductance / critical_inductance),
        "c_tes_under_critical_factor": float(critical_c_tes / c_tes),
    }
    return poles, zeros, stability


def complex_frequency(value: complex):
    return abs(value) / (2.0 * np.pi)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    base_parameters = json.loads(BASE_INPUT_PATH.read_text(encoding="utf-8"))
    parameters = dict(base_parameters)
    parameters.update(BEST_CHANGES)
    parameters["samples"] = 10_000
    (args.work_dir / "input.json").write_text(
        json.dumps(parameters, indent=2), encoding="utf-8"
    )

    rate = float(parameters["rate"])
    sample = 100_000
    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    measured_filtered, records = experimental_asd(
        rate, sample, float(parameters["cutoff"])
    )
    digital = analysis_magnitude(frequency, rate, float(parameters["cutoff"]))
    experimental = normalize_at_1khz(measured_filtered / digital, frequency)
    simulated = normalize_at_1khz(
        simulated_pre_analysis_asd(args.work_dir, frequency, rate, parameters),
        frequency,
    )

    fit_mask = (frequency >= 1_000.0) & (frequency <= 40_000.0)
    baseline_log = np.log(simulated[fit_mask] / experimental[fit_mask])
    baseline_factor = float(np.exp(np.sqrt(np.mean(baseline_log**2))))
    unconstrained_fits = [
        fit_candidate(name, frequency, experimental, simulated, fit_mask)
        for name in ("white", "lowpass", "bandpass")
    ]
    unconstrained_fits.sort(key=lambda item: item["rms_log"])

    # The independently normalized spectra do not determine their relative
    # absolute scale.  Anchor the physical-model contribution where the two
    # shapes re-approach, then infer the smallest extra PSD compatible with it.
    anchor_mask = (frequency >= 28_000.0) & (frequency <= 32_000.0)
    anchor_scale = float(np.median(experimental[anchor_mask] / simulated[anchor_mask]))
    residual_fit_mask = (frequency >= 1_000.0) & (frequency <= 30_000.0)
    best = fit_anchored_lowpass(
        frequency,
        experimental,
        simulated,
        anchor_scale,
        residual_fit_mask,
    )
    biquad = fit_multiplicative_biquad(
        frequency,
        experimental,
        simulated,
        residual_fit_mask,
    )

    stable_parameters = dict(base_parameters)
    stable_parameters.update(STABLE_INTERNAL_CHANGES)
    stable_parameters["samples"] = 10_000
    stable_work_dir = args.work_dir / "stable_internal"
    stable_work_dir.mkdir(parents=True, exist_ok=True)
    (stable_work_dir / "input.json").write_text(
        json.dumps(stable_parameters, indent=2), encoding="utf-8"
    )
    stable_simulated = normalize_at_1khz(
        simulated_pre_analysis_asd(
            stable_work_dir, frequency, rate, stable_parameters
        ),
        frequency,
    )

    fixed_residual_psd = experimental**2 - simulated**2
    scaled_residual_psd = experimental**2 - (
        best["detector_scale"] * simulated
    ) ** 2

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "frequency_hz,experiment_preanalysis_normalized,"
        "five_state_preanalysis_normalized,fixed_scale_residual_psd,"
        "best_scaled_residual_psd,best_source_asd,best_combined_asd,"
        "biquad_corrected_asd,stable_internal_five_state_asd\n"
    )
    data = np.column_stack(
        (
            frequency,
            experimental,
            simulated,
            fixed_residual_psd,
            scaled_residual_psd,
            best["source"],
            best["combined"],
            simulated * biquad["response"],
            stable_simulated,
        )
    )
    np.savetxt(args.csv, data, delimiter=",", header=header.rstrip(), comments="")

    # Validation-facing plot: show exactly the spectra after the same digital
    # analysis filter.  The stable internal-only curve is still a fitted
    # parameter candidate; agreement is not independent model validation.
    experiment_final = normalize_at_1khz(measured_filtered, frequency)
    five_state_final = normalize_at_1khz(simulated * digital, frequency)
    stable_final = normalize_at_1khz(stable_simulated * digital, frequency)
    plot_frequency, plot_arrays = binned_median(
        frequency,
        (experiment_final, five_state_final, stable_final),
    )
    plot_experiment, plot_five_state, plot_stable = plot_arrays
    current_ratio = plot_five_state / plot_experiment
    stable_ratio = plot_stable / plot_experiment
    comparison_band = (plot_frequency >= 1_000.0) & (plot_frequency <= 30_000.0)
    current_inside = 100.0 * np.mean(
        (current_ratio[comparison_band] >= 0.9)
        & (current_ratio[comparison_band] <= 1.1)
    )
    stable_inside = 100.0 * np.mean(
        (stable_ratio[comparison_band] >= 0.9)
        & (stable_ratio[comparison_band] <= 1.1)
    )

    args.plot.parent.mkdir(parents=True, exist_ok=True)
    figure, (spectrum_axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
        layout="constrained",
        gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.06},
    )
    spectrum_axis.loglog(
        plot_frequency,
        plot_experiment,
        color="black",
        linewidth=2.0,
        label="Experiment: final analyzed noise",
    )
    spectrum_axis.loglog(
        plot_frequency,
        plot_five_state,
        color="tab:orange",
        linewidth=1.7,
        linestyle="--",
        label="Best tested near-current-parameter model",
    )
    spectrum_axis.loglog(
        plot_frequency,
        plot_stable,
        color="tab:blue",
        linewidth=2.0,
        label="Stable internal-source-only five-state candidate",
    )
    spectrum_axis.set_ylabel("Normalized final ASD")
    spectrum_axis.set_title("CH0 final-noise agreement after identical analysis")
    spectrum_axis.grid(True, which="both", alpha=0.2)
    spectrum_axis.legend(loc="lower left")
    spectrum_axis.text(
        0.98,
        0.96,
        "Internal-only candidate was fitted here; independent parameter checks remain",
        transform=spectrum_axis.transAxes,
        horizontalalignment="right",
        verticalalignment="top",
        fontsize=9,
    )

    ratio_axis.axhspan(0.9, 1.1, color="tab:green", alpha=0.14, label="within +/-10%")
    ratio_axis.axhline(1.0, color="black", linewidth=1.0)
    ratio_axis.axvline(30_000.0, color="0.45", linewidth=1.0, linestyle=":")
    ratio_axis.semilogx(
        plot_frequency,
        current_ratio,
        color="tab:orange",
        linewidth=1.5,
        linestyle="--",
        label=f"Near-current model ({current_inside:.0f}% inside band)",
    )
    ratio_axis.semilogx(
        plot_frequency,
        stable_ratio,
        color="tab:blue",
        linewidth=1.8,
        label=f"Stable internal-only candidate ({stable_inside:.0f}% inside band)",
    )
    ratio_axis.set_xlim(1_000.0, 40_000.0)
    ratio_axis.set_ylim(0.45, 1.55)
    ratio_axis.set_xlabel("Frequency [Hz]")
    ratio_axis.set_ylabel("Simulation / experiment")
    ratio_axis.grid(True, which="both", alpha=0.2)
    ratio_axis.legend(loc="lower left", fontsize=9, ncol=2)
    plt.savefig(args.plot, dpi=220)
    plt.close()

    print(f"experimental records: {records}")
    print(f"baseline 1--40 kHz RMS factor: {baseline_factor:.6f}")
    print("unconstrained additive fits (shows normalization degeneracy):")
    for item in unconstrained_fits:
        summary = {
            key: value
            for key, value in item.items()
            if key not in {"combined", "source", "rms_log"}
        }
        print(json.dumps(summary, sort_keys=True))
    print("30 kHz anchored additive fit:", json.dumps({
        key: value for key, value in best.items()
        if key not in {"combined", "source"}
    }, sort_keys=True))
    print("multiplicative pole/zero fit:", json.dumps({
        key: value for key, value in biquad.items() if key != "response"
    }, sort_keys=True))
    print("stable internal-source-only candidate:", json.dumps(
        STABLE_INTERNAL_CHANGES, sort_keys=True
    ))
    print(f"fixed-scale negative PSD bins in 1--40 kHz: {np.count_nonzero(fixed_residual_psd[fit_mask] < 0)} / {np.count_nonzero(fit_mask)}")
    print("frequency_khz exp sim fixed_residual_asd best_source_asd best_combined")
    for target in (1_000, 3_000, 5_000, 7_000, 10_000, 15_000, 20_000, 30_000, 40_000):
        index = np.abs(frequency - target).argmin()
        fixed_asd = np.sqrt(max(fixed_residual_psd[index], 0.0))
        print(
            f"{frequency[index] / 1000:8.1f} {experimental[index]:.6f} "
            f"{simulated[index]:.6f} {fixed_asd:.6f} "
            f"{best['source'][index]:.6f} {best['combined'][index]:.6f}"
        )

    poles, zeros, stability = five_state_poles_and_zeros(parameters)
    print("five-state s-plane poles:", [str(value) for value in poles])
    print("five-state pole |f| (Hz):", [complex_frequency(value) for value in poles])
    finite_zeros = [value for value in zeros if complex_frequency(value) < 1e12]
    print("TES1 heat-to-CH0-current finite zeros (Hz):", [complex_frequency(value) for value in finite_zeros])
    print("stability:", json.dumps(stability, sort_keys=True))
    print(f"saved: {args.csv}")
    print(f"saved: {args.plot}")


if __name__ == "__main__":
    main()

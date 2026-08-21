"""Coarse L/R sweep against the CH0 experimental noise spectrum.

This is an isolated diagnostic: it writes only to the output directory given
by --work-dir, never to the nominal simulation output directory.
"""

import argparse
import json
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))
import PoST_Simulation as post  # noqa: E402


EXPERIMENT_PATH = Path(
    r"G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_"
    r"difftrig5e-5_rate500k_samples100k_gain5_day2"
)
BASE_INPUT_PATH = Path(r"H:\hata2025\new\input.json")


def experimental_asd(rate, sample, cutoff):
    """Rebuild the current CH0 analysis spectrum without saving it."""
    window = np.hanning(sample)
    power_gain = np.sqrt(np.mean(window**2))
    b, a = signal.bessel(2, cutoff / (rate / 2), "low")
    amplitude = np.zeros(sample // 2 + 1)
    count = 0
    for path in (EXPERIMENT_PATH / "CH0_noise" / "rawdata").glob("CH0_*.dat"):
        values = np.frombuffer(path.read_bytes()[4:], dtype=np.float64).copy()
        if len(values) != sample or values.max() - values.min() > 0.04:
            continue
        values -= values.mean()
        values = signal.filtfilt(b, a, values)
        if values.max() - values.min() > 0.04:
            continue
        amplitude += np.abs(np.fft.rfft(values * window) / power_gain)
        count += 1
    if count == 0:
        raise RuntimeError("No experimental records passed selection")
    return amplitude / count, count


def analysis_magnitude(frequency, rate, cutoff):
    b, a = signal.bessel(2, cutoff / (rate / 2), "low")
    _, response = signal.freqz(b, a, worN=2 * np.pi * frequency / rate)
    return np.abs(response) ** 2  # filtfilt


def simulated_asd(work_dir, frequency, rate, parameters):
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
    folded = post.general.AnalogBesselMagnitudeResponse(
        rate - frequency, 100_000.0, order=order
    )
    detector_main = np.interp(frequency, source_frequency, detector) * main
    detector_fold = np.interp(rate - frequency, source_frequency, detector) * folded
    white = float(parameters.get("post_filter_white_asd_A_rtHz", 0.0))
    return np.sqrt(detector_main**2 + detector_fold**2 + white**2)


def normalized(values, frequency):
    return values / values[np.abs(frequency - 1_000.0).argmin()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--sweep", choices=("lr", "etf", "thermal", "bath", "bath-lr"), default="lr")
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    parameters = json.loads(BASE_INPUT_PATH.read_text(encoding="utf-8"))
    rate = float(parameters["rate"])
    sample = int(parameters["samples"])
    frequency = np.fft.rfftfreq(sample, d=1 / rate)
    measured, records = experimental_asd(rate, sample, float(parameters["cutoff"]))
    measured = normalized(measured, frequency)
    filter_response = analysis_magnitude(frequency, rate, float(parameters["cutoff"]))
    band = (frequency >= 1_000) & (frequency <= 30_000)

    base_l = float(parameters["L"])
    base_r = float(parameters["R"])
    if args.sweep == "lr":
        trials = [
            {"L": inductance, "R": base_r * scale}
            for inductance in [base_l, 1.0e-6, 7.5e-7, 5.0e-7, 3.0e-7, 1.0e-7]
            for scale in [0.5, 0.75, 1.0, 1.25, 1.5]
        ]
    elif args.sweep == "etf":
        # Hold L/R at the coarse-sweep optimum while testing whether the
        # electrothermal-feedback parameters can form the observed shoulder.
        trials = [
            {
                "L": 1.0e-6,
                "R": base_r * 1.25,
                "alpha": parameters["alpha"] * alpha_scale,
                "beta": parameters["beta"] * beta_scale,
            }
            for alpha_scale in [0.5, 0.75, 1.0, 1.25, 1.5]
            for beta_scale in [0.5, 0.75, 1.0, 1.25, 1.5]
        ]
    elif args.sweep == "thermal":
        # One-at-a-time thermal scan.  The range is deliberately broad; a
        # candidate must improve the noise score before pulse work is run.
        thermal_keys = ("C_tes", "G_abs-tes", "G_tes-bath")
        trials = []
        for key in thermal_keys:
            for scale in [0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 10.0]:
                trials.append(
                    {
                        "L": 1.0e-6,
                        "R": base_r * 1.25,
                        key: parameters[key] * scale,
                        "_varied_key": key,
                        "_scale": scale,
                    }
                )
    elif args.sweep == "bath":
        # T_bath must remain below T_c.  This scan spans the plausible
        # sub-transition operating range without approaching the singular
        # T_bath == T_c limit.
        trials = [
            {
                "L": 1.0e-6,
                "R": base_r * 1.25,
                "T_bath": temperature,
                "_varied_key": "T_bath",
                "_scale": temperature / parameters["T_bath"],
            }
            for temperature in [0.195, 0.1975, 0.200, 0.2025, 0.205, 0.2075, 0.210]
        ]
    else:
        # Re-optimise the electrical knee at two physically plausible bath
        # temperatures; a fixed-L/R bath scan is not sufficient here.
        trials = [
            {
                "L": inductance,
                "R": base_r * r_scale,
                "T_bath": temperature,
                "_varied_key": "T_bath+L+R",
                "_scale": temperature / parameters["T_bath"],
            }
            for temperature in [0.195, 0.200]
            for inductance in [1.0e-6, 7.5e-7, 5.0e-7, 3.0e-7]
            for r_scale in [0.5, 0.75, 1.0, 1.25]
        ]
    rows = []
    for changes in trials:
        trial = dict(parameters)
        varied_key = changes.pop("_varied_key", "multiple")
        varied_scale = changes.pop("_scale", 1.0)
        trial.update(changes)
        # MakeNoise only needs a smooth source transfer function here.
        # 10k points gives 50 Hz resolution, while the comparison grid
        # remains the original 100k-point (5 Hz) analysis frequency grid.
        trial["samples"] = 10_000
        (args.work_dir / "input.json").write_text(
            json.dumps(trial, indent=2), encoding="utf-8"
        )
        incoming = simulated_asd(args.work_dir, frequency, rate, trial)
        simulated = normalized(incoming * filter_response, frequency)
        log_error = np.log10(simulated[band] / measured[band])
        score = float(np.sqrt(np.mean(log_error**2)))
        i10 = np.abs(frequency - 10_000).argmin()
        i30 = np.abs(frequency - 30_000).argmin()
        rows.append((score, varied_key, varied_scale, trial["L"], trial["R"], trial["C_tes"], trial["G_abs-tes"], trial["G_tes-bath"], trial["alpha"], trial["beta"], simulated[i10] / measured[i10], simulated[i30] / measured[i30]))
        print(
            f"L={trial['L']:.3e} H R={trial['R']:.5f} ohm "
            f"{varied_key}×{varied_scale:.2g} "
            f"alpha={trial['alpha']:.2f} beta={trial['beta']:.2f} "
            f"score={score:.5f} ratio10k={simulated[i10]/measured[i10]:.3f} "
            f"ratio30k={simulated[i30]/measured[i30]:.3f}",
            flush=True,
        )
    rows.sort()
    args.result.parent.mkdir(parents=True, exist_ok=True)
    with args.result.open("w", encoding="utf-8") as file:
        file.write("score_log_rms,varied_key,varied_scale,L_henry,R_ohm,C_tes,G_abs_tes,G_tes_bath,alpha,beta,sim_over_exp_10khz,sim_over_exp_30khz\n")
        for row in rows:
            file.write(",".join(value if isinstance(value, str) else f"{value:.12g}" for value in row) + "\n")
    print(f"experimental records: {records}")
    print(f"saved: {args.result}")
    print("best:", rows[0])


if __name__ == "__main__":
    main()

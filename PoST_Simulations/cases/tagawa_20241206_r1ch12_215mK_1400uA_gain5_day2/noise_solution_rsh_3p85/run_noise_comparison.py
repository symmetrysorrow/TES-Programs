"""Regenerate the archived PoST noise solution and compare it with CH0 data."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ARCHIVE_DIR = Path(__file__).resolve().parent
SIM_DIR = ARCHIVE_DIR.parents[2]
POST_SCRIPT = SIM_DIR / "PoST_Simulation.py"
INPUT_PATH = ARCHIVE_DIR / "input.json"
DEFAULT_MEASURED_NOISE = Path(
    r"G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_"
    r"difftrig5e-5_rate500k_samples100k_gain5_day2\CH0_noise\modelnoise.txt"
)
OUTPUT_NAME = "noise_total-bessel100k.dat"


def positive_spectrum(values: np.ndarray, rate_hz: float):
    frequency = np.arange(len(values), dtype=float) * (rate_hz / 2.0) / len(values)
    values = np.asarray(values, dtype=float)
    mask = (frequency > 0.0) & np.isfinite(values) & (values > 0.0)
    return frequency[mask], values[mask]


def normalize_at(frequency: np.ndarray, values: np.ndarray, reference_hz: float):
    index = int(np.argmin(np.abs(frequency - reference_hz)))
    if values[index] <= 0.0 or not np.isfinite(values[index]):
        raise ValueError(f"Invalid normalization point near {reference_hz:g} Hz")
    return values / values[index]


def generate_noise(output_dir: Path):
    with tempfile.TemporaryDirectory(prefix=".noise_run_", dir=ARCHIVE_DIR) as temp_name:
        temp_dir = Path(temp_name)
        shutil.copy2(INPUT_PATH, temp_dir / "input.json")
        subprocess.run(
            [
                sys.executable,
                str(POST_SCRIPT),
                "--noise-only",
                "--output",
                str(temp_dir),
            ],
            cwd=str(SIM_DIR),
            check=True,
        )
        generated_path = temp_dir / OUTPUT_NAME
        if not generated_path.is_file():
            raise FileNotFoundError(f"PoST did not generate {generated_path}")
        output_dir.mkdir(parents=True, exist_ok=True)
        archived_path = output_dir / OUTPUT_NAME
        shutil.copy2(generated_path, archived_path)
        generated_plot = temp_dir / "noise_total-bessel100k.png"
        if generated_plot.is_file():
            shutil.copy2(generated_plot, output_dir / generated_plot.name)
        return archived_path


def compare(measured_path: Path, generated_path: Path, output_dir: Path,
            measured_rate_hz: float, reference_hz: float,
            fit_min_hz: float, fit_max_hz: float):
    parameters = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
    measured = np.loadtxt(measured_path)
    generated = np.loadtxt(generated_path)
    measured_frequency, measured_asd = positive_spectrum(measured, measured_rate_hz)
    generated_frequency, generated_asd = positive_spectrum(
        generated,
        float(parameters["rate"]),
    )
    measured_normalized = normalize_at(measured_frequency, measured_asd, reference_hz)
    generated_normalized = normalize_at(generated_frequency, generated_asd, reference_hz)
    measured_on_model = np.interp(
        generated_frequency,
        measured_frequency,
        measured_normalized,
    )
    ratio = generated_normalized / measured_on_model
    fit_mask = (generated_frequency >= fit_min_hz) & (generated_frequency <= fit_max_hz)
    if not np.any(fit_mask):
        raise ValueError("The requested fit band contains no generated frequency bins")

    figure, (spectrum_axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    spectrum_axis.loglog(
        measured_frequency,
        measured_normalized,
        color="black",
        linewidth=2.0,
        label="Measured CH0",
    )
    spectrum_axis.loglog(
        generated_frequency,
        generated_normalized,
        color="tab:blue",
        linewidth=2.2,
        label="Archived PoST solution",
    )
    spectrum_axis.axvspan(
        fit_min_hz,
        fit_max_hz,
        color="gray",
        alpha=0.08,
        label="fit band",
    )
    spectrum_axis.set_ylabel("Normalized ASD")
    spectrum_axis.set_title("Archived PoST noise solution vs measured CH0 noise")
    spectrum_axis.grid(True, which="both", alpha=0.25)
    spectrum_axis.legend(fontsize=9)

    ratio_axis.semilogx(
        generated_frequency,
        ratio,
        color="tab:blue",
        linewidth=1.5,
        label="Model / measured",
    )
    ratio_axis.axhline(1.0, color="black", linewidth=1.0)
    ratio_axis.fill_between(
        [fit_min_hz, fit_max_hz],
        [0.9, 0.9],
        [1.1, 1.1],
        color="gray",
        alpha=0.15,
        label="±10%",
    )
    ratio_axis.set_yscale("log")
    ratio_axis.set_xlabel("Frequency [Hz]")
    ratio_axis.set_ylabel("Model / measured")
    ratio_axis.set_xlim(measured_frequency[0], measured_frequency[-1])
    ratio_axis.grid(True, which="both", alpha=0.25)
    ratio_axis.legend(fontsize=8)
    figure.tight_layout()

    plot_path = output_dir / "noise_comparison.png"
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)

    log_residual = np.log10(ratio[fit_mask])
    result = {
        "input": str(INPUT_PATH),
        "measured_noise": str(measured_path),
        "generated_noise": str(generated_path),
        "comparison_plot": str(plot_path),
        "R_SH_ohm": parameters["R_SH"],
        "R_TES_ohm": parameters["R"],
        "alpha": parameters["alpha"],
        "beta": parameters["beta"],
        "L_H": parameters["L"],
        "frequency_range_Hz": [
            float(measured_frequency[0]),
            float(measured_frequency[-1]),
        ],
        "fit_band_Hz": [fit_min_hz, fit_max_hz],
        "fit_score_mean_squared_log10_residual": float(np.mean(log_residual ** 2)),
        "fit_ratio_median": float(np.median(ratio[fit_mask])),
        "fit_ratio_p05_p95": [
            float(value) for value in np.percentile(ratio[fit_mask], [5, 95])
        ],
    }
    (output_dir / "comparison_result.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--measured-noise", type=Path, default=DEFAULT_MEASURED_NOISE)
    parser.add_argument("--output-dir", type=Path, default=ARCHIVE_DIR)
    parser.add_argument("--measured-rate-hz", type=float, default=500_000.0)
    parser.add_argument("--reference-hz", type=float, default=1_000.0)
    parser.add_argument("--fit-min-hz", type=float, default=1_000.0)
    parser.add_argument("--fit-max-hz", type=float, default=200_000.0)
    args = parser.parse_args()
    if not POST_SCRIPT.is_file():
        raise FileNotFoundError(f"Missing PoST simulator: {POST_SCRIPT}")
    if not INPUT_PATH.is_file():
        raise FileNotFoundError(f"Missing archived input: {INPUT_PATH}")
    if not args.measured_noise.is_file():
        raise FileNotFoundError(
            f"Missing measured noise file: {args.measured_noise}\n"
            "Pass its location with --measured-noise."
        )
    generated = generate_noise(args.output_dir)
    result = compare(
        args.measured_noise,
        generated,
        args.output_dir,
        args.measured_rate_hz,
        args.reference_hz,
        args.fit_min_hz,
        args.fit_max_hz,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

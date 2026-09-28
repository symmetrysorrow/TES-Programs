"""Compare the modeled CH0 noise spectrum for two TES inductances."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib.tes_noise_model import noise_components


L_VALUES_H = (12.3e-9, 1.0e-10)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("input.json"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("figures/noise_spectrum_L_comparison.png"),
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("figures/noise_spectrum_L_comparison.csv"),
    )
    args = parser.parse_args()

    parameters = json.loads(args.input.read_text(encoding="utf-8"))
    rate = float(parameters["rate"])
    frequencies = np.geomspace(1.0, rate / 2.0, 2400)
    asd_by_inductance = {}
    for inductance in L_VALUES_H:
        case_parameters = dict(parameters)
        case_parameters["L"] = inductance
        asd_by_inductance[inductance] = noise_components(
            case_parameters, frequencies
        )["total_ch0"] * 1e6  # A/sqrt(Hz) -> uA/sqrt(Hz)

    ratio = asd_by_inductance[L_VALUES_H[0]] / asd_by_inductance[L_VALUES_H[1]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.csv,
        np.column_stack(
            (frequencies, asd_by_inductance[L_VALUES_H[0]],
             asd_by_inductance[L_VALUES_H[1]], ratio)
        ),
        delimiter=",",
        header="frequency_Hz,ASD_CH0_L_12.3nH_uA_rtHz,ASD_CH0_L_0.1nH_uA_rtHz,ratio_12.3nH_over_0.1nH",
        comments="",
    )

    fig, (ax, ratio_ax) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
        constrained_layout=True,
    )
    ax.loglog(
        frequencies, asd_by_inductance[L_VALUES_H[0]],
        linewidth=2, label="L = 12.3 nH",
    )
    ax.loglog(
        frequencies, asd_by_inductance[L_VALUES_H[1]],
        linewidth=2, label="L = 0.1 nH (1e-10 H)",
    )
    ax.set_ylabel(r"CH0 noise ASD [$\mu$A/$\sqrt{\mathrm{Hz}}$]")
    ax.set_title("Noise spectrum comparison (only L varied)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()

    ratio_ax.semilogx(frequencies, ratio, color="tab:green", linewidth=1.7)
    ratio_ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    ratio_ax.set_xlabel("Frequency [Hz]")
    ratio_ax.set_ylabel("ASD ratio\n(12.3 nH / 0.1 nH)")
    ratio_ax.grid(True, which="both", alpha=0.25)

    fig.savefig(args.output, dpi=250)
    print(f"Saved figure: {args.output}")
    print(f"Saved spectrum data: {args.csv}")


if __name__ == "__main__":
    main()

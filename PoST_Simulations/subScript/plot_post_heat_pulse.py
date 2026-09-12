# -*- coding: utf-8 -*-

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]

TES_CPP_PYTHON = ROOT / "tes_cpp" / "python"
if str(TES_CPP_PYTHON) not in sys.path:
    sys.path.insert(0, str(TES_CPP_PYTHON))

from tes_cpp import posi2pulse

INPUT_PATH = ROOT / "PoST_Simulations" / "input.json"

HEAT_DISTRIBUTION = {
    17: 0.70,
    16: 0.20,
    15: 0.10,
}

CHANNEL = "ch0"

XMIN_MS = 0.0
XMAX_MS = 10.0

CURRENT_SCALE = 1e6
CURRENT_UNIT = r"$\mu$A"

SAVE_FIGURE = True
OUTPUT_FIGURE = (
    ROOT / "PoST_Simulations" / "subScript" / "split_heat_pulse.png"
)

# 合計は黒、部分は濃い灰色系
BLOCK_STYLES = [
    {"color": "dimgray",   "linestyle": "--"},
    {"color": "gray",      "linestyle": "-."},
    {"color": "gray", "linestyle": ":"},
    {"color": "gray",  "linestyle": "--"},
]

TOTAL_COLOR = "black"


def main():
    total_fraction = sum(HEAT_DISTRIBUTION.values())
    if not np.isclose(total_fraction, 1.0):
        raise ValueError(
            f"Heat fractions must sum to 1.0 (current sum = {total_fraction:.6f})"
        )

    positions = list(HEAT_DISTRIBUTION.keys())

    pulses = posi2pulse(INPUT_PATH, positions)
    if len(pulses) == 0:
        raise RuntimeError("No pulses were generated.")

    time_s = np.asarray(pulses[0].time, dtype=float)
    time_ms = time_s * 1e3

    weighted_pulses = {}

    for pulse in pulses:
        position = int(pulse.position)

        if CHANNEL == "ch0":
            waveform = np.asarray(pulse.ch0, dtype=float)
        elif CHANNEL == "ch1":
            waveform = np.asarray(pulse.ch1, dtype=float)
        else:
            raise ValueError("CHANNEL must be 'ch0' or 'ch1'")

        fraction = HEAT_DISTRIBUTION[position]
        weighted_pulses[position] = fraction * waveform

    total_waveform = np.zeros_like(next(iter(weighted_pulses.values())))
    for waveform in weighted_pulses.values():
        total_waveform += waveform

    total_plot = total_waveform * CURRENT_SCALE
    block_plot = {
        position: waveform * CURRENT_SCALE
        for position, waveform in weighted_pulses.items()
    }

    fig, ax = plt.subplots(figsize=(8.0, 5.5))

    # 部分波形
    for i, (position, waveform) in enumerate(block_plot.items()):
        fraction = HEAT_DISTRIBUTION[position]
        style = BLOCK_STYLES[i % len(BLOCK_STYLES)]

        ax.plot(
            time_ms,
            waveform,
            color=style["color"],
            linewidth=2.4,
            linestyle=style["linestyle"],
            alpha=0.95,
            label=f"block {position} ",
        )

    # 合計波形
    ax.plot(
        time_ms,
        total_plot,
        color=TOTAL_COLOR,
        linewidth=3.2,
        linestyle="-",
        label="Total",
        zorder=10,
    )

    ax.set_xlabel("Time [ms]", fontsize=25)
    ax.set_ylabel(f"Current [{CURRENT_UNIT}]", fontsize=25)
    ax.set_xlim(XMIN_MS, XMAX_MS)
    ax.tick_params(axis="both", labelsize=20)
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=20, frameon=True)

    fig.tight_layout()

    if SAVE_FIGURE:
        fig.savefig(OUTPUT_FIGURE, dpi=300, bbox_inches="tight",transparent=True)
        print(f"Saved figure: {OUTPUT_FIGURE}")

    plt.show()


if __name__ == "__main__":
    main()
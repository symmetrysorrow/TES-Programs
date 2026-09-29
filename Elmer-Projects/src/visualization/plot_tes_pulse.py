from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT = REPO_ROOT / "artifacts" / "series" / "tes_pulse_series.csv"
OUTPUT = REPO_ROOT / "artifacts" / "plots" / "tes_pulse_response.png"


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    data = np.genfromtxt(INPUT, delimiter=",", names=True)
    time_ms = data["time_s"] * 1e3
    temperature_mk = data["tes_temperature_K"] * 1e3
    current_ua = data["tes_current_A"] * 1e6

    peak_index = int(np.argmax(temperature_mk))
    minimum_current_index = int(np.argmin(current_ua))

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex="col")
    fig.suptitle("TES response to a 1332 keV absorber pulse", fontsize=14)

    series = (
        (temperature_mk, "TES temperature (mK)", "#c43c35"),
        (current_ua, "TES current (uA)", "#176b87"),
    )
    for row, (values, ylabel, color) in enumerate(series):
        for col in range(2):
            ax = axes[row, col]
            ax.plot(time_ms, values, color=color, linewidth=1.7)
            ax.set_ylabel(ylabel)
            ax.grid(True, color="#d9d9d9", linewidth=0.7)
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)

    axes[0, 0].set_title("Full response (0-100 ms)")
    axes[0, 1].set_title("Early response (0-2 ms)")
    axes[0, 1].set_xlim(0, 2)
    axes[1, 1].set_xlim(0, 2)
    axes[1, 0].set_xlabel("Time (ms)")
    axes[1, 1].set_xlabel("Time (ms)")

    axes[0, 1].scatter(
        time_ms[peak_index], temperature_mk[peak_index], color="#8f1d18", zorder=3
    )
    axes[0, 1].annotate(
        f"Peak {temperature_mk[peak_index]:.6f} mK\n"
        f"at {time_ms[peak_index]:.3f} ms",
        (time_ms[peak_index], temperature_mk[peak_index]),
        xytext=(0.62, 0.82),
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=9,
    )

    axes[1, 1].scatter(
        time_ms[minimum_current_index],
        current_ua[minimum_current_index],
        color="#104b60",
        zorder=3,
    )
    axes[1, 1].annotate(
        f"Minimum {current_ua[minimum_current_index]:.6f} uA\n"
        f"at {time_ms[minimum_current_index]:.3f} ms",
        (time_ms[minimum_current_index], current_ua[minimum_current_index]),
        xytext=(0.62, 0.12),
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=9,
    )

    fig.tight_layout()
    fig.savefig(OUTPUT, dpi=200, bbox_inches="tight")
    print(OUTPUT)


if __name__ == "__main__":
    main()

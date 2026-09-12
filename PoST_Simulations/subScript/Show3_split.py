from pathlib import Path

from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt

from show_data import load_position_resolution

jp_font = FontProperties(family="Yu Gothic", size=25)

BASE_DIR = Path(r"h:\hata")

DATASETS = {
    "662 keV": {
        "folder": BASE_DIR / "662_142_136_300split",
    },
    "1332 keV": {
        "folder": BASE_DIR / "1332_142_136_300split",
    },
}

CONDITIONS = {
    "MS": {
        "target": "Pulse_ms",
        "linestyle": "-",
        "marker": "o",
        "color": "#1f77b4",
    },
    "Noise": {
        "target": "Pulse_noise",
        "linestyle": "--",
        "marker": "^",
        "color": "#d62728",
    },
    "MS+Noise": {
        "target": "Pulse_ms_noise",
        "linestyle": "-.",
        "marker": "s",
        "color": "#2ca02c",
    },
}

DISPLAY_LABELS = {
    "MS": "多重散乱",
    "Noise": "ノイズ",
    "MS+Noise": "多重散乱 + ノイズ",
}


def plot_energy(energy_label: str, info: dict) -> None:
    folder = info["folder"]

    fig, ax = plt.subplots(figsize=(8, 8))
    plotted = False

    for condition, cond_info in CONDITIONS.items():
        try:
            positions, fwhm = load_position_resolution(
                folder,
                cond_info["target"],
            )
        except FileNotFoundError as error:
            print(f"Skipping missing result: {error}")
            continue

        ax.plot(
            positions,
            fwhm,
            label=DISPLAY_LABELS[condition],
            color=cond_info["color"],
            linestyle=cond_info["linestyle"],
            marker=cond_info["marker"],
            linewidth=2.0,
            markersize=5.5,
            markerfacecolor="white",
            markeredgewidth=1.4,
        )

        plotted = True

    title_map = {
        "662 keV": "(a) 662 keV",
        "1332 keV": "(b) 1332 keV",
    }

    ax.set_title(
        title_map.get(energy_label, energy_label),
        fontsize=50,
    )
    ax.set_xlabel("Position [mm]", fontsize=35)
    ax.set_ylabel("Position Resolution[mm]", fontsize=35)
    ax.set_ylim(0, 0.4)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=30,
    )

    ax.legend(
        ncol=1,
        frameon=True,
        prop=jp_font,
    )

    plt.tight_layout()

    if plotted:
        out_path = (
            BASE_DIR
            / f"position_resolution_{energy_label.replace(' ', '_')}_ms_noise_split.png"
        )

        plt.savefig(
            out_path,
            dpi=300,
            transparent=True,
        )
        plt.show()

        print(f"Saved: {out_path}")
    else:
        print(f"No data plotted for {energy_label}")


def main():
    for energy_label, info in DATASETS.items():
        if not info["folder"].exists():
            print(f"Skipping missing folder: {info['folder']}")
            continue

        plot_energy(energy_label, info)


if __name__ == "__main__":
    main()
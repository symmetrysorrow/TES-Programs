from pathlib import Path

import matplotlib.pyplot as plt

from show_data import load_position_resolution


BASE_DIR = Path(r"h:\hata")
DATASETS = {
    "662 keV": {
        "folder": BASE_DIR / "662_142_136_300split",
        "color": "#1f77b4",
    },
    "1332 keV": {
        "folder": BASE_DIR / "1332_142_136_300split",
        "color": "#d62728",
    },
}
CONDITIONS = {
    "MS": {"target": "Pulse_ms", "linestyle": "-", "marker": "o"},
    "Noise": {"target": "Pulse_noise", "linestyle": "--", "marker": "^"},
    "MS+Noise": {"target": "Pulse_ms_noise", "linestyle": "-.", "marker": "s"},
}
OUT_PATH = BASE_DIR / "position_resolution_662_1332_ms_noise_compare.png"


def main():
    fig, ax = plt.subplots(figsize=(11, 6))

    plotted = False

    for energy_label, info in DATASETS.items():
        folder = info["folder"]
        if not folder.exists():
            print(f"Skipping missing folder: {folder}")
            continue

        for condition, cond_info in CONDITIONS.items():
            try:
                positions, fwhm = load_position_resolution(folder, cond_info["target"])
            except FileNotFoundError as error:
                print(f"Skipping missing result: {error}")
                continue
            style = {
                "color": info["color"],
                "linestyle": cond_info["linestyle"],
                "marker": cond_info["marker"],
            }

            ax.plot(
                positions,
                fwhm,
                label=f"{energy_label} {condition}",
                linewidth=2.0,
                markersize=5.5,
                markerfacecolor="white",
                markeredgewidth=1.4,
                **style,
            )
            plotted = True

    ax.set_xlabel("Position [mm]", fontsize=20)
    ax.set_ylabel("Position Resolution", fontsize=20)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.legend(ncol=2, frameon=True, fontsize=16)
    plt.tight_layout()

    if plotted:
        plt.savefig(OUT_PATH, dpi=300)
        plt.show()
    else:
        print("No data plotted.")


if __name__ == "__main__":
    main()

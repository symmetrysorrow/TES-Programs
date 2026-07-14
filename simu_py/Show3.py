from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np


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
    "MS": {"file": "fwhms_Pulse_ms.txt", "linestyle": "-", "marker": "o"},
    "Noise": {"file": "fwhms_Pulse_noise.txt", "linestyle": "--", "marker": "^"},
    "MS+Noise": {"file": "fwhms_Pulse_ms_noise.txt", "linestyle": "-.", "marker": "s"},
}
OUT_PATH = BASE_DIR / "position_resolution_662_1332_ms_noise_compare.png"


def load_positions(folder: Path) -> np.ndarray:
    with open(folder / "input.json", "r", encoding="utf-8") as f:
        para = json.load(f)
    return np.asarray(para["position"], dtype=float)


def load_fwhm(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")
    return np.loadtxt(path, dtype=float)


def main():
    fig, ax = plt.subplots(figsize=(11, 6))

    plotted = False

    for energy_label, info in DATASETS.items():
        folder = info["folder"]
        if not folder.exists():
            print(f"Skipping missing folder: {folder}")
            continue

        positions = load_positions(folder)

        for condition, cond_info in CONDITIONS.items():
            fwhm_path = folder / cond_info["file"]
            if not fwhm_path.is_file():
                print(f"Skipping missing file: {fwhm_path}")
                continue

            fwhm = load_fwhm(fwhm_path)
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

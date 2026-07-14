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
DISPLAY_LABELS = {
    "MS": "multiple interactions",
    "Noise": "noise",
    "MS+Noise": "multiple interactions + noise",
}


def load_positions(folder: Path) -> np.ndarray:
    with open(folder / "input.json", "r", encoding="utf-8") as f:
        para = json.load(f)
    return np.asarray(para["position"], dtype=float)


def load_fwhm(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")
    return np.loadtxt(path, dtype=float)


def plot_energy(energy_label: str, info: dict) -> None:
    folder = info["folder"]
    positions = load_positions(folder)

    fig, ax = plt.subplots(figsize=(11, 6))
    plotted = False

    for condition, cond_info in CONDITIONS.items():
        fwhm_path = folder / cond_info["file"]
        if not fwhm_path.is_file():
            print(f"Skipping missing file: {fwhm_path}")
            continue

        fwhm = load_fwhm(fwhm_path)
        ax.plot(
            positions,
            fwhm,
            label=DISPLAY_LABELS[condition],
            color=info["color"],
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
    ax.set_title(title_map.get(energy_label, energy_label), fontsize=20)
    ax.set_xlabel("Position [mm]", fontsize=20)
    ax.set_ylabel("Position Resolution", fontsize=20)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.legend(ncol=2, frameon=True, fontsize=16)
    plt.tight_layout()

    if plotted:
        out_path = BASE_DIR / f"position_resolution_{energy_label.replace(' ', '_')}_ms_noise_split.png"
        plt.savefig(out_path, dpi=300)
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

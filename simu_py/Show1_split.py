from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(r"h:\hata")
ENERGY_FOLDERS = {
    "662 keV": BASE_DIR / "662_142_136_300split",
    "1332 keV": BASE_DIR / "1332_142_136_300split",
}
CSV_CANDIDATE = "ene_resos_Pulse_ms_noise.csv"


def load_dataset(folder: Path) -> pd.DataFrame:
    csv_path = folder / CSV_CANDIDATE
    if not csv_path.is_file():
        raise FileNotFoundError(f"No energy-resolution CSV found in {folder}: {CSV_CANDIDATE}")
    return pd.read_csv(csv_path, index_col=0)


def plot_energy(energy_label: str, folder: Path) -> None:
    df = load_dataset(folder)
    positions = df.index.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(11, 6))
    style_map = {
        "PS": {"color": "#1f77b4", "linestyle": "-", "marker": "o"},
        "ST": {"color": "#1f77b4", "linestyle": "--", "marker": "^"},
    }
    display_labels = {
        "PS": "Method 1",
        "ST": "Method 2",
    }
    column_map = {
        "PS": "Sum",
        "ST": "ST",
    }

    for metric in ("PS", "ST"):
        column = column_map[metric]
        if column not in df.columns:
            print(f"Skipping missing column '{column}' in {folder}")
            continue

        ax.plot(
            positions,
            df[column].to_numpy(dtype=float),
            label=display_labels[metric],
            linewidth=2.0,
            markersize=5.5,
            markerfacecolor="white",
            markeredgewidth=1.4,
            **style_map[metric],
        )

    title_map = {
        "662 keV": "(a) 662 keV, Method 1",
        "1332 keV": "(b) 1332 keV, Method 2",
    }
    ax.set_title(title_map.get(energy_label, energy_label), fontsize=20)
    ax.set_xlabel("Position[mm]", fontsize=20)
    ax.set_ylabel("Energy Resolution[keV]", fontsize=20)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.legend(ncol=2, frameon=True, fontsize=20)
    plt.tight_layout()

    out_path = BASE_DIR / f"energy_resolution_vs_position_{energy_label.replace(' ', '_').replace('keV', 'keV')}_split.png"
    plt.savefig(out_path, dpi=300)
    plt.show()
    print(f"Saved: {out_path}")


def main():
    for energy_label, folder in ENERGY_FOLDERS.items():
        if not folder.exists():
            print(f"Skipping missing folder: {folder}")
            continue
        plot_energy(energy_label, folder)


if __name__ == "__main__":
    main()

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(r"h:\hata")
ENERGY_FOLDERS = {
    "662 keV": BASE_DIR / "662_142_136_300split",
    "1332 keV": BASE_DIR / "1332_142_136_300split",
}
CSV_CANDIDATES = ("ene_resos_Pulse_ms_noise.csv",)
OUT_PATH = BASE_DIR / "energy_resolution_vs_position_662_1332_single.png"


def find_csv(folder: Path) -> Path:
    for name in CSV_CANDIDATES:
        path = folder / name
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"No energy-resolution CSV found in {folder}. "
        f"Checked: {', '.join(CSV_CANDIDATES)}"
    )


def load_dataset(folder: Path) -> pd.DataFrame:
    csv_path = find_csv(folder)
    return pd.read_csv(csv_path, index_col=0)


def main():
    fig, ax = plt.subplots(figsize=(11, 6))

    style_map = {
        ("662 keV", "Sum"): {"color": "#1f77b4", "linestyle": "-", "marker": "o"},
        ("662 keV", "ST"): {"color": "#1f77b4", "linestyle": "--", "marker": "^"},
        ("1332 keV", "Sum"): {"color": "#d62728", "linestyle": "-", "marker": "o"},
        ("1332 keV", "ST"): {"color": "#d62728", "linestyle": "--", "marker": "^"},
    }

    plotted = False

    for energy_label, folder in ENERGY_FOLDERS.items():
        if not folder.exists():
            print(f"Skipping missing folder: {folder}")
            continue

        df = load_dataset(folder)
        positions = df.index.to_numpy(dtype=float)

        for metric in ("Sum", "ST"):
            if metric not in df.columns:
                print(f"Skipping missing column '{metric}' in {folder}")
                continue

            style = style_map[(energy_label, metric)]
            ax.plot(
                positions,
                df[metric].to_numpy(dtype=float),
                label=f"{energy_label} {metric}",
                linewidth=2.0,
                markersize=5.5,
                markerfacecolor="white",
                markeredgewidth=1.4,
                **style,
            )
            plotted = True

    ax.set_xlabel("Position[mm]", fontsize=20)
    ax.set_ylabel("Energy Resolution[keV]", fontsize=20)
    #ax.set_title("Energy Resolution vs Position")
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.legend(ncol=2, frameon=True, fontsize=22)
    plt.tight_layout()

    if plotted:
        plt.savefig(OUT_PATH, dpi=300)
        plt.show()
    else:
        print("No data plotted.")


if __name__ == "__main__":
    main()

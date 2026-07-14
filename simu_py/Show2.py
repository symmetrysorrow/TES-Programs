from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(r"h:\hata\1332_142_136_300split")
CSV_FILES = {
    "MS": BASE_DIR / "ene_resos_Pulse_ms.csv",
    "noise": BASE_DIR / "ene_resos_Pulse_noise.csv",
    "ms_noise": BASE_DIR / "ene_resos_Pulse_ms_noise.csv",
}
OUT_PATH = Path(r"h:\hata\energy_resolution_1332_sum_st_ms_noise_compare.png")


def load_dataset(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing CSV: {path}")
    return pd.read_csv(path, index_col=0)


def main():
    fig, ax = plt.subplots(figsize=(11, 6))

    style_map = {
        "MS": {"color": "#1f77b4"},
        "noise": {"color": "#2ca02c"},
        "ms_noise": {"color": "#d62728"},
    }
    metric_style_map = {
        "Sum": {"linestyle": "-", "marker": "o"},
        "ST": {"linestyle": "--", "marker": "^"},
    }

    plotted = False

    for condition, csv_path in CSV_FILES.items():
        df = load_dataset(csv_path)
        positions = df.index.to_numpy(dtype=float)

        for metric in ("Sum", "ST"):
            if metric not in df.columns:
                print(f"Skipping missing column '{metric}' in {csv_path.name}")
                continue

            style = {}
            style.update(style_map[condition])
            style.update(metric_style_map[metric])

            ax.plot(
                positions,
                df[metric].to_numpy(dtype=float),
                label=f"{condition} {metric}",
                linewidth=2.0,
                markersize=5.5,
                markerfacecolor="white",
                markeredgewidth=1.4,
                **style,
            )
            plotted = True

    ax.set_xlabel("Position[mm]", fontsize=20)
    ax.set_ylabel("Energy Resolution[keV]", fontsize=20)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.legend(ncol=3, frameon=True, fontsize=16)
    plt.tight_layout()

    if plotted:
        plt.savefig(OUT_PATH, dpi=300)
        plt.show()
    else:
        print("No data plotted.")


if __name__ == "__main__":
    main()

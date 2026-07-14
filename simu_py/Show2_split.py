from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(r"h:\hata\1332_142_136_300split")
CSV_FILES = {
    "MS": BASE_DIR / "ene_resos_Pulse_ms.csv",
    "noise": BASE_DIR / "ene_resos_Pulse_noise.csv",
    "ms_noise": BASE_DIR / "ene_resos_Pulse_ms_noise.csv",
}
DISPLAY_LABELS = {
    "MS": "multiple interactions",
    "noise": "noise",
    "ms_noise": "multiple interactions + noise",
}


def load_dataset(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing CSV: {path}")
    return pd.read_csv(path, index_col=0)


def plot_metric(metric: str) -> None:
    fig, ax = plt.subplots(figsize=(11, 6))
    style_map = {
        "MS": {"color": "#1f77b4", "linestyle": "-", "marker": "o"},
        "noise": {"color": "#2ca02c", "linestyle": "--", "marker": "^"},
        "ms_noise": {"color": "#d62728", "linestyle": "-.", "marker": "s"},
    }
    metric_title_map = {
        "Sum": "(a) Method 1",
        "ST": "(b) Method 2",
    }

    plotted = False

    for condition, csv_path in CSV_FILES.items():
        df = load_dataset(csv_path)
        if metric not in df.columns:
            print(f"Skipping missing column '{metric}' in {csv_path.name}")
            continue

        positions = df.index.to_numpy(dtype=float)
        ax.plot(
            positions,
            df[metric].to_numpy(dtype=float),
            label=DISPLAY_LABELS[condition],
            linewidth=2.0,
            markersize=5.5,
            markerfacecolor="white",
            markeredgewidth=1.4,
            **style_map[condition],
        )
        plotted = True

    ax.set_title(metric_title_map.get(metric, metric), fontsize=20)
    ax.set_xlabel("Position[mm]", fontsize=20)
    ax.set_ylabel("Energy Resolution[keV]", fontsize=20)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.legend(ncol=2, frameon=True, fontsize=20)
    plt.tight_layout()

    out_path = BASE_DIR / f"energy_resolution_1332_{metric.lower()}_split.png"
    if plotted:
        plt.savefig(out_path, dpi=300)
        plt.show()
        print(f"Saved: {out_path}")
    else:
        print(f"No data plotted for {metric}.")


def main():
    for metric in ("Sum", "ST"):
        plot_metric(metric)


if __name__ == "__main__":
    main()

from pathlib import Path

import matplotlib.pyplot as plt
from show_data import load_energy_resolution


BASE_DIR = Path(r"h:\hata\new")
TARGETS = {
    "MS": "Pulse_ms",
    "noise": "Pulse_noise",
    "ms_noise": "Pulse_ms_noise",
}
DISPLAY_LABELS = {
    "MS": "multiple interactions",
    "noise": "noise",
    "ms_noise": "multiple interactions + noise",
}


def plot_metric(metric: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_box_aspect(1)
    style_map = {
        "MS": {"color": "#1f77b4", "linestyle": "-", "marker": "o"},
        "noise": {"color": "#2ca02c", "linestyle": "--", "marker": "^"},
        "ms_noise": {"color": "#d62728", "linestyle": "-.", "marker": "s"},
    }
    metric_title_map = {
        "Sum": "Sum",
        "ST": "ST",
    }

    plotted = False

    for condition, target in TARGETS.items():
        try:
            df = load_energy_resolution(BASE_DIR, target)
        except FileNotFoundError as error:
            print(f"Skipping missing result: {error}")
            continue
        if metric not in df.columns:
            print(f"Skipping missing column '{metric}' for {condition}")
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

    ax.set_title(metric_title_map.get(metric, metric), fontsize=35)
    ax.set_xlabel("Position[mm]", fontsize=30)
    ax.set_ylabel("Energy Resolution[keV]", fontsize=30)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=24)
    ax.legend(
        ncol=1,
        frameon=True,
        fontsize=16.5,
        markerscale=1.125,
        handlelength=2.25,
        labelspacing=0.525,
        borderpad=0.525,
    )
    plt.tight_layout()

    out_path = BASE_DIR / f"energy_resolution_1332_{metric.lower()}_split.png"
    if plotted:
        fig.savefig(out_path, dpi=300, transparent=True)
        plt.show()
        print(f"Saved: {out_path}")
    else:
        print(f"No data plotted for {metric}.")


def main():
    for metric in ("Sum", "ST"):
        plot_metric(metric)


if __name__ == "__main__":
    main()

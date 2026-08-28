from pathlib import Path

import matplotlib.pyplot as plt
import questionary
from show_data import load_energy_resolution


BASE_DIR = Path(r"h:\hata\new_noise_test")
TARGETS = {
    "MS": "Pulse_ms",
    "noise": "Pulse_noise",
    "ms_noise": "Pulse_ms_noise",
}
MODE_ALL = "MS・Noise・MS+Noiseを比較"
MODE_MS_NOISE_ONLY = "MS+Noiseのみ"


def ask_mode() -> str | None:
    """Choose whether to compare all conditions or plot MS+Noise alone."""
    return questionary.select(
        "表示モード:",
        choices=(MODE_ALL, MODE_MS_NOISE_ONLY),
        default=MODE_ALL,
    ).ask()


def ask_title() -> str | None:
    """Request the title displayed above the graph."""
    return questionary.text("グラフタイトル:", default="").ask()


def main():
    mode = ask_mode()
    if mode is None:
        return
    title = ask_title()
    if title is None:
        return

    targets = TARGETS if mode == MODE_ALL else {"ms_noise": TARGETS["ms_noise"]}
    out_path = BASE_DIR / (
        "energy_resolution_1332_sum_st_ms_noise_compare.png"
        if mode == MODE_ALL
        else "energy_resolution_1332_sum_st_ms_noise.png"
    )
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_box_aspect(1)

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

    for condition, target in targets.items():
        try:
            df = load_energy_resolution(BASE_DIR, target)
        except FileNotFoundError as error:
            print(f"Skipping missing result: {error}")
            continue
        positions = df.index.to_numpy(dtype=float)

        for metric in ("Sum", "ST"):
            if metric not in df.columns:
                print(f"Skipping missing column '{metric}' for {condition}")
                continue

            style = {}
            style.update(style_map[condition])
            style.update(metric_style_map[metric])
            label_metric = "Max" if metric == "Sum" else metric

            ax.plot(
                positions,
                df[metric].to_numpy(dtype=float),
                label=(
                    label_metric
                    if mode == MODE_MS_NOISE_ONLY
                    else f"{condition} {label_metric}"
                ),
                linewidth=2.0,
                markersize=5.5,
                markerfacecolor="white",
                markeredgewidth=1.4,
                **style,
            )
            plotted = True

    ax.set_title(title, fontsize=35)
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

    if plotted:
        fig.savefig(out_path, dpi=300, transparent=True)
        plt.show()
        print(f"Saved: {out_path}")
    else:
        print("No data plotted.")


if __name__ == "__main__":
    main()

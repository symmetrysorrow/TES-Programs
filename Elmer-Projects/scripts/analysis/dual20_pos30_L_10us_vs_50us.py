"""pos30 (inject x=-8mm) TES_L current: 10us vs 50us recovery-timestep comparison.

Both resolutions were confirmed converged against each other (docs/dual_tes_plan.md,
"Timestep convergence" record: 50us->10us RMS ~2nA over 21-60ms) -- this plot is a
visual sanity check of that agreement, not a resolution study.

Reads:
  results/case_dual20_pos30/tes_dual20_pos30_L_series.csv        (10us, current)
  <scratchpad>/pos30_50us/tes_dual20_pos30_L_series.csv           (50us snapshot)

Writes: generated/dual20_pos30_L_10us_vs_50us.png

Read-only. Usage: python scripts/analysis/dual20_pos30_L_10us_vs_50us.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT_50US = Path(
    r"C:\Users\symme\AppData\Local\Temp\claude\d--github-Elmer-Projects"
    r"\f7251838-64b7-4165-8a19-2a0c15ab4ef2\scratchpad\pos30_50us"
    r"\tes_dual20_pos30_L_series.csv"
)
CURRENT_10US = ROOT / "results" / "case_dual20_pos30" / "tes_dual20_pos30_L_series.csv"

COLUMNS = ["time_s", "tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W"]

PULSE_START_MS = 20.02

# dataviz chrome (light surface)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

COLOR_50US = "#86b6ef"   # light blue
COLOR_10US = "#0d366b"   # dark blue (finer)


def load(path: Path):
    import pandas as pd
    df = pd.read_csv(path, names=COLUMNS, header=0)
    df = df.sort_values("time_s").reset_index(drop=True)
    return df["time_s"].to_numpy() * 1e3, df["tes_current_A"].to_numpy() * 1e6


def main() -> None:
    if not SNAPSHOT_50US.exists():
        raise SystemExit(f"50us snapshot not found: {SNAPSHOT_50US}")
    if not CURRENT_10US.exists():
        raise SystemExit(f"10us series not found: {CURRENT_10US}")

    t50, i50 = load(SNAPSHOT_50US)
    t10, i10 = load(CURRENT_10US)

    fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.plot(t50, i50, color=COLOR_50US, lw=2.4, marker="o", ms=3.2, mew=0,
             label="50 µs recovery step", zorder=2)
    ax.plot(t10, i10, color=COLOR_10US, lw=1.2, marker="o", ms=1.6, mew=0,
             label="10 µs recovery step", zorder=3)

    ax.axvline(PULSE_START_MS, color=MUTED, ls=":", lw=1, zorder=1)
    ax.annotate("pulse", xy=(PULSE_START_MS, 0.98), xycoords=("data", "axes fraction"),
                xytext=(3, -2), textcoords="offset points", color=MUTED,
                fontsize=8, va="top")

    ax.set_xlim(19.8, 42)
    ax.set_xlabel("time [ms]", color=MUTED, fontsize=10)
    ax.set_ylabel("TES_L current [µA]", color=MUTED, fontsize=10)
    ax.set_title("TES_L current, pos30 (inject x = -8 mm): 10 µs vs 50 µs recovery step",
                 color=INK, fontsize=11, pad=10)

    ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9)

    leg = ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="lower right")

    fig.tight_layout()
    out = ROOT / "generated" / "dual20_pos30_L_10us_vs_50us.png"
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"plot written: {out}")


if __name__ == "__main__":
    main()

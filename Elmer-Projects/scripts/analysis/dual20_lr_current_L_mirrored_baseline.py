"""TES_L current (baseline-subtracted, dI) across the FULL injection range
x = -8..+8 mm, using the LOCAL-REFINE mesh and mirror symmetry (see
dual20_lr_current_L_mirrored.py for the full symmetry-mapping rationale).

Each curve is shown as dI = I(t) - I_baseline, where I_baseline is the mean
pre-pulse current (1-20ms window). This removes the L/R steady-state offset
(TES_L and TES_R are physically distinct circuits with a small residual
mesh-asymmetry-driven baseline gap -- ~0.41% for this local-refine mesh, see
docs/dual_tes_plan.md) so curves are compared on dip DEPTH/TIMING alone.

A case is only plotted if its series has data reaching at least
MIN_DURATION_MS -- guards against picking up a still-running case's partial,
mid-write CSV (skipped with a warning instead).

Read-only. Writes: generated/dual20_lr_current_L_mirrored_baseline.png
Usage: python scripts/analysis/dual20_lr_current_L_mirrored_baseline.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PULSE_START_MS = 20.02
MIN_DURATION_MS = 34.0  # guards against a still-running case's partial CSV
COLUMNS = ["time_s", "tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W"]

SOURCES: list[tuple[float, str, str]] = [
    (-8.0, "case_dual20_lr_pos30", "L"),
    (-6.0, "case_dual20_lr_pos60", "L"),
    (-4.0, "case_dual20_lr_pos90", "L"),
    (-2.0, "case_dual20_lr_pos120", "L"),
    (0.0, "case_dual20_lr_pos150", "L"),
    (2.0, "case_dual20_lr_pos120", "R"),
    (4.0, "case_dual20_lr_pos90", "R"),
    (6.0, "case_dual20_lr_pos60", "R"),
    (8.0, "case_dual20_lr_pos30", "R"),
]

DIV_BLUE = "#0d366b"
DIV_MID = "#f0efec"
DIV_RED = "#b23433"

SURFACE, INK, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#898781", "#e1e0d9", "#c3c2b7"


def find_series(case_dir: str, side: str) -> Path | None:
    for parent in (ROOT / "results" / case_dir, ROOT):
        for side_tag in (side, side.lower()):
            for c in parent.glob(f"*_{side_tag}_series.csv"):
                return c
    return None


def load(case_dir: str, side: str) -> pd.DataFrame | None:
    path = find_series(case_dir, side)
    if path is None:
        return None
    df = pd.read_csv(path, names=COLUMNS, header=0)
    return df.sort_values("time_s").reset_index(drop=True)


def main() -> None:
    cmap = LinearSegmentedColormap.from_list("dual20_pos_diverging", [DIV_BLUE, DIV_MID, DIV_RED])
    norm = Normalize(vmin=-8.0, vmax=8.0)

    fig, ax = plt.subplots(figsize=(8.6, 5.0), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    plotted = 0
    for x_mm, case_dir, side in SOURCES:
        df = load(case_dir, side)
        if df is None:
            print(f"WARNING: {case_dir}: TES_{side} series CSV not found -- skipping x={x_mm:+.0f}mm.")
            continue
        t_ms = df["time_s"].to_numpy() * 1e3
        if t_ms.max() < MIN_DURATION_MS:
            print(f"WARNING: {case_dir}: TES_{side} series only reaches {t_ms.max():.2f}ms "
                  f"(< {MIN_DURATION_MS}ms) -- likely still running, skipping x={x_mm:+.0f}mm.")
            continue
        cur_uA = df["tes_current_A"].to_numpy() * 1e6
        pre_mask = (t_ms > 1.0) & (t_ms < 20.0)
        baseline_uA = cur_uA[pre_mask].mean() if pre_mask.any() else cur_uA[0]
        d_cur_uA = cur_uA - baseline_uA
        mask = t_ms <= 35.0
        color = cmap(norm(x_mm))
        ax.plot(t_ms[mask], d_cur_uA[mask], color=color, lw=1.6, marker="o", ms=2.2,
                mew=0, zorder=3)
        plotted += 1

    if plotted == 0:
        print("WARNING: no data available -- not writing output.")
        plt.close(fig)
        return

    ax.axhline(0.0, color=AXIS, lw=0.8, zorder=1)
    ax.axvline(PULSE_START_MS, color=MUTED, ls=":", lw=1, zorder=1)
    ax.annotate("pulse", xy=(PULSE_START_MS, 0.98), xycoords=("data", "axes fraction"),
                xytext=(3, -2), textcoords="offset points", color=MUTED,
                fontsize=8, va="top")

    ax.set_xlim(19.8, 35)
    ax.set_xlabel("time [ms]", color=MUTED, fontsize=10)
    ax.set_ylabel("TES_L current change from baseline, dI [µA]", color=MUTED, fontsize=10)
    ax.set_title(f"TES_L baseline-corrected current (dI) via mirror symmetry\n"
                 f"(local-refine mesh, {plotted}/9 positions)",
                 color=INK, fontsize=11, pad=10)

    ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("injection position x [mm]", color=MUTED, fontsize=10)
    cbar.ax.tick_params(colors=MUTED, labelsize=9)
    cbar.outline.set_edgecolor(AXIS)

    fig.tight_layout()
    out = ROOT / "generated" / "dual20_lr_current_L_mirrored_baseline.png"
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"plot written: {out}  ({plotted}/9 positions)")


if __name__ == "__main__":
    main()

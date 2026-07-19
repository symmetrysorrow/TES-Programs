"""TES_L current across the FULL injection range x = -8..+8 mm, using the
LOCAL-REFINE mesh (mesh_dual_20mm_localrefine -- TES-column refinement only,
~65.6k nodes vs base's 56.8k / full-2x's 230k; see docs/dual_tes_plan.md
"TES local mesh refinement" record) and exploiting the mirror symmetry of
mesh_dual_20mm about x=0 (TES_L at x=-9.5mm, TES_R at x=+9.5mm): injecting at
+p mm and reading TES_L is physically equivalent to injecting at -p mm and
reading TES_R.

So L(+p) = R(-p). The five scan cases (x = -8,-6,-4,-2,0 mm) give L(x)
directly for x <= 0, and R(x) at the same five cases gives L(-x) for x >= 0:

    x [mm]   source
    -8       case_dual20_lr_pos30  TES_L series (direct)
    -6       case_dual20_lr_pos60  TES_L series (direct)
    -4       case_dual20_lr_pos90  TES_L series (direct)
    -2       case_dual20_lr_pos120 TES_L series (direct)
     0       case_dual20_lr_pos150 TES_L series (direct, L=R by symmetry at centre)
    +2       case_dual20_lr_pos120 TES_R series (mirror of x=-2)
    +4       case_dual20_lr_pos90  TES_R series (mirror of x=-4)
    +6       case_dual20_lr_pos60  TES_R series (mirror of x=-6)
    +8       case_dual20_lr_pos30  TES_R series (mirror of x=-8)

Cases whose series CSV isn't available yet (e.g. pos150 still running) are
skipped with a warning rather than raising, so this can be re-run as results
land. Position is encoded with the same diverging blue<->red colormap +
continuous colorbar as dual20_current_L_mirrored.py.

Read-only. Writes: generated/dual20_lr_current_L_mirrored.png
Usage: python scripts/analysis/dual20_lr_current_L_mirrored.py
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
COLUMNS = ["time_s", "tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W"]

# (display x [mm], case dir, side to read) -- side reflects the symmetry mapping above
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

# dataviz diverging pair: blue <-> red, neutral gray midpoint (palette.md)
DIV_BLUE = "#0d366b"    # dark blue pole (x = -8mm)
DIV_MID = "#f0efec"     # neutral gray midpoint (x = 0)
DIV_RED = "#b23433"     # dark red pole (x = +8mm), stepped from categorical red #e34948

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
            print(f"WARNING: {case_dir}: TES_{side} series CSV not found yet -- skipping x={x_mm:+.0f}mm.")
            continue
        t_ms = df["time_s"].to_numpy() * 1e3
        cur_uA = df["tes_current_A"].to_numpy() * 1e6
        mask = t_ms <= 35.0
        color = cmap(norm(x_mm))
        ax.plot(t_ms[mask], cur_uA[mask], color=color, lw=1.6, marker="o", ms=2.2,
                mew=0, zorder=3)
        plotted += 1

    if plotted == 0:
        print("WARNING: no data available -- not writing output.")
        plt.close(fig)
        return

    ax.axvline(PULSE_START_MS, color=MUTED, ls=":", lw=1, zorder=1)
    ax.annotate("pulse", xy=(PULSE_START_MS, 0.98), xycoords=("data", "axes fraction"),
                xytext=(3, -2), textcoords="offset points", color=MUTED,
                fontsize=8, va="top")

    ax.set_xlim(19.8, 35)
    ax.set_xlabel("time [ms]", color=MUTED, fontsize=10)
    ax.set_ylabel("TES_L current [µA]", color=MUTED, fontsize=10)
    ax.set_title(f"TES_L current across the injection range via mirror symmetry\n"
                 f"(local-refine mesh, {plotted}/9 positions available)",
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
    out = ROOT / "generated" / "dual20_lr_current_L_mirrored.png"
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"plot written: {out}  ({plotted}/9 positions)")


if __name__ == "__main__":
    main()

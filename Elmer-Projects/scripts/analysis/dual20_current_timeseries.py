"""20mm-absorber / 5-point-scan TES current TIME SERIES (docs/dual_tes_plan.md,
20mm dual-TES phase).

Produces TWO separate figures, one per readout TES:

  generated/dual20_current_L_timeseries.png   -- TES_L current vs time, 5 positions
  generated/dual20_current_R_timeseries.png   -- TES_R current vs time, 5 positions

Each figure overlays the five injection positions (pulse center x = -8/-6/-4/
-2/0 mm, cases case_dual20_pos30..pos150). Position is an *ordered* quantity, so
it is encoded with a single-hue blue ramp (light = far left end, dark = centre)
rather than categorical hues -- the progression is then readable at a glance and
is colour-vision-safe by lightness alone (dataviz skill: sequential/ordinal
encoding for magnitude).

CSV convention matches dual20_current_scan.py / dual_series_analysis.py:
columns time_s, tes_temperature_K, tes_current_A, tes_resistance_ohm,
tes_power_W. Cases whose series CSVs are not present yet are skipped with a
warning (the five pulse runs are long and finish at different times), so the
script can be run repeatedly as results arrive.

Read-only. Usage:
  python scripts/analysis/dual20_current_timeseries.py            # default window
  python scripts/analysis/dual20_current_timeseries.py --tmin 19.5 --tmax 60
  python scripts/analysis/dual20_current_timeseries.py --full     # whole series
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PULSE_START_MS = 20.02

COLUMNS = ["time_s", "tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W"]

# (case, series stem, pulse center x [mm], ordinal blue step light->dark).
# Steps come from the dataviz blue ramp; the lightest (250) clears the 2:1
# ordinal floor on the light surface, so no position is lost against the plane.
SCAN_POINTS: list[tuple[str, str, float, str]] = [
    ("case_dual20_pos30",  "tes_dual20_pos30",  -8.0, "#86b6ef"),
    ("case_dual20_pos60",  "tes_dual20_pos60",  -6.0, "#5598e7"),
    ("case_dual20_pos90",  "tes_dual20_pos90",  -4.0, "#2a78d6"),
    ("case_dual20_pos120", "tes_dual20_pos120", -2.0, "#1c5cab"),
    ("case_dual20_pos150", "tes_dual20_pos150",  0.0, "#0d366b"),
]

# dataviz chrome (light surface)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

TES_X = {"L": -9.5, "R": 9.5}


def find_series(case_dir: str, csv_base: str, side: str) -> Path | None:
    """Locate a case's per-side series CSV, or None if not run yet. Checks
    results/<case>/ and the repo root, both 'L'/'R' and 'l'/'r' casing."""
    for parent in (ROOT / "results" / case_dir, ROOT):
        for side_tag in (side, side.lower()):
            c = parent / f"{csv_base}_{side_tag}_series.csv"
            if c.exists():
                return c
    return None


def load(case_dir: str, csv_base: str, side: str) -> pd.DataFrame | None:
    path = find_series(case_dir, csv_base, side)
    if path is None:
        return None
    df = pd.read_csv(path, names=COLUMNS, header=0)
    return df.sort_values("time_s").reset_index(drop=True)


def plot_side(side: str, tmin: float, tmax: float | None, out: Path) -> Path | None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    plotted = 0
    for case_dir, csv_base, x_mm, color in SCAN_POINTS:
        df = load(case_dir, csv_base, side)
        if df is None:
            print(f"WARNING: {case_dir}: TES_{side} series CSV not found -- skipping.")
            continue
        t_ms = df["time_s"].to_numpy() * 1e3
        cur_uA = df["tes_current_A"].to_numpy() * 1e6
        label = f"x = {x_mm:+.0f} mm" + (" (centre)" if x_mm == 0 else "")
        ax.plot(t_ms, cur_uA, color=color, lw=1.8, marker="o", ms=2.6,
                mew=0, label=label, zorder=3)
        plotted += 1

    if plotted == 0:
        print(f"WARNING: no TES_{side} data available -- not writing {out}")
        plt.close(fig)
        return None

    ax.axvline(PULSE_START_MS, color=MUTED, ls=":", lw=1, zorder=1)
    ax.annotate("pulse", xy=(PULSE_START_MS, 0.98), xycoords=("data", "axes fraction"),
                xytext=(3, -2), textcoords="offset points", color=MUTED,
                fontsize=8, va="top")

    if tmax is not None:
        ax.set_xlim(tmin, tmax)
    else:
        ax.set_xlim(left=tmin)

    ax.set_xlabel("time [ms]", color=MUTED, fontsize=10)
    ax.set_ylabel("TES current [µA]", color=MUTED, fontsize=10)
    ax.set_title(f"TES_{side} current (readout at x = {TES_X[side]:+.1f} mm) "
                 f"— injection-position scan, mesh_dual_20mm",
                 color=INK, fontsize=11, pad=10)

    ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9)

    leg = ax.legend(frameon=False, fontsize=9, title="injection position",
                    labelcolor=INK, loc="best")
    leg.get_title().set_color(MUTED)
    leg.get_title().set_fontsize(9)

    fig.tight_layout()
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"plot written: {out}  ({plotted}/5 positions)")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tmin", type=float, default=19.5, help="x-axis min [ms] (default 19.5)")
    ap.add_argument("--tmax", type=float, default=60.0, help="x-axis max [ms] (default 60)")
    ap.add_argument("--full", action="store_true", help="show the whole series (ignore --tmax)")
    args = ap.parse_args()
    tmax = None if args.full else args.tmax

    outdir = ROOT / "generated"
    plot_side("L", args.tmin, tmax, outdir / "dual20_current_L_timeseries.png")
    plot_side("R", args.tmin, tmax, outdir / "dual20_current_R_timeseries.png")


if __name__ == "__main__":
    main()

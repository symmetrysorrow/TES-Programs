"""20mm-absorber / 5-point-scan current comparison (docs/dual_tes_plan.md,
20mm dual-TES phase).

For each of the five position-scan pulse cases (case_dual20_pos30 ..
case_dual20_pos150, pulse center x = -8/-6/-4/-2/0 mm respectively) this
reads the per-side series CSVs (columns: time_s, tes_temperature_K,
tes_current_A, tes_resistance_ohm, tes_power_W; same convention as
dual_series_analysis.py) and reports, for TES_L and TES_R:

  (a) the peak current dip after the pulse (max drop from the pre-pulse
      baseline current) and the time at which it occurs,
  (b) a plot of L/R dip depth vs. pulse position (generated/dual20_current_scan.png),
  (c) a numeric table on stdout.

Only case_dual20_steady + the five pos* pulse cases need to exist for this
to run; at the time this script was written only case_dual20_steady had been
executed (the five pulse cases are long-running and are launched separately
by the operator). Cases whose series CSVs are not found yet are skipped with
a warning rather than raising, so the script's I/O plumbing can be exercised
before all runs are complete: run it again once more cases have finished.

Read-only: does not touch the SIF/results pipeline.

Usage: python scripts/analysis/dual20_current_scan.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PULSE_START_S = 20.02e-3

COLOR_L = "#2a78d6"  # dataviz categorical slot 1 (blue)
COLOR_R = "#1baf7a"  # dataviz categorical slot 2 (aqua)

COLUMNS = ["time_s", "tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W"]

# (case name, series_file stem from elmer_project.json, pulse center x [mm])
SCAN_POINTS: list[tuple[str, str, float]] = [
    ("case_dual20_pos30", "tes_dual20_pos30", -8.0),
    ("case_dual20_pos60", "tes_dual20_pos60", -6.0),
    ("case_dual20_pos90", "tes_dual20_pos90", -4.0),
    ("case_dual20_pos120", "tes_dual20_pos120", -2.0),
    ("case_dual20_pos150", "tes_dual20_pos150", 0.0),
]


def find_series(case_dir: str, csv_base: str, side: str) -> Path | None:
    """Locate a case's per-side series CSV, or None if the case has not run
    yet. Checks both results/<case>/ and the repo root, and both 'L'/'R' and
    'l'/'r' casing (NTFS case-insensitive lookup; see dual_series_analysis.py)."""
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


def current_dip(df: pd.DataFrame) -> dict:
    """Peak post-pulse current dip relative to the pre-pulse baseline."""
    t = df["time_s"].to_numpy()
    cur = df["tes_current_A"].to_numpy()
    pre = t < PULSE_START_S
    post = t >= PULSE_START_S
    if not post.any():
        raise RuntimeError("no samples at/after pulse start")
    i_baseline = cur[pre][-1] if pre.any() else cur[0]
    idx_dip = int(np.argmin(cur[post]))
    i_dip = float(cur[post][idx_dip])
    t_dip = float(t[post][idx_dip])
    return {
        "i_baseline_A": float(i_baseline),
        "i_dip_A": i_dip,
        "dip_depth_A": i_baseline - i_dip,
        "dip_depth_pct": (i_dip - i_baseline) / i_baseline * 100.0,
        "t_dip_ms": (t_dip - PULSE_START_S) * 1e3,
    }


def collect() -> pd.DataFrame:
    rows = []
    for case_dir, csv_base, x_mm in SCAN_POINTS:
        df_l = load(case_dir, csv_base, "L")
        df_r = load(case_dir, csv_base, "R")
        if df_l is None or df_r is None:
            missing = [s for s, d in (("L", df_l), ("R", df_r)) if d is None]
            print(f"WARNING: {case_dir}: series CSV missing for side(s) {missing} "
                  f"(expected {csv_base}_<L|R>_series.csv under results/{case_dir}/ "
                  "or repo root) -- case not run yet, skipping.")
            continue
        dip_l = current_dip(df_l)
        dip_r = current_dip(df_r)
        rows.append({
            "case": case_dir,
            "pulse_x_mm": x_mm,
            "L_i_baseline_uA": dip_l["i_baseline_A"] * 1e6,
            "L_dip_depth_uA": dip_l["dip_depth_A"] * 1e6,
            "L_dip_depth_pct": dip_l["dip_depth_pct"],
            "L_t_dip_ms": dip_l["t_dip_ms"],
            "R_i_baseline_uA": dip_r["i_baseline_A"] * 1e6,
            "R_dip_depth_uA": dip_r["dip_depth_A"] * 1e6,
            "R_dip_depth_pct": dip_r["dip_depth_pct"],
            "R_t_dip_ms": dip_r["t_dip_ms"],
        })
    return pd.DataFrame(rows)


def report(df: pd.DataFrame) -> None:
    print("=" * 100)
    print("20mm dual-TES position scan: peak current dip vs. pulse position")
    print("=" * 100)
    if df.empty:
        print("(no cases have run yet -- nothing to report)")
        return
    pd.set_option("display.width", 200)
    pd.set_option("display.float_format", lambda v: f"{v:.4f}")
    print(df.to_string(index=False))


def plot(df: pd.DataFrame, out: Path) -> Path | None:
    if df.empty:
        print(f"WARNING: no data available -- not writing {out}")
        return None
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    ax.plot(df["pulse_x_mm"], df["L_dip_depth_uA"], marker="o", color=COLOR_L, lw=2, label="TES L dip")
    ax.plot(df["pulse_x_mm"], df["R_dip_depth_uA"], marker="o", color=COLOR_R, lw=2, label="TES R dip")
    ax.axvline(-9.5, color="0.6", ls=":", lw=1)
    ax.axvline(9.5, color="0.6", ls=":", lw=1)
    ax.set_xlabel("pulse center x [mm]  (TES_L at -9.5mm, TES_R at +9.5mm)")
    ax.set_ylabel("peak current dip [uA]  (baseline - post-pulse minimum)")
    ax.set_title("mesh_dual_20mm: L/R current-dip response vs. pulse position")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> None:
    df = collect()
    report(df)
    out_path = plot(df, ROOT / "generated" / "dual20_current_scan.png")
    if out_path is not None:
        print(f"\nplot written: {out_path}")


if __name__ == "__main__":
    main()

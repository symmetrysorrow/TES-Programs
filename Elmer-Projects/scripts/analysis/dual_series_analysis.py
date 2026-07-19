"""Phase D2 physical validation for the dual-TES cases (docs/dual_tes_plan.md
section 5, gates D-3.2..D-3.4).

Reads the six per-side series CSVs (columns: time_s, tes_temperature_K,
tes_current_A, tes_resistance_ohm, tes_power_W) for:
  - case_dual_shunt_transient (baseline flatness, no pulse)
  - case_dual_pulse_center     (x=0 pulse, L/R symmetry check)
  - case_dual_pulse_offset     (x=+1.5mm pulse, L/R asymmetry check)

and prints the gate-2..4 metrics, plus writes three comparison plots to
generated/. Read-only: does not touch the SIF/results pipeline.

Usage: python scripts/analysis/dual_series_analysis.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

# Steady operating point (docs/dual_tes_plan.md Phase D1 record).
T_STEADY_L = 0.1672581
T_STEADY_R = 0.1673430

PULSE_START_S = 20.02e-3

COLOR_L = "#2a78d6"  # dataviz categorical slot 1 (blue)
COLOR_R = "#1baf7a"  # dataviz categorical slot 2 (aqua)

COLUMNS = ["time_s", "tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W"]


def find_series(case_dir: str, csv_base: str, side: str) -> Path:
    """Locate a case's per-side series CSV. *csv_base* is the series_file
    stem from elmer_project.json (e.g. 'tes_dual_shunt_transient'; note it
    is NOT the same as the case name, which is 'case_dual_shunt_transient').
    Actual on-disk filenames may carry either 'L'/'R' or 'l'/'r' (Windows/
    NTFS is case-preserving on create but case-insensitive on lookup, and
    STATUS='REPLACE' in the UDF reuses whichever case an earlier run
    created), so search both, in both the repo root and results/<case>/."""
    for parent in (ROOT / "results" / case_dir, ROOT):
        for side_tag in (side, side.lower()):
            c = parent / f"{csv_base}_{side_tag}_series.csv"
            if c.exists():
                return c
    raise FileNotFoundError(f"no series CSV for case_dir={case_dir!r} csv_base={csv_base!r} side={side!r}")


def load(case_dir: str, csv_base: str, side: str) -> pd.DataFrame:
    path = find_series(case_dir, csv_base, side)
    df = pd.read_csv(path, names=COLUMNS, header=0)
    return df.sort_values("time_s").reset_index(drop=True)


def report_baseline(df_l: pd.DataFrame, df_r: pd.DataFrame) -> None:
    print("=" * 70)
    print("GATE 2: case_dual_shunt_transient baseline flatness")
    print("=" * 70)
    settle_s = 0.01  # first ~10ms is the restart-from-steady relaxation transient
    for label, df, t_ref in (("L", df_l, T_STEADY_L), ("R", df_r, T_STEADY_R)):
        t = df["tes_temperature_K"]
        i = df["tes_current_A"]
        r = df["tes_resistance_ohm"]
        span_uk = (t.max() - t.min()) * 1e6
        settled = df[df["time_s"] >= settle_s]
        settled_span_uk = (settled["tes_temperature_K"].max() - settled["tes_temperature_K"].min()) * 1e6
        drift_uk = (t.mean() - t_ref) * 1e6
        i_drift_pct = (i.mean() - i.iloc[0]) / i.iloc[0] * 100.0
        r_drift_pct = (r.mean() - r.iloc[0]) / r.iloc[0] * 100.0
        verdict = "PASS" if settled_span_uk < 10.0 else "FAIL"
        print(
            f"[{label}] T span (full) = {span_uk:.3f} uK  |  T span (t>={settle_s * 1e3:.0f}ms, "
            f"post restart-transient) = {settled_span_uk:.3f} uK  (mean-steady = {drift_uk:+.3f} uK)  "
            f"I drift = {i_drift_pct:+.4f}%  R drift = {r_drift_pct:+.4f}%   -> {verdict} (<10uK gate on settled span)"
        )


def rise_and_peak(df: pd.DataFrame, t_ref: float) -> dict:
    t = df["time_s"].to_numpy()
    temp = df["tes_temperature_K"].to_numpy()
    cur = df["tes_current_A"].to_numpy()
    post = t >= PULSE_START_S
    if not post.any():
        raise RuntimeError("no samples at/after pulse start")
    idx_peak = np.argmax(temp[post])
    t_peak = t[post][idx_peak]
    t_peak_max = temp[post][idx_peak]
    dT_peak = t_peak_max - t_ref
    rise_time = t_peak - PULSE_START_S
    pre = t < PULSE_START_S
    i_baseline = cur[pre][-1] if pre.any() else cur[0]
    idx_dip = np.argmin(cur[post])
    i_dip = cur[post][idx_dip]
    dip_depth_pct = (i_dip - i_baseline) / i_baseline * 100.0
    return {
        "dT_peak_K": dT_peak,
        "t_peak_s": t_peak,
        "rise_time_s": rise_time,
        "i_baseline_A": i_baseline,
        "i_dip_A": i_dip,
        "dip_depth_pct": dip_depth_pct,
    }


def report_pulse_center(df_l: pd.DataFrame, df_r: pd.DataFrame) -> None:
    print("=" * 70)
    print("GATE 3: case_dual_pulse_center L/R symmetry")
    print("=" * 70)
    n = min(len(df_l), len(df_r))
    tl, tr = df_l["tes_temperature_K"].to_numpy()[:n], df_r["tes_temperature_K"].to_numpy()[:n]
    il, ir = df_l["tes_current_A"].to_numpy()[:n], df_r["tes_current_A"].to_numpy()[:n]
    rl, rr = df_l["tes_resistance_ohm"].to_numpy()[:n], df_r["tes_resistance_ohm"].to_numpy()[:n]
    pl, pr = df_l["tes_power_W"].to_numpy()[:n], df_r["tes_power_W"].to_numpy()[:n]

    def max_rel_diff(a, b):
        denom = np.maximum(np.abs(a), np.abs(b))
        denom = np.where(denom == 0, 1.0, denom)
        return float(np.max(np.abs(a - b) / denom))

    t_reldiff = max_rel_diff(tl, tr)
    i_reldiff = max_rel_diff(il, ir)
    r_reldiff = max_rel_diff(rl, rr)
    p_reldiff = max_rel_diff(pl, pr)
    t_verdict = "PASS" if t_reldiff < 1e-3 else "FAIL"
    print(f"max relative diff (all samples): T={t_reldiff:.3e} ({t_verdict}, <1e-3 gate)  "
          f"I={i_reldiff:.3e}  R={r_reldiff:.3e}  P={p_reldiff:.3e}  (I/R/P: % order expected)")

    # Peak response is measured after the pulse (t >= PULSE_START_S); the
    # global max/min over the whole series is dominated by the restart
    # relaxation transient in the first couple of ms (see report_baseline),
    # not the pulse response.
    stats_l = rise_and_peak(df_l, T_STEADY_L)
    stats_r = rise_and_peak(df_r, T_STEADY_R)
    print(f"pulse response (post t>={PULSE_START_S * 1e3:.2f}ms): "
          f"L dT_peak={stats_l['dT_peak_K'] * 1e6:.3f} uK @ t_peak={stats_l['t_peak_s'] * 1e3:.4f}ms  "
          f"R dT_peak={stats_r['dT_peak_K'] * 1e6:.3f} uK @ t_peak={stats_r['t_peak_s'] * 1e3:.4f}ms")
    print(f"current dip (post-pulse, vs pre-pulse baseline): "
          f"L={stats_l['dip_depth_pct']:+.4f}%  R={stats_r['dip_depth_pct']:+.4f}%")


def report_pulse_offset(df_l: pd.DataFrame, df_r: pd.DataFrame) -> None:
    print("=" * 70)
    print("GATE 4: case_dual_pulse_offset asymmetry (R side is 0.5mm from")
    print("        the impact point, L side is 3.5mm)")
    print("=" * 70)
    stats_l = rise_and_peak(df_l, T_STEADY_L)
    stats_r = rise_and_peak(df_r, T_STEADY_R)
    for label, s in (("L", stats_l), ("R", stats_r)):
        print(
            f"[{label}] dT_peak={s['dT_peak_K'] * 1e6:.3f} uK  t_peak={s['t_peak_s'] * 1e3:.6f} ms  "
            f"rise_time={s['rise_time_s'] * 1e3:.6f} ms  "
            f"I_dip depth={s['dip_depth_pct']:+.4f}%"
        )
    amp_ratio = stats_r["dT_peak_K"] / stats_l["dT_peak_K"] if stats_l["dT_peak_K"] else float("nan")
    timing_diff_ms = (stats_l["t_peak_s"] - stats_r["t_peak_s"]) * 1e3
    verdict = "PASS" if (stats_r["dT_peak_K"] > stats_l["dT_peak_K"] and stats_r["t_peak_s"] <= stats_l["t_peak_s"]) else "FAIL"
    print(f"amplitude ratio dT_R/dT_L = {amp_ratio:.4f}")
    print(f"timing: R peaks {timing_diff_ms:+.6f} ms before L (positive = R leads)")
    print(f"qualitative gate (R leads & larger): {verdict}")


def plot_baseline(df_l: pd.DataFrame, df_r: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150)
    for ax in axes:
        ax.plot(df_l["time_s"] * 1e3, (df_l["tes_temperature_K"] - T_STEADY_L) * 1e6,
                color=COLOR_L, lw=2, label="TES L")
        ax.plot(df_r["time_s"] * 1e3, (df_r["tes_temperature_K"] - T_STEADY_R) * 1e6,
                color=COLOR_R, lw=2, label="TES R")
        ax.set_xlabel("time [ms]")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("T - T_steady [uK]")
    axes[0].set_title("full series (restart-from-steady transient visible <10ms)")
    axes[0].legend(frameon=False)
    axes[1].set_xlim(10, 99)
    settled = pd.concat([
        df_l.loc[df_l["time_s"] >= 0.01, "tes_temperature_K"] - T_STEADY_L,
        df_r.loc[df_r["time_s"] >= 0.01, "tes_temperature_K"] - T_STEADY_R,
    ]) * 1e6
    pad = max(1.0, (settled.max() - settled.min()) * 0.2)
    axes[1].set_ylim(settled.min() - pad, settled.max() + pad)
    axes[1].set_title("settled window (t >= 10ms)")
    fig.suptitle("case_dual_shunt_transient: baseline flatness (no pulse)")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_pulse(df_l: pd.DataFrame, df_r: pd.DataFrame, out: Path, title: str,
                overlay_l: pd.DataFrame | None = None, overlay_r: pd.DataFrame | None = None,
                overlay_label: str = "") -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150, sharey=False)
    t_ms_l = (df_l["time_s"] - PULSE_START_S) * 1e3
    t_ms_r = (df_r["time_s"] - PULSE_START_S) * 1e3
    axes[0].plot(t_ms_l, (df_l["tes_temperature_K"] - T_STEADY_L) * 1e6, color=COLOR_L, lw=2, label="TES L")
    axes[0].plot(t_ms_r, (df_r["tes_temperature_K"] - T_STEADY_R) * 1e6, color=COLOR_R, lw=2, label="TES R")
    axes[1].plot(t_ms_l, df_l["tes_current_A"] * 1e6, color=COLOR_L, lw=2, label="TES L")
    axes[1].plot(t_ms_r, df_r["tes_current_A"] * 1e6, color=COLOR_R, lw=2, label="TES R")
    if overlay_l is not None and overlay_r is not None:
        ot_l = (overlay_l["time_s"] - PULSE_START_S) * 1e3
        ot_r = (overlay_r["time_s"] - PULSE_START_S) * 1e3
        axes[0].plot(ot_l, (overlay_l["tes_temperature_K"] - T_STEADY_L) * 1e6,
                     color=COLOR_L, lw=1, ls="--", alpha=0.6, label=f"TES L ({overlay_label})")
        axes[0].plot(ot_r, (overlay_r["tes_temperature_K"] - T_STEADY_R) * 1e6,
                     color=COLOR_R, lw=1, ls="--", alpha=0.6, label=f"TES R ({overlay_label})")
        axes[1].plot(ot_l, overlay_l["tes_current_A"] * 1e6, color=COLOR_L, lw=1, ls="--", alpha=0.6)
        axes[1].plot(ot_r, overlay_r["tes_current_A"] * 1e6, color=COLOR_R, lw=1, ls="--", alpha=0.6)
    axes[0].set_xlabel("time - pulse start [ms]")
    axes[0].set_ylabel("T - T_steady [uK]")
    axes[0].set_xlim(-0.5, 5)
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel("time - pulse start [ms]")
    axes[1].set_ylabel("I_TES [uA]")
    axes[1].set_xlim(-0.5, 5)
    axes[1].grid(alpha=0.25)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main() -> None:
    df_shunt_l = load("case_dual_shunt_transient", "tes_dual_shunt_transient", "L")
    df_shunt_r = load("case_dual_shunt_transient", "tes_dual_shunt_transient", "R")
    df_center_l = load("case_dual_pulse_center", "tes_dual_pulse_center", "L")
    df_center_r = load("case_dual_pulse_center", "tes_dual_pulse_center", "R")
    df_offset_l = load("case_dual_pulse_offset", "tes_dual_pulse_offset", "L")
    df_offset_r = load("case_dual_pulse_offset", "tes_dual_pulse_offset", "R")

    report_baseline(df_shunt_l, df_shunt_r)
    report_pulse_center(df_center_l, df_center_r)
    report_pulse_offset(df_offset_l, df_offset_r)

    out_dir = ROOT / "generated"
    plot_baseline(df_shunt_l, df_shunt_r, out_dir / "dual_baseline_series.png")
    plot_pulse(df_center_l, df_center_r, out_dir / "dual_pulse_center_LR.png",
               "case_dual_pulse_center: L/R response (x=0, symmetric)")
    plot_pulse(df_offset_l, df_offset_r, out_dir / "dual_pulse_offset_LR.png",
               "case_dual_pulse_offset: L/R response (x=+1.5mm, R side closer)",
               overlay_l=df_center_l, overlay_r=df_center_r, overlay_label="center")
    print()
    print("plots written:")
    for name in ("dual_baseline_series.png", "dual_pulse_center_LR.png", "dual_pulse_offset_LR.png"):
        print(f"  {out_dir / name}")


if __name__ == "__main__":
    main()

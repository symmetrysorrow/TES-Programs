"""Plot the currently available part of the checkpointed Elmer run vs COMSOL."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
COMSOL = ROOT / "docs" / "Single-Pixel.txt"
ELMER = ROOT / "tes_pulse_20ms_3x_comsol_grid_series.csv"
OUT = ROOT / "artifacts" / "comparison" / "single_pixel_current" / "partial_current_timeseries_comparison.png"
PULSE_MS = 20.02


def main() -> None:
    c = np.loadtxt(COMSOL, comments="%", encoding="utf-8")
    e = np.genfromtxt(ELMER, delimiter=",", names=True)
    t_e = e["time_s"] * 1e3
    i_e = e["tes_current_A"] * 1e6
    # Use the last pre-pulse Elmer sample, so the comparison remains meaningful
    # while the pulse response is incomplete.
    baseline = float(i_e[t_e <= PULSE_MS][-1])
    end = float(t_e[-1])
    start = max(PULSE_MS - 0.002, float(t_e.min()))
    c_mask = (c[:, 0] >= start) & (c[:, 0] <= end)
    e_mask = (t_e >= start) & (t_e <= end)

    fig, axes = plt.subplots(2, 1, figsize=(8, 6.3), sharex=True, layout="constrained")
    for ax in axes:
        ax.axvline(PULSE_MS, color="0.45", ls=":", lw=1, label="pulse start")
        ax.grid(color="0.88")
    axes[0].plot(c[c_mask, 0], c[c_mask, 4], lw=1.7, color="#2468a2", label="COMSOL")
    axes[0].plot(t_e[e_mask], i_e[e_mask], lw=1.2, color="#d86438", label="Elmer (available)")
    axes[0].set_ylabel("TES current [µA]")
    axes[0].legend(frameon=False)
    axes[0].set_title(f"Partial comparison through {end:.9f} ms")
    axes[1].plot(c[c_mask, 0], c[c_mask, 4] - c[c_mask, 4][0], lw=1.7, color="#2468a2", label="COMSOL ΔI")
    axes[1].plot(t_e[e_mask], i_e[e_mask] - baseline, lw=1.2, color="#d86438", label="Elmer ΔI")
    axes[1].set_xlabel("time [ms]")
    axes[1].set_ylabel("current change [µA]")
    axes[1].legend(frameon=False)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(f"Wrote {OUT}; Elmer samples: {len(e)}, end={end:.12f} ms")


if __name__ == "__main__":
    main()

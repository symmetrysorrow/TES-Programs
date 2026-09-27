"""Compare live (possibly partial) MUMPS/HYPRE series with COMSOL.

Prints baseline currents and the baseline-subtracted current drop at common
post-pulse times.  Used while long transient runs are still in progress.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.analysis import compare_singlepixel_amgx_comsol as cmp  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("series", nargs="+", type=Path, help="label=path")
    parser.add_argument("--historical", type=Path, default=None)
    args = parser.parse_args()
    comsol = cmp.read_comsol(ROOT / "docs/Single-Pixel.txt")
    runs = {"COMSOL": comsol}
    for item in args.series:
        label, _, path = str(item).partition("=")
        runs[label] = cmp.read_elmer(Path(path))
    end = min(float(s.time_us[-1]) for s in runs.values())
    print(f"common end: {end:.3f} us")
    print("baseline [uA]: " + ", ".join(f"{k}={s.baseline_uA:.6f}" for k, s in runs.items()))
    times = [t for t in (0.5, 1, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100) if t <= end]
    print("t_us  " + "  ".join(f"{k:>12s}" for k in runs))
    for t in times:
        print(f"{t:5g} " + "  ".join(f"{float(np.interp(t, s.time_us, s.drop_uA)):12.6f}" for s in runs.values()))
    grid = np.linspace(0.0, end, max(2, int(end / 0.05) + 1))
    base = np.interp(grid, comsol.time_us, comsol.drop_uA)
    for k, s in runs.items():
        if k == "COMSOL":
            continue
        d = np.interp(grid, s.time_us, s.drop_uA) - base
        print(f"{k}: max|diff| vs COMSOL = {np.max(np.abs(d)):.6f} uA, RMSE = {np.sqrt(np.mean(d * d)):.6f} uA")
    labels = [k for k in runs if k != "COMSOL"]
    if len(labels) >= 2:
        a, b = runs[labels[0]], runs[labels[1]]
        d = np.interp(grid, a.time_us, a.drop_uA) - np.interp(grid, b.time_us, b.drop_uA)
        print(f"{labels[0]} - {labels[1]}: max|diff| = {np.max(np.abs(d)):.3e} uA; baseline diff = {a.baseline_uA - b.baseline_uA:.3e} uA")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

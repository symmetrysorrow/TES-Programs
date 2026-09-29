"""Plot TES current against simulation time from a TES series CSV."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--time-origin-ms", type=float, default=0.0)
    args = parser.parse_args()

    data = np.genfromtxt(args.input, delimiter=",", names=True)
    required = {"time_s", "tes_current_A"}
    missing = required.difference(data.dtype.names or ())
    if missing:
        raise SystemExit(f"missing columns in {args.input}: {', '.join(sorted(missing))}")

    time_ms = np.asarray(data["time_s"], dtype=float) * 1.0e3
    time_ms -= args.time_origin_ms
    current_uA = np.asarray(data["tes_current_A"], dtype=float) * 1.0e6

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.plot(time_ms, current_uA, color="#1565c0", linewidth=1.4)
    xlabel = "Time since pulse (ms)" if args.time_origin_ms else "Simulation time (ms)"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("TES current (µA)")
    ax.set_title("TES current response: 1 ms tail, 5 µs timestep")
    ax.grid(True, alpha=0.25)
    ax.set_xlim(float(np.nanmin(time_ms)), float(np.nanmax(time_ms)))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)

    print(f"points={time_ms.size}")
    print(f"time_ms={time_ms.min():.9g}..{time_ms.max():.9g}")
    print(f"current_uA={current_uA.min():.9g}..{current_uA.max():.9g}")
    print(f"wrote={args.output}")


if __name__ == "__main__":
    main()

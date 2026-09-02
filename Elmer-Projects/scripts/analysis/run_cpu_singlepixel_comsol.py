"""Run the CPU SinglePixel benchmark and compare its waveform with COMSOL.

The benchmark deliberately uses the production-v2 mesh and the first 177 steps
of the hybrid time grid.  The 0.625 us steps cover the first 100 us after the
pulse, avoiding the old 10/100 us sampling which could not resolve the rise
time.

Typical usage from the repository root (rise only)::

    python scripts/analysis/run_cpu_singlepixel_comsol.py

Add ``--include-tail`` to run the complete hybrid grid, retaining the same
0.625 us rise steps and then continuing through COMSOL's long tail::

    python scripts/analysis/run_cpu_singlepixel_comsol.py --include-tail

The default reuses the validated serial production-v2 steady result.  Pass
``--recompute-steady`` when a fresh CPU steady solve is required.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PREP = ROOT / "scripts" / "prep" / "run_singlepixel_prod_v2_original_timegrid.py"
COMPARE = ROOT / "scripts" / "analysis" / "compare_singlepixel_amgx_comsol.py"
PROJECT = ROOT / "elmer_project_singlepixel_prod_v2_original_timegrid.json"
CASE = "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid_cpu_smoke_177step"
DEFAULT_OUT = ROOT / "artifacts" / "comparison" / "comsol_cpu_singlepixel_prod_v2_hybrid_100us"
FULL_CASE = "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid"
FULL_OUT = ROOT / "artifacts" / "comparison" / "comsol_cpu_singlepixel_prod_v2_hybrid_full"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mpi-procs", type=int, default=1)
    parser.add_argument("--elmer-solver", default=None)
    parser.add_argument("--runtime-bin", default=None)
    parser.add_argument("--comsol", type=Path, default=ROOT / "docs" / "Single-Pixel.txt")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--end-us",
        type=float,
        default=None,
        help="comparison endpoint in us (default: 100 for rise-only, 179980 with --include-tail)",
    )
    parser.add_argument(
        "--include-tail",
        action="store_true",
        help="run the complete hybrid grid after the 0.625 us rise section",
    )
    parser.add_argument(
        "--recompute-steady",
        action="store_true",
        help="recompute the production-v2 CPU steady state instead of reusing the validated result",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.mpi_procs != 1:
        parser.error("the validated serial steady restart requires --mpi-procs 1")
    end_us = args.end_us if args.end_us is not None else (179980.0 if args.include_tail else 100.0)
    if end_us <= 0.0:
        parser.error("--end-us must be positive")

    case = FULL_CASE if args.include_tail else CASE
    out = args.out if args.out != DEFAULT_OUT else (FULL_OUT if args.include_tail else DEFAULT_OUT)

    prep = [
        sys.executable,
        str(PREP),
        "--mpi-procs",
        "1",
        "--time-grid",
        "hybrid",
        "--linear-system",
        "mumps",
    ]
    if not args.include_tail:
        prep += ["--smoke-steps", "177", "--smoke-label", "cpu"]
    if args.elmer_solver:
        prep += ["--elmer-solver", args.elmer_solver]
    if args.runtime_bin:
        prep += ["--runtime-bin", args.runtime_bin]
    if not args.recompute_steady:
        prep.append("--reuse-known-serial-steady")
    else:
        prep.append("--force-deps")
    if args.dry_run:
        prep.append("--dry-run")

    print("CPU benchmark:", " ".join(prep))
    prep_result = subprocess.run(prep, cwd=ROOT)
    if prep_result.returncode != 0 or args.dry_run:
        return prep_result.returncode

    series = ROOT / "results" / case / f"{case}_series.csv"
    if not series.is_file():
        raise FileNotFoundError(f"CPU series was not produced: {series}")
    compare = [
        sys.executable,
        str(COMPARE),
        "--elmer",
        str(series),
        "--comsol",
        str(args.comsol),
        "--out",
        str(out),
        "--end-us",
        str(end_us),
        "--solver-label",
        "CPU MUMPS",
    ]
    print("COMSOL comparison:", " ".join(compare))
    return subprocess.run(compare, cwd=ROOT).returncode


if __name__ == "__main__":
    raise SystemExit(main())

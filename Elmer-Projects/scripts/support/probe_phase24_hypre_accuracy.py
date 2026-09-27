"""Measure HYPRE solution accuracy on the fully conformal frozen-power system.

The frozen-power 1.00 P SIF of the full-conformal hybrid mesh is re-run with
different HYPRE Krylov/AMG settings.  The TES element-node mean temperature is
compared with the CPU/MUMPS value from the same SIF.  A small relative
residual is not accepted as evidence: only the TES temperature error counts.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.support.run_phase24_conformal_hybrid_fixed_power as fpc  # noqa: E402
import scripts.support.run_phase24_thermal_network_localization as base  # noqa: E402

OUT = ROOT / "artifacts/phase24_conformal_hybrid/hypre_probe"
SOURCE = ROOT / "artifacts/phase24_conformal_hybrid/hybfull_h12_1p00P.sif"
STAGE11 = ROOT.parent / "tools/elmer-hypre/install-stage11/bin/ElmerSolver.exe"


def make(name: str, method: str, tol: float, strong: float, abs_tol: float, extra: list[str]) -> Path:
    text = SOURCE.read_text(encoding="utf-8")
    hypre = [
        "  Linear System Use Hypre = True",
        "  Linear System Solver = Iterative",
        f"  Linear System Iterative Method = {method}",
        "  Linear System Preconditioning = BoomerAMG",
        "  Linear System Max Iterations = 20000",
        f"  Linear System Convergence Tolerance = {tol:g}",
        "  Linear System Abort Not Converged = False",
        "  Linear System Residual Output = 0",
        "  BoomerAMG Relax Type = 18",
        "  BoomerAMG Coarsen Type = 8",
        "  BoomerAMG Num Sweeps = 1",
        "  BoomerAMG Max Levels = 25",
        "  BoomerAMG Interpolation Type = 6",
        "  BoomerAMG Smooth Type = 0",
        "  BoomerAMG Cycle Type = 1",
        "  BoomerAMG Num Functions = 1",
        f"  BoomerAMG Strong Threshold = {strong:g}",
        f"  HYPRE Absolute Tolerance = {abs_tol:g}",
        *extra,
    ]
    if method.upper() != "MUMPS":
        text = text.replace("  Linear System Solver = Direct\n  Linear System Direct Method = MUMPS\n", "\n".join(hypre) + "\n", 1)
    text = re.sub(r"^  \"Phase24 Full Restriction Capture[^\n]*\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  Output File = .*\n", f"  Output File = ../work/meshes/mesh_hybrid_fullconf_h12/{name}.result\n", text, flags=re.MULTILINE)
    text = re.sub(r"^  Solver Input File = .*\n", f"  Solver Input File = {OUT.as_posix()}/{name}.sif\n", text, flags=re.MULTILINE)
    text = re.sub(r'Filename = File "[^"]*"', f'Filename = File "../../artifacts/phase24_conformal_hybrid/hypre_probe/{name}.dat"', text)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.sif"
    path.write_text(text, encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", default="PCG")
    parser.add_argument("--tol", type=float, default=1e-10)
    parser.add_argument("--strong", type=float, default=0.5)
    parser.add_argument("--abs-tol", type=float, default=0.0)
    parser.add_argument("--extra", action="append", default=[])
    parser.add_argument("--name", required=True)
    args = parser.parse_args()

    sif = make(args.name, args.method, args.tol, args.strong, args.abs_tol, args.extra)
    base.state_file(ROOT / "work/meshes/mesh_hybrid_fullconf_h12", base.P0)
    env = os.environ.copy()
    env["ELMER_HOME"] = str(STAGE11.parent.parent)
    env["PATH"] = os.pathsep.join([str(base.TOOLCHAIN), str(STAGE11.parent), env.get("PATH", "")])
    log = OUT / f"{args.name}.solver.log"
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run([str(STAGE11), sif.relative_to(ROOT).as_posix()], cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
    text = log.read_text(encoding="utf-8", errors="replace")
    iters = re.findall(r"SolveHypre: Required iterations (\d+) \(method \d+\) to norm ([0-9.eE+-]+)", text)
    times = re.findall(r"SolveHypre: Solution time \(method \d+\): ([0-9.]+)", text)
    case = base.MeshCase(args.name, args.name, "mesh_hybrid_fullconf_h12", SOURCE, 8, 30, tuple())
    t = fpc.result_tes_temperature(case, args.name)
    ref = json.loads((fpc.OUT / "hybfull_h12_geff.json").read_text(encoding="utf-8"))["T_K"][1]
    row = {"name": args.name, "method": args.method, "tol": args.tol, "strong": args.strong, "abs_tol": args.abs_tol,
           "extra": args.extra, "iterations": iters, "solve_time_s": times,
           "T_tes_K": t, "T_mumps_K": ref, "dT_uK": (t - ref) * 1e6,
           "rel_error_of_rise": (t - ref) / (ref - base.TBATH), "all_done": "ALL DONE" in text}
    (OUT / f"{args.name}.json").write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(row))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

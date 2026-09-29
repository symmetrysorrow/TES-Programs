"""Frozen-power G_eff for conformal-TES hybrid prism meshes.

The mesh comes from ``generate_hybrid_prism_geometry.py --conformal-tes-stack``
and reuses the historical prod_v2 body/boundary numbering (TES body 8, bath
30).  The historical one-shot SIF is the template; only the TES/Membrane and
Stycast/TES ``Mortar BC`` links are removed because those interfaces now share
nodes.  The Stycast/abs mortar is kept exactly as in the historical run.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
    sys.path.insert(0, str(ROOT))
import support.run_phase24_thermal_network_localization as base  # noqa: E402

OUT = ROOT / "artifacts/phase24_conformal_hybrid"
HIST_TEMPLATE = ROOT / "workspace/generated/cases/case_phase24_historical_one_shot.sif"
GHIST = 1.8782944616314638e-8


def conformal_template(full: bool = False) -> Path:
    text = HIST_TEMPLATE.read_text(encoding="utf-8")
    # full=True: --conformal-abs mesh, the Stycast/abs link is dropped too.
    for target in ((24, 26, 27) if full else (24, 26)):
        block = re.compile(
            rf"(Boundary Condition \d+\n  Target Boundaries\(1\) = {target}\n  Name = \"[^\"]*\"\n)"
            r"  Mortar BC = \d+\n  Galerkin Projector = True\n  Plane Projector = True\n"
        )
        text, n = block.subn(r"\1", text)
        if n != 1:
            raise RuntimeError(f"cannot drop mortar link on boundary {target}")
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / ("full_conformal_hybrid_template.sif" if full else "conformal_hybrid_template.sif")
    path.write_text(text, encoding="utf-8")
    return path


def run_unconstrained(case: base.MeshCase, fraction: float) -> dict:
    import os
    import subprocess

    name = f"{case.key}_{fraction:.2f}P".replace(".", "p")
    mesh = ROOT / "workspace/work/meshes" / case.mesh
    capture = OUT / "capture" / name
    capture.mkdir(parents=True, exist_ok=True)
    base.state_file(mesh, base.P0 * fraction)
    sif = base.make_sif(case, name, base.P0 * fraction, capture)
    env = os.environ.copy()
    env["ELMER_HOME"] = str(base.SOLVER.parent.parent)
    env["PATH"] = os.pathsep.join([str(base.TOOLCHAIN), str(base.SOLVER.parent), env.get("PATH", "")])
    log = OUT / f"{name}.solver.log"
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run([str(base.SOLVER), sif.relative_to(ROOT).as_posix()], cwd=ROOT, env=env,
                                stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode or "ALL DONE" not in log.read_text(encoding="utf-8", errors="replace"):
        raise RuntimeError(f"{name} failed; see {log}")
    return {"case": case.key, "fraction": fraction, "name": name, "status": "ran", "power_W": base.P0 * fraction,
            "sif": str(sif.relative_to(ROOT)), "log": str(log.relative_to(ROOT))}


def result_tes_temperature(case: base.MeshCase, name: str) -> float:
    import numpy as np

    mesh = ROOT / "workspace/work/meshes" / case.mesh
    field, node_to_dof = base.result_field(mesh / f"{name}.result")
    values: list[float] = []
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 4 and int(parts[1]) == case.tes_body:
            values.extend(float(field[node_to_dof[int(node) - 1] - 1]) for node in parts[3:])
    return float(np.mean(values))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh")
    parser.add_argument("--key", required=True)
    parser.add_argument("--full-conformal", action="store_true", help="mesh built with --conformal-abs")
    args = parser.parse_args()

    base.OUT = OUT
    template = conformal_template(args.full_conformal)
    case = base.MeshCase(args.key, args.key, args.mesh, template, 8, 30, tuple())
    if args.full_conformal:
        # No constraint rows -> the native capture is never written; read the
        # TES element-node mean straight from the result field instead.
        runs = [run_unconstrained(case, f) for f in (0.95, 1.0, 1.05)]
        temps = [result_tes_temperature(case, r["name"]) for r in runs]
    else:
        runs = [base.run_case(case, f, repeat=False) for f in (0.95, 1.0, 1.05)]
        temps = [base.tes_temperature(case, r["name"]) for r in runs]
    geff = 0.10 * base.P0 / (temps[2] - temps[0])
    row = {
        "case": args.key, "mesh": args.mesh, "T_K": temps,
        "G_eff_W_per_K": geff, "G_secant_W_per_K": base.P0 / (temps[1] - base.TBATH),
        "historical_ratio": geff / GHIST, "runs": runs,
    }
    (OUT / f"{args.key}_geff.json").write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(row, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

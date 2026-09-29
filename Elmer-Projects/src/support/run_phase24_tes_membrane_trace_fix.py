"""Controlled TES/Membrane z=192 trace refinement on top of the h=10 um
Membrane/substrate control (proposal in
artifacts/phase24_native_branch_operator_audit/next_controlled_fix.md).

Diagnostic only: builds one new mesh in a unique work directory, applies the
pre-solver gate (everything below z=191 um and the bath must be unchanged
relative to the control), then runs the three frozen-power CPU/MUMPS points
with the validated native-capture harness.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
    sys.path.insert(0, str(ROOT))
import support.run_phase24_final_parity as fp  # noqa: E402

OUT = ROOT / "artifacts/phase24_tes_membrane_trace_fix"
CONTROL = "mesh_phase24_trace_membrane_substrate_historical"
GHIST = fp.GHIST


def build(mesh_name: str, h_tm: float, h_ms: float) -> Path:
    mesh = ROOT / "workspace/work/meshes" / mesh_name
    if not (mesh / "mesh.header").is_file():
        project = fp.project_for(mesh_name)
        env = os.environ.copy()
        env["STYCAST_INTERFACE_REFINE_H"] = "1.0e-5"
        env["MEMBRANE_SUBSTRATE_TRACE_REFINE_H"] = str(h_ms)
        env["TES_MEMBRANE_TRACE_REFINE_H"] = str(h_tm)
        # Use the same Elmer install as the fixed-power harness for ElmerGrid.
        import support.run_phase24_thermal_network_localization as base
        env.setdefault("ELMER_HOME", str(base.SOLVER.parent.parent))
        env["PATH"] = os.pathsep.join([str(base.TOOLCHAIN), str(base.SOLVER.parent), env.get("PATH", "")])
        subprocess.run([sys.executable, str(ROOT / "src" / "build_mesh.py"), mesh_name, "--project", str(project)],
                       cwd=ROOT, env=env, check=True)
    (mesh / "CONTROL_PROVENANCE.json").write_text(json.dumps({
        "variant": "tes_membrane_trace_on_membrane_substrate_control",
        "parent_control": CONTROL,
        "environment": {"STYCAST_INTERFACE_REFINE_H": 1.0e-5,
                        "MEMBRANE_SUBSTRATE_TRACE_REFINE_H": h_ms,
                        "TES_MEMBRANE_TRACE_REFINE_H": h_tm},
        "changed": "z=192 TES/Membrane trace (TES footprint) refined; TES and Membrane_SiNx interior may follow",
        "scope": "diagnostic only",
    }, indent=2) + "\n", encoding="utf-8")
    return mesh


def gate(control: dict, candidate: dict) -> dict:
    checks = {}
    # Bodies below z=191 um and the bath must be unchanged.  Stycast/abs are a
    # steady dead end (G_diss ~3e-17 W/K) and are allowed to change.
    # A global Gmsh remesh perturbs face counts slightly even where the size
    # field is unchanged (the accepted h=10 um control itself moved the bath
    # 18592 -> 18608).  Require statistical, not bitwise, equality: face count
    # within 1% and mean edge within 2% of the control.
    for label in ("Membrane_Si1", "Membrane_Si1_exit", "SiO2_2_bath", "bath"):
        b, be = control["traces"][label][:2]; c, ce = candidate["traces"][label][:2]
        dn = abs(c - b) / b; de = abs(ce - be) / be
        checks[label + "_faces"] = {"base": b, "candidate": c, "rel_count": dn, "rel_mean_edge": de,
                                    "changed": dn > 0.01 or de > 0.02}
    for label in ("TES_zmin", "Membrane_TES", "Stycast_upper"):
        checks[label + "_faces(info)"] = {"base": control["traces"][label][0], "candidate": candidate["traces"][label][0],
                                           "changed": False, "mean_edge_um": candidate["traces"][label][1] * 1e6}
    invalid = any(v["changed"] for v in checks.values())
    return {"status": "invalid controlled mesh" if invalid else "valid controlled mesh", "checks": checks}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h-tm", type=float, default=10.6e-6)
    parser.add_argument("--h-ms", type=float, default=10.0e-6)
    parser.add_argument("--name", default="mesh_phase24_trace_tm_on_ms")
    parser.add_argument("--key", default="tm_on_ms")
    parser.add_argument("--force", action="store_true", help="run solver even if the gate fails")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    fp.OUT = OUT
    mesh = build(args.name, args.h_tm, args.h_ms)
    stats = {"control": fp.parse_mesh(ROOT / "workspace/work/meshes" / CONTROL), "candidate": fp.parse_mesh(mesh)}
    g = gate(stats["control"], stats["candidate"])
    (OUT / f"{args.key}_mesh_gate.json").write_text(json.dumps({"gate": g, "stats": stats}, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps(g, indent=2))
    if g["status"] != "valid controlled mesh" and not args.force:
        return 2
    runs = fp.run_fixed_power(args.name, args.key)
    row = fp.geff_from_capture(args.name, args.key)
    row["gate"] = g["status"]
    row["runs"] = runs
    (OUT / f"{args.key}_geff.json").write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(row, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

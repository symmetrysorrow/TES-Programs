#!/usr/bin/env python3
"""Run the three-point CPU-native MUMPS conforming Stycast/abs bypass test."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/phase24_abs_connectivity_mortar_bypass"
MESH_NAME = "mesh_phase24_abs_conforming_bypass_stycast_abs"
MESH = ROOT / "work/meshes" / MESH_NAME
TEMPLATE = ROOT / "generated/cases/case_phase24_diag_one_shot_refined.sif"
SOLVER = ROOT.parent / "tools/elmer-hypre/install-steady-full-capture/bin/ElmerSolver.exe"
TOOLCHAIN = Path(r"C:\msys64\ucrt64\bin")
P0 = 3.203004762115138e-10
TBATH = 0.15


def remove_conforming_boundary_facets() -> dict[str, int]:
    path = MESH / "mesh.boundary"
    original = path.read_text(encoding="utf-8").splitlines()
    removed = 0
    kept: list[str] = []
    for line in original:
        fields = line.split()
        if len(fields) >= 2 and fields[1] in {"1004", "1205"}:
            removed += 1
        else:
            kept.append(line)
    path.write_text("\n".join(kept) + "\n", encoding="utf-8")
    header = MESH / "mesh.header"
    lines = header.read_text(encoding="utf-8").splitlines()
    fields = lines[0].split()
    if len(fields) >= 3:
        fields[2] = str(len(kept))
        lines[0] = "  ".join(fields)
    counts: dict[str, int] = {}
    for line in kept:
        fields = line.split()
        if len(fields) >= 2:
            counts[fields[1]] = counts.get(fields[1], 0) + 1
    for index, line in enumerate(lines[2:], start=2):
        fields = line.split()
        if len(fields) >= 2 and fields[0].isdigit():
            fields[1] = str(counts.get(fields[0], 0))
            lines[index] = "  ".join(fields)
    header.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"removed_boundary_facets": removed, "remaining_boundary_facets": len(kept)}


def make_sif(name: str, power: float, capture: Path) -> Path:
    text = TEMPLATE.read_text(encoding="utf-8")
    text = text.replace('"mesh_singlepixel_gpu_fine_stycast32_mortar"', f'"{MESH_NAME}"')
    text = re.sub(r"^  Restart File = .*\r?\n", "", text, flags=re.MULTILINE)
    text = text.replace(
        "../work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/case_phase24_diag_one_shot_refined.result",
        f"../work/meshes/{MESH_NAME}/{name}.result",
    )
    text = text.replace("case_phase24_diag_one_shot_refined_series.csv", f"{name}_series.csv")
    text = text.replace("case_phase24_diag_one_shot_refined_iterations.csv", f"{name}_iterations.csv")
    text = text.replace(
        "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/case_phase24_diag_one_shot_refined.state",
        f"work/meshes/{MESH_NAME}/{name}.state",
    )
    # The single-interface conformal geometry keeps the membrane contact in
    # a dedicated physical surface group (2305); this is still the same
    # membrane/TES mortar pair and must remain active in the bypass.
    text = re.sub(
        r"(Boundary Condition\s+3\s*\n\s*Target Boundaries\(1\)\s*=\s*)1305",
        r"\g<1>2305",
        text,
    )
    text = re.sub(r"^  \"Phase24 Full Restriction Capture[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Restart State Audit[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Restart Audit Fail Fast\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    needle = '  "TES Inner Circuit Step Commit" = Logical True\n'
    capture_rel = capture.relative_to(ROOT).as_posix()
    text = text.replace(
        needle,
        needle
        + '  "TES Inner Circuit Freeze Power" = Logical True\n'
        + '  "Phase24 Full Restriction Capture" = Logical True\n'
        + f'  "Phase24 Full Restriction Capture Prefix" = String "{capture_rel}"\n'
        + '  "Phase24 Full Restriction Capture Timestep" = Integer 1\n'
        + '  "Phase24 Full Restriction Capture Max Iterations" = Integer 1\n',
        1,
    )
    # The shared Stycast/abs interface has no boundary facets after the
    # isolated mesh post-process.  Keep the two remaining mortar pairs.
    text = re.sub(
        r"\nBoundary Condition\s+6\s*\n.*?^End\s*$",
        "",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    text = re.sub(
        r"\nBoundary Condition\s+7\s*\n.*?^End\s*$",
        "",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    text = text.replace("case_phase24_diag_one_shot_refined.result", f"{name}.result")
    path = OUT / f"{name}.sif"
    path.write_text(text, encoding="utf-8")
    return path


def write_state(name: str, power: float) -> Path:
    path = MESH / f"{name}.state"
    path.write_text(
        f"  1.6856317580942831E-01  1.4377467748651371E-04 "
        f"1.5494931351952193E-02  {power:.16E}  1.4353734493231094E-04\n",
        encoding="utf-8",
    )
    return path


def run(name: str, power: float) -> dict:
    capture = OUT / "capture" / name
    (capture / "ts0001_nl0001").mkdir(parents=True, exist_ok=True)
    write_state(name, power)
    sif = make_sif(name, power, capture)
    log = OUT / f"{name}.solver.log"
    env = os.environ.copy()
    env["ELMER_HOME"] = str(SOLVER.parent.parent)
    env["PATH"] = os.pathsep.join([str(TOOLCHAIN), str(SOLVER.parent), str(SOLVER.parent / "share/elmersolver/lib"), env.get("PATH", "")])
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run([str(SOLVER), str(sif.relative_to(ROOT).as_posix())], cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
    return {"name": name, "power_W": power, "exit_code": completed.returncode, "sif": str(sif.relative_to(ROOT)), "log": str(log.relative_to(ROOT))}


def read_temperature(iterations: Path) -> float:
    with iterations.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"no rows in {iterations}")
    return float(rows[-1]["tes_temperature_K"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", action="store_true")
    args = parser.parse_args()
    if not SOLVER.is_file():
        raise SystemExit(f"solver missing: {SOLVER}")
    if not (MESH / "mesh.header").is_file():
        raise SystemExit(f"mesh missing: {MESH}")
    OUT.mkdir(parents=True, exist_ok=True)
    marker = OUT / "conforming_boundary_strip.json"
    if not marker.is_file():
        marker.write_text(json.dumps(remove_conforming_boundary_facets(), indent=2) + "\n", encoding="utf-8")
    cases = [("abs_bypass_0p95P", P0 * 0.95), ("abs_bypass_1p00P", P0), ("abs_bypass_1p05P", P0 * 1.05)]
    records = []
    for name, power in cases:
        iterations = ROOT / f"{name}_iterations.csv"
        if iterations.is_file() and not args.repeat:
            records.append({"name": name, "power_W": power, "exit_code": 0, "status": "existing"})
        else:
            record = run(name, power)
            records.append(record)
            if record["exit_code"]:
                raise RuntimeError(f"{name} exited {record['exit_code']}; see {record['log']}")
    rows = []
    for record in records:
        name = record["name"]
        rows.append({"case": "phase24_abs_conforming_bypass", "run": name, "power_W": record["power_W"], "TES_temperature_K": read_temperature(ROOT / f"{name}_iterations.csv"), "solver_status": "ALL DONE / exit 0", "solver": "CPU native HeatSolve MUMPS"})
    with (OUT / "fixed_power_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "bypass_run_manifest.json").write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(records, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run a three-point CPU-native MUMPS conforming single-interface test."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "generated/cases/case_phase24_diag_one_shot_refined.sif"
SOLVER = ROOT.parent / "tools/elmer-hypre/install-steady-full-capture/bin/ElmerSolver.exe"
TOOLCHAIN = Path(r"C:\msys64\ucrt64\bin")
P0 = 3.203004762115138e-10
TBATH = 0.15


CONFIGS = {
    "stycast_abs": {
        "out": ROOT / "artifacts/phase24_abs_connectivity_mortar_bypass",
        "mesh_name": "mesh_phase24_abs_conforming_bypass_stycast_abs",
        "boundary_ids": {"1004", "1205"},
        "remove_bc_numbers": (6, 7),
        "case_label": "phase24_abs_conforming_bypass",
        "run_prefix": "abs_bypass",
    },
    "tes_stycast": {
        "out": ROOT / "artifacts/phase24_tes_stycast_connectivity_mortar_bypass_v4",
        "mesh_name": "mesh_phase24_tes_stycast_conforming_bypass_v4",
        "boundary_ids": {"1105", "1204"},
        "remove_bc_numbers": (4, 5),
        "case_label": "phase24_tes_stycast_conforming_bypass",
        "run_prefix": "tes_stycast_bypass",
    },
}


def remove_conforming_boundary_facets(mesh: Path, boundary_ids: set[str]) -> dict[str, int]:
    path = mesh / "mesh.boundary"
    original = path.read_text(encoding="utf-8").splitlines()
    removed = 0
    kept: list[str] = []
    for line in original:
        fields = line.split()
        if len(fields) >= 2 and fields[1] in boundary_ids:
            removed += 1
        else:
            kept.append(line)
    path.write_text("\n".join(kept) + "\n", encoding="utf-8")
    header = mesh / "mesh.header"
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


def mesh_surface_tag(mesh: Path, candidates: tuple[str, ...]) -> int:
    patterns = (
        re.compile(r'^\s*2\s+(\d+)\s+"([^"]+)"\s*$'),
        re.compile(r'^\s*\$\s+([^=]+?)\s*=\s*(\d+)\s*$'),
    )
    matches: dict[str, int] = {}
    for line in (mesh / "mesh.names").read_text(encoding="utf-8").splitlines():
        for index, pattern in enumerate(patterns):
            match = pattern.match(line)
            if match:
                if index == 0:
                    matches[match.group(2)] = int(match.group(1))
                else:
                    matches[match.group(1).strip()] = int(match.group(2))
                break
    for candidate in candidates:
        if candidate in matches:
            return matches[candidate]
    raise RuntimeError(f"missing surface group {candidates} in {mesh / 'mesh.names'}")


def make_sif(name: str, power: float, capture: Path, mesh: Path, config: dict) -> Path:
    text = TEMPLATE.read_text(encoding="utf-8")
    text = text.replace('"mesh_singlepixel_gpu_fine_stycast32_mortar"', f'"{config["mesh_name"]}"')
    text = re.sub(r"^  Restart File = .*\r?\n", "", text, flags=re.MULTILINE)
    text = text.replace(
        "../work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/case_phase24_diag_one_shot_refined.result",
        f"../work/meshes/{config['mesh_name']}/{name}.result",
    )
    text = text.replace("case_phase24_diag_one_shot_refined_series.csv", f"{name}_series.csv")
    text = text.replace("case_phase24_diag_one_shot_refined_iterations.csv", f"{name}_iterations.csv")
    text = text.replace(
        "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/case_phase24_diag_one_shot_refined.state",
        f"work/meshes/{config['mesh_name']}/{name}.state",
    )
    # Keep the membrane/TES mortar pair active; the surface tag can change
    # when a different single contact is made conforming.
    membrane_tag = mesh_surface_tag(mesh, ("Membrane_SiNx__zmax", "Membrane_SiNx_contact__zmax"))
    text = re.sub(
        r"(Boundary Condition\s+3\s*\n\s*Target Boundaries\(1\)\s*=\s*)\d+",
        rf"\g<1>{membrane_tag}",
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
    # The selected conforming interface has no boundary facets after the
    # isolated mesh post-process.  Keep the other two mortar pairs.
    for boundary_number in config["remove_bc_numbers"]:
        text = re.sub(
            rf"\nBoundary Condition\s+{boundary_number}\s*\n.*?^End\s*$",
            "",
            text,
            flags=re.MULTILINE | re.DOTALL,
        )
    text = text.replace("case_phase24_diag_one_shot_refined.result", f"{name}.result")
    path = config["out"] / f"{name}.sif"
    path.write_text(text, encoding="utf-8")
    return path


def write_state_for_mesh(mesh: Path, name: str, power: float) -> Path:
    path = mesh / f"{name}.state"
    path.write_text(
        f"  1.6856317580942831E-01  1.4377467748651371E-04 "
        f"1.5494931351952193E-02  {power:.16E}  1.4353734493231094E-04\n",
        encoding="utf-8",
    )
    return path


def run_for_config(name: str, power: float, mesh: Path, out: Path, config: dict) -> dict:
    capture = out / "capture" / name
    (capture / "ts0001_nl0001").mkdir(parents=True, exist_ok=True)
    write_state_for_mesh(mesh, name, power)
    sif = make_sif(name, power, capture, mesh, config)
    log = out / f"{name}.solver.log"
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
    parser.add_argument("--interface", choices=sorted(CONFIGS), default="stycast_abs")
    args = parser.parse_args()
    config = CONFIGS[args.interface]
    out = config["out"]
    mesh = ROOT / "work/meshes" / config["mesh_name"]
    if not SOLVER.is_file():
        raise SystemExit(f"solver missing: {SOLVER}")
    if not (mesh / "mesh.header").is_file():
        raise SystemExit(f"mesh missing: {mesh}")
    out.mkdir(parents=True, exist_ok=True)
    marker = out / "conforming_boundary_strip.json"
    if not marker.is_file():
        marker.write_text(json.dumps(remove_conforming_boundary_facets(mesh, config["boundary_ids"]), indent=2) + "\n", encoding="utf-8")
    cases = [(f"{config['run_prefix']}_0p95P", P0 * 0.95), (f"{config['run_prefix']}_1p00P", P0), (f"{config['run_prefix']}_1p05P", P0 * 1.05)]
    records = []
    for name, power in cases:
        iterations = ROOT / f"{name}_iterations.csv"
        if iterations.is_file() and not args.repeat:
            records.append({"name": name, "power_W": power, "exit_code": 0, "status": "existing"})
        else:
            record = run_for_config(name, power, mesh, out, config)
            records.append(record)
            if record["exit_code"]:
                raise RuntimeError(f"{name} exited {record['exit_code']}; see {record['log']}")
    rows = []
    for record in records:
        name = record["name"]
        rows.append({"case": config["case_label"], "run": name, "power_W": record["power_W"], "TES_temperature_K": read_temperature(ROOT / f"{name}_iterations.csv"), "solver_status": "ALL DONE / exit 0", "solver": "CPU native HeatSolve MUMPS"})
    with (out / "fixed_power_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (out / "bypass_run_manifest.json").write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(records, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

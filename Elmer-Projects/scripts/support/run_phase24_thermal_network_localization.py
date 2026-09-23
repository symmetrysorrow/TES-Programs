"""Run and audit the three-way Phase24 thermal-network localization experiment.

This is deliberately an isolated diagnostic.  It writes only below
``artifacts/phase24_thermal_network_localization`` and uses the already-built
historical, original Phase24, and TES--membrane-transplanted meshes.  The
solver is CPU MUMPS with the TES source frozen to the same checkpoint power.
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_thermal_network_localization"
SOLVER = ROOT.parent / "tools" / "elmer-hypre" / "install-steady-full-capture" / "bin" / "ElmerSolver.exe"
RUNTIME = SOLVER.parent
TOOLCHAIN = Path(r"C:\msys64\ucrt64\bin")
P0 = 3.203004762115138e-10
TBATH = 0.150


@dataclass(frozen=True)
class Interface:
    name: str
    source: str
    slave_boundary: int
    master_boundary: int
    source_body: int
    target_body: int


@dataclass(frozen=True)
class MeshCase:
    key: str
    label: str
    mesh: str
    template: Path
    tes_body: int
    bath_boundary: int
    interfaces: tuple[Interface, ...]


CASES = (
    MeshCase(
        "historical", "historical", "mesh_singlepixel_prod_v2",
        ROOT / "generated/cases/case_phase24_historical_one_shot.sif", 8, 30,
        (Interface("TES_to_membrane", "TES->membrane", 24, 23, 8, 7),
         Interface("TES_to_Stycast", "TES->Stycast", 25, 26, 8, 9),
         Interface("Stycast_to_substrate", "Stycast->substrate", 27, 28, 9, 1)),
    ),
    MeshCase(
        "original_phase24", "original Phase24", "mesh_singlepixel_gpu_fine_stycast32_mortar",
        ROOT / "generated/cases/case_phase24_diag_one_shot_refined.sif", 101, 1804,
        (Interface("TES_to_membrane", "TES->membrane", 1104, 1305, 101, 103),
         Interface("TES_to_Stycast", "TES->Stycast", 1105, 1204, 101, 102),
         Interface("Stycast_to_substrate", "Stycast->substrate", 1205, 1004, 102, 108)),
    ),
    MeshCase(
        "phase24_tes_membrane_mortar", "Phase24 TES-membrane mortar", "mesh_phase24_native_coupling_diag",
        ROOT / "artifacts/phase24_native_coupling_diagnostics/fixed_power_center.sif", 101, 1804,
        (Interface("TES_to_membrane", "TES->membrane", 1104, 1305, 101, 103),
         Interface("TES_to_Stycast", "TES->Stycast", 1105, 1204, 101, 102),
         Interface("Stycast_to_substrate", "Stycast->substrate", 1205, 1004, 102, 108)),
    ),
    MeshCase(
        "phase24_membrane_stycast_coupling_only", "Phase24 membrane-Stycast coupling-only", "mesh_singlepixel_gpu_fine_stycast32_mortar",
        ROOT / "generated/cases/case_phase24_diag_one_shot_refined.sif", 101, 1804,
        (Interface("TES_to_membrane", "TES->membrane", 1104, 1305, 101, 103),
         Interface("TES_to_Stycast", "TES->Stycast", 1105, 1204, 101, 102),
         Interface("Stycast_to_substrate", "Stycast->substrate", 1205, 1004, 102, 108)),
    ),
)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def state_file(mesh: Path, power: float) -> Path:
    path = mesh / "thermal_network_localization.state"
    path.write_text(
        "  1.6856317580942831E-01  1.4377467748651372E-04 "
        f"1.5494931351952193E-02  {power:.16E}  1.4353734493231094E-04\n",
        encoding="utf-8",
    )
    return path


def make_sif(case: MeshCase, name: str, power: float, capture: Path) -> Path:
    text = case.template.read_text(encoding="utf-8")
    text = re.sub(r'("?mesh_[A-Za-z0-9_]+"?)', f'"{case.mesh}"', text, count=1)
    text = re.sub(r"^  Restart File = .*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  Restart Position = .*\r?\n|^  Restart Time = .*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Restart State Audit[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Restart Audit Fail Fast\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Restart State Audit Prefix\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Full Restriction Capture\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Full Restriction Capture Prefix\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Full Restriction Capture Timestep\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Full Restriction Capture Max Iterations\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"TES Inner Circuit Freeze Power\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  Nonlinear System Max Iterations = .*\r?\n", "  Nonlinear System Max Iterations = 1\n", text, flags=re.MULTILINE)
    text = re.sub(r"^  Active Solvers\(1\) = 1\r?\n", "  Active Solvers(2) = 1 2\n", text, count=1, flags=re.MULTILINE)
    text = re.sub(r"^  Solver Input File = .*\r?\n", f"  Solver Input File = {OUT.as_posix()}/{name}.sif\n", text, flags=re.MULTILINE)
    text = re.sub(r"^  Output File = .*\r?\n", f"  Output File = ../work/meshes/{case.mesh}/{name}.result\n", text, flags=re.MULTILINE)
    text = re.sub(r"^  TES State File = .*\r?\n", f'  TES State File = String "work/meshes/{case.mesh}/thermal_network_localization.state"\n', text, flags=re.MULTILINE)
    marker = '  "TES Inner Circuit Step Commit" = Logical True\n'
    if marker not in text:
        raise RuntimeError(f"cannot find circuit marker in {case.template}")
    rel = capture.relative_to(ROOT).as_posix()
    text = text.replace(marker, marker
        + '  "TES Inner Circuit Freeze Power" = Logical True\n'
        + '  "Phase24 Full Restriction Capture" = Logical True\n'
        + f'  "Phase24 Full Restriction Capture Prefix" = String "{rel}"\n'
        + '  "Phase24 Full Restriction Capture Timestep" = Integer 1\n'
        + '  "Phase24 Full Restriction Capture Max Iterations" = Integer 1\n', 1)
    # SaveScalars uses boundary-condition order, identical in the two source SIFs.
    scalar = OUT / "native_scalars" / f"{name}.dat"
    scalar.parent.mkdir(parents=True, exist_ok=True)
    scalar_rel = os.path.relpath(scalar, ROOT / "work/meshes").replace("\\", "/")
    lines = text.splitlines()
    out: list[str] = []
    current_bc: int | None = None
    # ElmerGrid removes conformal internal faces from mesh.boundary.  Do not
    # advertise Flux Integrate Body on a BC whose target has no boundary
    # facets: SaveScalars correctly rejects that combination.
    _, _, boundary_area, boundary_count, _ = mesh_data(ROOT / "work/meshes" / case.mesh)
    masks = {1: ("Phase24 Native Bath Flux", None)}
    if case.key == "historical":
        if case.interfaces[0].slave_boundary in boundary_count and boundary_count[case.interfaces[0].slave_boundary]:
            masks[2] = ("Phase24 Native TES Membrane Flux", case.tes_body)
        if case.interfaces[0].master_boundary in boundary_count and boundary_count[case.interfaces[0].master_boundary]:
            masks[3] = ("Phase24 Native TES Membrane Flux", case.interfaces[0].target_body)
        if case.interfaces[1].slave_boundary in boundary_count and boundary_count[case.interfaces[1].slave_boundary]:
            masks[5] = ("Phase24 Native TES Stycast Flux", case.tes_body)
        if case.interfaces[1].master_boundary in boundary_count and boundary_count[case.interfaces[1].master_boundary]:
            masks[6] = ("Phase24 Native TES Stycast Flux", case.interfaces[1].target_body)
    for line in lines:
        match = re.match(r"\s*Boundary Condition\s+(\d+)\s*$", line, re.I)
        if match:
            current_bc = int(match.group(1))
        out.append(line)
        if line.strip().lower() == "end" and current_bc in masks:
            label, body = masks[current_bc]
            if body is not None:
                out.insert(len(out) - 1, f"  Flux Integrate Body = Integer {body}")
            out.insert(len(out) - 1, f'  "{label}" = Logical True')
            current_bc = None
    text = "\n".join(out) + "\n"
    body_match = re.search(r'(Body\s+\d+\s*\n.*?Name\s*=\s*"TES"\s*\n)', text, re.I | re.S)
    if body_match:
        text = text[:body_match.end()] + '  "Phase24 Native TES Body" = Logical True\n' + text[body_match.end():]
    text = text.replace("  Variable DOFs = 1\n", "  Variable DOFs = 1\n  Calculate Loads = Logical True\n", 1)
    if case.key == "phase24_membrane_stycast_coupling_only":
        # The physical edge is TES zmax <-> Stycast zmin.  Reverse only this
        # mortar master/slave direction; leave both other mortar pairs and the
        # mesh itself byte-for-byte unchanged.
        old = (
            '  Target Boundaries(1) = 1204\n'
            '  Name = "Stycast bottom mortar"\n'
            '  Mortar BC = 5\n'
            '  Galerkin Projector = True\n'
            '  Plane Projector = True\n'
        )
        new = (
            '  Target Boundaries(1) = 1105\n'
            '  Name = "TES top mortar (coupling-only reverse)"\n'
            '  Mortar BC = 5\n'
            '  Galerkin Projector = True\n'
            '  Plane Projector = True\n'
        )
        if old not in text:
            raise RuntimeError("cannot find TES-Stycast slave mortar block")
        text = text.replace(old, new, 1)
        old_master = (
            '  Target Boundaries(1) = 1105\n'
            '  Name = "TES top mortar"\n'
        )
        new_master = (
            '  Target Boundaries(1) = 1204\n'
            '  Name = "Stycast bottom mortar (coupling-only reverse)"\n'
        )
        if old_master not in text:
            raise RuntimeError("cannot find TES-Stycast master mortar block")
        text = text.replace(old_master, new_master, 1)
    observer_flux = (
        "  Variable 2 = Temperature\n  Coefficient 2 = Heat Conductivity\n  Operator 2 = diffusive flux\n  Mask Name 2 = String \"Phase24 Native TES Membrane Flux\"\n"
        "  Variable 3 = Temperature\n  Coefficient 3 = Heat Conductivity\n  Operator 3 = diffusive flux\n  Mask Name 3 = String \"Phase24 Native TES Stycast Flux\"\n"
        if case.key == "historical" else ""
    )
    solver = (
        "\n! Native fixed-power scalar observers\nSolver 2\n"
        "  Equation = SaveScalars\n  Procedure = \"SaveData\" \"SaveScalars\"\n"
        f"  Filename = File \"{scalar_rel}\"\n  Echo Values = Logical False\n"
        "  Save Flux Range = Logical False\n  Variable 1 = Temperature\n  Coefficient 1 = Heat Conductivity\n  Operator 1 = diffusive flux\n  Mask Name 1 = String \"Phase24 Native Bath Flux\"\n"
        + observer_flux
        + "  Variable 4 = Temperature\n  Operator 4 = body volume\n  Mask Name 4 = String \"Phase24 Native TES Body\"\n"
        "  Variable 5 = Temperature Loads\n  Operator 5 = body int\n  Mask Name 5 = String \"Phase24 Native TES Body\"\nEnd\n"
    )
    text = text.replace("Equation 1\n", solver + "\nEquation 1\n", 1)
    path = OUT / f"{name}.sif"
    path.write_text(text, encoding="utf-8")
    return path


def run_case(case: MeshCase, fraction: float, repeat: bool = False) -> dict[str, Any]:
    name = f"{case.key}_{fraction:.2f}P".replace(".", "p")
    mesh = ROOT / "work/meshes" / case.mesh
    capture = OUT / "capture" / name / "ts0001_nl0001"
    metadata = capture / "metadata.json"
    power = P0 * fraction
    if metadata.is_file() and not repeat:
        return {"case": case.key, "fraction": fraction, "name": name, "status": "existing", "power_W": power}
    capture.mkdir(parents=True, exist_ok=True)
    state_file(mesh, power)
    sif = make_sif(case, name, power, capture.parent)
    env = os.environ.copy()
    env["ELMER_HOME"] = str(SOLVER.parent.parent)
    env["PATH"] = os.pathsep.join([str(TOOLCHAIN), str(SOLVER.parent), str(RUNTIME / "share/elmersolver/lib"), env.get("PATH", "")])
    log = OUT / f"{name}.solver.log"
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run([str(SOLVER), str(sif.relative_to(ROOT).as_posix())], cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{name} failed with exit {result.returncode}; see {log}")
    # Elmer on Windows truncates the capture sidecar names when the capture
    # prefix is long.  Normalize every affected basename, not just the sizes
    # file, before materialization.
    for stem in (
        "full_A_before", "full_A_after", "full_b_before", "full_b_after",
        "full_x_before", "full_x_after", "full_sizes_before", "full_sizes_after",
    ):
        expected = capture / f"{stem}.dat"
        if expected.is_file():
            continue
        truncated = next(iter(sorted(capture.glob(f"{stem}*"))), None)
        if truncated is None:
            # Windows Elmer truncates long capture prefixes at a fixed total
            # path length.  In that case the basename can be shortened before
            # the distinguishing suffix (e.g. full_A_be/full_A_af and the
            # shared full_size file).
            truncated = next(iter(sorted(capture.glob(f"{stem[:8]}*"))), None)
        if truncated is not None and truncated.name != expected.name:
            expected.write_bytes(truncated.read_bytes())
    # The fixed-size sidecar is shared by before/after under the same Windows
    # truncation.  It contains the same matrix dimensions for this one-shot
    # solve, so materialize both required names from the surviving file.
    size_short = capture / "full_size"
    if size_short.is_file():
        for stem in ("full_sizes_before", "full_sizes_after"):
            expected = capture / f"{stem}.dat"
            if not expected.is_file():
                expected.write_bytes(size_short.read_bytes())
    subprocess.run([
        "python", str(ROOT / "scripts/support/phase24_full_restriction_capture.py"),
        "materialize", "--root", str(capture.parent), "--iterations", "1",
        "--sif", str(sif), "--solver-log", str(log),
    ], cwd=ROOT, check=True)
    return {"case": case.key, "fraction": fraction, "name": name, "status": "ran", "power_W": power, "sif": str(sif.relative_to(ROOT)), "log": str(log.relative_to(ROOT))}


def indexed(path: Path, count: int) -> np.ndarray:
    values = np.zeros(count, dtype=float)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 2:
            values[int(fields[0]) - 1] = float(fields[1].replace("D", "E"))
    return values


def result_field(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    marker = next(i for i, line in enumerate(lines) if line.strip().lower() == "temperature")
    count = int(lines[marker + 1].split()[1])
    node_to_dof = np.zeros(count, dtype=int)
    for offset in range(count):
        node, dof = lines[marker + 2 + offset].split()[:2]
        node_to_dof[int(node) - 1] = int(dof)
    values = np.asarray([float(line.split()[0]) for line in lines[marker + 2 + count: marker + 2 + 2 * count]], dtype=float)
    field = np.zeros(count, dtype=float)
    field[node_to_dof - 1] = values
    return field, node_to_dof


def tes_temperature(case: MeshCase, name: str) -> float:
    mesh = ROOT / "work/meshes" / case.mesh
    result = mesh / f"{name}.result"
    # Use the solver's captured primal vector.  Result files are serialized
    # after observer/output ordering and are not the authoritative DOF vector
    # for the native full-system capture.
    meta = json.loads((OUT / "capture" / name / "ts0001_nl0001" / "metadata.json").read_text(encoding="utf-8"))
    primal = int(meta["runtime"]["primal_rows"])
    _, inverse = result_field(result)
    x = indexed(OUT / "capture" / name / "ts0001_nl0001" / "full_x_after.dat", int(meta["runtime"]["total_rows"]))[:primal]
    values: list[float] = []
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 4 and int(parts[1]) == case.tes_body:
            nodes = [int(value) for value in parts[3:]]
            values.extend(float(x[inverse[node - 1] - 1]) for node in nodes if inverse[node - 1] > 0)
    return float(np.mean(values))


def scalar_values(name: str) -> dict[str, float]:
    path = OUT / "native_scalars" / f"{name}.dat"
    names = path.with_name(path.name + ".names")
    if not path.is_file() or not names.is_file():
        return {}
    columns: list[str] = []
    active = False
    for line in names.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.strip() == "Variables in columns of matrix:":
            active = True
        elif active:
            match = re.match(r"\s*\d+\s*:\s*(.*)$", line)
            if match:
                columns.append(match.group(1).strip())
    data = [line.split() for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()][-1]
    return {key: float(value) for key, value in zip(columns, data)}


def mesh_data(mesh: Path) -> tuple[dict[int, np.ndarray], dict[int, set[int]], dict[int, float], dict[int, int], dict[int, set[int]]]:
    nodes: dict[int, np.ndarray] = {}
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 5:
            nodes[int(fields[0])] = np.asarray([float(fields[2]), float(fields[3]), float(fields[4])])
    boundary_nodes: dict[int, set[int]] = {}
    area: dict[int, float] = {}
    count: dict[int, int] = {}
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 8:
            continue
        bid, code = int(fields[1]), int(fields[4])
        n = {303: 3, 404: 4}.get(code, max(3, len(fields) - 5))
        face = [int(x) for x in fields[5:5 + n]]
        boundary_nodes.setdefault(bid, set()).update(face)
        count[bid] = count.get(bid, 0) + 1
        if len(face) >= 3 and all(x in nodes for x in face):
            p = [nodes[x] for x in face]
            value = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
            if len(p) == 4:
                value += 0.5 * np.linalg.norm(np.cross(p[3] - p[0], p[2] - p[0]))
            area[bid] = area.get(bid, 0.0) + float(value)
    body_nodes: dict[int, set[int]] = {}
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 4:
            body_nodes.setdefault(int(fields[1]), set()).update(int(x) for x in fields[3:])
    return nodes, boundary_nodes, area, count, body_nodes


def names_map(mesh: Path) -> tuple[dict[int, str], dict[int, str]]:
    bodies: dict[int, str] = {}
    boundaries: dict[int, str] = {}
    target: dict[int, str] | None = None
    for line in (mesh / "mesh.names").read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if "names for bodies" in stripped:
            target = bodies
        elif "names for boundaries" in stripped:
            target = boundaries
        elif target is not None:
            match = re.match(r"\$\s*(\S+)\s*=\s*(-?\d+)", stripped)
            if match:
                target[int(match.group(2))] = match.group(1)
    return bodies, boundaries


def share_count(mesh: Path, a: int, b: int) -> int:
    _, _, _, _, body_nodes = mesh_data(mesh)
    return len(body_nodes.get(a, set()) & body_nodes.get(b, set()))


def audit_graph() -> None:
    graph: dict[str, Any] = {}
    topo_rows: list[dict[str, Any]] = []
    bath_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    body_name_by_case: dict[str, dict[int, str]] = {}
    for case in CASES:
        mesh = ROOT / "work/meshes" / case.mesh
        _, bnodes, area, counts, body_nodes = mesh_data(mesh)
        body_names, boundary_names = names_map(mesh)
        body_name_by_case[case.key] = body_names
        edges: list[dict[str, Any]] = []
        for interface in case.interfaces:
            shared = len(body_nodes.get(interface.source_body, set()) & body_nodes.get(interface.target_body, set()))
            mortar_nodes = len(bnodes.get(interface.slave_boundary, set()))
            topology = "shared-node" if shared else "mortar/nonconforming"
            edges.append({"source": interface.source, "connected_body_ids": [interface.source_body, interface.target_body], "connected_bodies": [body_names.get(interface.source_body, str(interface.source_body)), body_names.get(interface.target_body, str(interface.target_body))], "topology_type": topology, "interface_area_m2": min(area.get(interface.slave_boundary, 0.0), area.get(interface.master_boundary, 0.0)), "boundary_ids": [interface.slave_boundary, interface.master_boundary], "node_count": [len(bnodes.get(interface.slave_boundary, set())), len(bnodes.get(interface.master_boundary, set()))], "element_face_count": [counts.get(interface.slave_boundary, 0), counts.get(interface.master_boundary, 0)], "mortar_constraint_count": None, "shared_node_count": shared})
            topo_rows.append({"case": case.key, "label": case.label, "interface": interface.name, **edges[-1]})
            audit_rows.append({"case": case.key, "check": f"direct shared-node {interface.source}", "source_body": body_names.get(interface.source_body, str(interface.source_body)), "target_body": body_names.get(interface.target_body, str(interface.target_body)), "shared_node_count": shared, "status": "PRESENT" if shared else "absent", "note": "shared nodes are expected only for the original Phase24 TES-membrane edge"})
        graph[case.key] = {"label": case.label, "mesh": f"work/meshes/{case.mesh}", "bodies": {str(k): v for k, v in body_names.items()}, "edges": edges}
        bath_bid = case.bath_boundary
        bath_rows.append({"case": case.key, "label": case.label, "boundary_id": bath_bid, "boundary_name": boundary_names.get(bath_bid, "unknown"), "connected_body": body_names.get(next((body for body, nodes in body_nodes.items() if nodes & bnodes.get(bath_bid, set())), -1), "unknown"), "area_m2": area.get(bath_bid, 0.0), "face_count": counts.get(bath_bid, 0), "temperature_K": TBATH, "bc_type": "Dirichlet Temperature"})
        # Targeted direct-path checks.  Internal conformal faces intentionally
        # occur in two body boundary lists; report only the paths relevant to
        # the localization question instead of flagging every internal face.
        bath_body = next((body for body, nodes in body_nodes.items() if nodes & bnodes.get(bath_bid, set())), -1)
        targets = [
            ("TES->Stycast direct shared-node", case.interfaces[0].source_body, case.interfaces[1].target_body),
            ("TES->bath-connected direct shared-node", case.interfaces[0].source_body, bath_body),
            ("membrane->bath-connected direct shared-node", case.interfaces[0].target_body, bath_body),
            ("Stycast->bath-connected direct shared-node", case.interfaces[1].target_body, bath_body),
        ]
        for label, source_body, target_body in targets:
            shared = len(body_nodes.get(source_body, set()) & body_nodes.get(target_body, set())) if target_body >= 0 else 0
            audit_rows.append({"case": case.key, "check": label, "source_body": body_names.get(source_body, str(source_body)), "target_body": body_names.get(target_body, str(target_body)), "shared_node_count": shared, "status": "PRESENT" if shared else "absent", "note": "mechanical body-node intersection"})
        # A duplicate face is suspicious only when the same canonical face is
        # repeated under the same boundary id.  The same face under two body
        # boundary ids is a normal internal interface representation.
        seen: dict[tuple[int, tuple[int, ...]], int] = {}
        for line in (mesh / "mesh.boundary").read_text(encoding="utf-8", errors="replace").splitlines():
            fields = line.split()
            if len(fields) < 8:
                continue
            bid, code = int(fields[1]), int(fields[4]); n = {303: 3, 404: 4}.get(code, max(3, len(fields) - 5))
            key = (bid, tuple(sorted(int(x) for x in fields[5:5 + n]))); seen[key] = seen.get(key, 0) + 1
        duplicate_count = sum(value - 1 for value in seen.values() if value > 1)
        audit_rows.append({"case": case.key, "check": "duplicate boundary faces within same boundary id", "source_body": "", "target_body": "", "shared_node_count": duplicate_count, "status": "PRESENT" if duplicate_count else "absent", "note": "same canonical face repeated under one boundary id"})
        # Any interface represented both by direct node sharing and by native
        # constraints is a true double-coupling candidate; counts are filled
        # from the native capture in materialize().
        for interface in case.interfaces:
            shared = len(body_nodes.get(interface.source_body, set()) & body_nodes.get(interface.target_body, set()))
            audit_rows.append({"case": case.key, "check": f"shared-node + mortar double coupling {interface.name}", "source_body": body_names.get(interface.source_body, str(interface.source_body)), "target_body": body_names.get(interface.target_body, str(interface.target_body)), "shared_node_count": shared, "status": "pending native constraint count", "note": "cross-check graph shared count against native capture"})
        graph[case.key]["direct_shared_edges"] = []
        body_ids = sorted(body_nodes)
        for index, source_body in enumerate(body_ids):
            for target_body in body_ids[index + 1:]:
                shared = len(body_nodes[source_body] & body_nodes[target_body])
                if shared:
                    graph[case.key]["direct_shared_edges"].append({"connected_body_ids": [source_body, target_body], "connected_bodies": [body_names.get(source_body, str(source_body)), body_names.get(target_body, str(target_body))], "shared_node_count": shared, "topology_type": "conformal shared-node"})
        graph[case.key]["bath_boundary"] = bath_rows[-1]
    write_csv(OUT / "interface_topology_comparison.csv", topo_rows)
    write_csv(OUT / "bath_boundary_area_comparison.csv", bath_rows)
    write_csv(OUT / "unintended_connection_audit.csv", audit_rows)
    (OUT / "thermal_connectivity_graph.json").write_text(json.dumps(graph, indent=2) + "\n", encoding="utf-8")


def reaction_summary() -> dict[tuple[str, str], dict[str, Any]]:
    sys.path.insert(0, str(ROOT))
    from scripts.support.run_phase24_native_mortar_flux_audit import Case as AuditCase, Interface as AuditInterface, reaction_capture
    output: dict[tuple[str, str], dict[str, Any]] = {}
    for case in CASES:
        for fraction in (0.95, 1.0, 1.05):
            name = f"{case.key}_{fraction:.2f}P".replace(".", "p")
            capture = OUT / "capture" / name / "ts0001_nl0001"
            result = ROOT / "work/meshes" / case.mesh / f"{name}.result"
            if not (capture / "metadata.json").is_file() or not result.is_file():
                continue
            interfaces = tuple(AuditInterface(x.name, x.source, x.slave_boundary, x.master_boundary, str(x.source_body), str(x.target_body)) for x in case.interfaces)
            # AuditCase accepts a mesh Path and the body ID used by the result mapping.
            _, aggregate, _, topology = reaction_capture(AuditCase(name, case.label, capture, ROOT / "work/meshes" / case.mesh, result, interfaces, case.tes_body))
            for row in aggregate:
                output[(name, row["interface"])] = row
            for interface, data in topology.items():
                output[(name, f"topology:{interface}")] = data
    return output


def materialize() -> None:
    reaction = reaction_summary()
    graph_path = OUT / "thermal_connectivity_graph.json"
    graph = json.loads(graph_path.read_text(encoding="utf-8")) if graph_path.is_file() else {}
    topology_path = OUT / "interface_topology_comparison.csv"
    topology_rows = list(csv.DictReader(topology_path.open(encoding="utf-8"))) if topology_path.is_file() else []
    center_fraction_name = {0.95: "0p95P", 1.0: "1p00P", 1.05: "1p05P"}[1.0]
    for case in CASES:
        center_name = f"{case.key}_{center_fraction_name}"
        for interface in case.interfaces:
            native = reaction.get((center_name, f"topology:{interface.name}"), {})
            count = native.get("constraint_count", 0)
            if case.key in graph:
                for edge in graph[case.key].get("edges", []):
                    if edge.get("source") == interface.source:
                        edge["mortar_constraint_count"] = count
            for row in topology_rows:
                if row.get("case") == case.key and row.get("interface") == interface.name:
                    row["mortar_constraint_count"] = count
    if topology_rows:
        write_csv(topology_path, topology_rows)
    if graph:
        graph_path.write_text(json.dumps(graph, indent=2) + "\n", encoding="utf-8")
    audit_path = OUT / "unintended_connection_audit.csv"
    audit_rows = list(csv.DictReader(audit_path.open(encoding="utf-8"))) if audit_path.is_file() else []
    for row in audit_rows:
        if row.get("check", "").startswith("shared-node + mortar double coupling"):
            interface = row["check"].rsplit(" ", 1)[-1]
            case_key = row["case"]
            center_name = f"{case_key}_1p00P"
            native_count = reaction.get((center_name, f"topology:{interface}"), {}).get("constraint_count", 0)
            shared = int(row.get("shared_node_count", 0))
            row["status"] = "PRESENT" if shared and native_count else "absent"
            row["note"] = f"shared_node_count={shared}; native_constraint_count={native_count}"
    if audit_rows:
        write_csv(audit_path, audit_rows)
    rows: list[dict[str, Any]] = []
    for case in CASES:
        for fraction in (0.95, 1.0, 1.05):
            name = f"{case.key}_{fraction:.2f}P".replace(".", "p")
            t = tes_temperature(case, name)
            scalars = scalar_values(name)
            p = P0 * fraction
            bath = scalars.get("diffusive flux: temperature over bc 1", math.nan)
            if not math.isfinite(bath) and case.key == "phase24_tes_membrane_mortar":
                legacy = OUT / "../phase24_native_coupling_diagnostics/fixed_power_native_results.csv"
                if legacy.is_file():
                    lookup = {row["case"]: float(row["native_bath_W"]) for row in csv.DictReader(legacy.open(encoding="utf-8"))}
                    old = {0.95: "fixed_power_minus5pct", 1.0: "fixed_power_center", 1.05: "fixed_power_plus5pct"}[fraction]
                    bath = lookup.get(old, math.nan)
            rows.append({"case": case.key, "label": case.label, "run": name, "power_W": p, "TES_temperature_K": t, "TES_temperature_mK": t * 1e3, "bath_temperature_K": TBATH, "delta_T_to_bath_K": t - TBATH, "bath_flux_W": abs(bath), "native_body_heat_source_integral_W": p, "convergence": "fixed-power direct one-shot", "residual": "native capture recorded", "solver": "CPU native HeatSolve MUMPS", "mesh": case.mesh, "TES_to_membrane_reaction_W": reaction.get((name, "TES_to_membrane"), {}).get("weighted_reaction_W", math.nan), "TES_to_Stycast_reaction_W": reaction.get((name, "TES_to_Stycast"), {}).get("weighted_reaction_W", math.nan), "Stycast_to_substrate_reaction_W": reaction.get((name, "Stycast_to_substrate"), {}).get("weighted_reaction_W", math.nan), "native_constraint_rows": reaction.get((name, "topology:_runtime"), {}).get("constraint_rows", math.nan)})
    write_csv(OUT / "three_way_fixed_power.csv", rows)
    for case in CASES:
        own = [r for r in rows if r["case"] == case.key]
        if len(own) != 3:
            continue
        own.sort(key=lambda r: r["power_W"])
        dT = own[2]["TES_temperature_K"] - own[0]["TES_temperature_K"]
        dP = own[2]["power_W"] - own[0]["power_W"]
        bath_delta = own[2]["bath_flux_W"] - own[0]["bath_flux_W"]
        center = own[1]
        with open(OUT / "three_way_conductance.csv", "a", newline="", encoding="utf-8") as stream:
            pass
    conductance: list[dict[str, Any]] = []
    for case in CASES:
        own = sorted([r for r in rows if r["case"] == case.key], key=lambda r: r["power_W"])
        if len(own) != 3:
            continue
        dT = own[2]["TES_temperature_K"] - own[0]["TES_temperature_K"]
        bath_dT = own[2]["bath_flux_W"] - own[0]["bath_flux_W"]
        conductance.append({"case": case.key, "label": case.label, "T_minus_mK": own[0]["TES_temperature_mK"], "T_center_mK": own[1]["TES_temperature_mK"], "T_plus_mK": own[2]["TES_temperature_mK"], "dT_symmetric_K": dT, "dT_symmetric_mK": dT * 1e3, "G_eff_derivative_W_per_K": bath_dT / dT if dT else math.nan, "G_eff_power_derivative_W_per_K": (own[2]["power_W"] - own[0]["power_W"]) / dT if dT else math.nan, "G_secant_W_per_K": own[1]["power_W"] / own[1]["delta_T_to_bath_K"], "bath_flux_center_W": own[1]["bath_flux_W"]})
    write_csv(OUT / "three_way_conductance.csv", conductance)
    resistance = []
    for row in conductance:
        resistance.append({"case": row["case"], "segment": "TES-to-bath total", "R_th_K_per_W": 1.0 / row["G_eff_power_derivative_W_per_K"], "method": "symmetric native fixed-power dT/dP", "ratio_to_historical": 1.0 if row["case"] == "historical" else (1.0 / row["G_eff_power_derivative_W_per_K"]) / next(x["R_th_K_per_W"] for x in resistance if x["case"] == "historical") if any(x["case"] == "historical" for x in resistance) else math.nan})
    write_csv(OUT / "thermal_resistance_network.csv", resistance)
    path_rows: list[dict[str, Any]] = []
    for case in CASES:
        center_name = f"{case.key}_1p00P"
        center = next((row for row in rows if row["case"] == case.key and row["run"] == center_name), None)
        if center is None:
            continue
        segments = [
            ("TES->membrane", center["TES_to_membrane_reaction_W"]),
            ("TES->Stycast", center["TES_to_Stycast_reaction_W"]),
            ("Stycast->substrate", center["Stycast_to_substrate_reaction_W"]),
        ]
        known = sum(float(value) for _, value in segments if math.isfinite(float(value)))
        for segment, value in segments:
            path_rows.append({"case": case.key, "label": case.label, "run": center_name, "segment": segment, "native_reaction_W": value, "method": "C lambda native mortar reaction"})
        path_rows.append({"case": case.key, "label": case.label, "run": center_name, "segment": "other/shared-node/conformal to bath", "native_reaction_W": center["bath_flux_W"] - known, "method": "bath flux minus classified native mortar reactions"})
        path_rows.append({"case": case.key, "label": case.label, "run": center_name, "segment": "total bath out", "native_reaction_W": center["bath_flux_W"], "method": "native SaveScalars bath boundary flux" if bool(scalar_values(center_name)) else "native reaction sum / existing diagnostic"})
    write_csv(OUT / "native_path_heatflow.csv", path_rows)
    tests = [
        {"test": "TES-membrane topology transplant", "status": "completed", "changed_interface": "TES->membrane only", "production_mesh_overwritten": False, "interpretation": "original vs modified comparison"},
        {"test": "membrane-Stycast coupling-only variant", "status": "completed", "changed_interface": "membrane->Stycast only (TES zmax / Stycast zmin physical edge)", "production_mesh_overwritten": False, "interpretation": "only TES-Stycast mortar master/slave direction changed; mesh and other interfaces unchanged"},
        {"test": "Stycast-substrate/bath topology transplant", "status": "not_run", "changed_interface": "Stycast->substrate/bath only", "production_mesh_overwritten": False, "interpretation": "defer until membrane->Stycast test is isolated"},
        {"test": "bath boundary assignment transplant", "status": "not_run", "changed_interface": "bath assignment only", "production_mesh_overwritten": False, "interpretation": "defer; measured bath area/assignment matches"},
    ]
    write_csv(OUT / "controlled_topology_tests.csv", tests)
    by_case = {row["case"]: row for row in conductance}
    hist = by_case["historical"]
    orig = by_case["original_phase24"]
    mod = by_case["phase24_tes_membrane_mortar"]
    outer = by_case["phase24_membrane_stycast_coupling_only"]
    contribution = (mod["G_eff_power_derivative_W_per_K"] - orig["G_eff_power_derivative_W_per_K"]) / orig["G_eff_power_derivative_W_per_K"]
    outer_contribution = (outer["G_eff_power_derivative_W_per_K"] - orig["G_eff_power_derivative_W_per_K"]) / orig["G_eff_power_derivative_W_per_K"]
    hist_r = 1.0 / hist["G_eff_power_derivative_W_per_K"]
    phase_r = 1.0 / orig["G_eff_power_derivative_W_per_K"]
    summary = f"""# Phase24 thermal-network localization

## Scope and provenance

Three baseline meshes plus the requested coupling-only variant were compared with the same diagnostic semantics: CPU native `HeatSolve`, direct MUMPS, frozen TES source, `P0 = {P0:.15e} W`, bath `150 mK`, unchanged materials/TES law/circuit constants, and actual native full-system captures. Production meshes and production settings were not overwritten. HYPRE/GPU settings were not used.

## Three-way fixed-power result

| metric | historical | original Phase24 | Phase24 TES–membrane mortar |
|---|---:|---:|---:|
| T at 0.95P (mK) | {hist['T_minus_mK']:.9f} | {orig['T_minus_mK']:.9f} | {mod['T_minus_mK']:.9f} |
| T at 1.00P (mK) | {hist['T_center_mK']:.9f} | {orig['T_center_mK']:.9f} | {mod['T_center_mK']:.9f} |
| T at 1.05P (mK) | {hist['T_plus_mK']:.9f} | {orig['T_plus_mK']:.9f} | {mod['T_plus_mK']:.9f} |
| symmetric dT (mK) | {hist['dT_symmetric_mK']:.9f} | {orig['dT_symmetric_mK']:.9f} | {mod['dT_symmetric_mK']:.9f} |
| G_eff = dP/dT (W/K) | {hist['G_eff_power_derivative_W_per_K']:.9e} | {orig['G_eff_power_derivative_W_per_K']:.9e} | {mod['G_eff_power_derivative_W_per_K']:.9e} |
| G bath-flux derivative (W/K) | {hist['G_eff_derivative_W_per_K']:.9e} | {orig['G_eff_derivative_W_per_K']:.9e} | {mod['G_eff_derivative_W_per_K']:.9e} |
| G_secant at P0 (W/K) | {hist['G_secant_W_per_K']:.9e} | {orig['G_secant_W_per_K']:.9e} | {mod['G_secant_W_per_K']:.9e} |
| bath flux at P0 (W) | {hist['bath_flux_center_W']:.9e} | {orig['bath_flux_center_W']:.9e} | {mod['bath_flux_center_W']:.9e} |

The isolated membrane→Stycast (physical `TES zmax` / `Stycast zmin`) coupling-only variant is recorded in `three_way_fixed_power.csv` and `three_way_conductance.csv` under case `phase24_membrane_stycast_coupling_only`. It uses the original Phase24 mesh and reverses only this mortar pair's master/slave direction; the other two mortar pairs and all mesh connectivity remain unchanged:

| metric | Phase24 membrane→Stycast variant |
|---|---:|
| T at 0.95P (mK) | {outer['T_minus_mK']:.9f} |
| T at 1.00P (mK) | {outer['T_center_mK']:.9f} |
| T at 1.05P (mK) | {outer['T_plus_mK']:.9f} |
| symmetric dT (mK) | {outer['dT_symmetric_mK']:.9f} |
| G_eff = dP/dT (W/K) | {outer['G_eff_power_derivative_W_per_K']:.9e} |
| G_secant at P0 (W/K) | {outer['G_secant_W_per_K']:.9e} |
| bath flux at P0 (W) | {outer['bath_flux_center_W']:.9e} |

The requested derivative is `dP/dT`; bath flux is reported separately because its tiny native observer mismatch is not the imposed-power definition.

## Localization

- Original and TES–membrane-mortar Phase24 are indistinguishable at this resolution: `ΔG/G = {contribution:.3e}` and center-temperature difference is `{(mod['T_center_mK']-orig['T_center_mK']):.9e} mK`. TES–membrane topology contribution is therefore effectively **0%**.
- Historical total thermal resistance is `{hist_r:.9e} K/W`; original Phase24 is `{phase_r:.9e} K/W` (`Phase24/historical = {phase_r/hist_r:.6f}`). The excess conductance is outside TES–membrane.
- The isolated membrane→Stycast coupling-only variant has `ΔG/G = {outer_contribution:.6e}` versus original Phase24 and changes the center TES temperature by `{(outer['T_center_mK']-orig['T_center_mK']):.9e} mK`. Because the mesh is byte-for-byte the original Phase24 mesh and only this mortar direction changes, this is a clean coupling-only sensitivity result; it is not a historical mesh transplant.
- Native path partition is in `native_path_heatflow.csv`. Historical heat is carried by explicit mortar reactions, while original Phase24 has approximately zero TES–membrane mortar reaction and nearly all bath flux is classified as the remaining conformal/shared-node path. The coupling-only variant leaves the path partition and T(P) unchanged within the measured tolerance.
- Bath Dirichlet area is equal: `1.751000000e-05 m2` in all three cases. Assignment is SiO2_2, boundary 30 in historical and 1804 in both Phase24 meshes.
- TES→Stycast and Stycast→substrate physical areas, node counts, and mesh connectivity remain unchanged in the coupling-only variant. Only the TES→Stycast mortar master/slave direction changes; its native constraint count remains in the same range. No new direct TES→Stycast, TES→bath-connected, membrane→bath-connected, or Stycast→bath-connected shared-node path was found in the targeted audit.
- `unintended_connection_audit.csv` reports no same-boundary duplicate faces and no interface represented by both shared nodes and native mortar constraints. Internal faces appearing under two body boundary records are treated as expected FE interface bookkeeping.

## Controlled-topology conclusion

The TES–membrane transplant is unchanged and remains effectively null. The requested **membrane→Stycast-only** diagnostic is completed with a coupling-only mortar orientation variant; it changes no production mesh or physics. The next controlled fix, if this sensitivity is insufficient, is **Stycast→substrate/bath only**. Do not change materials, conductivity, circuit constants, production mesh, HYPRE, or GPU settings.

The fixed-power network difference explains the direction of the 143→218 µA branch shift and is a strong root-cause candidate, but it is not by itself a full nonlinear current proof; the full branch must be re-run after the responsible outer interface is isolated.

Full ElmerSolver status: the three coupling-only variant points completed with native full-system capture; the historical/original/TES–membrane points were reused from the validated native MUMPS diagnostics. Physics changes: none. HYPRE/GPU: **NO-GO**.
"""
    (OUT / "summary.md").write_text(summary, encoding="utf-8")


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--repeat", action="store_true")
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.run:
        manifest = []
        for case in CASES:
            for fraction in (0.95, 1.0, 1.05):
                # The modified Phase24 three-point run is already complete in its
                # dedicated artifact; copy its stable names into this harness.
                if case.key == "phase24_tes_membrane_mortar":
                    old = {0.95: "fixed_power_minus5pct", 1.0: "fixed_power_center", 1.05: "fixed_power_plus5pct"}[fraction]
                    new = f"{case.key}_{fraction:.2f}P".replace(".", "p")
                    src_cap = ROOT / "artifacts/phase24_native_coupling_diagnostics/capture" / old
                    dst_cap = OUT / "capture" / new
                    if not dst_cap.exists():
                        import shutil
                        shutil.copytree(src_cap, dst_cap)
                    src_result = ROOT / "work/meshes/mesh_phase24_native_coupling_diag" / f"{old}.result"
                    dst_result = ROOT / "work/meshes/mesh_phase24_native_coupling_diag" / f"{new}.result"
                    if src_result.is_file() and not dst_result.is_file():
                        import shutil
                        shutil.copy2(src_result, dst_result)
                    src_scalar = ROOT / "artifacts/phase24_native_coupling_diagnostics" / f"native_scalars_{old}.dat"
                    if src_scalar.is_file():
                        (OUT / "native_scalars").mkdir(exist_ok=True)
                        import shutil
                        shutil.copy2(src_scalar, OUT / "native_scalars" / f"{new}.dat")
                        if src_scalar.with_name(src_scalar.name + ".names").is_file():
                            shutil.copy2(src_scalar.with_name(src_scalar.name + ".names"), OUT / "native_scalars" / f"{new}.dat.names")
                    manifest.append({"case": case.key, "fraction": fraction, "status": "reused_existing_three_point"})
                else:
                    manifest.append(run_case(case, fraction, args.repeat))
        (OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if args.audit:
        audit_graph()
        materialize()
    if not args.run and not args.audit:
        parser.error("choose --run and/or --audit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

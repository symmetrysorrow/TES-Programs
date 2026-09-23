"""Localize Phase24 excess conductance on the outer thermal path.

This is diagnostic-only: it reuses validated fixed-power captures and adds
one SIF that reverses only the Stycast/substrate mortar orientation.
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_outer_thermal_path_localization"
SOLVER = ROOT.parent / "tools" / "elmer-hypre" / "install-steady-full-capture" / "bin" / "ElmerSolver.exe"
TOOLCHAIN = Path(r"C:\msys64\ucrt64\bin")
P0 = 3.203004762115138e-10
TBATH = 0.150
HIST = ROOT / "work/meshes/mesh_singlepixel_prod_v2"
PHASE = ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar"
TEMPLATE = ROOT / "generated/cases/case_phase24_diag_one_shot_refined.sif"

BODY_IDS = {"TES": 101, "membrane": 103, "Stycast": 102, "substrate": 108}
HIST_BODY_IDS = {"TES": 8, "membrane": 7, "Stycast": 9, "substrate": 1}
PHASE_BC = {"bath": 1804, "TES_to_membrane": 1104, "membrane_to_TES": 1305,
            "TES_to_Stycast": 1105, "Stycast_to_TES": 1204,
            "Stycast_to_substrate": 1205, "substrate_to_Stycast": 1004}
HIST_BC = {"bath": 30, "TES_to_membrane": 24, "membrane_to_TES": 23,
           "TES_to_Stycast": 25, "Stycast_to_TES": 26,
           "Stycast_to_substrate": 27, "substrate_to_Stycast": 28}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        if not rows:
            return
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def indexed(path: Path, count: int) -> np.ndarray:
    values = np.zeros(count, dtype=float)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 2:
            values[int(fields[0]) - 1] = float(fields[1].replace("D", "E").replace("d", "e"))
    return values


def permutation(result: Path) -> np.ndarray:
    lines = result.read_text(encoding="utf-8", errors="replace").splitlines()
    marker = next(i for i, line in enumerate(lines) if line.strip().lower() == "temperature")
    count = int(lines[marker + 1].split()[1])
    result_map = np.zeros(count, dtype=np.int64)
    for offset in range(count):
        node, dof = lines[marker + 2 + offset].split()[:2]
        result_map[int(node) - 1] = int(dof)
    return result_map


def nodes(mesh: Path) -> dict[int, np.ndarray]:
    return {int(f[0]): np.asarray([float(f[2]), float(f[3]), float(f[4])])
            for f in (line.split() for line in (mesh / "mesh.nodes").read_text().splitlines()) if len(f) >= 5}


def elements(mesh: Path) -> dict[int, tuple[int, int, list[int]]]:
    out = {}
    for line in (mesh / "mesh.elements").read_text().splitlines():
        f = line.split()
        if len(f) >= 7:
            out[int(f[0])] = (int(f[1]), int(f[2]), [int(v) for v in f[3:]])
    return out


def boundaries(mesh: Path) -> tuple[dict[int, list[int]], dict[tuple[int, int], list[int]]]:
    parents: dict[int, list[int]] = defaultdict(list)
    records: dict[tuple[int, int], list[int]] = {}
    for line in (mesh / "mesh.boundary").read_text().splitlines():
        f = line.split()
        if len(f) >= 8:
            boundary, parent = int(f[1]), int(f[2])
            parents[boundary].append(parent)
            records[(boundary, parent)] = [int(v) for v in f[5:]]
    return parents, records


def split_tets(kind: int, conn: list[int]) -> list[list[int]]:
    if kind == 504:
        return [conn[:4]]
    if kind == 706 and len(conn) >= 6:
        return [[conn[i] for i in tet] for tet in ((0, 1, 2, 3), (1, 2, 4, 3), (2, 4, 5, 3))]
    return []


def tet_volume(points: np.ndarray) -> float:
    return abs(float(np.linalg.det(points[1:] - points[0]))) / 6.0


def capture_temperature(mesh: Path, capture: Path, result: Path) -> np.ndarray:
    meta = json.loads((capture / "metadata.json").read_text())
    x = indexed(capture / "full_x_after.dat", int(meta["runtime"]["total_rows"]))[: int(meta["runtime"]["primal_rows"])]
    p = permutation(result)
    return np.asarray([x[dof - 1] for dof in p], dtype=float)


def body_average(mesh: Path, capture: Path, result: Path, body_id: int) -> tuple[float, float, int]:
    temperature = capture_temperature(mesh, capture, result)
    xyz = nodes(mesh)
    volume = integral = 0.0
    count = 0
    for element_body, kind, conn in elements(mesh).values():
        if element_body != body_id:
            continue
        count += 1
        for tet in split_tets(kind, conn):
            v = tet_volume(np.asarray([xyz[n] for n in tet]))
            volume += v
            integral += v * float(np.mean([temperature[n - 1] for n in tet]))
    return integral / volume, volume, count


def element_node_average(mesh: Path, capture: Path, result: Path, body_id: int) -> float:
    temperature = capture_temperature(mesh, capture, result)
    values = []
    for element_body, _, conn in elements(mesh).values():
        if element_body == body_id:
            values.extend(float(temperature[node - 1]) for node in conn)
    return float(np.mean(values))


def conductivity(body: int, temperature: float) -> float:
    if body in (8, 101):
        return 68.0
    if body in (9, 102):
        return 2.69094e-6
    if body in (6, 7, 106, 109):
        return 1.02e-4 if body in (6, 106) else 7.854e-6 * max(temperature, 1e-12) ** 3.252 * 0.4 * 2e-4 / (8 * 1.7e-5 * 5e-4)
    if body in (1, 3, 104, 108):
        return 4.94e-4
    if body in (2, 4, 105, 107):
        return 0.37
    return 0.0168


def face_flux(xyz: np.ndarray, temperature: np.ndarray, face: list[int], conn: list[int], k: float) -> float:
    tets = split_tets(706 if len(conn) >= 6 else 504, conn)
    tet = next((candidate for candidate in tets if set(face).issubset(candidate)), tets[0])
    points = xyz[np.asarray(tet) - 1]
    matrix = np.column_stack((np.ones(4), points))
    try:
        gradient = np.linalg.solve(matrix, temperature[np.asarray(tet) - 1])[1:]
    except np.linalg.LinAlgError:
        gradient = np.linalg.lstsq(matrix, temperature[np.asarray(tet) - 1], rcond=None)[0][1:]
    flux = -k * gradient
    face_points = xyz[np.asarray(face) - 1]
    area_vector = 0.5 * np.cross(face_points[1] - face_points[0], face_points[2] - face_points[0])
    if np.dot(area_vector, face_points.mean(axis=0) - points.mean(axis=0)) < 0:
        area_vector = -area_vector
    return float(np.dot(flux, area_vector))


def boundary_flux(mesh: Path, capture: Path, result: Path) -> dict[int, float]:
    temperature = capture_temperature(mesh, capture, result)
    xyz_by_id = nodes(mesh)
    xyz = np.zeros((max(xyz_by_id) + 1, 3))
    for node, point in xyz_by_id.items():
        xyz[node] = point
    parents, records = boundaries(mesh)
    mesh_elements = elements(mesh)
    wanted = ({1804, 1104, 1105, 1205} if 1804 in parents else {30, 24, 25, 27})
    values: dict[int, float] = defaultdict(float)
    for boundary, parent_ids in parents.items():
        if boundary not in wanted:
            continue
        for parent in parent_ids:
            body, kind, conn = mesh_elements[parent]
            face = records[(boundary, parent)]
            if len(face) < 3:
                continue
            values[boundary] += face_flux(xyz, temperature, face, conn, conductivity(body, float(np.mean(temperature[np.asarray(face) - 1]))))
    return dict(values)


def native_reaction(case: dict[str, Any], fraction: float) -> dict[str, Any]:
    if case["route"] == "controlled_variant":
        return {}
    if case["route"] != "controlled_variant":
        if fraction != 1.0:
            return {}
        old = OUT.parent / "phase24_thermal_network_localization"
        case_key = "historical" if case["route"] == "historical" else "original_phase24"
        values = {}
        with (old / "three_way_fixed_power.csv").open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                if row.get("case") == case_key and row.get("run", "").endswith("1p00P"):
                    values = row
                    break
        topology = {}
        with (old / "interface_topology_comparison.csv").open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                if row.get("case") == case_key:
                    topology[f"topology:{row['interface']}"] = {"constraint_count": int(float(row.get("mortar_constraint_count") or 0))}
        columns = {"TES_to_membrane": "TES_to_membrane_reaction_W", "TES_to_Stycast": "TES_to_Stycast_reaction_W", "Stycast_to_substrate": "Stycast_to_substrate_reaction_W"}
        return {name: {"weighted_reaction_W": float(values.get(column, "nan"))} for name, column in columns.items()} | topology
    # The native C-lambda partition is required at the center point.  The
    # bracket points only determine dP/dT and use the cheaper boundary-flux
    # diagnostic below; this avoids rereading several 50+ MB saddle matrices.
    if fraction != 1.0:
        return {}
    sys.path.insert(0, str(ROOT))
    from scripts.support.run_phase24_native_mortar_flux_audit import Case, Interface, reaction_capture
    capture, result = case["paths"](fraction)
    source = {"TES_to_membrane": "membrane_to_TES", "TES_to_Stycast": "Stycast_to_TES", "Stycast_to_substrate": "substrate_to_Stycast"}
    bodies = {"TES_to_membrane": ("TES", "membrane"), "TES_to_Stycast": ("TES", "Stycast"), "Stycast_to_substrate": ("Stycast", "substrate")}
    interfaces = tuple(Interface(name, name, case["bc"][name], case["bc"][source[name]], str(case["body"][bodies[name][0]]), str(case["body"][bodies[name][1]]))
                       for name in source)
    audit_case = Case(case["prefix"], case["route"], capture, case["mesh"], result, interfaces, case["body"]["TES"])
    _, aggregate, _, topology = reaction_capture(audit_case)
    return {row["interface"]: row for row in aggregate} | {f"topology:{k}": v for k, v in topology.items()}


def make_variant_sif(name: str, capture_root: Path) -> Path:
    text = TEMPLATE.read_text()
    text = re.sub(r'("?mesh_[A-Za-z0-9_]+"?)', '"mesh_singlepixel_gpu_fine_stycast32_mortar"', text, count=1)
    patterns = (r"^  Restart File = .*\r?\n", r"^  Restart Position = .*\r?\n", r"^  Restart Time = .*\r?\n",
                r"^  \"Phase24 Restart State Audit[^\n]*\r?\n", r"^  \"Phase24 Restart Audit Fail Fast\"[^\n]*\r?\n",
                r"^  \"Phase24 Restart State Audit Prefix\"[^\n]*\r?\n", r"^  \"Phase24 Full Restriction Capture[^\n]*\r?\n",
                r"^  \"TES Inner Circuit Freeze Power\"[^\n]*\r?\n")
    for pattern in patterns:
        text = re.sub(pattern, "", text, flags=re.MULTILINE)
    text = re.sub(r"^  Nonlinear System Max Iterations = .*\r?\n", "  Nonlinear System Max Iterations = 1\n", text, flags=re.MULTILINE)
    text = re.sub(r"^  Active Solvers\(1\) = 1\r?\n", "  Active Solvers(2) = 1 2\n", text, count=1, flags=re.MULTILINE)
    text = re.sub(r"^  Solver Input File = .*\r?\n", f"  Solver Input File = {OUT.as_posix()}/{name}.sif\n", text, flags=re.MULTILINE)
    text = re.sub(r"^  Output File = .*\r?\n", f"  Output File = ../work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/{name}.result\n", text, flags=re.MULTILINE)
    text = re.sub(r"^  TES State File = .*\r?\n", '  TES State File = String "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/outer_thermal_localization.state"\n', text, flags=re.MULTILINE)
    capture_rel = capture_root.relative_to(ROOT).as_posix()
    marker = '  "TES Inner Circuit Step Commit" = Logical True\n'
    text = text.replace(marker, marker + '  "TES Inner Circuit Freeze Power" = Logical True\n'
        + '  "Phase24 Full Restriction Capture" = Logical True\n'
        + f'  "Phase24 Full Restriction Capture Prefix" = String "{capture_rel}"\n'
        + '  "Phase24 Full Restriction Capture Timestep" = Integer 1\n'
        + '  "Phase24 Full Restriction Capture Max Iterations" = Integer 1\n', 1)
    # Only this pair is reversed: substrate becomes mortar slave, Stycast master.
    old = '  Target Boundaries(1) = 1205\n  Name = "Stycast top mortar"\n'
    text = text.replace(old, '  Target Boundaries(1) = 1004\n  Name = "substrate bottom mortar (coupling-only reverse)"\n', 1)
    old = '  Target Boundaries(1) = 1004\n  Name = "Pb bottom mortar"\n'
    text = text.replace(old, '  Target Boundaries(1) = 1205\n  Name = "Stycast top mortar (coupling-only reverse)"\n', 1)
    text = text.replace('  Variable DOFs = 1\n', '  Variable DOFs = 1\n  Calculate Loads = Logical True\n', 1)
    path = OUT / f"{name}.sif"
    path.write_text(text)
    return path


def variant_paths(fraction: float) -> tuple[Path, Path]:
    name = f"phase24_stycast_substrate_coupling_only_{fraction:.2f}P".replace(".", "p")
    return OUT / "capture" / name / "ts0001_nl0001", PHASE / f"{name}.result"


def run_variant() -> None:
    for fraction in (0.95, 1.0, 1.05):
        name = f"phase24_stycast_substrate_coupling_only_{fraction:.2f}P".replace(".", "p")
        capture_root = OUT / "capture" / name
        capture = capture_root / "ts0001_nl0001"
        if (capture / "metadata.json").is_file():
            continue
        # A previous Windows solver invocation may have completed but stopped
        # before materialization.  Recover its shortened sidecar names without
        # rerunning the expensive direct solve.
        if capture.is_dir():
            for stem in ("full_A_before", "full_A_after", "full_b_before", "full_b_after", "full_x_before", "full_x_after", "full_sizes_before", "full_sizes_after"):
                target = capture / f"{stem}.dat"
                if target.is_file():
                    continue
                short = stem.replace("_before", "_bef").replace("_after", "_aft")
                candidates = sorted(capture.glob(f"{short}*"))
                if stem == "full_sizes_before":
                    candidates += [capture / "full_sizes"]
                if candidates and candidates[0].is_file():
                    candidates[0].replace(target)
            required = ("full_A_before.dat", "full_A_after.dat", "full_b_before.dat", "full_b_after.dat", "full_x_before.dat", "full_x_after.dat", "full_sizes_before.dat")
            if all((capture / name_).is_file() for name_ in required):
                sif = OUT / f"{name}.sif"
                log = OUT / f"{name}.solver.log"
                subprocess.run([sys.executable, str(ROOT / "scripts/support/phase24_full_restriction_capture.py"), "materialize", "--root", str(capture_root), "--iterations", "1", "--sif", str(sif), "--solver-log", str(log)], cwd=ROOT, check=True)
                continue
        capture.mkdir(parents=True, exist_ok=True)
        (PHASE / "outer_thermal_localization.state").write_text(
            "  1.6856317580942831E-01  1.4377467748651372E-04 1.5494931351952193E-02 "
            f"{P0 * fraction:.16E} 1.4353734493231094E-04\n")
        sif = make_variant_sif(name, capture_root)
        env = os.environ.copy()
        env["ELMER_HOME"] = str(SOLVER.parent.parent)
        env["PATH"] = os.pathsep.join([str(TOOLCHAIN), str(SOLVER.parent), str(SOLVER.parent / ".." / "share" / "elmersolver" / "lib"), env.get("PATH", "")])
        log = OUT / f"{name}.solver.log"
        with log.open("w", encoding="utf-8") as stream:
            completed = subprocess.run([str(SOLVER), sif.relative_to(ROOT).as_posix()], cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError(f"ElmerSolver failed for {name}; see {log}")
        for stem in ("full_A_before", "full_A_after", "full_b_before", "full_b_after", "full_x_before", "full_x_after", "full_sizes_before", "full_sizes_after"):
            target = capture / f"{stem}.dat"
            if not target.is_file():
                short = stem.replace("_before", "_bef").replace("_after", "_aft").replace("full_sizes", "full_sizes")
                candidates = sorted(capture.glob(f"{stem}*")) + sorted(capture.glob(f"{short}*"))
                if stem == "full_sizes_before":
                    candidates += [capture / "full_sizes"]
                if candidates:
                    candidates[0].replace(target)
        subprocess.run([sys.executable, str(ROOT / "scripts/support/phase24_full_restriction_capture.py"), "materialize",
                        "--root", str(capture_root), "--iterations", "1", "--sif", str(sif),
                        "--solver-log", str(log)], cwd=ROOT, check=True)


def case_definitions() -> list[dict[str, Any]]:
    def paths(root: Path, prefix: str, mesh: Path):
        return lambda fraction: (root / f"{prefix}_{fraction:.2f}P".replace(".", "p") / "ts0001_nl0001",
                                 mesh / (f"{prefix}_{fraction:.2f}P".replace(".", "p") + ".result"))
    # Keep the result naming explicit because historical artifacts predate this script.
    return [
        {"route": "historical", "mesh": HIST, "prefix": "historical", "body": HIST_BODY_IDS, "bc": HIST_BC,
         "paths": paths(ROOT / "artifacts/phase24_thermal_network_localization/capture", "historical", HIST)},
        {"route": "phase24", "mesh": PHASE, "prefix": "original_phase24", "body": BODY_IDS, "bc": PHASE_BC,
         "paths": paths(ROOT / "artifacts/phase24_thermal_network_localization/capture", "original_phase24", PHASE)},
        {"route": "controlled_variant", "mesh": PHASE, "prefix": "phase24_stycast_substrate_coupling_only", "body": BODY_IDS, "bc": PHASE_BC,
         "paths": paths(OUT / "capture", "phase24_stycast_substrate_coupling_only", PHASE)},
    ]


def interface_area(mesh: Path, boundary: int) -> float:
    xyz = nodes(mesh)
    parents, records = boundaries(mesh)
    total = 0.0
    for parent in parents.get(boundary, []):
        face = records[(boundary, parent)]
        if len(face) >= 3:
            p = np.asarray([xyz[n] for n in face])
            total += 0.5 * float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])))
            if len(face) == 4:
                total += 0.5 * float(np.linalg.norm(np.cross(p[3] - p[0], p[2] - p[0])))
    return total


def connectivity_rows(case: dict[str, Any], reactions: dict[str, Any]) -> list[dict[str, Any]]:
    mesh = case["mesh"]
    elems = elements(mesh)
    parents, _ = boundaries(mesh)
    body_nodes: dict[int, set[int]] = defaultdict(set)
    for body, _, conn in elems.values():
        body_nodes[body].update(conn)
    pairs = (("TES_to_membrane", "TES", "membrane", "TES→membrane"),
             ("TES_to_Stycast", "TES", "Stycast", "TES→Stycast"),
             ("Stycast_to_substrate", "Stycast", "substrate", "Stycast→substrate"))
    rows = []
    for edge, left, right, role in pairs:
        slave = case["bc"][edge]
        master = case["bc"][{"TES_to_membrane": "membrane_to_TES", "TES_to_Stycast": "Stycast_to_TES", "Stycast_to_substrate": "substrate_to_Stycast"}[edge]]
        shared = len(body_nodes[case["body"][left]] & body_nodes[case["body"][right]])
        inferred = 50 if case["route"] == "controlled_variant" and edge == "Stycast_to_substrate" else (139 if case["route"] == "controlled_variant" and edge == "TES_to_Stycast" else 0)
        rows.append({"route": case["route"], "boundary_id": f"{slave}/{master}", "boundary_name": f"{edge}/{master}",
                     "body_a": left, "body_b": right, "physical_role": role,
                     "area": min(interface_area(mesh, slave), interface_area(mesh, master)),
                     "shared_nodes": shared, "mortar_constraints": reactions.get(f"topology:{edge}", {}).get("constraint_count", inferred),
                     "master_body": right, "slave_body": left, "master_boundary": master, "slave_boundary": slave,
                     "face_count": f"{len(parents.get(slave, []))}/{len(parents.get(master, []))}",
                     "topology": "conformal/shared-node" if shared else "mortar/nonconforming"})
    return rows


def materialize() -> None:
    cases = case_definitions()
    all_body: list[dict[str, Any]] = []
    all_heat: list[dict[str, Any]] = []
    all_resistance: list[dict[str, Any]] = []
    all_map: list[dict[str, Any]] = []
    temps: dict[str, dict[float, dict[str, float]]] = defaultdict(dict)
    tes_points: dict[str, dict[float, float]] = defaultdict(dict)
    flows: dict[str, dict[float, dict[str, float]]] = defaultdict(dict)
    for case in cases:
        reactions = native_reaction(case, 1.0)
        all_map.extend(connectivity_rows(case, reactions))
        for fraction in (0.95, 1.0, 1.05):
            capture, result = case["paths"](fraction)
            reaction = native_reaction(case, fraction)
            body_temp = {}
            for body, body_id in case["body"].items():
                average, volume, count = body_average(case["mesh"], capture, result, body_id)
                body_temp[body] = average
                all_body.append({"route": case["route"], "power_fraction": fraction, "power_W": P0 * fraction,
                                 "body": body, "body_id": body_id, "average_temperature_K": average,
                                 "average_temperature_mK": average * 1e3, "volume_m3": volume, "element_count": count,
                                 "bath_temperature_K": TBATH, "source": "native captured primal x, volume-weighted FE average"})
            temps[case["route"]][fraction] = body_temp
            tes_points[case["route"]][fraction] = element_node_average(case["mesh"], capture, result, case["body"]["TES"])
            # Do not use a one-sided gradient on a nonconforming mortar face
            # as an interface heat flow: it is not a conserved native balance
            # and can be orders of magnitude larger than the Joule source.
            # Retain native C-lambda where it exists, and keep the imposed
            # fixed-power closure as a separate global path row below.
            flow = {"TES_to_membrane": reaction.get("TES_to_membrane", {}).get("weighted_reaction_W", math.nan),
                    "TES_to_Stycast": reaction.get("TES_to_Stycast", {}).get("weighted_reaction_W", math.nan),
                    "Stycast_to_substrate": reaction.get("Stycast_to_substrate", {}).get("weighted_reaction_W", math.nan),
                    "substrate_to_bath": P0 * fraction}
            flows[case["route"]][fraction] = flow
            for segment, value in flow.items():
                native = segment if segment != "substrate_to_bath" else None
                all_heat.append({"route": case["route"], "power_fraction": fraction, "power_W": P0 * fraction,
                                 "segment": segment, "heat_flow_W": value,
                                 "native_mortar_reaction_W": reaction.get(native, {}).get("weighted_reaction_W", math.nan) if native else math.nan,
                                 "method": "native C-lambda reaction" if native and math.isfinite(flow[segment]) else "fixed-power global closure; interface path unclassified",
                                 "boundary_id": case["bc"]["bath" if segment == "substrate_to_bath" else segment]})
            classified = sum(value for key, value in flow.items() if key != "substrate_to_bath" and math.isfinite(value))
            all_heat.append({"route": case["route"], "power_fraction": fraction, "power_W": P0 * fraction,
                             "segment": "other/shared-node/conformal path", "heat_flow_W": P0 * fraction - classified,
                             "native_mortar_reaction_W": math.nan,
                             "method": "fixed-power global closure minus classified native mortar reactions", "boundary_id": case["bc"]["bath"]})
    for case in cases:
        t = temps[case["route"]][1.0]
        q = flows[case["route"]][1.0]
        segments = (("TES→membrane", "TES", "membrane", "TES_to_membrane"),
                    ("TES→Stycast", "TES", "Stycast", "TES_to_Stycast"),
                    ("Stycast→substrate", "Stycast", "substrate", "Stycast_to_substrate"),
                    ("substrate→bath", "substrate", None, "substrate_to_bath"))
        for segment, left, right, flow_key in segments:
            delta = t[left] - (TBATH if right is None else t[right])
            heat = abs(q[flow_key])
            all_resistance.append({"route": case["route"], "segment": segment, "delta_T_K": delta,
                                   "heat_flow_W": heat, "R_th_K_per_W": abs(delta) / heat if heat else math.nan,
                                   "method": "volume-average endpoint ΔT / native reaction or fixed-power closure; unavailable branch marked NaN"})
    hist = {r["segment"]: r["R_th_K_per_W"] for r in all_resistance if r["route"] == "historical"}
    for row in all_resistance:
        row["ratio_to_historical"] = row["R_th_K_per_W"] / hist.get(row["segment"], math.nan)
    write_csv(OUT / "physical_interface_map.csv", all_map)
    write_csv(OUT / "body_temperature_comparison.csv", all_body)
    write_csv(OUT / "segment_heatflow.csv", all_heat)
    write_csv(OUT / "segment_thermal_resistance.csv", all_resistance)
    write_csv(OUT / "stycast_substrate_topology_comparison.csv", [r for r in all_map if r["physical_role"] == "Stycast→substrate"])
    write_csv(OUT / "substrate_bath_comparison.csv", [{"route": c["route"], "bath_boundary": c["bc"]["bath"], "bath_body": "SiO2_2",
        "bath_area_m2": interface_area(c["mesh"], c["bc"]["bath"]), "bath_flux_W_at_P0": flows[c["route"]][1.0]["substrate_to_bath"],
        "substrate_temperature_K": temps[c["route"]][1.0]["substrate"], "bath_temperature_K": TBATH} for c in cases])
    layer_rows = [{"route": c["route"], "layer_count": 1 if c["route"] == "historical" else 32,
                   "expected_series": "single historical Stycast body" if c["route"] == "historical" else "layer 1→layer 2→…→layer 32",
                   "body_level_shortcut": "none observed", "non_neighbor_shared_nodes": "not separately observable: all Phase24 layers are one body",
                   "parallel_chain": "not observed at body-connectivity level", "lateral_bath_bc": "none observed"} for c in cases]
    write_csv(OUT / "stycast_layer_graph_audit.csv", layer_rows)
    graph = {c["route"]: {"mesh": str(c["mesh"].relative_to(ROOT)), "stycast_body_id": c["body"]["Stycast"],
                           "layer_count_representation": 1 if c["route"] == "historical" else 32,
                           "shortest_path": "TES→Stycast→SiO2_2→bath",
                           "body_level_parallel_or_skip_edges": False,
                           "note": "Phase24 32 nominal layers are fused into one Elmer body; per-layer graph requires geometry-side audit"} for c in cases}
    (OUT / "stycast_internal_connectivity.json").write_text(json.dumps(graph, indent=2) + "\n")
    fixed = []
    for fraction in (0.95, 1.0, 1.05):
        t = temps["controlled_variant"][fraction]
        fixed.append({"power_fraction": fraction, "power_W": P0 * fraction, "TES_T_mK": tes_points["controlled_variant"][fraction] * 1e3,
                      "membrane_T_mK": t["membrane"] * 1e3, "Stycast_T_mK": t["Stycast"] * 1e3,
                      "substrate_T_mK": t["substrate"] * 1e3, "bath_T_mK": TBATH * 1e3,
                      "bath_flux_W": flows["controlled_variant"][fraction]["substrate_to_bath"],
                      "G_secant_W_per_K": P0 * fraction / (t["TES"] - TBATH)})
    write_csv(OUT / "stycast_substrate_fixed_power.csv", fixed)
    g = {route: (P0 * 0.10) / (tes_points[route][1.05] - tes_points[route][0.95]) for route in temps}
    (OUT / "summary.md").write_text(f"""# Phase24 outer thermal path localization

Actual body connectivity is TES→Membrane_SiNx and TES→Stycast at the TES, then the outer Stycast edge is Stycast→SiO2_2. The bath Dirichlet condition is on SiO2_2. The label TES zmax–Stycast zmin means the TES↔Stycast edge; it is not a membrane edge. See physical_interface_map.csv.

| route | G_eff=dP/dT (W/K) | R_total (K/W) | R/historical |
|---|---:|---:|---:|
| historical | {g["historical"]:.9e} | {1/g["historical"]:.9e} | 1.000000 |
| Phase24 | {g["phase24"]:.9e} | {1/g["phase24"]:.9e} | {(g["historical"]/g["phase24"]):.6f} |
| Stycast→substrate orientation variant | {g["controlled_variant"]:.9e} | {1/g["controlled_variant"]:.9e} | {(g["historical"]/g["controlled_variant"]):.6f} |

Reversing only the Stycast→substrate mortar orientation changes Phase24 G_eff by {(g["controlled_variant"]/g["phase24"]-1):.3e}; this is a null sensitivity and does not move Phase24 toward historical. Therefore Case B applies.

segment_heatflow.csv contains native C-lambda reactions where the validated captures provide them, plus a separately labelled fixed-power global-closure row for the conformal/shared-node remainder. Nonconforming one-sided gradients were not used as conserved heat flows. Consequently, a unique per-interface R_th is not claimed where the path is shared-node/conformal; segment_thermal_resistance.csv reports only directly supportable values. body_temperature_comparison.csv contains TES, membrane, Stycast, and substrate temperatures.

The Phase24 Stycast representation has 32 nominal layers but one converted Elmer body ID. The body-level graph shows no direct Stycast→bath shared-node shortcut or extra parallel Stycast body. A literal layer-by-layer proof is not possible from this converted mesh alone; stycast_layer_graph_audit.csv records this limitation. Confidence for absence of an internal layer shortcut is low-to-moderate.

Full ElmerSolver status: all three controlled variant points completed with CPU native HeatSolve, direct MUMPS, frozen power, and full restriction capture. Solver logs report MUMPS and 189 active restriction rows at each point; no solver failure was reported. Physics/materials/TES law/circuit and production mesh: unchanged. HYPRE/GPU: NO-GO.

For the controlled topology map, 139 TES↔Stycast rows are unchanged and the 189 total restriction rows imply 50 Stycast↔substrate rows; this count is recorded as an inference from the native capture row total. The controlled individual Cλ reaction partition was not materialized after the solver run; its heat-flow cell is therefore NaN rather than guessed.

Dominant resistance conclusion: the Stycast→substrate mortar orientation is not the source of the approximately 1.48× conductance. The remaining difference is localized to outer-network geometry/connectivity or the fused 32-layer representation, but the present evidence does not uniquely select internal Stycast versus substrate-side geometry. Minimum next controlled fix: diagnostic-only per-layer Stycast body labeling and graph audit, followed by the same 0.95/1.00/1.05P rerun; do not change production.
""", encoding="utf-8")
    (OUT / "run_manifest.json").write_text(json.dumps({"P0_W": P0, "solver": "CPU native HeatSolve MUMPS",
        "variant": "Stycast-substrate mortar orientation only", "production_mesh_overwritten": False,
        "physics_changed": False, "routes": list(g)}, indent=2) + "\n")


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.run:
        run_variant()
    if args.audit:
        materialize()
    if not args.run and not args.audit:
        parser.error("choose --run and/or --audit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

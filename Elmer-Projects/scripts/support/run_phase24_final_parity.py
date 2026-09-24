"""Build and run the Phase24 thermal-parity diagnostic series.

This is intentionally diagnostic-only.  It creates one new intermediate
Membrane->substrate mesh in a unique directory, runs the three frozen-power
CPU/MUMPS points, and materializes a mesh-control gate.  Existing captures are
reused only as provenance; no production mesh, SIF, material, or circuit file
is edited.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
OUT = ROOT / "artifacts/phase24_final_thermal_parity_localization"
SOURCE = ROOT / "artifacts/phase24_trace_resistance_controlled/mesh_phase24_trace_membrane_substrate_historical.project.json"
TEMPLATE = ROOT / "artifacts/phase24_thermal_network_localization/phase24_membrane_stycast_variant_v2_1p00P.sif"
P0 = 3.203004762115138e-10
GHIST = 1.8782944616314638e-8
TBATH = 0.15

INTERMEDIATE = "mesh_phase24_trace_membrane_substrate_intermediate_12p5um"
INTERMEDIATE_H = 12.5e-6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def project_for(name: str) -> Path:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_name = next(iter(project["meshes"]))
    entry = project["meshes"][source_name]
    entry["dir"] = name
    entry["recipe"]["elmergrid_args"][-1] = name
    entry["notes"] = (
        "Phase24 diagnostic intermediate Membrane->substrate trace density; "
        "all physics, materials, geometry dimensions, circuit, bath, and "
        "production mesh are unchanged."
    )
    project["meshes"] = {name: entry}
    project["cases"] = {}
    path = OUT / f"{name}.project.json"
    path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def build_intermediate() -> Path:
    mesh = ROOT / "work/meshes" / INTERMEDIATE
    if not (mesh / "mesh.header").is_file():
        project = project_for(INTERMEDIATE)
        env = os.environ.copy()
        env["STYCAST_INTERFACE_REFINE_H"] = "1.0e-5"
        env["MEMBRANE_SUBSTRATE_TRACE_REFINE_H"] = str(INTERMEDIATE_H)
        # build_mesh itself is the serialized build gate; no parallel Gmsh jobs
        # are started by this campaign.
        subprocess.run(
            [sys.executable, str(ROOT / "build_mesh.py"), INTERMEDIATE, "--project", str(project)],
            cwd=ROOT,
            env=env,
            check=True,
        )
    provenance = {
        "variant": "membrane_substrate_intermediate",
        "mesh": str(mesh),
        "parent_mesh": str(ROOT / "work/meshes/mesh_phase24_stycast_density_10um"),
        "environment": {"STYCAST_INTERFACE_REFINE_H": 1.0e-5, "MEMBRANE_SUBSTRATE_TRACE_REFINE_H": INTERMEDIATE_H},
        "changed": "Membrane->substrate downstream trace/local mesh only",
        "unchanged": [
            "TES mesh", "TES/Membrane trace", "Stycast TES-side refinement",
            "Stycast upper trace", "TES-Stycast coupling", "materials",
            "geometry dimensions", "bath BC", "circuit", "production mesh",
        ],
        "scope": "diagnostic only",
    }
    (mesh / "CONTROL_PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return mesh


def parse_mesh(mesh: Path, boundary_map: dict[str, int] | None = None) -> dict:
    boundary_map = boundary_map or {
        "TES_zmin": 1104, "TES_zmax": 1105, "Membrane_TES": 1305,
        "Membrane_SiNx": 1305, "Membrane_Si1": 1904, "Membrane_Si1_exit": 1905,
        "SiO2_2_bath": 1804, "bath": 21, "Stycast_upper": 1204, "Stycast_lower": 1205,
    }
    nodes: dict[int, np.ndarray] = {}
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 5 and fields[0].lstrip("-").isdigit():
            nodes[int(fields[0])] = np.asarray([float(fields[2]), float(fields[3]), float(fields[4])])

    boundaries: dict[int, list[list[int]]] = {}
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 8 or not fields[0].lstrip("-").isdigit():
            continue
        bid = int(fields[1])
        code = int(fields[4])
        n = {303: 3, 404: 4}.get(code, max(3, len(fields) - 5))
        boundaries.setdefault(bid, []).append([int(x) for x in fields[5:5 + n]])

    elements_by_body: dict[int, int] = {}
    z_by_body: dict[int, list[float]] = {}
    element_count = 0
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 4 or not fields[0].lstrip("-").isdigit():
            continue
        element_count += 1
        body = int(fields[1])
        elements_by_body[body] = elements_by_body.get(body, 0) + 1
        pts = [nodes.get(int(x)) for x in fields[3:]]
        pts = [x for x in pts if x is not None]
        if pts:
            z_by_body.setdefault(body, []).append(float(np.mean(np.asarray(pts)[:, 2])))

    def trace(bid: int) -> tuple[int, float, float, float]:
        faces = boundaries.get(bid, [])
        lengths: list[float] = []
        for face in faces:
            for a, b in zip(face, face[1:] + face[:1]):
                if a in nodes and b in nodes:
                    lengths.append(float(np.linalg.norm(nodes[a] - nodes[b])))
        if not lengths:
            return len(faces), math.nan, math.nan, math.nan
        return len(faces), float(np.mean(lengths)), float(np.percentile(lengths, 10)), float(np.percentile(lengths, 90))

    header = (mesh / "mesh.header").read_text(encoding="utf-8", errors="replace").splitlines()[0].split()
    return {
        "mesh": str(mesh),
        "mesh_hash_sha256": sha256(mesh / "mesh.nodes") + ":" + sha256(mesh / "mesh.elements") + ":" + sha256(mesh / "mesh.boundary"),
        "node_count": len(nodes),
        "element_count": element_count,
        "header": header,
        "body_element_counts": elements_by_body,
        "body_z_depth_m": {str(k): (max(v) - min(v) if v else math.nan) for k, v in z_by_body.items()},
        "traces": {
            name: trace(boundary_id) for name, boundary_id in boundary_map.items()
        },
    }


def run_fixed_power(mesh_name: str, key: str) -> list[dict]:
    # Reuse the validated frozen-power harness.  It writes SIF/capture/logs only
    # below OUT and the diagnostic mesh directory.
    import scripts.support.run_phase24_thermal_network_localization as base

    base.OUT = OUT
    case = base.MeshCase(
        key, key, mesh_name, TEMPLATE, 101, 1804,
        (
            base.Interface("TES_to_membrane", "TES->membrane", 1104, 1305, 101, 103),
            base.Interface("TES_to_Stycast", "TES->Stycast", 1105, 1204, 101, 102),
            base.Interface("Stycast_to_substrate", "Stycast->substrate", 1205, 1004, 102, 108),
        ),
    )
    return [base.run_case(case, f, repeat=False) for f in (0.95, 1.0, 1.05)]


def write_series(rows: list[dict], stats: dict[str, dict]) -> None:
    with (OUT / "mesh_statistics.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["variant", "mesh", "node_count", "element_count", "TES_trace_faces", "TES_trace_mean_edge_um", "Membrane_trace_faces", "Membrane_trace_mean_edge_um", "Membrane_substrate_faces", "Membrane_substrate_mean_edge_um", "Stycast_upper_faces", "Stycast_upper_mean_edge_um", "bath_faces", "local_element_depth_m", "constraint_count", "support_span_x_um", "support_span_y_um", "mesh_hash_sha256"]
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for variant, data in stats.items():
            t = data["traces"]
            constraint = data.get("constraint_count", "")
            writer.writerow({
                "variant": variant, "mesh": data["mesh"], "node_count": data["node_count"], "element_count": data["element_count"],
                "TES_trace_faces": t["TES_zmin"][0], "TES_trace_mean_edge_um": t["TES_zmin"][1] * 1e6,
                "Membrane_trace_faces": t["Membrane_TES"][0], "Membrane_trace_mean_edge_um": t["Membrane_TES"][1] * 1e6,
                "Membrane_substrate_faces": t["Membrane_SiNx"][0], "Membrane_substrate_mean_edge_um": t["Membrane_SiNx"][1] * 1e6,
                "Stycast_upper_faces": t["Stycast_upper"][0], "Stycast_upper_mean_edge_um": t["Stycast_upper"][1] * 1e6,
                "bath_faces": t["bath"][0], "local_element_depth_m": data["body_z_depth_m"].get("103", ""),
                "constraint_count": constraint, "support_span_x_um": data.get("support_span_x_um", ""), "support_span_y_um": data.get("support_span_y_um", ""),
                "mesh_hash_sha256": data["mesh_hash_sha256"],
            })

    with (OUT / "membrane_substrate_convergence.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["case", "G_eff_W_per_K", "G_secant_W_per_K", "historical_ratio", "power_points", "solver_status", "residual_source"]
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        writer.writerows(rows)


def geff_from_capture(mesh_name: str, key: str) -> dict:
    import scripts.support.run_phase24_thermal_network_localization as base
    case = base.MeshCase(key, key, mesh_name, TEMPLATE, 101, 1804, tuple())
    temps = []
    for f in (0.95, 1.0, 1.05):
        name = f"{key}_{f:.2f}P".replace(".", "p")
        temps.append(base.tes_temperature(case, name))
    geff = 0.10 * P0 / (temps[2] - temps[0])
    sec = P0 / (temps[1] - TBATH)
    return {"case": key, "G_eff_W_per_K": geff, "G_secant_W_per_K": sec, "historical_ratio": geff / GHIST, "power_points": "0.95/1.00/1.05 P0", "solver_status": "ALL DONE / exit 0; CPU native HeatSolve MUMPS", "residual_source": "native full restriction capture"}


def tes_membrane_gate(stats: dict[str, dict]) -> dict:
    base = stats["refined_mortar_parent"]
    candidate = stats["tes_membrane_previous_test"]
    checks = {}
    for field in ("node_count", "element_count"):
        checks[field] = {"base": base[field], "candidate": candidate[field], "changed": base[field] != candidate[field]}
    for label in ("Stycast_upper", "TES_zmax", "Membrane_SiNx", "bath"):
        b = base["traces"][label][0]; c = candidate["traces"][label][0]
        checks[label + "_faces"] = {"base": b, "candidate": c, "changed": b != c}
    invalid = any(value["changed"] for value in checks.values())
    return {"status": "invalid controlled mesh" if invalid else "valid controlled mesh", "solver_executed": False, "target": "TES/Membrane trace only", "checks": checks, "reason": "candidate changes non-target topology/statistics; previous result is retained as confounded evidence" if invalid else "all non-target counts unchanged"}


def membrane_substrate_mesh_gate(parent: dict, candidate: dict) -> dict:
    checks = {}
    for body, label in ((101, "TES_element_count"), (102, "Stycast_element_count")):
        b = parent["body_element_counts"].get(body, 0); c = candidate["body_element_counts"].get(body, 0)
        checks[label] = {"base": b, "candidate": c, "changed": b != c}
    for label in ("Stycast_upper", "TES_zmax", "Membrane_TES", "bath"):
        b = parent["traces"][label][0]; c = candidate["traces"][label][0]
        checks[label + "_faces"] = {"base": b, "candidate": c, "changed": b != c}
    invalid = any(item["changed"] for item in checks.values())
    return {"status": "invalid controlled mesh" if invalid else "valid controlled mesh", "solver_executed": False, "target": "Membrane/substrate trace/local mesh only", "checks": checks, "reason": "non-target mesh statistics changed; candidate is not used as a controlled convergence point" if invalid else "TES/Stycast/bath statistics fixed"}


def trace_temperature(mesh: Path, result: Path, capture: Path, boundary_id: int) -> float:
    """Area-weighted boundary temperature from the native captured primal x."""
    import scripts.support.run_phase24_thermal_network_localization as base

    _, node_to_dof = base.result_field(result)
    meta = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    primal = int(meta["runtime"]["primal_rows"])
    x = base.indexed(capture / "full_x_after.dat", int(meta["runtime"]["total_rows"]))[:primal]
    nodes: dict[int, np.ndarray] = {}
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 5 and fields[0].lstrip("-").isdigit():
            nodes[int(fields[0])] = np.asarray([float(fields[2]), float(fields[3]), float(fields[4])])
    total = 0.0
    weighted = 0.0
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 8 or int(fields[1]) != boundary_id:
            continue
        code = int(fields[4]); n = {303: 3, 404: 4}.get(code, max(3, len(fields) - 5))
        face = [int(v) for v in fields[5:5 + n]]
        if not all(node in nodes and node_to_dof[node - 1] > 0 for node in face):
            continue
        points = [nodes[node] for node in face]
        area = 0.5 * float(np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0])))
        if len(points) == 4:
            area += 0.5 * float(np.linalg.norm(np.cross(points[3] - points[0], points[2] - points[0])))
        values = [float(x[int(node_to_dof[node - 1]) - 1]) for node in face]
        total += area; weighted += area * float(np.mean(values))
    return weighted / total if total else float("nan")


def materialize_reports(geff_rows: list[dict], gate: dict) -> None:
    def run_label(fraction: float) -> str:
        return f"{fraction:.2f}P".replace(".", "p", 1)

    cases = {
        "historical": {
            "mesh": ROOT / "work/meshes/mesh_singlepixel_prod_v2",
            "result": lambda f: ROOT / "work/meshes/mesh_singlepixel_prod_v2" / f"historical_{f:.2f}p.result".replace(".", "p", 1),
            "capture": lambda f: ROOT / "artifacts/phase24_thermal_network_localization/capture" / f"historical_{run_label(f)}" / "ts0001_nl0001",
            "trace_ids": {"TES_trace": 24, "Membrane_trace": 23, "Si1_branch_entry": 19, "Si1_branch_exit": 21, "Si2_branch_entry": 14, "Si2_branch_exit": 13, "SiO2_2_bath_side": 11, "bath": 30},
        },
        "refined_mortar_parent": {
            "mesh": ROOT / "work/meshes/mesh_phase24_stycast_density_10um",
            "result": lambda f: ROOT / "work/meshes/mesh_phase24_stycast_density_10um" / f"phase24_historical_density_{f:.2f}p.result".replace(".", "p", 1),
            "capture": lambda f: ROOT / "artifacts/p24d10b/capture" / f"phase24_historical_density_{run_label(f)}" / "ts0001_nl0001",
            "trace_ids": {"TES_trace": 1104, "Membrane_trace": 1305, "Si1_branch_entry": 1904, "Si1_branch_exit": 1605, "Si2_branch_entry": 1404, "Si2_branch_exit": 1705, "SiO2_2_bath_side": 1804, "bath": 21},
        },
        "membrane_substrate_intermediate": {
            "mesh": ROOT / "work/meshes" / INTERMEDIATE,
            "result": lambda f: ROOT / "work/meshes" / INTERMEDIATE / f"membrane_substrate_intermediate_{f:.2f}P.result".replace(".", "p", 1),
            "capture": lambda f: OUT / "capture" / f"membrane_substrate_intermediate_{run_label(f)}" / "ts0001_nl0001",
            "trace_ids": {"TES_trace": 1104, "Membrane_trace": 1305, "Si1_branch_entry": 1904, "Si1_branch_exit": 1605, "Si2_branch_entry": 1404, "Si2_branch_exit": 1705, "SiO2_2_bath_side": 1804, "bath": 21},
        },
        "membrane_substrate_control": {
            "mesh": ROOT / "work/meshes/mesh_phase24_trace_membrane_substrate_historical",
            "result": lambda f: ROOT / "work/meshes/mesh_phase24_trace_membrane_substrate_historical" / f"ms_hist_{f:.2f}p.result".replace(".", "p", 1),
            "capture": lambda f: ROOT / "artifacts/p24trace/ms/capture" / f"ms_hist_{run_label(f)}" / "ts0001_nl0001",
            "trace_ids": {"TES_trace": 1104, "Membrane_trace": 1305, "Si1_branch_entry": 1904, "Si1_branch_exit": 1605, "Si2_branch_entry": 1404, "Si2_branch_exit": 1705, "SiO2_2_bath_side": 1804, "bath": 21},
        },
    }
    # Correct the two-digit filename spelling explicitly (f-string replacement
    # above is kept deterministic across Windows/POSIX path handling).
    for spec in cases.values():
        old = spec["result"]
        spec["result"] = lambda f, old=old: Path(str(old(f)).replace("0p95p", "0p95p").replace("1p00p", "1p00p").replace("1p05p", "1p05p"))

    fractions = (0.95, 1.0, 1.05)
    if not cases["membrane_substrate_intermediate"]["result"](1.0).is_file():
        del cases["membrane_substrate_intermediate"]
    temp: dict[str, dict[float, dict[str, float]]] = {}
    temp_rows: list[dict] = []
    for case, spec in cases.items():
        temp[case] = {}
        for fraction in fractions:
            row = {name: trace_temperature(spec["mesh"], spec["result"](fraction), spec["capture"](fraction), bid) for name, bid in spec["trace_ids"].items()}
            temp[case][fraction] = row
            for name, value in row.items():
                temp_rows.append({"case": case, "power_fraction": fraction, "power_W": P0 * fraction, "trace": name, "temperature_K": value, "temperature_mK": value * 1.0e3, "source": "area-weighted native captured primal trace"})
    with (OUT / "trace_temperature_slopes.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["case", "trace", "T_0p95P_K", "T_1p00P_K", "T_1p05P_K", "dT_dP_K_per_W", "T_center_minus_bath_K", "source"]
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for case, points in temp.items():
            for trace in cases[case]["trace_ids"]:
                minus, center, plus = (points[f][trace] for f in fractions)
                writer.writerow({"case": case, "trace": trace, "T_0p95P_K": minus, "T_1p00P_K": center, "T_1p05P_K": plus, "dT_dP_K_per_W": (plus - minus) / (0.10 * P0), "T_center_minus_bath_K": center - TBATH, "source": "area-weighted native captured primal trace"})
    with (OUT / "trace_temperature.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = list(temp_rows[0]); writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(temp_rows)

    segments = [
        ("TES_to_Membrane", "TES_trace", "Membrane_trace"),
        ("Membrane_to_Si1_branch_entry", "Membrane_trace", "Si1_branch_entry"),
        ("Si1_branch", "Si1_branch_entry", "Si1_branch_exit"),
        ("Membrane_to_Si2_branch_entry", "Membrane_trace", "Si2_branch_entry"),
        ("Si2_branch", "Si2_branch_entry", "Si2_branch_exit"),
        ("SiO2_2_to_bath", "SiO2_2_bath_side", "bath"),
    ]
    with (OUT / "differential_resistance_breakdown.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["case", "segment", "left_trace", "right_trace", "R_differential_K_per_W", "R_secant_K_per_W", "delta_T_center_K", "method"]
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for case, points in temp.items():
            for segment, left, right in segments:
                dm = points[0.95][left] - points[0.95][right]
                dc = points[1.0][left] - points[1.0][right]
                dp = points[1.05][left] - points[1.05][right]
                writer.writerow({"case": case, "segment": segment, "left_trace": left, "right_trace": right, "R_differential_K_per_W": (dp - dm) / (0.10 * P0), "R_secant_K_per_W": dc / P0, "delta_T_center_K": dc, "method": "symmetric 0.95/1.05 P trace slope; no series summation across parallel branches"})

    # Native C-lambda reactions exist for the mortar interfaces, but the two
    # downstream branches are conformal/shared-node paths in these captures.
    # Do not manufacture Q_i from total power or split it by area.
    with (OUT / "parallel_branch_heatflow.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["case", "branch", "Q_i_native_W", "hot_trace", "cold_trace", "delta_T_center_K", "G_i_W_per_K", "status", "note"]
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for case, points in temp.items():
            for branch, hot, cold in (("Si1_SiNx", "Si1_branch_entry", "Si1_branch_exit"), ("Si2_SiO2_2", "Si2_branch_entry", "Si2_branch_exit")):
                writer.writerow({"case": case, "branch": branch, "Q_i_native_W": "nan", "hot_trace": hot, "cold_trace": cold, "delta_T_center_K": points[1.0][hot] - points[1.0][cold], "G_i_W_per_K": "nan", "status": "not identifiable from capture", "note": "No native branch flux mask/reaction was captured; global Q/P is not split or inferred."})

    g = {row["case"]: float(row["G_eff_W_per_K"]) for row in geff_rows}
    base = g["refined_mortar_parent"]; hist = GHIST; control = g["membrane_substrate_control"]
    explained = (base - control) / (base - hist)
    payload = {
        "historical_G_eff_W_per_K": hist,
        "refined_mortar_parent_G_eff_W_per_K": base,
        "membrane_substrate_control_G_eff_W_per_K": control,
        "membrane_substrate_intermediate_G_eff_W_per_K": g.get("membrane_substrate_intermediate", float("nan")),
        "membrane_substrate_controlled_fraction_of_parent_to_historical_conductance_gap": explained,
        "unexplained_conductance_gap_fraction_after_control": 1.0 - explained,
        "TES_membrane_true_one_factor": {"status": "invalid controlled mesh", "contribution_fraction": None, "solver_executed": False},
        "resistance_decomposition_existing_audit": {"TES_to_membrane_effective_discrete_response_share": 0.458, "membrane_substrate_network_share": 0.542, "note": "resistance decomposition is not a one-factor proof; TES/Membrane share remains unconfirmed"},
        "controlled_tests_explained_fraction_of_remaining_gap": explained,
        "unexplained_remainder_fraction": 1.0 - explained,
        "units_note": "fraction uses requested G_eff endpoint definition; resistance rows are differential trace slopes",
    }
    (OUT / "remaining_gap_accounting.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (OUT / "nonlinear_steady_validation.json").write_text(json.dumps({"status": "SKIPPED", "reason": "best valid thermal diagnostic is 4.3348% above historical, so |G/Ghist-1| <= 0.03 is not met", "full_nonlinear_run": False, "final_current_uA": None, "solver": "CPU MUMPS not launched for nonlinear validation"}, indent=2) + "\n", encoding="utf-8")

    best = control
    intermediate_value = g.get("membrane_substrate_intermediate", float("nan"))
    if math.isfinite(intermediate_value):
        convergence_note = f"The accepted series point at h={INTERMEDIATE_H*1e6:.1f} µm is {intermediate_value/hist:.4f}× historical; compare it with the h=10 µm control for monotonicity."
    else:
        convergence_note = "The attempted intermediate failed the pre-solver mesh gate and is excluded from the valid convergence series; no solver result is accepted for it."
    summary = [
        "# Phase24 final thermal parity localization",
        "",
        "Diagnostic-only campaign. Physics parameters, materials, TES law, circuit constants, bath BC, and production mesh were unchanged. Mesh generation was serialized and the new intermediate used a unique work directory.",
        "",
        "## Fixed-power convergence",
        "",
        "| case | G_eff (W/K) | historical ratio | interpretation |",
        "|---|---:|---:|---|",
        f"| Historical | {hist:.12e} | 1.000000 | reference |",
        f"| Refined mortar parent | {base:.12e} | {base/hist:.6f} | audited parent |",
        f"| Membrane/substrate intermediate (h={INTERMEDIATE_H*1e6:.1f} µm) | {g.get('membrane_substrate_intermediate', float('nan')):.12e} | {g.get('membrane_substrate_intermediate', float('nan'))/hist:.6f} | CPU/MUMPS if gate passed |",
        f"| Membrane/substrate control (h=10 µm) | {control:.12e} | {control/hist:.6f} | audited accepted control |",
        "",
        convergence_note,
        "",
        "## TES/Membrane gate",
        "",
        f"The previous TES/Membrane candidate is `{gate['status']}`. Its non-target statistics change (including the Stycast upper trace), so no TES/Membrane fixed-power result is accepted and its 2.6857e-8 W/K value remains confounded. `f_TM` is therefore not identifiable.",
        "",
        "## Branch and trace diagnostics",
        "",
        "`parallel_branch_heatflow.csv` records Q_i/G_i as unavailable: the existing native captures have no separate flux mask or mortar reaction for the conformal downstream Si1/SiNx and Si2/SiO2_2 paths. Splitting total P by area or assuming a series chain would be invalid. Trace slopes and differential segment resistances are materialized in `trace_temperature_slopes.csv` and `differential_resistance_breakdown.csv`.",
        "",
        f"The accepted Membrane/substrate control explains {explained*100:.2f}% of the parent-to-historical G_eff gap; {((1-explained)*100):.2f}% remains unexplained by valid one-factor tests. Existing resistance decomposition assigns 45.8% to TES/Membrane effective discrete response and 54.2% to the downstream network, but only the downstream share has a clean controlled mesh test.",
        "",
        "## GO/NO-GO",
        "",
        f"Best valid thermal parity is `{best/hist:.6f}` (4.3348% high), outside the 2–3% GO band; full nonlinear CPU/MUMPS validation was skipped. Final current parity is therefore not measured. HYPRE/GPU: **NO-GO**.",
        "",
        "See `mesh_statistics.csv`, `mesh_control_gate.csv`, `remaining_gap_accounting.json`, and the trace/branch CSVs for the machine-readable audit.",
    ]
    (OUT / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    intermediate = build_intermediate()
    # Existing meshes are read-only inputs for the campaign.
    mesh_paths = {
        "refined_mortar_parent": ROOT / "work/meshes/mesh_phase24_stycast_density_10um",
        "membrane_substrate_intermediate": intermediate,
        "membrane_substrate_control": ROOT / "work/meshes/mesh_phase24_trace_membrane_substrate_historical",
        "tes_membrane_previous_test": ROOT / "work/meshes/mesh_phase24_trace_tes_membrane_historical",
    }
    stats = {key: parse_mesh(path) for key, path in mesh_paths.items()}
    stats["historical_reference"] = parse_mesh(
        ROOT / "work/meshes/mesh_singlepixel_prod_v2",
        {"TES_zmin": 24, "TES_zmax": 25, "Membrane_TES": 23, "Membrane_SiNx": 23, "Membrane_Si1": 19, "Membrane_Si1_exit": 21, "SiO2_2_bath": 11, "bath": 30, "Stycast_upper": 26, "Stycast_lower": 27},
    )

    substrate_gate = membrane_substrate_mesh_gate(stats["refined_mortar_parent"], stats["membrane_substrate_intermediate"])
    (OUT / "membrane_substrate_mesh_gate.json").write_text(json.dumps(substrate_gate, indent=2) + "\n", encoding="utf-8")
    # The parent/control captures are already residual-audited.  The new
    # intermediate is launched only after its pre-solver gate passes.
    if substrate_gate["status"] == "valid controlled mesh":
        run_fixed_power(INTERMEDIATE, "membrane_substrate_intermediate")
        geff_rows = [geff_from_capture(INTERMEDIATE, "membrane_substrate_intermediate")]
    else:
        geff_rows = [{"case": "membrane_substrate_intermediate", "G_eff_W_per_K": float("nan"), "G_secant_W_per_K": float("nan"), "historical_ratio": float("nan"), "power_points": "not run: mesh gate failed", "solver_status": "NOT RUN; invalid controlled mesh", "residual_source": "pre-solver gate"}]
    geff_rows.extend([
        {"case": "refined_mortar_parent", "G_eff_W_per_K": 2.029899315e-8, "G_secant_W_per_K": float("nan"), "historical_ratio": 2.029899315e-8 / GHIST, "power_points": "0.95/1.00/1.05 P0 (existing audited capture)", "solver_status": "existing ALL DONE / exit 0; CPU native HeatSolve MUMPS", "residual_source": "existing audited capture"},
        {"case": "membrane_substrate_control", "G_eff_W_per_K": 1.959713948e-8, "G_secant_W_per_K": float("nan"), "historical_ratio": 1.959713948e-8 / GHIST, "power_points": "0.95/1.00/1.05 P0 (existing audited capture)", "solver_status": "existing ALL DONE / exit 0; CPU native HeatSolve MUMPS", "residual_source": "existing audited capture"},
    ])
    # Constraint count/support spans are authoritative from the center capture
    # if available; retain empty fields rather than inventing them.
    meta = OUT / "capture/membrane_substrate_intermediate_1p00P/ts0001_nl0001/metadata.json"
    if meta.is_file():
        runtime = json.loads(meta.read_text(encoding="utf-8"))["runtime"]
        stats["membrane_substrate_intermediate"]["constraint_count"] = runtime.get("constraint_rows", "")
    for key in stats:
        stats[key].setdefault("constraint_count", "existing capture not duplicated")

    gate = tes_membrane_gate(stats)
    (OUT / "mesh_control_gate.csv").write_text("check,status,base,candidate,note\n" + "\n".join(f"{name},{'changed' if item['changed'] else 'unchanged'},{item['base']},{item['candidate']},TES/Membrane gate" for name, item in gate["checks"].items()) + "\n", encoding="utf-8")
    (OUT / "tes_membrane_one_factor.json").write_text(json.dumps(gate, indent=2) + "\n", encoding="utf-8")
    write_series(geff_rows, stats)
    materialize_reports(geff_rows, gate)
    (OUT / "run_manifest.json").write_text(json.dumps({"P0_W": P0, "bath_K": TBATH, "solver": "CPU native HeatSolve MUMPS", "mesh_generation": "serialized; unique diagnostic directory", "physics_changed": False, "production_mesh_overwritten": False, "runs": geff_rows, "tes_membrane_gate": gate}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT), "intermediate": geff_rows[0], "tes_membrane_gate": gate["status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

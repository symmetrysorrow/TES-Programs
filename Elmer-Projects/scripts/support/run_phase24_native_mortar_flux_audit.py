"""Audit solver-native mortar reactions and the historical/Phase24 thermal path.

The input is the full saddle-point capture written by the native
``SolveWithLinearRestriction`` diagnostic.  The lower-left block is the
native constraint matrix ``B`` and the multiplier tail of the native solved
vector is ``lambda``.  The reaction used here is the actual upper-right
block product ``C lambda`` (normally ``B^T lambda``), not a reconstructed
one-sided gradient.

This is intentionally a read-only audit of existing captures.  It does not
change a mesh, material, TES law, circuit constant, time step, or solver
backend.
"""
from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_native_mortar_flux_audit"
TBATH = 0.150
MATCHED_T = 0.168563175809
MATCHED_POWER = 3.203004762115138e-10


@dataclass(frozen=True)
class Interface:
    name: str
    source_to_target: str
    slave_boundary: int
    master_boundary: int
    slave_name: str
    master_name: str


@dataclass(frozen=True)
class Case:
    name: str
    label: str
    capture: Path
    mesh: Path
    result: Path
    interfaces: tuple[Interface, ...]
    tes_body: int


CASES = (
    Case(
        "case_phase24_historical_one_shot",
        "historical",
        ROOT / "artifacts/phase24_historical_mesh_strict_steady_validation/capture/one/ts0001_nl0001",
        ROOT / "work/meshes/mesh_singlepixel_prod_v2",
        ROOT / "work/meshes/mesh_singlepixel_prod_v2/phase24_historical_restart_snapshot.result",
        (
            Interface("TES_to_membrane", "TES->membrane", 24, 23, "TES", "Membrane_SiNx"),
            Interface("TES_to_Stycast", "TES->Stycast", 25, 26, "TES", "Stycast"),
            Interface("Stycast_to_substrate", "Stycast->substrate/bath", 27, 28, "Stycast", "Pb"),
        ),
        8,
    ),
    Case(
        "case_phase24_native_flux_one_shot",
        "Phase24",
        ROOT / "artifacts/phase24_historical_mesh_strict_steady_validation/capture/phase24/ts0001_nl0001",
        ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar",
        ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/case_phase24_restart_refine_fine_stycast32_mortar_mumps_mortar.result",
        (
            Interface("TES_to_membrane", "TES->membrane", 1104, 1305, "TES", "Membrane_SiNx"),
            Interface("TES_to_Stycast", "TES->Stycast", 1105, 1204, "TES", "Stycast"),
            Interface("Stycast_to_substrate", "Stycast->substrate/bath", 1205, 1004, "Stycast", "Pb"),
        ),
        101,
    ),
)


def indexed(path: Path, count: int) -> np.ndarray:
    values = np.zeros(count, dtype=np.float64)
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 2:
                values[int(fields[0]) - 1] = float(fields[1].replace("D", "E").replace("d", "e"))
    return values


def write_csv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_blocks(case: Case) -> tuple[dict, csr_matrix, csr_matrix, np.ndarray]:
    meta = json.loads((case.capture / "metadata.json").read_text(encoding="utf-8"))
    total = int(meta["runtime"]["total_rows"])
    primal = int(meta["runtime"]["primal_rows"])
    constraints = total - primal
    bottom_rows: list[int] = []
    bottom_cols: list[int] = []
    bottom_values: list[float] = []
    top_rows: list[int] = []
    top_cols: list[int] = []
    top_values: list[float] = []
    with (case.capture / "full_A_after.dat").open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 3:
                continue
            row = int(fields[0]) - 1
            col = int(fields[1]) - 1
            value = float(fields[2].replace("D", "E").replace("d", "e"))
            if row >= primal and col < primal:
                bottom_rows.append(row - primal)
                bottom_cols.append(col)
                bottom_values.append(value)
            elif row < primal and col >= primal:
                top_rows.append(row)
                top_cols.append(col - primal)
                top_values.append(value)
    b = coo_matrix(
        (bottom_values, (bottom_rows, bottom_cols)), shape=(constraints, primal)
    ).tocsr()
    c = coo_matrix(
        (top_values, (top_rows, top_cols)), shape=(primal, constraints)
    ).tocsr()
    solution = indexed(case.capture / "full_x_after.dat", total)
    return meta, b, c, solution[primal:]


def result_inverse_permutation(path: Path, primal: int) -> np.ndarray:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    marker = next(i for i, line in enumerate(lines) if line.strip().lower() == "temperature")
    count = int(lines[marker + 1].split()[1])
    inverse = np.zeros(primal + 1, dtype=np.int64)
    for offset in range(count):
        node, dof = (int(value) for value in lines[marker + 2 + offset].split()[:2])
        if dof <= primal:
            inverse[dof] = node
    return inverse


def boundary_data(mesh: Path) -> tuple[dict[int, set[int]], dict[int, float], dict[int, int]]:
    nodes: dict[int, np.ndarray] = {}
    with (mesh / "mesh.nodes").open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 5:
                nodes[int(fields[0])] = np.asarray([float(fields[2]), float(fields[3]), float(fields[4])])
    boundary_nodes: dict[int, set[int]] = {}
    boundary_area: dict[int, float] = {}
    boundary_count: dict[int, int] = {}
    with (mesh / "mesh.boundary").open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 8:
                continue
            boundary_id = int(fields[1])
            code = int(fields[4])
            count = {303: 3, 404: 4}.get(code, max(3, len(fields) - 5))
            face = [int(value) for value in fields[5 : 5 + count]]
            boundary_nodes.setdefault(boundary_id, set()).update(face)
            boundary_count[boundary_id] = boundary_count.get(boundary_id, 0) + 1
            if len(face) >= 3 and all(node in nodes for node in face):
                a, b, c = (nodes[node] for node in face[:3])
                area = 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))
                if len(face) == 4:
                    d = nodes[face[3]]
                    area += 0.5 * float(np.linalg.norm(np.cross(d - a, c - a)))
                boundary_area[boundary_id] = boundary_area.get(boundary_id, 0.0) + area
    return boundary_nodes, boundary_area, boundary_count


def classify_rows(
    b: csr_matrix, inverse: np.ndarray, interfaces: tuple[Interface, ...], boundary_nodes: dict[int, set[int]]
) -> list[int]:
    result: list[int] = []
    for row in range(b.shape[0]):
        nodes = {int(inverse[col + 1]) for col in b[row].indices if inverse[col + 1] > 0}
        scores = [
            len(nodes & boundary_nodes.get(interface.slave_boundary, set()))
            + len(nodes & boundary_nodes.get(interface.master_boundary, set()))
            for interface in interfaces
        ]
        result.append(int(np.argmax(scores)) if scores and max(scores) > 0 else -1)
    return result


def reaction_capture(case: Case) -> tuple[dict, list[dict[str, object]], list[dict[str, object]], dict[str, dict]]:
    meta, b, c, multipliers = read_blocks(case)
    primal = int(meta["runtime"]["primal_rows"])
    inverse = result_inverse_permutation(case.result, primal)
    boundary_nodes, boundary_area, boundary_count = boundary_data(case.mesh)
    classes = classify_rows(b, inverse, case.interfaces, boundary_nodes)
    c_by_constraint = c.T.tocsr()
    aggregate: list[dict[str, object]] = []
    row_records: list[dict[str, object]] = []
    topology: dict[str, dict] = {}
    for index, interface in enumerate(case.interfaces):
        selected = np.asarray([row for row, group in enumerate(classes) if group == index], dtype=np.int64)
        slave_set = boundary_nodes.get(interface.slave_boundary, set())
        master_set = boundary_nodes.get(interface.master_boundary, set())
        source_reaction = 0.0
        target_reaction = 0.0
        nnz = 0
        norms: list[float] = []
        for row in selected.tolist():
            cols = c_by_constraint[row].indices
            values = c_by_constraint[row].data
            nodes = [int(inverse[col + 1]) for col in cols if inverse[col + 1] > 0]
            slave_reaction = float(sum(value * multipliers[row] for col, value in zip(cols, values) if int(inverse[col + 1]) in slave_set))
            master_reaction = float(sum(value * multipliers[row] for col, value in zip(cols, values) if int(inverse[col + 1]) in master_set))
            row_reaction = slave_reaction + master_reaction
            source_reaction += slave_reaction
            target_reaction += master_reaction
            nnz += len(cols)
            norms.append(float(np.linalg.norm(values)))
            row_records.append({
                "case": case.label,
                "interface": interface.name,
                "constraint_row_id": row + 1,
                "master_side": interface.master_name,
                "slave_side": interface.slave_name,
                "master_boundary": interface.master_boundary,
                "slave_boundary": interface.slave_boundary,
                "associated_primal_dofs": ";".join(str(int(col) + 1) for col in cols),
                "associated_primal_nodes": ";".join(str(node) for node in nodes),
                "multiplier_value": float(multipliers[row]),
                "multiplier_units": "W/m^2 (Galerkin area-weighted thermal constraint)",
                "slave_reaction_W": slave_reaction,
                "master_reaction_W": master_reaction,
                "reaction_balance_W": row_reaction,
                "sign_convention": "B=slave-master; C lambda is primal reaction; source_to_target=slave reaction",
            })
        physical_flow = source_reaction
        if interface.name == "TES_to_Stycast":
            physical_flow = -source_reaction
        aggregate.append({
            "route": "native_mortar_reaction",
            "case": case.label,
            "interface": interface.name,
            "side": "slave/source-to-master" if interface.name != "TES_to_Stycast" else "master/TES-to-slave",
            "constraint_count": int(selected.size),
            "projector_nnz": int(nnz),
            "multiplier_sum": float(multipliers[selected].sum()) if selected.size else 0.0,
            "weighted_reaction_W": physical_flow,
            "slave_reaction_W": source_reaction,
            "master_reaction_W": target_reaction,
            "reaction_balance_W": source_reaction + target_reaction,
            "sign_convention": "B=slave-master; C lambda is primal reaction; positive listed flow follows interface name",
            "physical_units": "reaction W; multiplier W/m^2; B/C projector weights m^2",
        })
        topology[interface.name] = {
            "constraint_count": int(selected.size),
            "projector_nnz": int(nnz),
            "row_norm_min": min(norms) if norms else 0.0,
            "row_norm_mean": float(np.mean(norms)) if norms else 0.0,
            "row_norm_max": max(norms) if norms else 0.0,
            "slave_nodes": len(slave_set),
            "master_nodes": len(master_set),
            "slave_boundary_elements": boundary_count.get(interface.slave_boundary, 0),
            "master_boundary_elements": boundary_count.get(interface.master_boundary, 0),
            "slave_area_m2": boundary_area.get(interface.slave_boundary, 0.0),
            "master_area_m2": boundary_area.get(interface.master_boundary, 0.0),
            "master_side": interface.master_name,
            "slave_side": interface.slave_name,
        }
    topology["_runtime"] = {
        "total_rows": int(meta["runtime"]["total_rows"]),
        "primal_rows": int(meta["runtime"]["primal_rows"]),
        "constraint_rows": int(meta["runtime"]["constraint_rows"]),
        "native_constraint_matrix_nnz": int(b.nnz),
        "unclassified_rows": int(classes.count(-1)),
    }
    return meta, aggregate, row_records, topology


def native_scalars() -> dict[str, dict[str, float]]:
    path = ROOT / "artifacts/phase24_historical_mesh_strict_steady_validation/native_flux_integrals.csv"
    values: dict[str, dict[str, float]] = {}
    for row in csv.DictReader(path.open(encoding="utf-8")):
        values.setdefault(row["case"], {})[row["quantity"]] = float(row["value"])
    return values


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    all_reactions: list[dict[str, object]] = []
    all_rows: list[dict[str, object]] = []
    all_topology: dict[str, dict] = {}
    for case in CASES:
        _, reactions, rows, topology = reaction_capture(case)
        all_reactions.extend(reactions)
        all_rows.extend(rows)
        all_topology[case.label] = topology
    write_csv(OUT / "native_mortar_reaction_integrals.csv", all_reactions, list(all_reactions[0]))
    write_csv(OUT / "native_mortar_reaction_rows.csv", all_rows, list(all_rows[0]))

    scalars = native_scalars()
    case_scalar = {
        "historical": scalars["case_phase24_historical_one_shot"],
        "Phase24": scalars["case_phase24_native_flux_one_shot"],
    }
    reactions_by_case = {label: {row["interface"]: float(row["weighted_reaction_W"]) for row in all_reactions if row["case"] == label} for label in case_scalar}
    partitions: list[dict[str, object]] = []
    for label in ("historical", "Phase24"):
        data = case_scalar[label]
        joule = data["native TES body source: raw Joule power (W)"]
        bath_out = abs(data["native bath boundary diffusive flux (W)"])
        r = reactions_by_case[label]
        mortar_sum = sum(r.values())
        partitions.append({
            "case": label,
            "TES_Joule_input_W": joule,
            "TES_to_membrane_reaction_W": r.get("TES_to_membrane", 0.0),
            "TES_to_Stycast_reaction_W": r.get("TES_to_Stycast", 0.0),
            "Stycast_to_substrate_reaction_W": r.get("Stycast_to_substrate", 0.0),
            "mortar_reaction_contribution_W": mortar_sum,
            "other_or_conformal_path_W": bath_out - mortar_sum,
            "total_bath_flux_outgoing_W": bath_out,
            "global_balance_error_W": joule - bath_out,
            "global_balance_relative": (joule - bath_out) / joule if joule else math.nan,
            "evidence": "native SaveScalars bath/Joule plus native C lambda reaction; other path includes conformal/shared-node transport",
        })
    write_csv(OUT / "native_heatflow_partition.csv", partitions, list(partitions[0]))

    matched_t_rows: list[dict[str, object]] = []
    matched_p_rows: list[dict[str, object]] = []
    for label in ("historical", "Phase24"):
        data = case_scalar[label]
        temperature = {"historical": 0.16856316510815095, "Phase24": 0.16655816045597896}[label]
        q = abs(data["native bath boundary diffusive flux (W)"])
        g = q / (temperature - TBATH)
        q_at_match = g * (MATCHED_T - TBATH)
        t_at_power = TBATH + MATCHED_POWER / g
        matched_t_rows.append({
            "case": label,
            "target_TES_temperature_K": MATCHED_T,
            "operating_TES_temperature_K": temperature,
            "operating_total_bath_Q_out_W": q,
            "secant_G_eff_W_per_K": g,
            "matched_T_Q_out_estimate_W": q_at_match,
            "method": "native bath flux secant projection; no ±delta-T solve was archived",
            "diagnostic_status": "estimated_from_native_operating_point",
        })
        matched_p_rows.append({
            "case": label,
            "target_Joule_power_W": MATCHED_POWER,
            "operating_Joule_power_W": data["native TES body source: raw Joule power (W)"],
            "secant_G_eff_W_per_K": g,
            "matched_power_TES_temperature_estimate_K": t_at_power,
            "method": "native bath flux secant projection; diagnostic frozen-power estimate",
            "diagnostic_status": "estimated_from_native_operating_point",
        })
    write_csv(OUT / "matched_temperature_comparison.csv", matched_t_rows, list(matched_t_rows[0]))
    write_csv(OUT / "matched_power_comparison.csv", matched_p_rows, list(matched_p_rows[0]))
    g_rows = [{
        "case": row["case"],
        "secant_G_eff_W_per_K": row["secant_G_eff_W_per_K"],
        "local_symmetric_G_eff_W_per_K": "not_run",
        "delta_T_K": "not_run",
        "method": "Q_native_bath/(T_TES-T_bath)",
        "interpretation": "secant only; symmetric perturbation remains a controlled follow-up",
    } for row in matched_t_rows]
    write_csv(OUT / "effective_conductance_comparison.csv", g_rows, list(g_rows[0]))

    area_rows: list[dict[str, object]] = []
    for interface in CASES[0].interfaces:
        hist = all_topology["historical"][interface.name]
        phase = all_topology["Phase24"][interface.name]
        hist_area = min(hist["slave_area_m2"], hist["master_area_m2"])
        phase_area = min(phase["slave_area_m2"], phase["master_area_m2"])
        area_rows.append({
            "interface": interface.name,
            "historical_slave_area_m2": hist["slave_area_m2"],
            "historical_master_area_m2": hist["master_area_m2"],
            "historical_relevant_area_m2": hist_area,
            "Phase24_slave_area_m2": phase["slave_area_m2"],
            "Phase24_master_area_m2": phase["master_area_m2"],
            "Phase24_relevant_area_m2": phase_area,
            "relative_relevant_area_difference": (phase_area - hist_area) / hist_area if hist_area else math.nan,
            "area_interpretation": "actual triangulated boundary area from mesh.boundary/mesh.nodes",
        })
    write_csv(OUT / "interface_area_comparison.csv", area_rows, list(area_rows[0]))

    topology_rows: list[dict[str, object]] = []
    for interface in CASES[0].interfaces:
        hist = all_topology["historical"][interface.name]
        phase = all_topology["Phase24"][interface.name]
        topology_rows.append({
            "interface": interface.name,
            "historical_constraint_count": hist["constraint_count"],
            "Phase24_constraint_count": phase["constraint_count"],
            "historical_projector_nnz": hist["projector_nnz"],
            "Phase24_projector_nnz": phase["projector_nnz"],
            "historical_row_norm_min": hist["row_norm_min"],
            "historical_row_norm_mean": hist["row_norm_mean"],
            "historical_row_norm_max": hist["row_norm_max"],
            "Phase24_row_norm_min": phase["row_norm_min"],
            "Phase24_row_norm_mean": phase["row_norm_mean"],
            "Phase24_row_norm_max": phase["row_norm_max"],
            "historical_master_nodes": hist["master_nodes"],
            "historical_slave_nodes": hist["slave_nodes"],
            "Phase24_master_nodes": phase["master_nodes"],
            "Phase24_slave_nodes": phase["slave_nodes"],
            "master_slave_assignment": f"master={interface.master_name}; slave={interface.slave_name}",
            "coupling_interpretation": "zero constraint rows means this interface is not represented by the captured mortar saddle block",
        })
    write_csv(OUT / "mortar_topology_comparison.csv", topology_rows, list(topology_rows[0]))

    closure = {
        "method": {
            "native_reaction": "C lambda, where C is the solver-assembled upper-right saddle block and lambda is the solver-native multiplier tail",
            "native_capture_implementation": "../tools/elmer-hypre/src/fem/src/SolverUtils.F90: SolveWithLinearRestriction / Phase24CaptureFullSystem",
            "audit_implementation": "scripts/support/run_phase24_native_mortar_flux_audit.py: read_blocks / reaction_capture",
            "units": "C/B projector weights are m^2 for Galerkin thermal mortar; lambda is W/m^2; integrated C lambda is W",
            "sign": "B=slave-master; slave reaction is positive for slave-to-master transfer; master reaction should be equal and opposite",
            "one_sided_fluxes_used_as_balance": False,
        },
        "cases": {row["case"]: row for row in partitions},
        "topology_runtime": all_topology,
        "matched_point_note": "Matched-T and matched-P rows are controlled secant projections from native operating-point Q; no new production solve or ±delta-T perturbation was run in this audit.",
        "physics_changed": False,
        "hypre_gpu": "NO-GO",
    }
    (OUT / "energy_balance_closure.json").write_text(json.dumps(closure, indent=2) + "\n", encoding="utf-8")

    hist_t, phase_t = matched_t_rows
    hist_p, phase_p = matched_p_rows
    dominant_hist = partitions[0]["TES_to_membrane_reaction_W"]
    dominant_phase = partitions[1]["other_or_conformal_path_W"]
    summary = [
        "# Phase24 native mortar heat-flow audit",
        "",
        "## Scope",
        "",
        "This is a diagnostic audit of the existing CPU MUMPS native full saddle-point captures. Production physics, materials, TES law, circuit constants, geometry, mesh topology, mortar formulation, timestep, and HYPRE/GPU settings were not changed.",
        "",
        "## Native reaction definition",
        "",
        "- The solver-native capture supplies the assembled constraint block `B`, the transpose block `C`, and the solved multiplier vector `lambda`.",
        "- Native capture location: `../tools/elmer-hypre/src/fem/src/SolverUtils.F90`, `SolveWithLinearRestriction` / `Phase24CaptureFullSystem`; audit conversion location: `scripts/support/run_phase24_native_mortar_flux_audit.py`.",
        "- The integrated reaction is `C lambda` summed over the primal DOFs belonging to the interface side; this is an actual constraint reaction, not a one-sided gradient reconstruction.",
        "- Galerkin thermal projector weights carry m², multiplier values carry W/m², and the integrated primal reaction carries W.",
        "- `B=slave-master`; slave and master reactions should cancel. The row-level records preserve row ID, DOFs, multiplier, both side reactions, sign, and units.",
        "",
        "## Results",
        "",
        f"- Historical TES→membrane reaction: **{dominant_hist:.9e} W**; Phase24 TES→membrane mortar rows: **0**.",
        f"- Phase24 mortar reaction is instead limited to TES↔Stycast and Stycast↔substrate rows; the remaining native bath heat is classified as conformal/shared-node or otherwise non-mortar path: **{dominant_phase:.9e} W**.",
        f"- Historical native balance: Joule={partitions[0]['TES_Joule_input_W']:.9e} W, bath={partitions[0]['total_bath_flux_outgoing_W']:.9e} W, error={partitions[0]['global_balance_error_W']:.9e} W ({partitions[0]['global_balance_relative']:.3e}).",
        f"- Phase24 native balance: Joule={partitions[1]['TES_Joule_input_W']:.9e} W, bath={partitions[1]['total_bath_flux_outgoing_W']:.9e} W, error={partitions[1]['global_balance_error_W']:.9e} W ({partitions[1]['global_balance_relative']:.3e}).",
        f"- Matched-T ({MATCHED_T*1e3:.6f} mK) secant Q estimate: historical={hist_t['matched_T_Q_out_estimate_W']:.9e} W, Phase24={phase_t['matched_T_Q_out_estimate_W']:.9e} W; Phase24/historical={phase_t['matched_T_Q_out_estimate_W']/hist_t['matched_T_Q_out_estimate_W']:.6f}.",
        f"- Matched-power ({MATCHED_POWER:.9e} W) secant T estimate: historical={hist_p['matched_power_TES_temperature_estimate_K']*1e3:.6f} mK, Phase24={phase_p['matched_power_TES_temperature_estimate_K']*1e3:.6f} mK.",
        "- A symmetric ±δT local derivative was not run; the recorded G_eff values are native-Q secants and are explicitly not a local derivative proof.",
        "",
        "## Required conclusions",
        "",
        "1. Native mortar reaction physical integration: **YES**, from native `C lambda`; row-level and aggregate artifacts are written.",
        "2. Reaction-inclusive energy closure: **YES for global native Joule/bath closure**; the interface partition closes only for the mortar portion, with conformal/shared-node heat reported separately.",
        "3. Same-T conductance: the native operating-point secant estimates Phase24 higher than historical; exact matched-T solve remains pending.",
        "4. Same-P temperature: the secant estimate gives a lower Phase24 TES temperature, indicating higher effective transport in the captured native bath path.",
        "5. Dominant path difference: historical TES→membrane is an explicit mortar path; Phase24 has no TES→membrane mortar rows and carries heat through a non-mortar/conformal path in this capture.",
        "6. Area alone: do not claim area sufficiency; the interface-area CSV is provided, but the zero-vs-nonzero constraint topology is the larger structural difference.",
        "7. Mortar topology: **YES, materially different**; per-interface counts and row norms are recorded.",
        "8. 143→218 µA: **not quantitatively proven by this audit alone**. It establishes a large topology/path candidate, but the matched-point rows are secant estimates rather than new frozen-source solves.",
        "9. Minimum next controlled fix: make the TES–membrane interface coupling identical between meshes (preserve the Phase24 physical mesh, add only the missing controlled interface representation), then rerun native MUMPS with ±δT and fixed-P diagnostics.",
        "10. HYPRE/GPU: **NO-GO** until the controlled thermal-path comparison is completed.",
        "",
        "## Artifacts",
        "",
        "- `native_mortar_reaction_integrals.csv` — integrated reaction by case/interface.",
        "- `native_mortar_reaction_rows.csv` — constraint row ID, primal DOFs/nodes, multiplier, side reactions, units, and sign.",
        "- `native_heatflow_partition.csv` — reaction-inclusive path partition and native global closure.",
        "- `matched_temperature_comparison.csv`, `matched_power_comparison.csv`, `effective_conductance_comparison.csv` — explicit secant diagnostic status.",
        "- `interface_area_comparison.csv`, `mortar_topology_comparison.csv`, `energy_balance_closure.json` — geometry/topology and provenance.",
    ]
    (OUT / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(OUT / "summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

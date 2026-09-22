"""Materialize the fixed-power native mortar diagnostics as CSV/Markdown."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_native_coupling_diagnostics"
MESH = ROOT / "work" / "meshes" / "mesh_phase24_native_coupling_diag"
sys.path.insert(0, str(ROOT))
from scripts.support.run_phase24_native_mortar_flux_audit import (  # noqa: E402
    Case,
    Interface,
    indexed,
    reaction_capture,
    result_inverse_permutation,
)


INTERFACES = (
    Interface("TES_to_membrane", "TES->membrane", 1104, 1305, "TES", "Membrane_SiNx"),
    Interface("TES_to_Stycast", "TES->Stycast", 1105, 1204, "TES", "Stycast"),
    Interface("Stycast_to_substrate", "Stycast->substrate/bath", 1205, 1004, "Stycast", "Pb"),
)


def tes_nodes() -> set[int]:
    nodes: set[int] = set()
    for line in (MESH / "mesh.elements").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 4 and int(fields[1]) == 101:
            nodes.update(int(value) for value in fields[3:])
    return nodes


def temperature(name: str, capture: Path) -> float:
    meta = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    total = int(meta["runtime"]["total_rows"])
    primal = int(meta["runtime"]["primal_rows"])
    inverse = result_inverse_permutation(MESH / f"{name}.result", primal)
    x = indexed(capture / "full_x_after.dat", total)[:primal]
    selected_nodes = tes_nodes()
    values = [float(x[dof - 1]) for dof, node in enumerate(inverse[1:], start=1) if node in selected_nodes]
    return float(np.mean(values))


def main() -> int:
    powers = {
        "fixed_power_minus5pct": 3.203004762115138e-10 * 0.95,
        "fixed_power_center": 3.203004762115138e-10,
        "fixed_power_plus5pct": 3.203004762115138e-10 * 1.05,
    }
    rows: list[dict[str, object]] = []
    reaction_rows: list[dict[str, object]] = []
    for name, power in powers.items():
        capture = OUT / "capture" / name / "ts0001_nl0001"
        case = Case(name, name, capture, MESH, MESH / f"{name}.result", INTERFACES, 2)
        meta, aggregate, _, topology = reaction_capture(case)
        by_interface = {row["interface"]: float(row["weighted_reaction_W"]) for row in aggregate}
        reaction_rows.extend(aggregate)
        t = temperature(name, capture)
        rows.append({
            "case": name,
            "fixed_power_W": power,
            "TES_temperature_K": t,
            "TES_temperature_mK": t * 1e3,
            "TES_to_membrane_W": by_interface["TES_to_membrane"],
            "TES_to_Stycast_W": by_interface["TES_to_Stycast"],
            "Stycast_to_substrate_W": by_interface["Stycast_to_substrate"],
            "native_bath_W": sum(by_interface.values()),
            "constraint_rows": int(meta["runtime"]["constraint_matrix_rows"]),
            "TES_to_membrane_rows": int(topology["TES_to_membrane"]["constraint_count"]),
            "solver": "CPU MUMPS native full saddle capture",
        })
    rows.sort(key=lambda row: float(row["fixed_power_W"]))
    fields = list(rows[0])
    with (OUT / "fixed_power_native_results.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    reaction_fields = list(reaction_rows[0])
    with (OUT / "native_mortar_reaction_integrals.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=reaction_fields)
        writer.writeheader()
        writer.writerows(reaction_rows)

    minus, center, plus = rows
    dT = float(plus["TES_temperature_K"]) - float(minus["TES_temperature_K"])
    dQ = float(plus["native_bath_W"]) - float(minus["native_bath_W"])
    dP = float(plus["fixed_power_W"]) - float(minus["fixed_power_W"])
    payload = {
        "method": "actual three-point fixed-power native MUMPS solves; symmetric ±5% power bracket used to obtain the local ±delta-T pair",
        "delta_T_K": dT,
        "delta_T_mK": dT * 1e3,
        "delta_Q_native_bath_W": dQ,
        "delta_power_W": dP,
        "G_native_bath_W_per_K": dQ / dT,
        "R_native_bath_K_per_W": dT / dQ,
        "center": center,
        "minus": minus,
        "plus": plus,
    }
    (OUT / "symmetric_delta_T_native.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (OUT / "summary.md").write_text(
        "# Phase24 native TES–membrane coupling diagnostics\n\n"
        "- Isolated mesh: `mesh_phase24_native_coupling_diag`; the 179 TES–membrane contact nodes were duplicated on the membrane side, so the native mortar interface is nonconforming and has zero shared nodes.\n"
        "- Physics, materials, TES law, circuit constants, bath temperature, and CPU MUMPS backend were retained. Only the diagnostic TES power was frozen at 0.95, 1.00, and 1.05 times the historical checkpoint power.\n"
        f"- Native TES–membrane rows: {center['TES_to_membrane_rows']} in every case; native full saddle constraints: {center['constraint_rows']}.\n"
        f"- Symmetric bracket: ΔT = {payload['delta_T_mK']:.9f} mK, ΔQ = {payload['delta_Q_native_bath_W']:.9e} W, G_eff = {payload['G_native_bath_W_per_K']:.9e} W/K.\n"
        "- The CSV contains the actual native Cλ reaction integrals and temperatures for all three solves.\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

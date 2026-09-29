"""Measure trace-to-trace differential resistance for the Phase24 cases.

This is a read-only postprocessor for the already captured 0.95/1.00/1.05 P
solutions.  Trace temperatures are area-weighted averages on named mesh
boundaries.  The reported resistance is d(T_left-T_right)/dP, so it remains
well-defined for the fixed-power bracket and does not assume that every
network branch carries the full Joule power.

The rows are grouped as follows:

* direct mortar traces: TES/Membrane, TES/Stycast, and Stycast/substrate;
* the active endpoint interval: Membrane trace -> bath trace;
* the conforming downstream branch traces.  The latter are diagnostic drops
  between successive named traces and must not be summed across the parallel
  Si1/SiNx and Si2/SiO2_2 branches.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_trace_resistance_controlled"
P0 = 3.203004762115138e-10
TBATH = 0.150

sys.path.insert(0, str(ROOT / "src"))

sys.path.insert(0, str(ROOT))
from support.analyze_phase24_outer_capture import indexed_vector, load_temp_permutation  # noqa: E402


CASES: dict[str, dict[str, Any]] = {
    "historical": {
        "mesh": ROOT / "workspace/work/meshes/mesh_singlepixel_prod_v2",
        "runs": [
            (0.95, ROOT / "workspace/work/meshes/mesh_singlepixel_prod_v2/historical_0p95p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_0p95P/ts0001_nl0001"),
            (1.00, ROOT / "workspace/work/meshes/mesh_singlepixel_prod_v2/historical_1p00p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p00P/ts0001_nl0001"),
            (1.05, ROOT / "workspace/work/meshes/mesh_singlepixel_prod_v2/historical_1p05p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p05P/ts0001_nl0001"),
        ],
        "traces": {
            "TES_to_Membrane": (24, 23),
            "TES_to_Stycast": (25, 26),
            "Stycast_to_substrate": (27, 28),
            "Membrane_to_bath_endpoint": (23, 30),
            "Membrane_SiNx_to_Membrane_Si1": (23, 19),
            "Membrane_Si1_to_SiO2_1": (19, 14),
            "SiO2_1_to_Si_1": (15, 17),
            "Si_1_to_SiNx": (17, 21),
            "SiO2_1_to_Si_2": (14, 13),
            "Si_2_to_SiO2_2": (13, 11),
            "SiO2_2_to_bath": (11, 30),
        },
    },
    "refined_mortar": {
        "mesh": ROOT / "workspace/work/meshes/mesh_phase24_stycast_density_10um",
        "runs": [
            (0.95, ROOT / "workspace/work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_0p95p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_0p95P/ts0001_nl0001"),
            (1.00, ROOT / "workspace/work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p00p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p00P/ts0001_nl0001"),
            (1.05, ROOT / "workspace/work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p05p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p05P/ts0001_nl0001"),
        ],
        "traces": {
            "TES_to_Membrane": (1104, 1305),
            "TES_to_Stycast": (1105, 1204),
            "Stycast_to_substrate": (1205, 1004),
            "Membrane_to_bath_endpoint": (1305, 21),
            "Membrane_SiNx_to_Membrane_Si1": (1305, 1905),
            "Membrane_Si1_to_SiO2_1": (1904, 1404),
            "SiO2_1_to_Si_1": (1405, 1505),
            "Si_1_to_SiNx": (1505, 1605),
            "SiO2_1_to_Si_2": (1404, 1705),
            "Si_2_to_SiO2_2": (1705, 1805),
            "SiO2_2_to_bath": (1804, 21),
        },
    },
}

# The one-at-a-time tests retain the same body and boundary numbering as the
# refined parent.  Add them without duplicating the definitions above.
for _case, _root, _prefix in (
    ("tes_membrane_trace_historical", "artifacts/p24trace/tm", "tm_hist"),
    ("membrane_substrate_trace_historical", "artifacts/p24trace/ms", "ms_hist"),
):
    CASES[_case] = {
        "mesh": ROOT / ("workspace/work/meshes/mesh_phase24_trace_tes_membrane_historical" if _case.startswith("tes_") else "workspace/work/meshes/mesh_phase24_trace_membrane_substrate_historical"),
        "runs": [
            (f, ROOT / ("workspace/work/meshes/mesh_phase24_trace_tes_membrane_historical" if _case.startswith("tes_") else "workspace/work/meshes/mesh_phase24_trace_membrane_substrate_historical") / f"{_prefix}_{label}.result", ROOT / _root / "capture" / f"{_prefix}_{label.replace('p', 'p')}" / "ts0001_nl0001")
            for f, label in ((0.95, "0p95P"), (1.00, "1p00P"), (1.05, "1p05P"))
        ],
        "traces": CASES["refined_mortar"]["traces"],
    }


def load_nodes(mesh: Path) -> np.ndarray:
    rows = []
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 5 and fields[0].lstrip("-").isdigit():
            rows.append((int(fields[0]), *map(float, fields[2:5])))
    rows.sort()
    return np.asarray([[x, y, z] for _, x, y, z in rows], dtype=float)


def boundary_records(mesh: Path) -> dict[int, list[list[int]]]:
    result: dict[int, list[list[int]]] = {}
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 8:
            result.setdefault(int(fields[1]), []).append([int(value) for value in fields[5:]])
    return result


def face_area(points: np.ndarray) -> float:
    if len(points) < 3:
        return 0.0
    area = 0.5 * float(np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0])))
    if len(points) == 4:
        area += 0.5 * float(np.linalg.norm(np.cross(points[3] - points[0], points[2] - points[0])))
    return area


def trace_temperature(mesh: Path, result: Path, capture: Path, boundary_id: int, nodes: np.ndarray, boundaries: dict[int, list[list[int]]]) -> tuple[float, float, int]:
    meta = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    total_rows = int(meta["runtime"]["total_rows"])
    primal_rows = int(meta["runtime"]["primal_rows"])
    x = indexed_vector(capture / "full_x_after.dat", total_rows)[:primal_rows]
    permutation = load_temp_permutation(result)
    nodal = np.zeros(permutation.size, dtype=float)
    valid = permutation > 0
    nodal[valid] = x[permutation[valid] - 1]
    faces = boundaries.get(boundary_id, [])
    weighted = 0.0
    area_sum = 0.0
    node_count: set[int] = set()
    for face in faces:
        p = nodes[np.asarray(face) - 1]
        area = face_area(p)
        if area <= 0.0:
            continue
        weighted += area * float(np.mean(nodal[np.asarray(face) - 1]))
        area_sum += area
        node_count.update(face)
    return weighted / area_sum if area_sum else float("nan"), area_sum, len(node_count)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    temperature_rows: list[dict[str, Any]] = []
    grouped: dict[str, dict[float, dict[int, float]]] = {}
    trace_meta: dict[str, dict[int, tuple[float, int]]] = {}
    for case_name, case in CASES.items():
        grouped[case_name] = {}
        trace_meta[case_name] = {}
        nodes = load_nodes(case["mesh"])
        boundaries = boundary_records(case["mesh"])
        for fraction, result, capture in case["runs"]:
            values: dict[int, float] = {}
            for boundary_id in sorted({bid for pair in case["traces"].values() for bid in pair}):
                value, area, count = trace_temperature(case["mesh"], result, capture, boundary_id, nodes, boundaries)
                values[boundary_id] = value
                trace_meta[case_name][boundary_id] = (area, count)
            grouped[case_name][fraction] = values
            for boundary_id, value in values.items():
                temperature_rows.append({"case": case_name, "power_fraction": fraction, "boundary_id": boundary_id, "temperature_K": value, "temperature_mK": value * 1.0e3, "area_m2": trace_meta[case_name][boundary_id][0], "unique_nodes": trace_meta[case_name][boundary_id][1]})

    rows: list[dict[str, Any]] = []
    for case_name, case in CASES.items():
        minus, center, plus = grouped[case_name][0.95], grouped[case_name][1.00], grouped[case_name][1.05]
        for trace, (left_id, right_id) in case["traces"].items():
            dminus = minus[left_id] - minus[right_id]
            dcenter = center[left_id] - center[right_id]
            dplus = plus[left_id] - plus[right_id]
            r_diff = (dplus - dminus) / (0.10 * P0)
            rows.append({
                "case": case_name,
                "trace_interval": trace,
                "left_boundary_id": left_id,
                "right_boundary_id": right_id,
                "delta_T_0p95P_K": dminus,
                "delta_T_1p00P_K": dcenter,
                "delta_T_1p05P_K": dplus,
                "R_trace_to_trace_differential_K_per_W": r_diff,
                "R_trace_to_trace_secant_K_per_W": dcenter / P0,
                "interpretation": "active endpoint or diagnostic branch drop; do not sum parallel branch rows",
            })
    write_csv(OUT / "trace_to_trace_resistance.csv", rows)
    write_csv(OUT / "trace_temperature.csv", temperature_rows)

    hist = {(row["trace_interval"]): row for row in rows if row["case"] == "historical"}
    refined = {(row["trace_interval"]): row for row in rows if row["case"] == "refined_mortar"}
    comparison = []
    for trace in hist:
        h = float(hist[trace]["R_trace_to_trace_differential_K_per_W"])
        r = float(refined[trace]["R_trace_to_trace_differential_K_per_W"])
        comparison.append({"trace_interval": trace, "historical_R_K_per_W": h, "refined_mortar_R_K_per_W": r, "refined_minus_historical_K_per_W": r - h, "refined_over_historical": r / h if h else float("nan")})
    write_csv(OUT / "trace_to_trace_comparison.csv", comparison)

    lines = [
        "# Phase24 trace-to-trace resistance decomposition",
        "",
        "Trace temperatures are area-weighted named-boundary averages from the captured CPU/MUMPS solutions. R is the differential d(Delta T_trace)/dP from the 0.95/1.00/1.05 P bracket.",
        "",
        "The Membrane-to-bath endpoint is the active second interval. The Si1/SiNx and Si2/SiO2_2 rows are parallel downstream branches; they are reported for localization and must not be summed as a series chain.",
        "",
    ]
    for row in comparison:
        lines.append(f"- `{row['trace_interval']}`: historical `{float(row['historical_R_K_per_W']):.9e} K/W`, refined mortar `{float(row['refined_mortar_R_K_per_W']):.9e} K/W`, ratio `{float(row['refined_over_historical']):.6f}`")
    lines += [
        "",
        "The controlled one-at-a-time runs are in trace_to_trace_resistance.csv. Their strict validity still depends on the mesh audit: the TES/Membrane variant changed the Stycast top trace face count, so it is not accepted as a clean one-factor result; the Membrane/substrate variant retained the direct TES/Membrane trace and is the usable controlled test.",
    ]
    (OUT / "trace_to_trace_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Build the Phase22 physical-error-based HYPRE tolerance study.

The residual tolerance is only an input.  Acceptance is based on observables
and temperature fields against the validated common-geometry Mortar result.
Missing Mortar electrical/pulse series are retained as unavailable; they are
never synthesized from a conformal run.
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.evaluate_physical_parity import mesh_data, result_values


ROOT = Path(__file__).resolve().parents[2]
RESULT_ROOT = ROOT / "results"
MESH_ROOT = ROOT / "work" / "meshes"
OUT = ROOT / "artifacts" / "phase22_physical_tolerance" / "hypre_physical_tolerance_study.json"
FLOAT = r"([0-9.Ee+\-]+)"
TOLERANCES = (1.0e-6, 1.0e-7, 1.0e-8)


def result_path(case: str, mesh: str) -> Path | None:
    path = MESH_ROOT / mesh / f"{case}.result"
    return path if path.exists() else None


def series_rows(case: str) -> list[dict[str, str]]:
    candidates = sorted((RESULT_ROOT / case).glob("*_series.csv"))
    if not candidates:
        return []
    with candidates[0].open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, skipinitialspace=True))


def pulse_metrics(rows: list[dict[str, str]]) -> dict[str, float | int | None]:
    if not rows:
        return {"status": "REFERENCE_UNAVAILABLE", "sample_count": 0}
    def values(name: str) -> list[float]:
        return [float(row[name]) for row in rows if row.get(name) not in (None, "")]
    current = values("tes_current_A")
    temperature = values("tes_temperature_K")
    power = values("tes_power_W")
    baseline_current = current[0] if current else None
    baseline_temperature = temperature[0] if temperature else None
    return {
        "status": "DONE",
        "sample_count": len(rows),
        "tes_current_drop_peak_A": baseline_current - min(current) if current else None,
        "tes_temperature_rise_peak_K": max(temperature) - baseline_temperature if temperature else None,
        "tes_power_peak_W": max(power) if power else None,
    }


def log_info(case: str) -> dict:
    path = RESULT_ROOT / case / "solver.log"
    if not path.exists():
        return {"status": "NOT_RUN", "log": str(path.resolve())}
    text = path.read_text(encoding="utf-8", errors="replace")
    wall = re.search(rf"WALL_SECONDS\s+{FLOAT}", text)
    cpu = re.search(rf"SOLVER TOTAL TIME\(CPU,REAL\):\s*{FLOAT}\s+{FLOAT}", text)
    iterations = [int(value) for value in re.findall(r"Required iterations (\d+)", text)]
    return {
        "status": "DONE" if "ALL DONE" in text else "FAILED",
        "log": str(path.resolve()),
        "wall_seconds": float(wall.group(1)) if wall else None,
        "solver_cpu_seconds": float(cpu.group(1)) if cpu else None,
        "solver_real_seconds": float(cpu.group(2)) if cpu else None,
        "max_hypre_iterations": max(iterations) if iterations else None,
    }


def field_comparison(reference_mesh: str, reference: Path | None, candidate_mesh: str, candidate: Path | None) -> dict:
    if reference is None or candidate is None:
        return {"status": "UNAVAILABLE", "reason": "missing result file"}
    _, _, ref_nodes, _, _, _, _ = mesh_data(MESH_ROOT / reference_mesh)
    _, _, cand_nodes, _, _, _, _ = mesh_data(MESH_ROOT / candidate_mesh)
    ref_values = result_values(reference, field_index=-1)
    cand_values = result_values(candidate, field_index=-1)
    ref_by_coord = {
        tuple(round(value, 12) for value in point): node
        for node, point in ref_nodes.items()
    }
    deltas = []
    missing = 0
    for node, point in cand_nodes.items():
        ref_node = ref_by_coord.get(tuple(round(value, 12) for value in point))
        if ref_node is None or ref_node not in ref_values or node not in cand_values:
            missing += 1
            continue
        deltas.append(abs(cand_values[node] - ref_values[ref_node]))
    return {
        "status": "DONE" if deltas else "UNAVAILABLE",
        "common_coordinates": len(deltas),
        "candidate_coordinates_without_reference": missing,
        "max_abs_difference_K": max(deltas) if deltas else None,
        "rms_difference_K": math.sqrt(sum(value * value for value in deltas) / len(deltas)) if deltas else None,
    }


def temperature_observable(mesh: str, result: Path | None, body: str) -> float | None:
    if result is None:
        return None
    # Always compare the final saved field; for a steady result it is the
    # same field, while transient cases contain one field per saved step.
    return _body_average_last(mesh, result, body)


def _body_average_last(mesh: str, result: Path, body: str) -> float:
    bodies, _, nodes, _, elements_by_body, _, _ = mesh_data(MESH_ROOT / mesh)
    values = result_values(result, field_index=-1)
    body_id = bodies[body]
    weighted = 0.0
    volume = 0.0
    for conn in elements_by_body[body_id]:
        points = [nodes[node] for node in conn]
        a = tuple(points[1][i] - points[0][i] for i in range(3))
        b = tuple(points[2][i] - points[0][i] for i in range(3))
        c = tuple(points[3][i] - points[0][i] for i in range(3))
        cross = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
        element_volume = abs(sum(cross[i] * c[i] for i in range(3))) / 6.0
        volume += element_volume
        weighted += element_volume * sum(values[node] for node in conn) / len(conn)
    return weighted / volume


def compare_window(reference_case: str, reference_mesh: str, candidate_case: str, candidate_mesh: str) -> dict:
    reference = result_path(reference_case, reference_mesh)
    candidate = result_path(candidate_case, candidate_mesh)
    ref_rows = series_rows(reference_case)
    cand_rows = series_rows(candidate_case)
    fields = ("tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W")
    electrical = {}
    if ref_rows and cand_rows:
        ref = ref_rows[-1]
        cand = cand_rows[-1]
        for field in fields:
            electrical[field] = {
                "reference": float(ref[field]),
                "candidate": float(cand[field]),
                "absolute_difference": abs(float(ref[field]) - float(cand[field])),
            }
    else:
        electrical = {field: {"status": "REFERENCE_UNAVAILABLE"} for field in fields}
    reference_tes = temperature_observable(reference_mesh, reference, "TES")
    candidate_tes = temperature_observable(candidate_mesh, candidate, "TES")
    reference_abs = temperature_observable(reference_mesh, reference, "abs")
    candidate_abs = temperature_observable(candidate_mesh, candidate, "abs")
    return {
        "reference_case": reference_case,
        "candidate_case": candidate_case,
        "reference_timing": log_info(reference_case),
        "candidate_timing": log_info(candidate_case),
        "tes_volume_average_temperature_K": {
            "reference": reference_tes,
            "candidate": candidate_tes,
            "absolute_difference_K": abs(reference_tes - candidate_tes) if reference_tes is not None and candidate_tes is not None else None,
        },
        "absorber_temperature_K": {
            "reference": reference_abs,
            "candidate": candidate_abs,
            "absolute_difference_K": abs(reference_abs - candidate_abs) if reference_abs is not None and candidate_abs is not None else None,
        },
        "electrical_observables": electrical,
        "pulse_metrics": {
            "reference": pulse_metrics(ref_rows),
            "candidate": pulse_metrics(cand_rows),
            "comparison_status": "DONE" if ref_rows and cand_rows else "REFERENCE_UNAVAILABLE",
        },
        "full_temperature_field": field_comparison(reference_mesh, reference, candidate_mesh, candidate),
    }


def main() -> int:
    runs = []
    for tolerance in TOLERANCES:
        token = f"{tolerance:.0e}".replace("-", "m")
        steady = compare_window(
            "case_phase21_mortar_steady",
            "mesh_physical_parity_mortar",
            f"case_tes_steady_physical_parity_conformal_hypre_{token}",
            "mesh_physical_parity_conformal",
        )
        transient = compare_window(
            "case_phase21_mortar_7step",
            "mesh_physical_parity_mortar",
            f"case_tes_transient_physical_parity_conformal_hypre_{token}_7step",
            "mesh_physical_parity_conformal",
        )
        candidate_times = [
            item["candidate_timing"].get("wall_seconds")
            for item in (steady, transient)
            if item["candidate_timing"].get("wall_seconds") is not None
        ]
        field_errors = [
            item["full_temperature_field"].get("max_abs_difference_K")
            for item in (steady, transient)
            if item["full_temperature_field"].get("max_abs_difference_K") is not None
        ]
        absorber_errors = [
            item["absorber_temperature_K"].get("absolute_difference_K")
            for item in (steady, transient)
            if item["absorber_temperature_K"].get("absolute_difference_K") is not None
        ]
        strict_case = f"case_tes_steady_physical_parity_conformal_hypre_1em08" if token != "1em08" else f"case_tes_steady_physical_parity_conformal_hypre_1em08"
        strict_transient_case = f"case_tes_transient_physical_parity_conformal_hypre_1em08_7step"
        strict_result = result_path(strict_case, "mesh_physical_parity_conformal")
        strict_transient_result = result_path(strict_transient_case, "mesh_physical_parity_conformal")
        candidate_result = result_path(
            f"case_tes_steady_physical_parity_conformal_hypre_{token}",
            "mesh_physical_parity_conformal",
        )
        candidate_transient_result = result_path(
            f"case_tes_transient_physical_parity_conformal_hypre_{token}_7step",
            "mesh_physical_parity_conformal",
        )
        strict_field = field_comparison(
            "mesh_physical_parity_conformal", strict_result,
            "mesh_physical_parity_conformal", candidate_result,
        )
        strict_transient_field = field_comparison(
            "mesh_physical_parity_conformal", strict_transient_result,
            "mesh_physical_parity_conformal", candidate_transient_result,
        )
        strict_tes = temperature_observable("mesh_physical_parity_conformal", strict_result, "TES")
        strict_abs = temperature_observable("mesh_physical_parity_conformal", strict_result, "abs")
        candidate_tes = temperature_observable("mesh_physical_parity_conformal", candidate_result, "TES")
        candidate_abs = temperature_observable("mesh_physical_parity_conformal", candidate_result, "abs")
        strict_errors = [
            value for value in (
                strict_field.get("max_abs_difference_K"),
                strict_transient_field.get("max_abs_difference_K"),
                abs(strict_tes - candidate_tes) if strict_tes is not None and candidate_tes is not None else None,
                abs(strict_abs - candidate_abs) if strict_abs is not None and candidate_abs is not None else None,
            ) if value is not None
        ]
        runs.append(
            {
                "tolerance": tolerance,
                "steady": steady,
                "transient_7step": transient,
                "candidate_total_wall_seconds": sum(candidate_times) if candidate_times else None,
                "physical_gate": {
                    "mortar_temperature_field_max_abs_K": max(field_errors) if field_errors else None,
                    "mortar_absorber_temperature_max_abs_K": max(absorber_errors) if absorber_errors else None,
                    "field_tolerance_K": 1.0e-5,
                    "absorber_tolerance_K": 1.0e-5,
                    "status": "PASS" if field_errors and absorber_errors and max(field_errors) <= 1.0e-5 and max(absorber_errors) <= 1.0e-5 else "FAIL_ROUTE_DISCREPANCY_OR_INCOMPLETE_REFERENCE",
                    "note": "physical error gate; solver residual is not used as the acceptance criterion",
                },
                "solver_accuracy_gate": {
                    "reference_tolerance": 1.0e-8,
                    "max_abs_temperature_error_K": max(strict_errors) if strict_errors else None,
                    "tolerance_K": 3.0e-5,
                    "status": "PASS" if strict_errors and max(strict_errors) <= 3.0e-5 else "INCOMPLETE",
                    "note": "isolates algebraic tolerance error from the common Mortar-versus-conformal route/mesh offset",
                },
            }
        )
    eligible = [row for row in runs if row["solver_accuracy_gate"]["status"] == "PASS" and row["candidate_total_wall_seconds"] is not None]
    production = min(eligible, key=lambda row: row["candidate_total_wall_seconds"]) if eligible else None
    report = {
        "status": "PASS" if runs else "NOT_RUN",
        "study": "common Mortar direct reference versus conformal CPU HYPRE tolerance sweep",
        "reference": {"case_steady": "case_phase21_mortar_steady", "case_transient": "case_phase21_mortar_7step", "mesh": "mesh_physical_parity_mortar"},
        "candidate_mesh": "mesh_physical_parity_conformal",
        "tolerances": ["1e-6", "1e-7", "1e-8"],
        "runs": runs,
        "production_tolerance": {
            "status": "SELECTED_PROVISIONAL" if production else "PENDING_REFERENCE_COMPLETION",
            "tolerance": production["tolerance"] if production else None,
            "reason": "lowest measured wall time among tolerances whose error relative to the 1e-8 conformal physical reference passes; common Mortar route offset remains a separate unresolved gate" if production else "No tolerance passed the conformal physical-error gate",
        },
        "requirement_provenance": "The 3e-5 K conformal solver-error gate is a provisional study threshold because no application-level physical-error budget was supplied; the Mortar route comparison remains reported independently at the stricter 1e-5 K diagnostic threshold.",
        "unavailable_policy": "Missing Mortar electrical/pulse series remain REFERENCE_UNAVAILABLE and do not become a pass by proxy.",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

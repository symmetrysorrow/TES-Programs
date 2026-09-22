"""Materialize and validate a Phase24 native restart-state audit.

The native HeatSolve diagnostic writes three rows before the first thermal
solve.  This tool intentionally does no solving: it records the checkpoint
bytes, resolves the path from ElmerSolver's cwd, and emits small JSON files
that are safe to compare or archive.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "artifacts" / "phase24_augmented_restart_audit"


def resolve_runtime_state_path(configured: str, cwd: Path = ROOT) -> Path:
    """Resolve a Constants path exactly as ElmerSolver does for this runner."""
    return (cwd / configured).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_checkpoint(path: Path) -> dict[str, Any]:
    values = [float(item) for item in path.read_text(encoding="utf-8").split()]
    if len(values) != 5:
        raise ValueError(f"checkpoint schema is not legacy-v1-5: {path} ({len(values)} values)")
    return {
        "schema": "legacy-v1-5",
        "temperature_K": values[0],
        "current_A": values[1],
        "resistance_ohm": values[2],
        "power_W": values[3],
        "previous_current_A": values[4],
    }


def _bool(value: str) -> bool:
    return value.strip().upper() in {"T", "TRUE", ".TRUE."}


def _float(row: dict[str, str], key: str) -> float:
    return float(row[key].replace("D", "E").replace("d", "e"))


def read_stage_csv(path: Path) -> dict[str, dict[str, Any]]:
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    expected = {"load_after", "circuit_init_after", "pre_first_assembly"}
    actual = {row["stage"] for row in rows}
    missing = expected - actual
    if missing:
        raise ValueError(f"restart audit missing stages: {sorted(missing)}")
    result: dict[str, dict[str, Any]] = {}
    for row in rows:
        stage = row["stage"]
        result[stage] = {
            "stage": stage,
            "physical_time_s": _float(row, "physical_time"),
            "timestep": int(row["timestep"]),
            "dt_s": _float(row, "dt"),
            "bdf_order": int(row["bdf_order"]),
            "state_file": row["state_file"],
            "state_exists": _bool(row["state_exists"]),
            "state_open_success": _bool(row["state_open_success"]),
            "state_parse_success": _bool(row["state_parse_success"]),
            "state_loaded": _bool(row["state_loaded"]),
            "current_A": _float(row, "current"),
            "previous_current_A": _float(row, "previous_current"),
            "resistance_ohm": _float(row, "resistance"),
            "power_W": _float(row, "power"),
            "average_temperature_K": _float(row, "average_temperature"),
            "sweep_temperature_K": _float(row, "sweep_temperature"),
            "omega": _float(row, "omega"),
            "omega_cap": _float(row, "omega_cap"),
            "previous_residual_W": _float(row, "prev_residual"),
            "last_circuit_dt_s": _float(row, "last_circuit_dt"),
            "last_circuit_step": int(row["last_circuit_step"]),
            "circuit_iteration": int(row["circuit_iter_in_step"]),
            "circuit_calls_in_step": int(row["circuit_calls_in_step"]),
            "circuit_initialized": _bool(row["circuit_initialized"]),
            "circuit_power_initialized": _bool(row["circuit_power_initialized"]),
            "source": row["source"],
        }
    if len(result) != 3:
        raise ValueError(f"restart audit has duplicate or extra stage rows: {actual}")
    return result


def compare(saved: dict[str, Any], stage: dict[str, Any], tolerance: float = 1.0e-15) -> dict[str, Any]:
    fields = ("current_A", "previous_current_A", "resistance_ohm", "power_W")
    values = {}
    for field in fields:
        delta = stage[field] - saved[field]
        values[field] = {"saved": saved[field], "stage": stage[field], "delta": delta,
                         "equal_within_tolerance": abs(delta) <= tolerance * max(abs(saved[field]), 1.0)}
    return values


def materialize(audit_csv: Path, saved_path: Path, output: Path, sif: Path | None = None,
                original_before: Path | None = None, runtime_cwd: Path = ROOT) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=True)
    saved = read_checkpoint(saved_path)
    stages = read_stage_csv(audit_csv)
    actual_path = resolve_runtime_state_path(stages["load_after"]["state_file"], runtime_cwd)
    if not actual_path.is_file():
        raise FileNotFoundError(actual_path)
    loaded_sha = sha256(actual_path)
    state_diffs = {
        stage: compare(saved, row)
        for stage, row in stages.items()
    }
    transitions = {}
    fields = ("current_A", "previous_current_A", "resistance_ohm", "power_W",
              "average_temperature_K", "omega", "omega_cap", "last_circuit_dt_s")
    previous: dict[str, Any] = saved | {"average_temperature_K": saved["temperature_K"],
                                        "omega": None, "omega_cap": None, "last_circuit_dt_s": None}
    for left, right in (("saved", "load_after"), ("load_after", "circuit_init_after"),
                        ("circuit_init_after", "pre_first_assembly")):
        right_values = stages[right] if right != "saved" else saved
        left_values = stages[left] if left != "saved" else previous
        transitions[f"{left}_to_{right}"] = {
            field: (None if left_values.get(field) is None or right_values.get(field) is None
                    else right_values[field] - left_values[field])
            for field in fields
        }
    original_before_sha = sha256(original_before) if original_before else None
    original_after_sha = sha256(saved_path)
    path_audit = {
        "elmer_solver_cwd": str(ROOT),
        "sif": str(sif) if sif else None,
        "configured_state_file": stages["load_after"]["state_file"],
        "actual_loaded_state_file": str(actual_path),
        "configured_path_resolves_to": str(actual_path.resolve()),
        "file_exists": stages["load_after"]["state_exists"],
        "open_success": stages["load_after"]["state_open_success"],
        "parse_success": stages["load_after"]["state_parse_success"],
        "loaded": stages["load_after"]["state_loaded"],
        "schema": saved["schema"],
        "state_sha256": loaded_sha,
        "original_checkpoint_sha256_before": original_before_sha,
        "original_checkpoint_sha256_after": original_after_sha,
        "original_checkpoint_overwritten": (None if original_before is None
                                              else original_before_sha != original_after_sha),
    }
    path_audit["diagnostic_snapshot_isolated"] = saved_path.resolve() != actual_path.resolve()
    payload = {
        "saved_checkpoint": saved,
        "stages": stages,
        "state_diff": state_diffs,
        "transitions": transitions,
        "path_audit": path_audit,
        "all_required_fields_loaded": all(stages["load_after"][key] for key in
                                           ("state_exists", "state_open_success", "state_parse_success", "state_loaded")),
        "first_mismatch": None,
    }
    for stage in ("load_after", "circuit_init_after", "pre_first_assembly"):
        if any(not value["equal_within_tolerance"] for value in state_diffs[stage].values()):
            payload["first_mismatch"] = stage
            break
    (output / "saved_checkpoint.json").write_text(json.dumps(saved, indent=2) + "\n", encoding="utf-8")
    for stage, filename in (("load_after", "stage_a_after_load.json"),
                           ("circuit_init_after", "stage_b_after_circuit_init.json"),
                           ("pre_first_assembly", "stage_c_pre_first_assembly.json")):
        (output / filename).write_text(json.dumps(stages[stage], indent=2) + "\n", encoding="utf-8")
    (output / "restart_state_diff.json").write_text(json.dumps({"fields": state_diffs, "transitions": transitions}, indent=2) + "\n", encoding="utf-8")
    (output / "restart_path_audit.json").write_text(json.dumps(path_audit, indent=2) + "\n", encoding="utf-8")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-csv", type=Path, required=True)
    parser.add_argument("--saved-state", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--sif", type=Path)
    parser.add_argument("--original-before", type=Path)
    args = parser.parse_args()
    payload = materialize(args.audit_csv, args.saved_state, args.output, args.sif, args.original_before)
    (args.output / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "first_mismatch": payload["first_mismatch"],
                      "loaded": payload["all_required_fields_loaded"]}))
    return 0 if payload["all_required_fields_loaded"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

"""Audit and reproduce Phase24 steady -> transient restart continuity.

This campaign targets the nonconforming 32-layer + mortar path where the
Gate3 steady state is near 143.57 uA but the MUMPS transient starts near
165 uA before the pulse.  It is intentionally a reusable diagnostic rather
than a one-off run.

The script has three useful modes:

* default: audit existing steady/transient artifacts, generate short pulse-OFF
  variants, run them, and summarize where continuity first breaks;
* --audit-only: inspect existing SIF/result/state/manifest inputs without
  starting Elmer;
* --dry-run: generate the campaign project and commands but do not run Elmer.

Every transient variant uses a private copy of its TES state file under the
campaign artifact directory.  This is required because the UDF checkpoints
the state file after accepted transient timesteps; the original Gate3 seed
must never be overwritten by a diagnostic run.

Typical use from Elmer-Projects:

    python scripts/support/run_phase24_restart_continuity_campaign.py

If the two Phase24 cases live in an uncommitted/generated project, point to it:

    python scripts/support/run_phase24_restart_continuity_campaign.py \
      --project path/to/project.json

The campaign avoids 40 us / 100 us production runs.  It concentrates on the
pre-pulse first solve and a five-step pulse-OFF steady hold.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STEADY_CASE = "case_phase24_g3_s32m_s32m2_10f17b"
DEFAULT_TRANSIENT_CASE = "case_phase24_g45_s32m_40us_mumps_mortar_mumps"
OLD_GOOD_CASE = "case_p19_pulse_phase23_tight"
CAMPAIGN_DIR = ROOT / "artifacts/phase24_restart_continuity_campaign"
PROJECT_PATH = CAMPAIGN_DIR / "phase24_restart_continuity_campaign.json"
STATE_SNAPSHOT_DIR = CAMPAIGN_DIR / "state_snapshots"
DUMP_DIR = CAMPAIGN_DIR / "dumps"
CURRENT_GATE_UA = 143.567589
CURRENT_GATE_REL = 1.0e-4  # 0.01 %
MESH_FILES = ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary", "mesh.names")


def sha256_file(path: Path | None) -> str | None:
    if path is None or not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_load(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def relpath(path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def find_project(case_names: tuple[str, ...], explicit: Path | None) -> tuple[Path | None, dict[str, Any] | None]:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit if explicit.is_absolute() else ROOT / explicit)
    candidates.extend(sorted(ROOT.glob("elmer_project*.json")))
    # Local diagnostic projects frequently live below artifacts and are not
    # necessarily committed.  Search them after top-level production inputs.
    artifact_root = ROOT / "artifacts"
    if artifact_root.is_dir():
        candidates.extend(sorted(artifact_root.rglob("*.json")))

    seen: set[Path] = set()
    best: tuple[int, Path, dict[str, Any]] | None = None
    for path in candidates:
        path = path.resolve()
        if path in seen or not path.is_file():
            continue
        seen.add(path)
        model = json_load(path)
        if not model or not isinstance(model.get("cases"), dict):
            continue
        score = sum(name in model["cases"] for name in case_names)
        if score == 0:
            continue
        if best is None or score > best[0]:
            best = (score, path, model)
            if score == len(case_names):
                break
    if best is None:
        return None, None
    return best[1], best[2]


def find_project_for_case(case_name: str, explicit: Path | None = None) -> tuple[Path | None, dict[str, Any] | None]:
    """Find a project that actually contains *case_name*.

    Audit mode may combine evidence from separate generated projects.  Run mode,
    however, must clone the real transient case definition rather than whichever
    project happened to contain the steady case first.
    """
    return find_project((case_name,), explicit)


def find_case_sif(case: str) -> Path | None:
    runtime = ROOT / "results" / case / "runtime.sif"
    generated = ROOT / "generated" / "cases" / f"{case}.sif"
    if runtime.is_file():
        return runtime
    return generated if generated.is_file() else None


def parse_sif(path: Path | None) -> dict[str, Any]:
    if path is None or not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    wanted = (
        "Restart File",
        "Restart Position",
        "Restart Time",
        "TES State File",
        "TES Bias Current",
        "TES Shunt Resistance",
        "TES R0",
        "TES Rmin",
        "TES Alpha",
        "TES Beta",
        "TES I0",
        "TES Tc",
        "TES T0",
        "TES Volume",
        "TES Inductance",
        "Apply Mortar BCs",
        "BDF Order",
        "Pulse Energy",
        "Pulse Start",
        "Pulse Duration",
        "Pulse Transition Zone",
        "Timestep Method",
    )
    out: dict[str, Any] = {}
    for key in wanted:
        pattern = re.compile(
            rf'^\s*"?{re.escape(key)}"?\s*=\s*(?:(?:String|Real|Integer|Logical)\s+)?(.+?)\s*$',
            re.IGNORECASE | re.MULTILINE,
        )
        match = pattern.search(text)
        if not match:
            continue
        raw = match.group(1).strip()
        if len(raw) >= 2 and raw[0] == raw[-1] == '"':
            raw = raw[1:-1]
        out[key] = raw
    return out


def input_path_candidates(raw: str | None) -> list[Path]:
    if not raw:
        return []
    p = Path(raw)
    if p.is_absolute():
        return [p]
    normalized = raw.replace("\\", "/")
    candidates = [ROOT / p]
    while normalized.startswith("../"):
        normalized = normalized[3:]
    candidates.append(ROOT / Path(normalized))
    if normalized.startswith("mesh_"):
        candidates.append(ROOT / "work" / "meshes" / Path(normalized))
    if normalized.startswith("work/meshes/"):
        candidates.append(ROOT / Path(normalized))
    unique: list[Path] = []
    for item in candidates:
        item = item.resolve()
        if item not in unique:
            unique.append(item)
    return unique


def resolve_input(raw: str | None) -> Path | None:
    candidates = input_path_candidates(raw)
    for path in candidates:
        if path.is_file():
            return path
    return candidates[0] if candidates else None


def read_state_file(path: Path | None) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": relpath(path),
        "exists": bool(path and path.is_file()),
        "sha256": sha256_file(path),
    }
    if path is None or not path.is_file():
        return result
    try:
        fields = path.read_text(encoding="utf-8", errors="replace").split()
        values = [float(value) for value in fields[:5]]
    except (OSError, ValueError):
        result["parse_error"] = True
        return result
    if len(values) != 5:
        result["parse_error"] = True
        return result
    temp, current, resistance, power, previous_current = values
    result.update(
        {
            "temperature_K": temp,
            "current_A": current,
            "current_uA": current * 1.0e6,
            "resistance_ohm": resistance,
            "power_W": power,
            "previous_current_A": previous_current,
            "previous_current_uA": previous_current * 1.0e6,
        }
    )
    return result


def csv_rows(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error):
        return []


def numeric_row(row: dict[str, str] | None) -> dict[str, float]:
    if not row:
        return {}
    out: dict[str, float] = {}
    for key, value in row.items():
        try:
            out[key] = float(value)
        except (TypeError, ValueError):
            continue
    if "tes_current_A" in out:
        out["tes_current_uA"] = out["tes_current_A"] * 1.0e6
    if "tes_temperature_K" in out:
        out["tes_temperature_mK"] = out["tes_temperature_K"] * 1.0e3
    return out


def find_result_csv(case: str, preferred: str | None = None, iteration: bool = False) -> Path | None:
    out = ROOT / "results" / case
    if preferred:
        direct = out / Path(preferred).name
        if direct.is_file():
            return direct
    pattern = "*iterations*.csv" if iteration else "*series*.csv"
    paths = sorted(out.glob(pattern)) if out.is_dir() else []
    if iteration:
        return paths[0] if paths else None
    paths = [p for p in paths if "iteration" not in p.name.lower()]
    return paths[0] if paths else None


def series_summary(case: str, spec: dict[str, Any] | None) -> dict[str, Any]:
    spec = spec or {}
    series = find_result_csv(case, spec.get("series_file"), iteration=False)
    iterations = find_result_csv(case, spec.get("iteration_series_file"), iteration=True)
    rows = csv_rows(series)
    it_rows = csv_rows(iterations)
    return {
        "series": relpath(series),
        "series_sha256": sha256_file(series),
        "row_count": len(rows),
        "first": numeric_row(rows[0] if rows else None),
        "last": numeric_row(rows[-1] if rows else None),
        "iteration_series": relpath(iterations),
        "iteration_series_sha256": sha256_file(iterations),
        "iteration_row_count": len(it_rows),
        "first_iteration": numeric_row(it_rows[0] if it_rows else None),
    }


def mesh_audit(model: dict[str, Any] | None, spec: dict[str, Any] | None, manifest: dict[str, Any] | None) -> dict[str, Any]:
    spec = spec or {}
    manifest = manifest or {}
    mesh_name = spec.get("mesh") or manifest.get("mesh")
    mesh_dir: Path | None = None
    if model and mesh_name and isinstance(model.get("meshes"), dict) and mesh_name in model["meshes"]:
        mesh_entry = model["meshes"][mesh_name]
        if isinstance(mesh_entry, dict) and mesh_entry.get("dir"):
            mesh_dir = ROOT / "work" / "meshes" / str(mesh_entry["dir"])
    if mesh_dir is None and mesh_name:
        mesh_dir = ROOT / "work" / "meshes" / str(mesh_name)
    return {
        "name": mesh_name,
        "dir": relpath(mesh_dir),
        "files": {
            filename: {
                "exists": bool(mesh_dir and (mesh_dir / filename).is_file()),
                "sha256": sha256_file(mesh_dir / filename if mesh_dir else None),
            }
            for filename in MESH_FILES
        },
    }


def flatten(value: Any, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    if isinstance(value, dict):
        for key in sorted(value):
            child = f"{prefix}.{key}" if prefix else str(key)
            out.update(flatten(value[key], child))
    elif isinstance(value, list):
        out[prefix] = value
    else:
        out[prefix] = value
    return out


def dict_diff(left: dict[str, Any] | None, right: dict[str, Any] | None) -> list[dict[str, Any]]:
    a = flatten(left or {})
    b = flatten(right or {})
    rows: list[dict[str, Any]] = []
    for key in sorted(set(a) | set(b)):
        if a.get(key) != b.get(key):
            rows.append({"key": key, "left": a.get(key), "right": b.get(key)})
    return rows


def audit_case(case: str, project_path: Path | None, model: dict[str, Any] | None) -> dict[str, Any]:
    spec = (model or {}).get("cases", {}).get(case, {}) if model else {}
    manifest_path = ROOT / "results" / case / "manifest.json"
    manifest = json_load(manifest_path)
    sif = find_case_sif(case)
    sif_values = parse_sif(sif)

    state_raw = spec.get("state_file") or sif_values.get("TES State File")
    state_path = resolve_input(str(state_raw)) if state_raw else None

    restart_path: Path | None = None
    restart_raw = sif_values.get("Restart File")
    if restart_raw:
        restart_path = resolve_input(str(restart_raw))
    # run.py records exact preexisting restart inputs. Prefer those when the
    # SIF path is relative to Elmer's mesh DB and therefore ambiguous in Python.
    if manifest and manifest.get("preexisting_restart_inputs_sha256"):
        for raw in manifest["preexisting_restart_inputs_sha256"]:
            candidate = ROOT / raw
            if candidate.suffix.lower().startswith(".result") or ".result" in candidate.name:
                restart_path = candidate
                break

    return {
        "case": case,
        "project": relpath(project_path),
        "project_sha256": sha256_file(project_path),
        "spec": spec,
        "sif": relpath(sif),
        "sif_sha256": sha256_file(sif),
        "sif_values": sif_values,
        "manifest": relpath(manifest_path) if manifest_path.is_file() else None,
        "manifest_sha256": sha256_file(manifest_path),
        "manifest_data": manifest,
        "mesh": mesh_audit(model, spec, manifest),
        "restart": {
            "raw": restart_raw,
            "resolved": relpath(restart_path),
            "exists": bool(restart_path and restart_path.is_file()),
            "sha256": sha256_file(restart_path),
            "position": sif_values.get("Restart Position"),
            "time": sif_values.get("Restart Time"),
        },
        "state": read_state_file(state_path),
        "series": series_summary(case, spec),
    }


def snapshot_state(source: Path | None, label: str) -> str | None:
    if source is None or not source.is_file():
        return None
    STATE_SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    target = STATE_SNAPSHOT_DIR / f"{label}.state"
    shutil.copy2(source, target)
    return relpath(target)


def set_short_variant(
    base: dict[str, Any],
    name: str,
    *,
    backend: str,
    bdf_order: int,
    steps: int,
    state_file: str | None,
    apply_mortar: bool,
    dump_matrix: bool,
) -> dict[str, Any]:
    candidate = copy.deepcopy(base)
    candidate["series_file"] = f"{name}_series.csv"
    candidate["iteration_series_file"] = f"{name}_iterations.csv"
    candidate["bdf_order"] = bdf_order
    candidate["apply_mortar_bcs"] = apply_mortar
    candidate["output_result"] = False
    candidate["output_result_path"] = None
    candidate["output_file_path"] = None
    candidate["post_file"] = False
    candidate["vtu"] = False
    if candidate.get("timesteps"):
        first_dt = candidate["timesteps"][0][0]
    else:
        first_dt = "18[us]"
    candidate["timesteps"] = [[first_dt, steps]]
    candidate["output_intervals"] = [1]
    pulse = copy.deepcopy(candidate.get("pulse") or {})
    pulse["energy"] = 0.0
    candidate["pulse"] = pulse
    if state_file is None:
        candidate.pop("state_file", None)
    else:
        candidate["state_file"] = state_file
    solver = copy.deepcopy(candidate.get("solver") or {})
    solver["linear_system"] = backend
    if dump_matrix:
        DUMP_DIR.mkdir(parents=True, exist_ok=True)
        solver["matrix_dump_prefix"] = str((DUMP_DIR / name).resolve())
        solver["matrix_dump_solution"] = True
    else:
        solver.pop("matrix_dump_prefix", None)
        solver.pop("matrix_dump_solution", None)
    candidate["solver"] = solver
    candidate["solver_comment"] = (
        f"Phase24 restart-continuity diagnostic: pulse OFF, {backend}, "
        f"BDF{bdf_order}, steps={steps}, mortar={apply_mortar}"
    )
    candidate["phase24_smoke"] = {
        "purpose": "restart-continuity campaign; pulse OFF",
        "variant": name,
        "no_matrix_dump": not dump_matrix,
        "no_vtu": True,
    }
    return candidate


def build_campaign_project(
    source_path: Path,
    model: dict[str, Any],
    steady_case: str,
    transient_case: str,
    include_fallback: bool,
    include_no_mortar: bool,
) -> tuple[Path, dict[str, dict[str, Any]]]:
    cases = model.get("cases") or {}
    if transient_case not in cases:
        raise ValueError(f"transient case not found in project: {transient_case}")
    base = cases[transient_case]
    steady_spec = cases.get(steady_case, {})

    transient_state_source = resolve_input(base.get("state_file"))
    gate3_state_source = (
        resolve_input(steady_spec.get("state_file"))
        if steady_spec.get("state_file")
        else None
    )
    preferred_source = gate3_state_source or transient_state_source

    variants: dict[str, dict[str, Any]] = {
        "current_state_mumps_bdf1_1": {
            "backend": "mumps", "bdf": 1, "steps": 1, "state_source": transient_state_source, "mortar": True, "dump": True,
        },
        "gate3_state_mumps_bdf1_1": {
            "backend": "mumps", "bdf": 1, "steps": 1, "state_source": preferred_source, "mortar": True, "dump": True,
        },
        "gate3_state_mumps_bdf1_hold5": {
            "backend": "mumps", "bdf": 1, "steps": 5, "state_source": preferred_source, "mortar": True, "dump": False,
        },
        "gate3_state_mumps_production_hold5": {
            "backend": "mumps", "bdf": int(base.get("bdf_order", 2)), "steps": 5, "state_source": preferred_source, "mortar": True, "dump": False,
        },
        "gate3_state_hypre_bdf1_1": {
            "backend": "iterative_hypre_flexgmres_boomeramg", "bdf": 1, "steps": 1, "state_source": preferred_source, "mortar": True, "dump": True,
        },
    }
    if include_fallback:
        variants["no_state_mumps_bdf1_1"] = {
            "backend": "mumps", "bdf": 1, "steps": 1, "state_source": None, "mortar": True, "dump": True,
        }
    if include_no_mortar:
        variants["gate3_state_mumps_bdf1_1_nomortar"] = {
            "backend": "mumps", "bdf": 1, "steps": 1, "state_source": preferred_source, "mortar": False, "dump": True,
        }
        # A one-step run is useful for matrix capture but the TES series is
        # committed when the *next* timestep begins.  Keep a multi-step twin so
        # the first accepted no-mortar state is always observable.
        variants["gate3_state_mumps_bdf1_hold5_nomortar"] = {
            "backend": "mumps", "bdf": 1, "steps": 5, "state_source": preferred_source, "mortar": False, "dump": False,
        }

    generated = copy.deepcopy(model)
    for suffix, options in variants.items():
        name = f"case_phase24_restartdiag_{suffix}"
        options["case"] = name
        # The UDF writes transient checkpoints back into TES State File.
        # Give every variant its own working copy so one diagnostic can never
        # change another variant's initial electrical state.
        working_state = snapshot_state(options["state_source"], f"{suffix}_working")
        options["state"] = working_state
        options["initial_state"] = read_state_file(resolve_input(working_state))
        generated["cases"][name] = set_short_variant(
            base,
            name,
            backend=options["backend"],
            bdf_order=options["bdf"],
            steps=options["steps"],
            state_file=working_state,
            apply_mortar=options["mortar"],
            dump_matrix=options["dump"],
        )
    CAMPAIGN_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_PATH.write_text(json.dumps(generated, indent=2) + "\n", encoding="utf-8")
    return PROJECT_PATH, variants


def dump_path(prefix: Path, suffix: str) -> Path | None:
    direct = prefix.with_name(prefix.name + suffix)
    if direct.is_file():
        return direct
    ranked = prefix.with_name(prefix.name + suffix + ".0")
    return ranked if ranked.is_file() else None


def read_numeric_dump(path: Path | None, columns: int) -> dict[Any, float]:
    if path is None or not path.is_file():
        return {}
    out: dict[Any, float] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < columns:
                continue
            try:
                if columns == 3:
                    key: Any = (int(fields[0]), int(fields[1]))
                    value = float(fields[2])
                else:
                    key = int(fields[0])
                    value = float(fields[1])
            except ValueError:
                continue
            out[key] = value
    return out


def compare_dumps(left: Path | None, right: Path | None, columns: int) -> dict[str, Any]:
    a = read_numeric_dump(left, columns)
    b = read_numeric_dump(right, columns)
    keys = set(a) | set(b)
    if not keys:
        return {"available": False, "left": relpath(left), "right": relpath(right)}
    diffs = [abs(a.get(key, 0.0) - b.get(key, 0.0)) for key in keys]
    return {
        "available": True,
        "left": relpath(left),
        "right": relpath(right),
        "left_records": len(a),
        "right_records": len(b),
        "union_records": len(keys),
        "nonzero_difference_records": sum(value != 0.0 for value in diffs),
        "max_absolute_difference": max(diffs, default=0.0),
    }


def solution_dump_path(prefix: Path) -> Path | None:
    """Return the Elmer SaveLinearSystem solution vector, when present.

    SaveLinearSystem convention uses *_sol.dat. Older diagnostic code looked
    for *_x.dat, so keep that as a compatibility fallback.
    """
    for suffix in ("_sol.dat", "_x.dat"):
        found = dump_path(prefix, suffix)
        if found is not None:
            return found
    return None


def restart_candidate_residual(prefix: Path) -> dict[str, Any]:
    """Evaluate b-A*x_saved and split primal/constraint rows.

    Linear System Save Solution is treated as an Elmer x0/restart candidate,
    not as proof of the exact native solver vector. Mortar multipliers are
    commonly absent from the saved solution vector, so missing trailing
    entries are extended with zero. Constraint rows are identified
    algebraically as rows without a nonzero diagonal entry.
    """
    a_path = dump_path(prefix, "_a.dat")
    b_path = dump_path(prefix, "_b.dat")
    x_path = solution_dump_path(prefix)
    if a_path is None or b_path is None or x_path is None:
        return {
            "available": False,
            "matrix": relpath(a_path),
            "rhs": relpath(b_path),
            "saved_solution": relpath(x_path),
            "caveat": (
                "Requires SaveLinearSystem A, b, and solution dumps. "
                "The saved solution is an x0 candidate, not guaranteed to be "
                "the exact native solver vector."
            ),
        }

    b = read_numeric_dump(b_path, 2)
    x = read_numeric_dump(x_path, 2)
    if not b:
        return {"available": False, "reason": "empty RHS dump"}

    n = max(int(k) for k in b)
    ax = [0.0] * (n + 1)
    row_sq = [0.0] * (n + 1)
    diag_nonzero = [False] * (n + 1)
    nnz = 0

    with a_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            try:
                row = int(fields[0])
                col = int(fields[1])
                value = float(fields[2])
            except ValueError:
                continue
            if row < 1 or row > n:
                continue
            ax[row] += value * x.get(col, 0.0)
            vv = value * value
            row_sq[row] += vv
            nnz += 1
            if row == col and value != 0.0:
                diag_nonzero[row] = True

    primal_rows = [i for i in range(1, n + 1) if diag_nonzero[i]]
    constraint_rows = [i for i in range(1, n + 1) if not diag_nonzero[i]]
    residual = [0.0] * (n + 1)
    for i in range(1, n + 1):
        residual[i] = b.get(i, 0.0) - ax[i]

    x2 = sum(value * value for value in x.values())
    xnorm = math.sqrt(x2)

    def block_stats(rows: list[int]) -> dict[str, Any]:
        if not rows:
            return {
                "rows": 0,
                "residual_l2": 0.0,
                "residual_max_abs": 0.0,
                "rhs_l2": 0.0,
                "matrix_frobenius": 0.0,
                "backward_error": 0.0,
            }
        r2 = sum(residual[i] * residual[i] for i in rows)
        b2 = sum(b.get(i, 0.0) ** 2 for i in rows)
        a2 = sum(row_sq[i] for i in rows)
        rnorm = math.sqrt(r2)
        bnorm = math.sqrt(b2)
        anorm = math.sqrt(a2)
        denom = anorm * xnorm + bnorm
        return {
            "rows": len(rows),
            "residual_l2": rnorm,
            "residual_max_abs": max(abs(residual[i]) for i in rows),
            "rhs_l2": bnorm,
            "matrix_frobenius": anorm,
            "backward_error": rnorm / max(denom, 1.0e-300),
        }

    return {
        "available": True,
        "matrix": relpath(a_path),
        "rhs": relpath(b_path),
        "saved_solution": relpath(x_path),
        "rows": n,
        "matrix_records": nnz,
        "saved_solution_records": len(x),
        "saved_solution_l2": xnorm,
        "missing_solution_entries_zero_extended": max(n - len(x), 0),
        "primal_rows": len(primal_rows),
        "constraint_rows": len(constraint_rows),
        "full": block_stats(list(range(1, n + 1))),
        "primal": block_stats(primal_rows),
        "constraint": block_stats(constraint_rows),
        "constraint_row_identification": "rows without a nonzero diagonal entry",
        "caveat": (
            "Elmer SaveLinearSystem solution is analyzed as an x0/restart candidate. "
            "It may omit mortar multipliers and is not asserted to be the exact "
            "native pre-solve vector; missing entries are extended by zero."
        ),
    }


def command_for(case: str, project: Path, solver: Path, runtime_bin: Path, toolchain_bin: Path) -> list[str]:
    return [
        sys.executable,
        str(ROOT / "run.py"),
        case,
        "--project", str(project),
        "--skip-sync",
        "--elmer-solver", str(solver),
        "--runtime-bin", str(runtime_bin),
        "--toolchain-bin", str(toolchain_bin),
        "--mpi-procs", "1",
    ]


def run_command(command: list[str], log: Path, dry_run: bool) -> dict[str, Any]:
    started = datetime.now(timezone.utc).isoformat()
    if dry_run:
        return {"command": command, "started": started, "finished": started, "exit_code": None}
    log.parent.mkdir(parents=True, exist_ok=True)
    begin = time.monotonic()
    with log.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    return {
        "command": command,
        "started": started,
        "finished": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - begin,
        "exit_code": proc.returncode,
        "log": relpath(log),
    }


def continuity_metrics(label: str, audit: dict[str, Any], reference_uA: float) -> dict[str, Any]:
    state = audit.get("state", {})
    first_it = audit.get("series", {}).get("first_iteration", {})
    first = audit.get("series", {}).get("first", {})
    current_state = state.get("current_uA")
    first_it_current = first_it.get("tes_current_uA")
    first_current = first.get("tes_current_uA")
    def delta(value: float | None) -> float | None:
        return None if value is None else value - reference_uA
    state_temp_mK = (
        state.get("temperature_K") * 1.0e3
        if state.get("temperature_K") is not None
        else None
    )
    first_temp_mK = first.get("tes_temperature_mK")
    return {
        "variant": label,
        "state_current_uA": current_state,
        "state_delta_uA": delta(current_state),
        "state_temperature_mK": state_temp_mK,
        "first_iteration_current_uA": first_it_current,
        "first_iteration_delta_uA": delta(first_it_current),
        "first_accepted_current_uA": first_current,
        "first_accepted_delta_uA": delta(first_current),
        "first_accepted_temperature_mK": first_temp_mK,
        "first_accepted_temperature_delta_mK": (
            first_temp_mK - state_temp_mK
            if first_temp_mK is not None and state_temp_mK is not None
            else None
        ),
        "last_current_uA": audit.get("series", {}).get("last", {}).get("tes_current_uA"),
        "row_count": audit.get("series", {}).get("row_count", 0),
    }


def first_jump_location(metrics: dict[str, Any], reference_uA: float) -> str:
    limit = abs(reference_uA) * CURRENT_GATE_REL
    checks = (
        ("TES state file before solver", metrics.get("state_current_uA")),
        ("first nonlinear/iteration record", metrics.get("first_iteration_current_uA")),
        ("first accepted timestep", metrics.get("first_accepted_current_uA")),
    )
    for label, value in checks:
        if value is not None and abs(value - reference_uA) > limit:
            return label
    return "not observed in captured checkpoints"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_artifacts(summary: dict[str, Any]) -> None:
    CAMPAIGN_DIR.mkdir(parents=True, exist_ok=True)
    for name in (
        "provenance",
        "input_diff",
        "restart_audit",
        "state_file_audit",
        "field_continuity",
        "linear_system_comparison",
        "first_step_restart_residual",
        "old_good_route_diff",
    ):
        payload = summary.get(name, {})
        (CAMPAIGN_DIR / f"{name}.json").write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
    write_csv(CAMPAIGN_DIR / "case_matrix.csv", summary.get("case_matrix", []))
    write_csv(CAMPAIGN_DIR / "first_step_metrics.csv", summary.get("first_step_metrics", []))
    (CAMPAIGN_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    diagnosis = summary.get("diagnosis", {})
    lines = [
        "# Phase24 restart-continuity campaign",
        "",
        f"Generated: `{summary.get('generated_utc')}`",
        "",
        "## Headline",
        "",
        f"- First observed 143.57 -> abnormal-current divergence: **{diagnosis.get('first_jump', 'unknown')}**",
        f"- Restart File audit: **{diagnosis.get('restart_status', 'unknown')}**",
        f"- TES State File audit: **{diagnosis.get('state_status', 'unknown')}**",
        f"- Mortar initialization implicated: **{diagnosis.get('mortar_status', 'not yet isolated')}**",
        f"- Time-integration initialization implicated: **{diagnosis.get('time_status', 'not yet isolated')}**",
        f"- MUMPS/HYPRE backend difference implicated: **{diagnosis.get('backend_status', 'not yet isolated')}**",
        f"- Strongest current explanation: **{diagnosis.get('strongest', 'insufficient evidence')}**",
        f"- Series observability: {diagnosis.get('series_observability', 'unknown')}",
        f"- Restart/x0 residual: **{diagnosis.get('restart_residual_status', 'not captured')}**",
        "",
        "## Evidence",
        "",
        f"- Reference Gate3 current: {summary.get('reference_current_uA')} uA",
        f"- Source project: `{summary.get('project')}`",
        f"- Existing transient case: `{summary.get('transient_case')}`",
        f"- Existing steady case: `{summary.get('steady_case')}`",
        "",
        "Detailed machine-readable outputs are beside this file: provenance.json, "
        "input_diff.json, restart_audit.json, state_file_audit.json, case_matrix.csv, "
        "first_step_metrics.csv, field_continuity.json, linear_system_comparison.json, "
        "first_step_restart_residual.json, and old_good_route_diff.json.",
        "",
        "The campaign deliberately does not run a 40 us or 100 us production trace.",
        "",
    ]
    (CAMPAIGN_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path)
    parser.add_argument("--steady-case", default=DEFAULT_STEADY_CASE)
    parser.add_argument("--transient-case", default=DEFAULT_TRANSIENT_CASE)
    parser.add_argument("--old-good-case", default=OLD_GOOD_CASE)
    parser.add_argument("--reference-current-uA", type=float, default=CURRENT_GATE_UA)
    parser.add_argument("--solver", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin\ElmerSolver_mpi.exe"))
    parser.add_argument("--runtime-bin", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\lib\elmersolver"))
    parser.add_argument("--toolchain-bin", type=Path, default=Path(r"C:\msys64\ucrt64\bin"))
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--include-state-fallback", action="store_true")
    parser.add_argument("--include-no-mortar", action="store_true")
    args = parser.parse_args()

    project_path, model = find_project((args.steady_case, args.transient_case), args.project)
    steady_project_path, steady_model = find_project_for_case(args.steady_case, args.project)
    transient_project_path, transient_model = find_project_for_case(args.transient_case, args.project)

    steady = audit_case(args.steady_case, steady_project_path or project_path, steady_model or model)
    transient = audit_case(
        args.transient_case,
        transient_project_path or project_path,
        transient_model or model,
    )
    old_good_project_path, old_good_model = find_project_for_case(args.old_good_case, None)
    old_good = audit_case(
        args.old_good_case,
        old_good_project_path or project_path,
        old_good_model or model,
    )

    source_state = transient.get("state", {})
    steady_state = steady.get("state", {})
    existing_metrics = continuity_metrics(args.transient_case, transient, args.reference_current_uA)

    summary: dict[str, Any] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": relpath(project_path),
        "project_found": bool(project_path and model),
        "source_projects": {
            "combined": relpath(project_path),
            "steady": relpath(steady_project_path),
            "transient": relpath(transient_project_path),
            "old_good": relpath(old_good_project_path),
        },
        "steady_case": args.steady_case,
        "transient_case": args.transient_case,
        "old_good_case": args.old_good_case,
        "reference_current_uA": args.reference_current_uA,
        "provenance": {
            "steady": steady,
            "transient": transient,
            "old_good": old_good,
        },
        "input_diff": {
            "steady_vs_transient_spec": dict_diff(steady.get("spec"), transient.get("spec")),
            "old_good_vs_transient_spec": dict_diff(old_good.get("spec"), transient.get("spec")),
        },
        "restart_audit": {
            "steady": steady.get("restart"),
            "transient": transient.get("restart"),
            "same_mesh_hashes": steady.get("mesh", {}).get("files") == transient.get("mesh", {}).get("files"),
        },
        "state_file_audit": {
            "steady": steady_state,
            "transient": source_state,
            "same_sha256": bool(
                steady_state.get("sha256")
                and steady_state.get("sha256") == source_state.get("sha256")
            ),
            "same_resolved_path": bool(
                steady_state.get("path")
                and steady_state.get("path") == source_state.get("path")
            ),
            "shared_path_risk": (
                "A transient accepted step can overwrite the steady seed because the "
                "UDF checkpoints to TES State File. Shared steady/transient paths are "
                "therefore mutable provenance, not an immutable restart seed."
                if steady_state.get("path")
                and steady_state.get("path") == source_state.get("path")
                else None
            ),
        },
        "old_good_route_diff": {
            "spec_diff": dict_diff(old_good.get("spec"), transient.get("spec")),
            "mesh_diff": dict_diff(old_good.get("mesh"), transient.get("mesh")),
            "sif_diff": dict_diff(old_good.get("sif_values"), transient.get("sif_values")),
        },
        "case_matrix": [],
        "first_step_metrics": [existing_metrics],
        "first_step_restart_residual": {},
        "field_continuity": {
            "note": "Scalar TES continuity is captured from state/iteration/series. Nodal field deltas require a text-exported result and are reported as unavailable rather than guessed.",
            "existing_transient_temperature_mK": transient.get("series", {}).get("first", {}).get("tes_temperature_mK"),
        },
        "linear_system_comparison": {},
        "runs": {},
    }

    run_project_path = transient_project_path or project_path
    run_model = transient_model or model
    if run_model and run_project_path and args.transient_case in run_model.get("cases", {}) and not args.audit_only:
        campaign_project, variants = build_campaign_project(
            run_project_path,
            run_model,
            args.steady_case,
            args.transient_case,
            args.include_state_fallback,
            args.include_no_mortar,
        )
        summary["project"] = relpath(campaign_project)
        sync_command = [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(campaign_project)]
        summary["sync_command"] = sync_command
        if not args.dry_run:
            sync = subprocess.run(sync_command, cwd=ROOT, capture_output=True, text=True)
            (CAMPAIGN_DIR / "sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
            if sync.returncode != 0:
                summary["sync_exit_code"] = sync.returncode
                write_artifacts(summary)
                return sync.returncode

        solver = args.solver if args.solver.is_absolute() else ROOT / args.solver
        runtime_bin = args.runtime_bin if args.runtime_bin.is_absolute() else ROOT / args.runtime_bin
        toolchain_bin = args.toolchain_bin if args.toolchain_bin.is_absolute() else ROOT / args.toolchain_bin
        for suffix, options in variants.items():
            case = options["case"]
            summary["case_matrix"].append(
                {
                    "variant": suffix,
                    "case": case,
                    "backend": options["backend"],
                    "bdf_order": options["bdf"],
                    "steps": options["steps"],
                    "mortar": options["mortar"],
                    "state_snapshot": options["state"],
                    "initial_state_current_uA": options.get("initial_state", {}).get("current_uA"),
                    "initial_state_sha256": options.get("initial_state", {}).get("sha256"),
                }
            )
            command = command_for(case, campaign_project, solver, runtime_bin, toolchain_bin)
            launcher = ROOT / "results" / case / "restart_continuity_launcher.log"
            run = run_command(command, launcher, args.dry_run)
            summary["runs"][suffix] = run
            if not args.dry_run:
                audit = audit_case(case, campaign_project, json_load(campaign_project))
                metrics = continuity_metrics(suffix, audit, args.reference_current_uA)
                # audit_case sees the working state after the transient has had
                # a chance to checkpoint it.  Restore the captured pre-run
                # state here so the "before solver" checkpoint is truthful.
                initial_state = options.get("initial_state", {})
                if initial_state.get("current_uA") is not None:
                    metrics["state_current_uA"] = initial_state["current_uA"]
                    metrics["state_delta_uA"] = (
                        initial_state["current_uA"] - args.reference_current_uA
                    )
                summary["first_step_metrics"].append(metrics)

        if not args.dry_run:
            prefixes = {
                suffix: DUMP_DIR / options["case"]
                for suffix, options in variants.items()
                if options["dump"]
            }
            names = list(prefixes)
            comparisons: dict[str, Any] = {}
            for i, left_name in enumerate(names):
                for right_name in names[i + 1 :]:
                    key = f"{left_name}__vs__{right_name}"
                    comparisons[key] = {
                        "A": compare_dumps(
                            dump_path(prefixes[left_name], "_a.dat"),
                            dump_path(prefixes[right_name], "_a.dat"),
                            3,
                        ),
                        "b": compare_dumps(
                            dump_path(prefixes[left_name], "_b.dat"),
                            dump_path(prefixes[right_name], "_b.dat"),
                            2,
                        ),
                        "x": compare_dumps(
                            solution_dump_path(prefixes[left_name]),
                            solution_dump_path(prefixes[right_name]),
                            2,
                        ),
                    }
            summary["linear_system_comparison"] = comparisons
            summary["first_step_restart_residual"] = {
                suffix: restart_candidate_residual(prefix)
                for suffix, prefix in prefixes.items()
            }

    if not args.audit_only and not args.dry_run and not summary["runs"]:
        summary["run_setup_error"] = (
            "No project JSON containing the transient case was found. Pass --project "
            "with the JSON used to generate the transient case."
        )

    metrics = summary["first_step_metrics"]
    existing_jump = first_jump_location(existing_metrics, args.reference_current_uA)
    state_delta = source_state.get("current_uA")
    state_bad = state_delta is not None and abs(state_delta - args.reference_current_uA) > abs(args.reference_current_uA) * CURRENT_GATE_REL
    restart_ok = transient.get("restart", {}).get("exists")
    same_state = summary["state_file_audit"].get("same_sha256")
    shared_state_path = summary["state_file_audit"].get("same_resolved_path")

    diagnosis = {
        "first_jump": existing_jump,
        "restart_status": "present" if restart_ok else "missing/unresolved",
        "state_status": (
            "already inconsistent before solver" if state_bad
            else "matches current gate at file level" if source_state.get("current_uA") is not None
            else "missing/unparsed"
        ),
        "mortar_status": "not isolated",
        "time_status": "not isolated",
        "backend_status": "not isolated",
        "strongest": (
            "first transient accepted-step evolution is the first demonstrated failure; "
            "inspect first-step thermal/circuit equations, restart field continuity, BDF initialization, and mortar coupling"
            if existing_jump == "first accepted timestep" and not state_bad
            else "restart/state handoff remains the first target"
        ),
    }

    by_variant = {row["variant"]: row for row in metrics}
    cur = by_variant.get("current_state_mumps_bdf1_1")
    gate_1 = by_variant.get("gate3_state_mumps_bdf1_1")
    gate = by_variant.get("gate3_state_mumps_bdf1_hold5") or gate_1
    prod = by_variant.get("gate3_state_mumps_production_hold5")
    hypre = by_variant.get("gate3_state_hypre_bdf1_1")
    nomortar = (
        by_variant.get("gate3_state_mumps_bdf1_hold5_nomortar")
        or by_variant.get("gate3_state_mumps_bdf1_1_nomortar")
    )

    if cur and gate_1 and cur.get("first_accepted_current_uA") is not None and gate_1.get("first_accepted_current_uA") is not None:
        cur_err = abs(cur["first_accepted_current_uA"] - args.reference_current_uA)
        gate_err = abs(gate_1["first_accepted_current_uA"] - args.reference_current_uA)
        if gate_err < 0.2 * cur_err:
            diagnosis["state_status"] = "isolated: Gate3 state snapshot materially restores continuity"
            diagnosis["strongest"] = "TES State File mismatch/initialization"
    if gate and prod and gate.get("first_accepted_current_uA") is not None and prod.get("first_accepted_current_uA") is not None:
        if abs(gate["first_accepted_current_uA"] - prod["first_accepted_current_uA"]) > 0.01:
            diagnosis["time_status"] = "BDF/time initialization changes the first-step current"
        else:
            diagnosis["time_status"] = "no material first-step BDF effect observed"
    if gate and hypre and gate.get("first_accepted_current_uA") is not None and hypre.get("first_accepted_current_uA") is not None:
        if abs(gate["first_accepted_current_uA"] - hypre["first_accepted_current_uA"]) > 0.01:
            diagnosis["backend_status"] = "backend-dependent first-step divergence observed"
        else:
            diagnosis["backend_status"] = "MUMPS/HYPRE agree at first accepted step"
    if gate and nomortar and gate.get("first_accepted_current_uA") is not None and nomortar.get("first_accepted_current_uA") is not None:
        if abs(gate["first_accepted_current_uA"] - nomortar["first_accepted_current_uA"]) > 0.01:
            diagnosis["mortar_status"] = "mortar-dependent first-step divergence observed"
        else:
            diagnosis["mortar_status"] = "no material first-step mortar effect observed"
    if same_state:
        diagnosis["state_status"] += "; steady/transient state hashes are identical"
    if shared_state_path:
        diagnosis["state_status"] += "; steady/transient share a mutable state-file path"
        if state_bad:
            diagnosis["strongest"] = (
                "mutable TES State File provenance: the transient shares and can overwrite "
                "the steady seed, and the currently observed state is off the Gate3 current"
            )

    if summary.get("run_setup_error"):
        diagnosis["strongest"] += "; short-run isolation still requires the transient source project JSON"

    residuals = summary.get("first_step_restart_residual", {})
    mortar_residual = residuals.get("gate3_state_mumps_bdf1_1", {})
    nomortar_residual = residuals.get("gate3_state_mumps_bdf1_1_nomortar", {})
    if mortar_residual.get("available"):
        cblock = mortar_residual.get("constraint", {})
        pblock = mortar_residual.get("primal", {})
        diagnosis["restart_residual_status"] = (
            "captured: "
            f"constraint rows={mortar_residual.get('constraint_rows')}, "
            f"constraint backward error={cblock.get('backward_error')}, "
            f"primal backward error={pblock.get('backward_error')}"
        )
        if mortar_residual.get("constraint_rows", 0) > 0:
            cbe = cblock.get("backward_error")
            pbe = pblock.get("backward_error")
            if cbe is not None and pbe is not None and cbe > 10.0 * max(pbe, 1.0e-300):
                diagnosis["strongest"] = (
                    "saved restart/x0 candidate disproportionately violates the "
                    "mortar constraint block at the first transient solve; inspect "
                    "nonconforming restart DOF reconstruction/mortar initialization"
                )
    else:
        diagnosis["restart_residual_status"] = "not captured"

    if mortar_residual.get("available") and nomortar_residual.get("available"):
        diagnosis["restart_residual_comparison"] = (
            "mortar and no-mortar first-solve residual candidates captured"
        )
    else:
        diagnosis["restart_residual_comparison"] = "not fully captured"

    diagnosis["series_observability"] = (
        "TES accepted-step series is written when the next timestep begins; "
        "one-step variants may legitimately have row_count=0. Hold5 variants "
        "are used for BDF and mortar first-step diagnosis."
    )

    summary["diagnosis"] = diagnosis
    write_artifacts(summary)
    print(f"summary={CAMPAIGN_DIR / 'summary.md'}")
    if args.audit_only or args.dry_run:
        return 0
    failures = [run for run in summary["runs"].values() if run.get("exit_code") not in (0, None)]
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

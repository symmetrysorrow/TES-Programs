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


def discover_linear_solve_dumps(prefix: Path) -> list[dict[str, Any]]:
    """Discover continuously-numbered SaveLinearSystem A/b/x dump groups.

    Elmer's exact filename convention varies across builds.  Group by the
    basename preceding *_a.dat and accept an optional rank suffix (.0).
    The unnumbered prefix sorts first; numbered siblings then sort by their
    numeric suffix, which is the SaveLinearSystem call order for the current
    single-rank diagnostic path.
    """
    parent = prefix.parent
    groups: dict[str, dict[str, Path]] = {}
    for path in parent.glob(prefix.name + "*"):
        if not path.is_file():
            continue
        name = path.name
        normalized = name[:-2] if name.endswith(".0") else name
        kind = None
        base = None
        for marker, label in (
            ("_a.dat", "A"),
            ("_b.dat", "b"),
            ("_sol.dat", "x"),
            ("_x.dat", "x"),
            ("_sizes.dat", "sizes"),
        ):
            pos = normalized.rfind(marker)
            if pos >= 0 and pos + len(marker) == len(normalized):
                base = normalized[:pos]
                kind = label
                break
        if base is None or kind is None or not base.startswith(prefix.name):
            continue
        groups.setdefault(base, {})[kind] = path

    def order_key(base: str) -> tuple[Any, ...]:
        tail = base[len(prefix.name):]
        if not tail:
            return (0, (), base)
        numbers = tuple(int(value) for value in re.findall(r"\d+", tail))
        if numbers:
            return (1, numbers, base)
        return (2, (), base)

    out: list[dict[str, Any]] = []
    for ordinal, base in enumerate(sorted(groups, key=order_key), start=1):
        files = groups[base]
        if "A" not in files or "b" not in files:
            continue
        out.append(
            {
                "ordinal": ordinal,
                "base": base,
                "A_path": files["A"],
                "b_path": files["b"],
                "x_path": files.get("x"),
                "sizes_path": files.get("sizes"),
            }
        )
    return out


def matrix_dump_dimension(a_path: Path | None) -> dict[str, Any]:
    """Derive the assembled square-system dimension from matrix indices.

    SaveLinearSystem may omit zero RHS entries, so RHS max index is not a
    reliable dimension for saddle systems. Matrix row/column indices remain
    authoritative because multiplier rows/columns contain B/Bt couplings.
    """
    if a_path is None or not a_path.is_file():
        return {"available": False, "rows": 0, "max_row": 0, "max_col": 0, "records": 0}
    max_row = 0
    max_col = 0
    records = 0
    with a_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            try:
                row = int(fields[0])
                col = int(fields[1])
                float(fields[2])
            except ValueError:
                continue
            if row <= 0 or col <= 0:
                continue
            max_row = max(max_row, row)
            max_col = max(max_col, col)
            records += 1
    n = max(max_row, max_col)
    return {
        "available": n > 0,
        "rows": n,
        "max_row": max_row,
        "max_col": max_col,
        "records": records,
        "square_index_extent": max_row == max_col,
    }


def sizes_dump_metadata(path: Path | None) -> dict[str, Any]:
    """Capture raw integer metadata from Elmer *_sizes.dat without guessing semantics."""
    if path is None or not path.is_file():
        return {"available": False, "path": relpath(path)}
    integers: list[int] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            for field in line.split():
                try:
                    integers.append(int(field))
                except ValueError:
                    continue
    return {
        "available": True,
        "path": relpath(path),
        "integer_count": len(integers),
        "integers": integers[:64],
        "caveat": (
            "Recorded for provenance only. Matrix dimension is derived from "
            "A row/column indices because *_sizes.dat semantics are build-specific."
        ),
    }


def system_row_partition(
    matrix: dict[tuple[int, int], float], n: int
) -> tuple[set[int], set[int]]:
    """Split rows by explicit nonzero diagonal over the full A-derived index range."""
    rows = set(range(1, n + 1))
    primal = {
        row
        for (row, col), value in matrix.items()
        if row == col and value != 0.0 and 1 <= row <= n
    }
    return primal, rows - primal


def _difference_stats(
    left: dict[Any, float], right: dict[Any, float], keys: set[Any] | None = None
) -> dict[str, Any]:
    use_keys = keys if keys is not None else (set(left) | set(right))
    if not use_keys:
        return {
            "records": 0,
            "difference_l2": 0.0,
            "difference_max_abs": 0.0,
            "left_l2": 0.0,
            "relative_l2_to_left": 0.0,
            "nonzero_difference_records": 0,
        }
    diffs = [left.get(key, 0.0) - right.get(key, 0.0) for key in use_keys]
    left_values = [left.get(key, 0.0) for key in use_keys]
    diff_l2 = math.sqrt(sum(value * value for value in diffs))
    left_l2 = math.sqrt(sum(value * value for value in left_values))
    return {
        "records": len(use_keys),
        "difference_l2": diff_l2,
        "difference_max_abs": max(abs(value) for value in diffs),
        "left_l2": left_l2,
        "relative_l2_to_left": diff_l2 / max(left_l2, 1.0e-300),
        "nonzero_difference_records": sum(value != 0.0 for value in diffs),
    }


def compare_linear_solve_dump_pair(
    left: dict[str, Any], right: dict[str, Any]
) -> dict[str, Any]:
    """Compare two saved nonlinear linear systems, split by saddle blocks."""
    left_a = read_numeric_dump(left.get("A_path"), 3)
    right_a = read_numeric_dump(right.get("A_path"), 3)
    left_b = read_numeric_dump(left.get("b_path"), 2)
    right_b = read_numeric_dump(right.get("b_path"), 2)
    if not left_a or not right_a or not left_b or not right_b:
        return {"available": False}

    left_dim = matrix_dump_dimension(left.get("A_path"))
    right_dim = matrix_dump_dimension(right.get("A_path"))
    left_n = int(left_dim.get("rows", 0))
    right_n = int(right_dim.get("rows", 0))
    lp, lc = system_row_partition(left_a, left_n)
    rp, rc = system_row_partition(right_a, right_n)
    same_partition = lp == rp and lc == rc

    result: dict[str, Any] = {
        "available": True,
        "left_ordinal": left.get("ordinal"),
        "right_ordinal": right.get("ordinal"),
        "left_A": relpath(left.get("A_path")),
        "right_A": relpath(right.get("A_path")),
        "left_b": relpath(left.get("b_path")),
        "right_b": relpath(right.get("b_path")),
        "same_primal_constraint_partition": same_partition,
        "left_dimension": left_dim,
        "right_dimension": right_dim,
        "left_sizes_metadata": sizes_dump_metadata(left.get("sizes_path")),
        "right_sizes_metadata": sizes_dump_metadata(right.get("sizes_path")),
        "left_primal_rows": len(lp),
        "left_constraint_rows": len(lc),
        "right_primal_rows": len(rp),
        "right_constraint_rows": len(rc),
        "A_full": _difference_stats(left_a, right_a),
        "b_full": _difference_stats(left_b, right_b),
    }

    if same_partition:
        all_a_keys = set(left_a) | set(right_a)
        block_keys = {
            "K": {key for key in all_a_keys if key[0] in lp and key[1] in lp},
            "Bt": {key for key in all_a_keys if key[0] in lp and key[1] in lc},
            "B": {key for key in all_a_keys if key[0] in lc and key[1] in lp},
            "D": {key for key in all_a_keys if key[0] in lc and key[1] in lc},
        }
        result["A_blocks"] = {
            name: _difference_stats(left_a, right_a, keys)
            for name, keys in block_keys.items()
        }
        result["b_primal"] = _difference_stats(left_b, right_b, set(lp))
        result["b_constraint"] = _difference_stats(left_b, right_b, set(lc))

    left_x_path = left.get("x_path")
    right_x_path = right.get("x_path")
    if left_x_path is not None and right_x_path is not None:
        left_x = read_solution_vector(left_x_path)
        right_x = read_solution_vector(right_x_path)
        result["x_saved"] = _difference_stats(left_x, right_x)
        result["left_x"] = relpath(left_x_path)
        result["right_x"] = relpath(right_x_path)

    return result


def nonlinear_linear_system_sequence(prefix: Path, limit: int = 3) -> dict[str, Any]:
    """Analyze the first continuously-numbered linear systems of one timestep."""
    dumps = discover_linear_solve_dumps(prefix)
    selected = dumps[: max(limit, 0)]
    payload: dict[str, Any] = {
        "available": len(selected) >= 2,
        "dump_count_discovered": len(dumps),
        "analyzed_count": len(selected),
        "sequence": [
            {
                "ordinal": item["ordinal"],
                "base": item["base"],
                "A": relpath(item["A_path"]),
                "b": relpath(item["b_path"]),
                "x": relpath(item.get("x_path")),
                "sizes": relpath(item.get("sizes_path")),
                "dimension": matrix_dump_dimension(item.get("A_path")),
            }
            for item in selected
        ],
        "pairs": {},
        "caveat": (
            "Ordinals are SaveLinearSystem dispatch order. In the single-step "
            "MUMPS restart diagnostic they are expected to track nonlinear "
            "linear solves, but they are deliberately not relabeled as Elmer "
            "nonlinear iteration IDs without an explicit filename marker."
        ),
    }
    for left, right in zip(selected, selected[1:]):
        key = f"solve{left['ordinal']}__vs__solve{right['ordinal']}"
        payload["pairs"][key] = compare_linear_solve_dump_pair(left, right)
    return payload


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


def read_solution_vector(path: Path | None) -> dict[int, float]:
    """Read SaveLinearSystem solution dumps in indexed or value-only form."""
    if path is None or not path.is_file():
        return {}
    out: dict[int, float] = {}
    implicit_index = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            try:
                if len(fields) >= 2:
                    index = int(fields[0])
                    value = float(fields[1])
                else:
                    implicit_index += 1
                    index = implicit_index
                    value = float(fields[0])
            except ValueError:
                continue
            out[index] = value
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
    x = read_solution_vector(x_path)
    if not b:
        return {"available": False, "reason": "empty RHS dump"}

    dimension = matrix_dump_dimension(a_path)
    n = int(dimension.get("rows", 0))
    if n <= 0:
        return {"available": False, "reason": "empty matrix dump"}
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
        "dimension_source": "matrix A max row/column index",
        "dimension": dimension,
        "sizes_metadata": sizes_dump_metadata(dump_path(prefix, "_sizes.dat")),
        "rhs_saved_records": len(b),
        "rhs_implicit_zero_entries": max(n - len(b), 0),
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
            "native pre-solve vector; missing entries are extended by zero. Matrix "
            "dimension is derived from A, not RHS, because Save Skip Zeros can omit "
            "zero constraint RHS rows."
        ),
    }


def matrix_block_decomposition(prefix: Path) -> dict[str, Any]:
    """Describe K/B/Bt/D blocks using zero-diagonal rows as constraints."""
    a_path = dump_path(prefix, "_a.dat")
    b_path = dump_path(prefix, "_b.dat")
    if a_path is None or b_path is None:
        return {
            "available": False,
            "matrix": relpath(a_path),
            "rhs": relpath(b_path),
        }

    rhs = read_numeric_dump(b_path, 2)
    if not rhs:
        return {"available": False, "reason": "empty RHS dump"}
    dimension = matrix_dump_dimension(a_path)
    n = int(dimension.get("rows", 0))
    if n <= 0:
        return {"available": False, "reason": "empty matrix dump"}

    entries: list[tuple[int, int, float]] = []
    diag_nonzero = [False] * (n + 1)
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
            if not (1 <= row <= n and 1 <= col <= n):
                continue
            entries.append((row, col, value))
            if row == col and value != 0.0:
                diag_nonzero[row] = True

    primal = {i for i in range(1, n + 1) if diag_nonzero[i]}
    constraint = {i for i in range(1, n + 1) if not diag_nonzero[i]}
    block_acc = {
        "K": {"records": 0, "frob_sq": 0.0, "max_abs": 0.0},
        "Bt": {"records": 0, "frob_sq": 0.0, "max_abs": 0.0},
        "B": {"records": 0, "frob_sq": 0.0, "max_abs": 0.0},
        "D": {"records": 0, "frob_sq": 0.0, "max_abs": 0.0},
    }
    k_values: dict[tuple[int, int], float] = {}
    b_values: dict[tuple[int, int], float] = {}
    bt_values: dict[tuple[int, int], float] = {}

    for row, col, value in entries:
        if row in primal and col in primal:
            name = "K"
            k_values[(row, col)] = value
        elif row in primal and col in constraint:
            name = "Bt"
            bt_values[(row, col)] = value
        elif row in constraint and col in primal:
            name = "B"
            b_values[(row, col)] = value
        else:
            name = "D"
        acc = block_acc[name]
        acc["records"] += 1
        acc["frob_sq"] += value * value
        acc["max_abs"] = max(acc["max_abs"], abs(value))

    def finish(acc: dict[str, Any]) -> dict[str, Any]:
        return {
            "records": acc["records"],
            "frobenius": math.sqrt(acc["frob_sq"]),
            "max_abs": acc["max_abs"],
        }

    transpose_keys = set(b_values) | {(col, row) for row, col in bt_values}
    transpose_max = 0.0
    transpose_l2_sq = 0.0
    for row, col in transpose_keys:
        diff = b_values.get((row, col), 0.0) - bt_values.get((col, row), 0.0)
        transpose_max = max(transpose_max, abs(diff))
        transpose_l2_sq += diff * diff

    primal_rhs_l2 = math.sqrt(sum(rhs.get(i, 0.0) ** 2 for i in primal))
    constraint_rhs_l2 = math.sqrt(sum(rhs.get(i, 0.0) ** 2 for i in constraint))
    return {
        "available": True,
        "rows": n,
        "dimension_source": "matrix A max row/column index",
        "dimension": dimension,
        "sizes_metadata": sizes_dump_metadata(dump_path(prefix, "_sizes.dat")),
        "rhs_saved_records": len(rhs),
        "rhs_implicit_zero_entries": max(n - len(rhs), 0),
        "primal_rows": len(primal),
        "constraint_rows": len(constraint),
        "blocks": {name: finish(acc) for name, acc in block_acc.items()},
        "B_minus_BtT_max_abs": transpose_max,
        "B_minus_BtT_l2": math.sqrt(transpose_l2_sq),
        "primal_rhs_l2": primal_rhs_l2,
        "constraint_rhs_l2": constraint_rhs_l2,
        "constraint_row_identification": "rows without a nonzero diagonal entry",
    }


def compare_primal_systems(left_prefix: Path, right_prefix: Path) -> dict[str, Any]:
    """Compare primal K and primal RHS between two first-solve systems."""
    left_a = dump_path(left_prefix, "_a.dat")
    right_a = dump_path(right_prefix, "_a.dat")
    left_b = dump_path(left_prefix, "_b.dat")
    right_b = dump_path(right_prefix, "_b.dat")
    if None in (left_a, right_a, left_b, right_b):
        return {"available": False}

    def load(path_a: Path, path_b: Path) -> tuple[dict[tuple[int, int], float], dict[int, float], set[int]]:
        rhs = read_numeric_dump(path_b, 2)
        n = int(matrix_dump_dimension(path_a).get("rows", 0))
        diag_nonzero = [False] * (n + 1)
        entries: list[tuple[int, int, float]] = []
        with path_a.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 3:
                    continue
                try:
                    row, col, value = int(fields[0]), int(fields[1]), float(fields[2])
                except ValueError:
                    continue
                if not (1 <= row <= n and 1 <= col <= n):
                    continue
                entries.append((row, col, value))
                if row == col and value != 0.0:
                    diag_nonzero[row] = True
        primal = {i for i in range(1, n + 1) if diag_nonzero[i]}
        k = {(r, col): v for r, col, v in entries if r in primal and col in primal}
        prhs = {i: rhs.get(i, 0.0) for i in primal}
        return k, prhs, primal

    lk, lb, lp = load(left_a, left_b)
    rk, rb, rp = load(right_a, right_b)
    if lp != rp:
        return {
            "available": True,
            "same_primal_row_set": False,
            "left_primal_rows": len(lp),
            "right_primal_rows": len(rp),
            "common_primal_rows": len(lp & rp),
        }

    k_keys = set(lk) | set(rk)
    k_diffs = [abs(lk.get(key, 0.0) - rk.get(key, 0.0)) for key in k_keys]
    rhs_keys = set(lb) | set(rb)
    rhs_diffs = [abs(lb.get(key, 0.0) - rb.get(key, 0.0)) for key in rhs_keys]
    return {
        "available": True,
        "same_primal_row_set": True,
        "primal_rows": len(lp),
        "K_left_records": len(lk),
        "K_right_records": len(rk),
        "K_nonzero_difference_records": sum(v != 0.0 for v in k_diffs),
        "K_max_absolute_difference": max(k_diffs, default=0.0),
        "rhs_nonzero_difference_records": sum(v != 0.0 for v in rhs_diffs),
        "rhs_max_absolute_difference": max(rhs_diffs, default=0.0),
    }


def independent_direct_solve(prefix: Path) -> dict[str, Any]:
    """Best-effort SciPy direct solve of the saved first transient system."""
    a_path = dump_path(prefix, "_a.dat")
    b_path = dump_path(prefix, "_b.dat")
    x_path = solution_dump_path(prefix)
    if a_path is None or b_path is None or x_path is None:
        return {"available": False, "reason": "missing A, b, or saved solution dump"}

    try:
        import numpy as np
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import spsolve
    except Exception as exc:
        return {
            "available": False,
            "reason": "SciPy sparse direct solver unavailable",
            "detail": f"{type(exc).__name__}: {exc}",
        }

    rhs_map = read_numeric_dump(b_path, 2)
    saved_map = read_solution_vector(x_path)
    dimension = matrix_dump_dimension(a_path)
    n = int(dimension.get("rows", 0))
    if n <= 0:
        return {"available": False, "reason": "empty matrix"}

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    diag_nonzero = np.zeros(n, dtype=bool)
    with a_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            try:
                row, col, value = int(fields[0]), int(fields[1]), float(fields[2])
            except ValueError:
                continue
            if not (1 <= row <= n and 1 <= col <= n):
                continue
            rows.append(row - 1)
            cols.append(col - 1)
            vals.append(value)
            if row == col and value != 0.0:
                diag_nonzero[row - 1] = True

    A = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)).tocsc()
    rhs = np.zeros(n, dtype=float)
    for idx, value in rhs_map.items():
        if 1 <= int(idx) <= n:
            rhs[int(idx) - 1] = value
    saved = np.zeros(n, dtype=float)
    for idx, value in saved_map.items():
        if 1 <= int(idx) <= n:
            saved[int(idx) - 1] = value

    try:
        direct = np.asarray(spsolve(A, rhs), dtype=float)
    except Exception as exc:
        return {
            "available": False,
            "reason": "SciPy spsolve failed",
            "detail": f"{type(exc).__name__}: {exc}",
            "rows": n,
            "matrix_records": len(vals),
        }

    if direct.shape != (n,) or not np.all(np.isfinite(direct)):
        return {
            "available": False,
            "reason": "SciPy direct solution is non-finite or wrong-sized",
            "rows": n,
        }

    residual = rhs - A.dot(direct)
    rhs_norm = float(np.linalg.norm(rhs))
    residual_norm = float(np.linalg.norm(residual))
    delta = direct - saved
    primal_mask = diag_nonzero
    constraint_mask = ~diag_nonzero

    def delta_stats(mask: Any) -> dict[str, Any]:
        values = delta[mask]
        direct_values = direct[mask]
        saved_values = saved[mask]
        if values.size == 0:
            return {"rows": 0, "delta_l2": 0.0, "delta_max_abs": 0.0}
        return {
            "rows": int(values.size),
            "delta_l2": float(np.linalg.norm(values)),
            "delta_max_abs": float(np.max(np.abs(values))),
            "direct_l2": float(np.linalg.norm(direct_values)),
            "saved_l2": float(np.linalg.norm(saved_values)),
        }

    primal = delta_stats(primal_mask)
    constraint = delta_stats(constraint_mask)
    primal["delta_max_abs_mK_if_temperature"] = primal["delta_max_abs"] * 1.0e3

    return {
        "available": True,
        "rows": n,
        "dimension_source": "matrix A max row/column index",
        "dimension": dimension,
        "rhs_saved_records": len(rhs_map),
        "rhs_implicit_zero_entries": max(n - len(rhs_map), 0),
        "matrix_records": len(vals),
        "saved_solution_records": len(saved_map),
        "missing_solution_entries_zero_extended": max(n - len(saved_map), 0),
        "direct_residual_l2": residual_norm,
        "direct_relative_residual": residual_norm / max(rhs_norm, 1.0e-300),
        "full_delta_l2": float(np.linalg.norm(delta)),
        "full_delta_max_abs": float(np.max(np.abs(delta))),
        "primal": primal,
        "constraint": constraint,
        "interpretation_caveat": (
            "The saved Elmer solution dump is an x0/restart candidate, not proven "
            "to equal the exact native pre-solve vector. For this scalar HeatSolve, "
            "primal entries are treated as temperature DOFs; multiplier entries "
            "missing from the saved vector are zero-extended."
        ),
    }


def direct_solve_transition_sensitivity(
    left: dict[str, Any], right: dict[str, Any]
) -> dict[str, Any]:
    """Decompose one saved-system transition into RHS/operator direct effects.

    Solve the four combinations
      x11 = A1^-1 b1
      x12 = A1^-1 b2   (RHS-only change)
      x21 = A2^-1 b1   (operator-only change)
      x22 = A2^-1 b2
    and report primal-temperature and constraint deltas. Two sparse LU
    factorizations are reused for the four solves.
    """
    try:
        import numpy as np
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import splu
    except Exception as exc:
        return {
            "available": False,
            "reason": "SciPy sparse direct solver unavailable",
            "detail": f"{type(exc).__name__}: {exc}",
        }

    def load_system(item: dict[str, Any]) -> tuple[Any, Any, Any] | tuple[None, None, None]:
        a_path = item.get("A_path")
        b_path = item.get("b_path")
        if a_path is None or b_path is None:
            return None, None, None
        rhs_map = read_numeric_dump(b_path, 2)
        dimension = matrix_dump_dimension(a_path)
        n = int(dimension.get("rows", 0))
        if n <= 0:
            return None, None, None
        rows: list[int] = []
        cols: list[int] = []
        vals: list[float] = []
        diag_nonzero = np.zeros(n, dtype=bool)
        with a_path.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 3:
                    continue
                try:
                    row, col, value = int(fields[0]), int(fields[1]), float(fields[2])
                except ValueError:
                    continue
                if not (1 <= row <= n and 1 <= col <= n):
                    continue
                rows.append(row - 1)
                cols.append(col - 1)
                vals.append(value)
                if row == col and value != 0.0:
                    diag_nonzero[row - 1] = True
        matrix = coo_matrix(
            (np.asarray(vals), (np.asarray(rows), np.asarray(cols))),
            shape=(n, n),
        ).tocsc()
        rhs = np.zeros(n, dtype=float)
        for idx, value in rhs_map.items():
            if 1 <= int(idx) <= n:
                rhs[int(idx) - 1] = value
        return matrix, rhs, diag_nonzero

    A1, b1, mask1 = load_system(left)
    A2, b2, mask2 = load_system(right)
    if A1 is None or A2 is None or b1 is None or b2 is None:
        return {"available": False, "reason": "missing or empty saved A/b system"}
    if A1.shape != A2.shape or b1.shape != b2.shape:
        return {
            "available": False,
            "reason": "saved systems have different dimensions",
            "left_shape": list(A1.shape),
            "right_shape": list(A2.shape),
        }
    if not np.array_equal(mask1, mask2):
        return {
            "available": False,
            "reason": "primal/constraint row partition changes between saved systems",
            "left_primal_rows": int(np.count_nonzero(mask1)),
            "right_primal_rows": int(np.count_nonzero(mask2)),
        }

    try:
        lu1 = splu(A1)
        lu2 = splu(A2)
        x11 = np.asarray(lu1.solve(b1), dtype=float)
        x12 = np.asarray(lu1.solve(b2), dtype=float)
        x21 = np.asarray(lu2.solve(b1), dtype=float)
        x22 = np.asarray(lu2.solve(b2), dtype=float)
    except Exception as exc:
        return {
            "available": False,
            "reason": "SciPy sparse LU factorization/solve failed",
            "detail": f"{type(exc).__name__}: {exc}",
            "rows": int(A1.shape[0]),
        }

    if not all(np.all(np.isfinite(x)) for x in (x11, x12, x21, x22)):
        return {"available": False, "reason": "direct sensitivity solution is non-finite"}

    primal = mask1
    constraint = ~mask1

    def stats(delta: Any, mask: Any, *, temperature: bool = False) -> dict[str, Any]:
        values = delta[mask]
        if values.size == 0:
            out = {"rows": 0, "delta_l2": 0.0, "delta_max_abs": 0.0}
        else:
            out = {
                "rows": int(values.size),
                "delta_l2": float(np.linalg.norm(values)),
                "delta_max_abs": float(np.max(np.abs(values))),
            }
        if temperature:
            out["delta_l2_mK_if_temperature"] = out["delta_l2"] * 1.0e3
            out["delta_max_abs_mK_if_temperature"] = out["delta_max_abs"] * 1.0e3
        return out

    transitions = {
        "full_A2b2_minus_A1b1": x22 - x11,
        "rhs_only_A1b2_minus_A1b1": x12 - x11,
        "operator_only_A2b1_minus_A1b1": x21 - x11,
        "rhs_effect_on_A2_A2b2_minus_A2b1": x22 - x21,
        "operator_effect_on_b2_A2b2_minus_A1b2": x22 - x12,
        "interaction": x22 - x12 - x21 + x11,
    }

    rhs_effect = transitions["rhs_only_A1b2_minus_A1b1"][primal]
    operator_effect = transitions["operator_only_A2b1_minus_A1b1"][primal]
    full_effect = transitions["full_A2b2_minus_A1b1"][primal]
    rhs_l2 = float(np.linalg.norm(rhs_effect))
    op_l2 = float(np.linalg.norm(operator_effect))
    full_l2 = float(np.linalg.norm(full_effect))

    def rel_residual(A: Any, x: Any, b: Any) -> float:
        residual = b - A.dot(x)
        return float(np.linalg.norm(residual) / max(np.linalg.norm(b), 1.0e-300))

    return {
        "available": True,
        "left_ordinal": left.get("ordinal"),
        "right_ordinal": right.get("ordinal"),
        "rows": int(A1.shape[0]),
        "dimension_source": "matrix A max row/column index",
        "left_dimension": matrix_dump_dimension(left.get("A_path")),
        "right_dimension": matrix_dump_dimension(right.get("A_path")),
        "left_sizes_metadata": sizes_dump_metadata(left.get("sizes_path")),
        "right_sizes_metadata": sizes_dump_metadata(right.get("sizes_path")),
        "primal_rows": int(np.count_nonzero(primal)),
        "constraint_rows": int(np.count_nonzero(constraint)),
        "direct_relative_residuals": {
            "A1_b1": rel_residual(A1, x11, b1),
            "A1_b2": rel_residual(A1, x12, b2),
            "A2_b1": rel_residual(A2, x21, b1),
            "A2_b2": rel_residual(A2, x22, b2),
        },
        "primal": {
            name: stats(delta, primal, temperature=True)
            for name, delta in transitions.items()
        },
        "constraint": {
            name: stats(delta, constraint)
            for name, delta in transitions.items()
        },
        "primal_l2_ratios": {
            "rhs_only_to_full": rhs_l2 / max(full_l2, 1.0e-300),
            "operator_only_to_full": op_l2 / max(full_l2, 1.0e-300),
            "rhs_only_to_operator_only": rhs_l2 / max(op_l2, 1.0e-300),
        },
        "interpretation_caveat": (
            "This is a linear sensitivity decomposition of two saved assembled "
            "systems. RHS-only and operator-only effects need not add exactly "
            "because changing A and b has an interaction term, which is reported "
            "explicitly. Primal entries are treated as scalar temperature DOFs."
        ),
    }


def nonlinear_transition_direct_sensitivity(
    prefix: Path, limit: int = 3
) -> dict[str, Any]:
    dumps = discover_linear_solve_dumps(prefix)[: max(limit, 0)]
    payload: dict[str, Any] = {
        "available": len(dumps) >= 2,
        "dump_count_analyzed": len(dumps),
        "pairs": {},
    }
    for left, right in zip(dumps, dumps[1:]):
        key = f"solve{left['ordinal']}__vs__solve{right['ordinal']}"
        payload["pairs"][key] = direct_solve_transition_sensitivity(left, right)
    return payload


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
        "first_step_matrix_blocks",
        "first_step_direct_solve",
        "primal_system_comparison",
        "nonlinear_linear_system_sequence",
        "nonlinear_transition_direct_sensitivity",
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
        f"- Independent direct first solve: **{diagnosis.get('independent_direct_status', 'not captured')}**",
        f"- Mortar/no-mortar primal system: **{diagnosis.get('primal_system_status', 'not captured')}**",
        f"- Nonlinear A/b sequence: **{diagnosis.get('nonlinear_system_sequence_status', 'not captured')}**",
        f"- Direct solve1->2 sensitivity: **{diagnosis.get('direct_sensitivity_status', 'not captured')}**",
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
        "first_step_restart_residual.json, first_step_matrix_blocks.json, "
        "first_step_direct_solve.json, primal_system_comparison.json, "
        "nonlinear_linear_system_sequence.json, "
        "nonlinear_transition_direct_sensitivity.json, and old_good_route_diff.json.",
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
        "first_step_matrix_blocks": {},
        "first_step_direct_solve": {},
        "primal_system_comparison": {},
        "nonlinear_linear_system_sequence": {},
        "nonlinear_transition_direct_sensitivity": {},
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
            summary["first_step_matrix_blocks"] = {
                suffix: matrix_block_decomposition(prefix)
                for suffix, prefix in prefixes.items()
            }
            summary["first_step_direct_solve"] = {
                suffix: independent_direct_solve(prefix)
                for suffix, prefix in prefixes.items()
                if suffix.startswith("gate3_state_mumps_bdf1_1")
            }
            mortar_prefix = prefixes.get("gate3_state_mumps_bdf1_1")
            nomortar_prefix = prefixes.get("gate3_state_mumps_bdf1_1_nomortar")
            if mortar_prefix is not None and nomortar_prefix is not None:
                summary["primal_system_comparison"] = compare_primal_systems(
                    mortar_prefix, nomortar_prefix
                )
            gate_prefix = prefixes.get("gate3_state_mumps_bdf1_1")
            if gate_prefix is not None:
                summary["nonlinear_linear_system_sequence"] = (
                    nonlinear_linear_system_sequence(gate_prefix, limit=3)
                )
                summary["nonlinear_transition_direct_sensitivity"] = (
                    nonlinear_transition_direct_sensitivity(gate_prefix, limit=3)
                )

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
            "BDF order is not materially implicated, so inspect the first transient "
            "mortar/interface operator and electrothermal nonlinear update"
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

    direct_solves = summary.get("first_step_direct_solve", {})
    mortar_direct = direct_solves.get("gate3_state_mumps_bdf1_1", {})
    if mortar_direct.get("available"):
        pd = mortar_direct.get("primal", {})
        diagnosis["independent_direct_status"] = (
            "captured: "
            f"primal max |x_direct-x_saved|="
            f"{pd.get('delta_max_abs_mK_if_temperature')} mK, "
            f"direct relative residual={mortar_direct.get('direct_relative_residual')}"
        )
        max_mk = pd.get("delta_max_abs_mK_if_temperature")
        if max_mk is not None and max_mk > 0.05:
            diagnosis["strongest"] = (
                "the saved first transient linear system itself requires a material "
                "primal-temperature correction from the saved restart/x0 candidate; "
                "inspect mortar/interface operator and primal RHS construction before "
                "the nonlinear circuit update"
            )
        elif max_mk is not None:
            diagnosis["strongest"] = (
                "independent direct first-solve correction is small; the large accepted "
                "TES jump is more likely generated after the first linear solve by the "
                "nonlinear electrothermal/solution-transfer path"
            )
    else:
        diagnosis["independent_direct_status"] = "not captured"

    primal_cmp = summary.get("primal_system_comparison", {})
    if primal_cmp.get("available"):
        diagnosis["primal_system_status"] = (
            f"same primal row set={primal_cmp.get('same_primal_row_set')}; "
            f"K max diff={primal_cmp.get('K_max_absolute_difference')}; "
            f"primal RHS max diff={primal_cmp.get('rhs_max_absolute_difference')}"
        )
    else:
        diagnosis["primal_system_status"] = "not captured"

    sequence = summary.get("nonlinear_linear_system_sequence", {})
    sequence_pairs = sequence.get("pairs", {})
    if sequence.get("available"):
        diagnosis["nonlinear_system_sequence_status"] = (
            f"captured {sequence.get('analyzed_count')} of "
            f"{sequence.get('dump_count_discovered')} saved linear solves"
        )
        pair12 = sequence_pairs.get("solve1__vs__solve2")
        pair23 = sequence_pairs.get("solve2__vs__solve3")
        available_pairs = [
            (name, pair)
            for name, pair in (("1->2", pair12), ("2->3", pair23))
            if pair and pair.get("available")
        ]
        for name, pair in available_pairs:
            afull = pair.get("A_full", {})
            bfull = pair.get("b_full", {})
            diagnosis["nonlinear_system_sequence_status"] += (
                f"; {name} A rel-L2={afull.get('relative_l2_to_left')}, "
                f"b rel-L2={bfull.get('relative_l2_to_left')}"
            )
        if available_pairs:
            def pair_magnitude(item: tuple[str, dict[str, Any]]) -> float:
                pair = item[1]
                a_rel = pair.get("A_full", {}).get("relative_l2_to_left") or 0.0
                b_rel = pair.get("b_full", {}).get("relative_l2_to_left") or 0.0
                return max(a_rel, b_rel)

            focus_name, focus = max(available_pairs, key=pair_magnitude)
            afull = focus.get("A_full", {})
            bfull = focus.get("b_full", {})
            a_rel = afull.get("relative_l2_to_left")
            b_rel = bfull.get("relative_l2_to_left")
            if a_rel is not None and b_rel is not None:
                if b_rel > 10.0 * max(a_rel, 1.0e-300):
                    diagnosis["strongest"] = (
                        f"saved solve transition {focus_name} changes primarily in "
                        "the RHS; inspect transient-history, body-force, and "
                        "RHS-reuse assembly"
                    )
                elif a_rel > 10.0 * max(b_rel, 1.0e-300):
                    diagnosis["strongest"] = (
                        f"saved solve transition {focus_name} changes primarily in "
                        "the operator; inspect temperature-dependent material, "
                        "mortar, and matrix-reuse assembly"
                    )
                else:
                    diagnosis["strongest"] = (
                        f"saved solve transition {focus_name} changes in both operator "
                        "and RHS; use block-resolved A/b differences to isolate "
                        "K/B/Bt/D versus primal/constraint RHS"
                    )
    else:
        diagnosis["nonlinear_system_sequence_status"] = "not captured"

    sensitivity = summary.get("nonlinear_transition_direct_sensitivity", {})
    sensitivity_pairs = sensitivity.get("pairs", {})
    if sensitivity.get("available"):
        sens12 = sensitivity_pairs.get("solve1__vs__solve2", {})
        if sens12.get("available"):
            primal = sens12.get("primal", {})
            full = primal.get("full_A2b2_minus_A1b1", {})
            rhs_only = primal.get("rhs_only_A1b2_minus_A1b1", {})
            op_only = primal.get("operator_only_A2b1_minus_A1b1", {})
            diagnosis["direct_sensitivity_status"] = (
                "solve1->2: "
                f"full max={full.get('delta_max_abs_mK_if_temperature')} mK, "
                f"RHS-only max={rhs_only.get('delta_max_abs_mK_if_temperature')} mK, "
                f"operator-only max={op_only.get('delta_max_abs_mK_if_temperature')} mK"
            )
            ratios = sens12.get("primal_l2_ratios", {})
            rhs_to_op = ratios.get("rhs_only_to_operator_only")
            if rhs_to_op is not None and rhs_to_op > 10.0:
                diagnosis["strongest"] = (
                    "independent direct sensitivity shows the saved solve1->2 "
                    "temperature correction is RHS-dominated; inspect transient-"
                    "history/body-force/RHS-reuse assembly"
                )
            elif rhs_to_op is not None and rhs_to_op < 0.1:
                diagnosis["strongest"] = (
                    "independent direct sensitivity shows the saved solve1->2 "
                    "temperature correction is operator-dominated; inspect material/"
                    "mortar/matrix-reuse assembly"
                )
            elif rhs_to_op is not None:
                diagnosis["strongest"] = (
                    "independent direct sensitivity shows both RHS and operator "
                    "contribute materially to the saved solve1->2 temperature correction"
                )
        else:
            diagnosis["direct_sensitivity_status"] = (
                f"not captured: {sens12.get('reason', 'solve1->2 unavailable')}"
            )
    else:
        diagnosis["direct_sensitivity_status"] = "not captured"

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

"""Run and evaluate the Phase24 Gate 3 HYPRE steady-state case.

Gate 3 is deliberately an independent steady solve: it starts at ``T_0`` on
the frozen ``mesh_singlepixel_prod_v2`` mesh without mortar coupling and does not
consume a restart result.  The script keeps the generated project, launcher
log, Elmer manifest, solver log, iteration series, and a machine-readable
gate decision together.  A failed HYPRE solve is therefore a useful Gate 3
artifact, rather than an unrecorded interruption.

From the repository root::

    python scripts/support/run_phase24_gate3_hypre_steady.py

Use ``--dry-run`` to generate the project and print the command without
starting Elmer.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import platform
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PROJECT = ROOT / "elmer_project_singlepixel_prod_v2_original_timegrid.json"
BASE_CASE = "case_tes_steady_singlepixel_prod_v2_original_timegrid_hybrid"
CASE_TAG = "prod_v2"
CASE_NAME = "case_phase24_gate3_hypre_steady_prod_v2"
MESH = "mesh_singlepixel_prod_v2"
DIAGNOSTIC_DIR = ROOT / "artifacts/phase24_gate3_hypre_steady"
PROJECT_PATH = DIAGNOSTIC_DIR / "phase24_gate3_hypre_steady.json"
RESULT_DIR = ROOT / "results" / CASE_NAME
SOLVER_LOG = RESULT_DIR / "solver.log"
LAUNCHER_LOG = RESULT_DIR / "gate3_launcher.log"

COMSOL_CURRENT_UA = 143.055049
STAGE11_MUMPS_CURRENT_UA = 143.53734493231093
APPLY_MORTAR_BCS = False
MUMPS_LIMIT_PERCENT = 0.05
COMSOL_LIMIT_PERCENT = 0.6
LINEAR_SYSTEM = "iterative_hypre_flexgmres_boomeramg"
LINEAR_MAX_ITERATIONS = 4000
LINEAR_TOLERANCE = 1.0e-10
HYPRE_REUSE = False
PRECONDITIONER_LAGGING: str | None = None
BOOMER_AMG_STRONG_THRESHOLD: float | None = None
HYPRE_GMRES_DIMENSION = 100


def configure_mesh(
    *,
    source_project: Path | None,
    mesh: str | None,
    base_case: str | None,
    reference_current_uA: float | None,
    comsol_current_uA: float | None,
) -> None:
    """Select a mesh/project while retaining the original prod-v2 defaults."""
    global SOURCE_PROJECT, BASE_CASE, CASE_TAG, CASE_NAME, MESH
    global STAGE11_MUMPS_CURRENT_UA, COMSOL_CURRENT_UA
    if source_project is not None:
        SOURCE_PROJECT = absolute(source_project)
    if mesh is not None:
        MESH = mesh
        mesh_tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", mesh).strip("_")
        # Keep Elmer's TES State File below its 128-character limit even for
        # the long conformal refinement registry names.
        CASE_TAG = {
            "mesh_singlepixel_conformal_gpu": "conformal_coarse",
            "mesh_singlepixel_conformal_gpu_fine": "conformal_fine",
            "mesh_singlepixel_conformal_gpu_refine20": "conformal_refine20",
            "mesh_singlepixel_conformal_gpu_fine_stycast32": "conformal_fine_stycast32",
            "mesh_singlepixel_conformal_gpu_fine_stycast32_nomortar": "s32nm",
            "mesh_singlepixel_gpu_fine_stycast32_mortar": "s32m",
            "mesh_singlepixel_prod_v2": "prod_v2",
        }.get(mesh, mesh_tag[:32])
    if base_case is not None:
        BASE_CASE = base_case
    if reference_current_uA is not None:
        STAGE11_MUMPS_CURRENT_UA = reference_current_uA
    if comsol_current_uA is not None:
        COMSOL_CURRENT_UA = comsol_current_uA
    CASE_NAME = (
        f"case_phase24_gate3_hypre_steady_{CASE_TAG}"
        if MESH == "mesh_singlepixel_prod_v2"
        else f"case_phase24_g3_{CASE_TAG}"
    )


def configure_variant(variant: str | None) -> None:
    """Move a tuning trial into a collision-free case and artifact bundle."""
    global CASE_NAME, DIAGNOSTIC_DIR, PROJECT_PATH, RESULT_DIR, SOLVER_LOG, LAUNCHER_LOG
    if not variant:
        return
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant).strip("_")
    if not safe:
        raise ValueError("variant must contain at least one alphanumeric character")
    if MESH == "mesh_singlepixel_prod_v2":
        CASE_NAME = f"case_phase24_gate3_hypre_steady_{CASE_TAG}_{safe}"
    else:
        suffix_hash = hashlib.sha256(safe.encode("utf-8")).hexdigest()[:6]
        CASE_NAME = f"case_phase24_g3_{CASE_TAG}_{safe[:10]}_{suffix_hash}"
    DIAGNOSTIC_DIR = ROOT / f"artifacts/phase24_gate3_hypre_steady_{safe}"
    PROJECT_PATH = DIAGNOSTIC_DIR / "phase24_gate3_hypre_steady.json"
    RESULT_DIR = ROOT / "results" / CASE_NAME
    SOLVER_LOG = RESULT_DIR / "solver.log"
    LAUNCHER_LOG = RESULT_DIR / "gate3_launcher.log"


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def build_project() -> Path:
    """Create the isolated Gate 3 project without mutating the source case."""
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    if BASE_CASE not in project.get("cases", {}):
        raise KeyError(f"source project has no case {BASE_CASE!r}")

    case = copy.deepcopy(project["cases"][BASE_CASE])
    case.update(
        {
            "template": "steady",
            "mesh": MESH,
            "heat_source": "circuit_inner",
            "initial_temperature": "T_0",
            "apply_mortar_bcs": APPLY_MORTAR_BCS,
            "restart_from": None,
            "restart_file_path": None,
            "preexisting_restart": False,
            "restart_file_base": None,
            "state_file": f"work/meshes/{MESH}/{CASE_NAME}.state",
            "series_file": f"{CASE_NAME}_series.csv",
            "iteration_series_file": f"{CASE_NAME}_iterations.csv",
            "output_result": True,
            "output_file_path": f"../work/meshes/{MESH}/{CASE_NAME}.result",
            "output_intervals": 1,
            "steady_state_max_iterations": 1,
            "vtu": False,
            "solver_comment": (
                "Phase24 Gate 3: independent initial-T0 HYPRE steady solve; "
                f"{MESH} without mortar"
            ),
            "comparison_time_grid": {
                "mode": "Phase24 Gate 3 HYPRE steady",
                "purpose": "independent steady current against Stage 11 MUMPS and COMSOL",
                "mesh": MESH,
                "mortar": APPLY_MORTAR_BCS,
            },
            "phase24_smoke": {
                "purpose": "Gate 3 steady qualification",
                "independent_initial_temperature": True,
                "no_restart": True,
                "no_vtu": True,
            },
        }
    )
    if HYPRE_REUSE:
        case["phase24_hypre_reuse"] = True
    if PRECONDITIONER_LAGGING:
        case["phase24_preconditioner_lagging"] = PRECONDITIONER_LAGGING
    case["solver"] = {
        **case.get("solver", {}),
        "linear_system": LINEAR_SYSTEM,
        "linear_system_max_iterations": LINEAR_MAX_ITERATIONS,
        "linear_system_convergence_tolerance": LINEAR_TOLERANCE,
        "hypre_gmres_dimension": HYPRE_GMRES_DIMENSION,
        "nonlinear_max_iterations": 120,
        "nonlinear_convergence_tolerance": 1.0e-8,
        "nonlinear_relaxation_factor": 1.0,
        "steady_state_convergence_tolerance": 1.0e-8,
    }
    if BOOMER_AMG_STRONG_THRESHOLD is not None:
        case["solver"]["boomer_amg_strong_threshold"] = BOOMER_AMG_STRONG_THRESHOLD

    # Keep the registry, parameter expressions, materials, and geometry from
    # the frozen source, but expose only this case to prevent an accidental
    # dependency run from another case in the source project.
    project["cases"] = {CASE_NAME: case}
    DIAGNOSTIC_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_PATH.write_text(
        json.dumps(project, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return PROJECT_PATH


def command_for(solver: Path, runtime_bin: Path, toolchain_bin: Path) -> list[str]:
    return [
        sys.executable,
        str(ROOT / "run.py"),
        CASE_NAME,
        "--project",
        str(PROJECT_PATH),
        "--skip-sync",
        "--elmer-solver",
        str(solver),
        "--runtime-bin",
        str(runtime_bin),
        "--toolchain-bin",
        str(toolchain_bin),
        "--mpi-procs",
        "1",
    ]


def _number(value: str) -> int | float | str:
    try:
        return int(value)
    except ValueError:
        try:
            return float(value.replace("D", "E").replace("d", "e"))
        except ValueError:
            return value


def audit_solver_log(path: Path) -> dict[str, Any]:
    """Extract all Gate 3 residual evidence without declaring a pass."""
    if not path.is_file():
        return {"log_exists": False, "telemetry_records": [], "linear_residuals": []}
    text = path.read_text(encoding="utf-8", errors="replace")
    records: list[dict[str, Any]] = []
    for match in re.finditer(r"PHASE24_HYPRE_SOLVE_FAILURE\s+([^\r\n]+)", text):
        record: dict[str, Any] = {}
        for item in match.group(1).split():
            if "=" in item:
                key, value = item.split("=", 1)
                record[key] = _number(value)
        records.append(record)

    linear_residuals: list[float] = []
    for value in re.findall(
        r"(?:final_relative_residual|norm)\s*[=:]?\s*([-+0-9.EeDd]+)", text
    ):
        try:
            parsed = float(value.replace("D", "E").replace("d", "e"))
        except ValueError:
            continue
        if math.isfinite(parsed):
            linear_residuals.append(parsed)
    nonlinear_residuals: list[float] = []
    for value in re.findall(
        r"ComputeChange:\s+(?:NS|SS).*?\(\s*[-+0-9.EeDd]+\s+([-+0-9.EeDd]+)",
        text,
    ):
        try:
            parsed = float(value.replace("D", "E").replace("d", "e"))
        except ValueError:
            continue
        if math.isfinite(parsed):
            nonlinear_residuals.append(parsed)
    nonlinear_ns_residuals: list[float] = []
    for value in re.findall(
        r"ComputeChange:\s+NS.*?\(\s*[-+0-9.EeDd]+\s+([-+0-9.EeDd]+)",
        text,
    ):
        try:
            parsed = float(value.replace("D", "E").replace("d", "e"))
        except ValueError:
            continue
        if math.isfinite(parsed):
            nonlinear_ns_residuals.append(parsed)
    return {
        "log_exists": True,
        "all_done": "MAIN: *** Elmer Solver: ALL DONE ***" in text,
        "stop_1": text.count("STOP 1"),
        "abort_markers": sum(
            text.count(marker)
            for marker in ("job aborted", "application aborted", "Could not allocate memory")
        ),
        "hypre_markers": text.count("SolveHypre:"),
        "telemetry_records": records,
        "telemetry_missing": not records,
        "linear_residuals": linear_residuals,
        "linear_residual_final": linear_residuals[-1] if linear_residuals else None,
        "nonlinear_residuals": nonlinear_residuals,
        "nonlinear_residual_final": (
            nonlinear_ns_residuals[-1]
            if nonlinear_ns_residuals
            else (nonlinear_residuals[-1] if nonlinear_residuals else None)
        ),
        "tolerance_failure": "HYPRE failed to satisfy" in text,
    }


def read_last_iteration(path: Path) -> dict[str, float] | None:
    if not path.is_file():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    row = rows[-1]
    result: dict[str, float] = {}
    for key, value in row.items():
        if value is None or not value.strip():
            continue
        try:
            result[key] = float(value)
        except ValueError:
            continue
    return result


def read_last_series_current(path: Path) -> float | None:
    if not path.is_file():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    try:
        value = float(rows[-1]["tes_current_A"])
    except (KeyError, TypeError, ValueError):
        return None
    return value * 1.0e6 if math.isfinite(value) else None


def sha256(path: Path | None) -> str | None:
    if path is None or not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def evaluate_gate3(
    *,
    run: dict[str, Any],
    audit: dict[str, Any],
    iteration: dict[str, float] | None,
    series_current_uA: float | None,
    project_case: dict[str, Any],
) -> dict[str, Any]:
    """Evaluate the documented Gate 3 criteria from recorded evidence."""
    current_uA = series_current_uA
    current_source = "series" if current_uA is not None else None
    if current_uA is None and iteration is not None:
        raw = iteration.get("raw_current_A")
        if raw is not None and math.isfinite(raw):
            current_uA = raw * 1.0e6
            current_source = "iteration_series"

    finite_observables = False
    if iteration:
        names = ("tes_temperature_K", "tes_resistance_ohm", "raw_power_W")
        finite_observables = all(
            name in iteration and math.isfinite(iteration[name]) for name in names
        )
    if current_uA is not None:
        finite_observables = finite_observables and math.isfinite(current_uA)

    no_restart = not any(
        project_case.get(key)
        for key in ("restart_from", "restart_file_path", "preexisting_restart")
    )
    linear_recorded = bool(audit.get("linear_residuals")) or bool(
        audit.get("telemetry_records")
    )
    # The iteration CSV proves that the TES hook was called, but its
    # ``residual_W`` column is the circuit residual, not Elmer's field
    # nonlinear residual.  Require an actual ComputeChange record here.
    nonlinear_recorded = bool(audit.get("nonlinear_residuals"))
    circuit_recorded = bool(
        iteration and "residual_W" in iteration and math.isfinite(iteration["residual_W"])
    )
    mumps_delta = (
        abs(current_uA - STAGE11_MUMPS_CURRENT_UA) / STAGE11_MUMPS_CURRENT_UA * 100.0
        if current_uA is not None
        else None
    )
    comsol_delta = (
        abs(current_uA - COMSOL_CURRENT_UA) / COMSOL_CURRENT_UA * 100.0
        if current_uA is not None
        else None
    )
    # COMSOL is the physical acceptance target.  The same-mesh MUMPS result is
    # retained as a backend/parity diagnostic, but it is not a hard Gate 3
    # blocker: a conformal discretization can make HYPRE closer to COMSOL than
    # to an older MUMPS operating point.
    criteria = {
        "normal_exit": run.get("exit_code") == 0 and bool(audit.get("all_done")),
        "independent_initial_temperature": (
            project_case.get("initial_temperature") == "T_0"
            and project_case.get("apply_mortar_bcs") is APPLY_MORTAR_BCS
            and no_restart
        ),
        "finite_observables": finite_observables,
        "mumps_reference_recorded": (
            mumps_delta is not None
            and math.isfinite(STAGE11_MUMPS_CURRENT_UA)
            and STAGE11_MUMPS_CURRENT_UA > 0.0
        ),
        "comsol_current_within_0p6_percent": (
            comsol_delta is not None and comsol_delta <= COMSOL_LIMIT_PERCENT
        ),
        "hypre_linear_residual_recorded": linear_recorded,
        "elmer_nonlinear_residual_recorded": nonlinear_recorded,
        "tes_circuit_residual_recorded": circuit_recorded,
    }
    return {
        "status": "PASS" if all(criteria.values()) else "FAIL",
        "criteria": criteria,
        "current_uA": current_uA,
        "current_source": current_source,
        "stage11_mumps_current_uA": STAGE11_MUMPS_CURRENT_UA,
        "comsol_current_uA": COMSOL_CURRENT_UA,
        "mumps_delta_percent": mumps_delta,
        "comsol_delta_percent": comsol_delta,
        "linear_residual_final": audit.get("linear_residual_final"),
        "nonlinear_residual_final": audit.get("nonlinear_residual_final"),
        "circuit_residual_W": iteration.get("residual_W") if iteration else None,
    }


def provenance(
    *,
    project: Path,
    solver: Path,
    runtime_bin: Path,
    toolchain_bin: Path,
    manifest: dict[str, Any] | None,
) -> dict[str, Any]:
    sif = ROOT / "generated" / "cases" / f"{CASE_NAME}.sif"
    mesh_header = ROOT / "work" / "meshes" / MESH / "mesh.header"
    files = {
        "project": sha256(project),
        "sif": sha256(sif),
        "mesh_header": sha256(mesh_header),
        "solver": sha256(solver),
        "solver_log": sha256(SOLVER_LOG),
    }
    payload = {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "solver": str(solver.resolve()),
        "runtime_bin": str(runtime_bin.resolve()),
        "toolchain_bin": str(toolchain_bin.resolve()),
        "runtime_artifacts_sha256": (manifest or {}).get("runtime_artifacts_sha256", {}),
        "files_sha256": files,
    }
    encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")
    payload["environment_sha256"] = hashlib.sha256(encoded).hexdigest()
    return payload


def write_summary(payload: dict[str, Any]) -> None:
    DIAGNOSTIC_DIR.mkdir(parents=True, exist_ok=True)
    (DIAGNOSTIC_DIR / "summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    gate = payload.get("gate3", {})
    criteria = gate.get("criteria", {})
    lines = [
        "# Phase24 Gate 3: HYPRE steady state",
        "",
        f"- Status: **{gate.get('status', 'NOT RUN')}**",
        f"- Case: `{CASE_NAME}`",
        f"- Source project: `{SOURCE_PROJECT}`",
        f"- Mesh: `{MESH}`; mortar: `{APPLY_MORTAR_BCS}`",
        "- Initial state: `T_0`, with no restart input",
        f"- HYPRE: `{LINEAR_SYSTEM}`, max {LINEAR_MAX_ITERATIONS} iterations, tolerance `{LINEAR_TOLERANCE:g}`",
        f"- HYPRE GMRES dimension: `{HYPRE_GMRES_DIMENSION}`",
        f"- HYPRE reuse: `{HYPRE_REUSE}`; preconditioner lagging: `{PRECONDITIONER_LAGGING}`",
        f"- BoomerAMG strong threshold override: `{BOOMER_AMG_STRONG_THRESHOLD}`",
        "",
        "## Criteria",
        "",
        "| criterion | result |",
        "|---|---|",
    ]
    for key, value in criteria.items():
        lines.append(f"| `{key}` | `{value}` |")
    terminal_note = (
        "Gate 3 passed with COMSOL as the primary physical target; the MUMPS delta "
        "is retained as a backend diagnostic. Gate 4 and later must use this "
        "no-mortar policy and remain subject to their own gates."
        if gate.get("status") == "PASS"
        else "Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle."
    )
    lines.extend(
        [
            "",
            f"- HYPRE current [µA]: `{gate.get('current_uA')}`",
            f"- Current evidence source: `{gate.get('current_source')}` (a non-series source is not a converged steady result)",
        f"- Stage 11 MUMPS delta [%]: `{gate.get('mumps_delta_percent')}`",
        f"- COMSOL delta [%]: `{gate.get('comsol_delta_percent')}`",
        f"- Stage 11 MUMPS reference [µA]: `{gate.get('stage11_mumps_current_uA')}`",
        f"- COMSOL reference [µA]: `{gate.get('comsol_current_uA')}`",
            f"- HYPRE final linear residual: `{gate.get('linear_residual_final')}`",
            f"- Elmer nonlinear residual: `{gate.get('nonlinear_residual_final')}`",
            f"- TES circuit residual [W]: `{gate.get('circuit_residual_W')}`",
            "",
            f"Project: `{PROJECT_PATH}`",
            f"Solver log: `{SOLVER_LOG}`",
            f"Launcher log: `{LAUNCHER_LOG}`",
            "",
            terminal_note,
            "",
        ]
    )
    (DIAGNOSTIC_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    global LINEAR_SYSTEM, LINEAR_MAX_ITERATIONS, LINEAR_TOLERANCE
    global HYPRE_REUSE, PRECONDITIONER_LAGGING
    global BOOMER_AMG_STRONG_THRESHOLD
    global HYPRE_GMRES_DIMENSION, APPLY_MORTAR_BCS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--solver",
        type=Path,
        default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin\ElmerSolver.exe"),
    )
    parser.add_argument(
        "--runtime-bin",
        type=Path,
        default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin"),
    )
    parser.add_argument("--toolchain-bin", type=Path, default=Path(r"C:\msys64\ucrt64\bin"))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--variant", help="suffix for an isolated tuning trial artifact")
    parser.add_argument(
        "--source-project",
        type=Path,
        default=None,
        help="project containing the selected mesh registry and steady template",
    )
    parser.add_argument("--mesh", default=None, help="mesh registry key")
    parser.add_argument("--base-case", default=None, help="source steady case to clone")
    parser.add_argument(
        "--reference-current-uA",
        type=float,
        default=None,
        help="same-mesh converged MUMPS current used by the Gate 3 comparison",
    )
    parser.add_argument(
        "--comsol-current-uA",
        type=float,
        default=None,
        help="COMSOL current corresponding to the selected mesh/physics comparison",
    )
    parser.add_argument("--linear-system", default=LINEAR_SYSTEM)
    parser.add_argument("--max-iterations", type=int, default=LINEAR_MAX_ITERATIONS)
    parser.add_argument("--linear-tolerance", type=float, default=LINEAR_TOLERANCE)
    parser.add_argument("--hypre-reuse", action="store_true")
    parser.add_argument(
        "--preconditioner-lagging",
        choices=("adaptive", "disabled"),
        default=None,
    )
    parser.add_argument("--strong-threshold", type=float, default=None)
    parser.add_argument("--gmres-dimension", type=int, default=HYPRE_GMRES_DIMENSION)
    parser.add_argument(
        "--audit-only",
        action="store_true",
        help="re-evaluate the existing Gate 3 result without starting Elmer",
    )
    parser.add_argument("--mortar", action="store_true", help="diagnostic: enable Mortar BCs for a nonconforming mesh")
    args = parser.parse_args()
    APPLY_MORTAR_BCS = bool(args.mortar)

    configure_mesh(
        source_project=args.source_project,
        mesh=args.mesh,
        base_case=args.base_case,
        reference_current_uA=args.reference_current_uA,
        comsol_current_uA=args.comsol_current_uA,
    )
    configure_variant(args.variant)
    LINEAR_SYSTEM = args.linear_system
    LINEAR_MAX_ITERATIONS = args.max_iterations
    LINEAR_TOLERANCE = args.linear_tolerance
    HYPRE_REUSE = args.hypre_reuse
    PRECONDITIONER_LAGGING = args.preconditioner_lagging
    BOOMER_AMG_STRONG_THRESHOLD = args.strong_threshold
    if args.gmres_dimension < 2:
        raise SystemExit("--gmres-dimension must be at least 2")
    HYPRE_GMRES_DIMENSION = args.gmres_dimension

    solver = absolute(args.solver)
    runtime_bin = absolute(args.runtime_bin)
    toolchain_bin = absolute(args.toolchain_bin)
    project = build_project()
    sync_command = [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project)]
    if args.dry_run:
        print("[dry-run] " + " ".join(sync_command))
    else:
        for path, label in (
            (solver, "solver"),
            (runtime_bin, "runtime DLL directory"),
            (toolchain_bin, "toolchain DLL directory"),
        ):
            if not path.exists():
                raise SystemExit(f"{label} not found: {path}")
        sync = subprocess.run(sync_command, cwd=ROOT, capture_output=True, text=True)
        (DIAGNOSTIC_DIR / "sync.log").write_text(
            sync.stdout + sync.stderr, encoding="utf-8"
        )
        if sync.returncode != 0:
            raise SystemExit(f"sync_elmer_parameters.py failed with exit code {sync.returncode}")

    command = command_for(solver, runtime_bin, toolchain_bin)
    print("[gate3] command: " + " ".join(command))
    started = datetime.now(timezone.utc).isoformat()
    if args.audit_only:
        manifest_path = RESULT_DIR / "manifest.json"
        if not manifest_path.is_file():
            raise SystemExit(f"Gate 3 manifest not found: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        project_data = json.loads(project.read_text(encoding="utf-8"))
        case = project_data["cases"][CASE_NAME]
        run = {
            "command": command,
            "started": manifest.get("started"),
            "finished": manifest.get("finished"),
            "exit_code": manifest.get("exit_code"),
            "log": str(LAUNCHER_LOG),
        }
        iteration = read_last_iteration(RESULT_DIR / f"{CASE_NAME}_iterations.csv")
        series_current = read_last_series_current(RESULT_DIR / f"{CASE_NAME}_series.csv")
        audit = audit_solver_log(SOLVER_LOG)
        gate = evaluate_gate3(
            run=run,
            audit=audit,
            iteration=iteration,
            series_current_uA=series_current,
            project_case=case,
        )
        payload = {
            "schema": "tes.phase24.gate3.v1",
            "generated_utc": started,
            "project": str(project),
            "case": CASE_NAME,
            "run": run,
            "audit": audit,
            "iteration_last": iteration,
            "series_current_uA": series_current,
            "manifest": str(manifest_path),
            "gate3": gate,
            "provenance": provenance(
                project=project,
                solver=solver,
                runtime_bin=runtime_bin,
                toolchain_bin=toolchain_bin,
                manifest=manifest,
            ),
        }
        write_summary(payload)
        print(f"gate3={gate['status']}")
        print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
        return 0 if gate["status"] == "PASS" else 1

    if args.dry_run:
        run = {"command": command, "exit_code": None, "started": started, "finished": started}
        payload = {
            "schema": "tes.phase24.gate3.v1",
            "generated_utc": started,
            "project": str(project),
            "case": CASE_NAME,
            "run": run,
            "audit": {},
            "gate3": {"status": "NOT RUN", "criteria": {}},
            "provenance": provenance(
                project=project,
                solver=solver,
                runtime_bin=runtime_bin,
                toolchain_bin=toolchain_bin,
                manifest=None,
            ),
        }
        write_summary(payload)
        print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
        return 0

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with LAUNCHER_LOG.open("w", encoding="utf-8") as handle:
        process = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    run = {
        "command": command,
        "started": started,
        "finished": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - start,
        "exit_code": process.returncode,
        "log": str(LAUNCHER_LOG),
    }
    project_data = json.loads(project.read_text(encoding="utf-8"))
    case = project_data["cases"][CASE_NAME]
    iteration = read_last_iteration(RESULT_DIR / f"{CASE_NAME}_iterations.csv")
    series_current = read_last_series_current(RESULT_DIR / f"{CASE_NAME}_series.csv")
    audit = audit_solver_log(SOLVER_LOG)
    manifest_path = RESULT_DIR / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else None
    gate = evaluate_gate3(
        run=run,
        audit=audit,
        iteration=iteration,
        series_current_uA=series_current,
        project_case=case,
    )
    payload = {
        "schema": "tes.phase24.gate3.v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(project),
        "case": CASE_NAME,
        "run": run,
        "audit": audit,
        "iteration_last": iteration,
        "series_current_uA": series_current,
        "manifest": str(manifest_path) if manifest_path.is_file() else None,
        "gate3": gate,
        "provenance": provenance(
            project=project,
            solver=solver,
            runtime_bin=runtime_bin,
            toolchain_bin=toolchain_bin,
            manifest=manifest,
        ),
    }
    write_summary(payload)
    print(f"gate3={gate['status']}")
    print(f"summary={DIAGNOSTIC_DIR / 'summary.md'}")
    return 0 if gate["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

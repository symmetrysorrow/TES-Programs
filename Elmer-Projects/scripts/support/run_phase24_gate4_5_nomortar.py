"""Run the Phase24 transient Gate 4/5 campaign.

The runner uses a selected mesh, a converged MUMPS restart, and the same HYPRE
policy that passed Gate 3.  It can run a short 0.9-us diagnostic, a targeted
40-us diagnosis, the full 100-us Gate 4 window, or a 1-ms Gate 5 stability
window.  Each run is isolated under
``artifacts/phase24_gate4_5_nomortar``.

For a mortar restart produced by an iterative Gate 3 solve, ``--refine-reference``
adds one full constrained MUMPS steady solve before the transient cases.  This
projects the approximate Gate 3 field back onto the exact primal equilibrium
without changing the physical model or the transient SIF formulation.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE_PROJECT = ROOT / "elmer_project_singlepixel_prod_v2_original_timegrid.json"
TIMEGRID_PROJECT = ROOT / "elmer_project_comsol_timegrid.json"
COMSOL = ROOT / "docs" / "Single-Pixel.txt"
MESH = "mesh_singlepixel_prod_v2"
MESH_TAG = "prod_v2"
REFERENCE_CASE = "case_tes_steady_prod_v2_nomortar"
REFERENCE_RESULT = ROOT / "work" / "meshes" / MESH / f"{REFERENCE_CASE}.result"
REFERENCE_STATE = ROOT / "work" / "meshes" / MESH / f"{REFERENCE_CASE}.state"
COMSOL_CURRENT_UA = 143.055049
MUMPS_CURRENT_UA = 143.53734493231093
MUMPS_LIMIT_PERCENT = 0.05
WAVEFORM_MAX_LIMIT_UA = 0.05
WAVEFORM_RMSE_LIMIT_UA = 0.03
TIME_CONSTANT_LIMIT_US = 1.0
HYPRE_SYSTEM = "iterative_hypre_flexgmres_boomeramg"
HYPRE_MAX_ITERATIONS = 4000
HYPRE_TOLERANCE = 2.0e-8
HYPRE_STRONG_THRESHOLD = 0.5
SHORT_STAGES = [
    ["18[us]", 1],
    ["1[us]", 2],
    ["1[ns]", 1],
    ["10[ns]", 10],
    ["100[ns]", 9],
]
PULSE_EVENT_S = 0.020020
OUT_ROOT = ROOT / "artifacts" / "phase24_gate4_5_nomortar"
APPLY_MORTAR_BCS = False
RUN_TAG = ""
NONLINEAR_TOLERANCE = 1.0e-8
NONLINEAR_MAX_ITERATIONS = 120


def configure_mesh(
    *,
    base_project: Path | None,
    mesh: str | None,
    mesh_tag: str | None,
    reference_case: str | None,
    reference_current_uA: float | None,
) -> None:
    global BASE_PROJECT, MESH, MESH_TAG, REFERENCE_CASE
    global REFERENCE_RESULT, REFERENCE_STATE, MUMPS_CURRENT_UA
    if base_project is not None:
        BASE_PROJECT = base_project if base_project.is_absolute() else ROOT / base_project
    if mesh is not None:
        MESH = mesh
    if mesh_tag is not None:
        MESH_TAG = mesh_tag
    elif mesh is not None:
        MESH_TAG = {
            "mesh_singlepixel_prod_v2": "prod_v2",
            "mesh_singlepixel_conformal_gpu": "conformal_coarse",
            "mesh_singlepixel_conformal_gpu_fine": "conformal_fine",
            "mesh_singlepixel_conformal_gpu_refine20": "conformal_refine20",
            "mesh_singlepixel_conformal_gpu_fine_stycast32": "conformal_fine_stycast32",
            "mesh_singlepixel_conformal_gpu_fine_stycast32_nomortar": "conformal_fine_stycast32_nomortar",
            "mesh_singlepixel_gpu_fine_stycast32_mortar": "fine_stycast32_mortar",
        }.get(mesh, "mesh")
    if reference_case is not None:
        REFERENCE_CASE = reference_case
    if reference_current_uA is not None:
        MUMPS_CURRENT_UA = reference_current_uA
    REFERENCE_RESULT = ROOT / "work" / "meshes" / MESH / f"{REFERENCE_CASE}.result"
    REFERENCE_STATE = ROOT / "work" / "meshes" / MESH / f"{REFERENCE_CASE}.state"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def seconds(token: str) -> float:
    match = re.fullmatch(r"([0-9.eE+-]+)\[([a-z]+)\]", token)
    if not match:
        raise ValueError(f"unsupported timestep token: {token}")
    return float(match.group(1)) * {
        "s": 1.0,
        "ms": 1.0e-3,
        "us": 1.0e-6,
        "ns": 1.0e-9,
    }[match.group(2)]


def hybrid_schedule() -> list[list[object]]:
    sys.path.insert(0, str(ROOT))
    from scripts.prep.run_singlepixel_prod_v2_original_timegrid import hybrid_timesteps

    source = json.loads(TIMEGRID_PROJECT.read_text(encoding="utf-8"))
    original = source["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]["timesteps"]
    return hybrid_timesteps(original)


def trim_schedule(schedule: list[list[object]], window: str) -> list[list[object]]:
    if window == "short":
        return copy.deepcopy(SHORT_STAGES)
    target_us = {"40us": 40.0, "100us": 100.0, "1ms": 1000.0}[window]
    target_s = PULSE_EVENT_S + target_us * 1.0e-6
    elapsed = 0.0
    result: list[list[object]] = []
    for token_obj, count_obj in schedule:
        token = str(token_obj)
        count = int(count_obj)
        dt = seconds(token)
        available = target_s - (0.020000 + elapsed)
        if available <= 1.0e-15:
            break
        ratio = available / dt
        nearest = round(ratio)
        # Avoid dropping a final timestep to binary floating-point roundoff
        # when the requested endpoint lies exactly on the timestep grid.
        if abs(ratio - nearest) < 1.0e-3:
            take = min(count, max(0, int(nearest)))
        else:
            take = min(count, max(0, int(math.floor(ratio))))
        if take:
            result.append([token, take])
            elapsed += dt * take
        if take < count:
            break
    if not result or elapsed < target_us * 1.0e-6 - 1.0e-9:
        raise ValueError(f"time grid does not cover {window}: elapsed={elapsed}")
    return result


def case_name(window: str, backend: str) -> str:
    suffix = "_mortar" if APPLY_MORTAR_BCS else ""
    tag = f"_{RUN_TAG}" if RUN_TAG else ""
    return f"case_phase24_g45_{MESH_TAG}_{window}_{backend}{suffix}{tag}"


def refinement_case_name() -> str:
    suffix = "_mortar" if APPLY_MORTAR_BCS else ""
    tag = f"_{RUN_TAG}" if RUN_TAG else ""
    return f"case_phase24_restart_refine_{MESH_TAG}_mumps{suffix}{tag}"


def refinement_spec(
    template: dict[str, Any],
    *,
    name: str,
    source_case: str,
) -> dict[str, Any]:
    """Build the direct steady projection case used before a transient.

    The source result/state remain read-only inputs.  The refinement has its
    own state/result names so a subsequent transient cannot overwrite the
    input checkpoint or accidentally feed a partially completed run back into
    the chain.
    """
    candidate = copy.deepcopy(template)
    candidate.update(
        {
            "template": "steady",
            "mesh": MESH,
            "heat_source": "circuit_inner",
            "apply_mortar_bcs": APPLY_MORTAR_BCS,
            "initial_temperature": "T_0",
            "restart_from": None,
            "restart_file_base": source_case,
            "restart_file_path": f"../work/meshes/{MESH}/{source_case}.result",
            "preexisting_restart": True,
            "restart_time": 0.020,
            "steady_state_max_iterations": 1,
            "output_intervals": 1,
            "series_file": f"{name}_series.csv",
            "iteration_series_file": f"{name}_iterations.csv",
            "state_file": f"work/meshes/{MESH}/{name}.state",
            "output_result": True,
            "output_file_path": f"../work/meshes/{MESH}/{name}.result",
            "post_file": False,
            "vtu": False,
            "inner_circuit_step_commit": True,
            "solver_comment": (
                "Phase24 restart refinement: full constrained steady MUMPS "
                f"projection of {source_case}"
            ),
            "phase24_smoke": {
                "purpose": "exact constrained restart projection before transient",
                "source_restart": source_case,
                "backend": "mumps",
                "mortar": APPLY_MORTAR_BCS,
            },
        }
    )
    candidate.pop("pulse", None)
    candidate.pop("timesteps", None)
    candidate["solver"] = {
        **candidate.get("solver", {}),
        "linear_system": "mumps",
        "nonlinear_max_iterations": NONLINEAR_MAX_ITERATIONS,
        "nonlinear_convergence_tolerance": NONLINEAR_TOLERANCE,
        "nonlinear_relaxation_factor": 1.0,
        "steady_state_convergence_tolerance": 1.0e-9,
    }
    return candidate


def prepare_refinement_project() -> tuple[Path, str]:
    """Create the isolated direct-refinement project and seed its state."""
    if not REFERENCE_RESULT.is_file():
        raise FileNotFoundError(
            f"restart refinement input result not found: {REFERENCE_RESULT}"
        )
    if not REFERENCE_STATE.is_file():
        raise FileNotFoundError(
            f"restart refinement input state not found: {REFERENCE_STATE}"
        )

    project = json.loads(BASE_PROJECT.read_text(encoding="utf-8"))
    template_project = json.loads(TIMEGRID_PROJECT.read_text(encoding="utf-8"))
    template = copy.deepcopy(
        template_project["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]
    )
    name = refinement_case_name()
    out = OUT_ROOT / "restart_refinement" / ("mortar" if APPLY_MORTAR_BCS else "nomortar")
    out.mkdir(parents=True, exist_ok=True)
    state_path = ROOT / "work" / "meshes" / MESH / f"{name}.state"
    state_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(REFERENCE_STATE, state_path)

    project["cases"] = {
        name: refinement_spec(template, name=name, source_case=REFERENCE_CASE)
    }
    project_path = out / "project.json"
    project_path.write_text(
        json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return project_path, name


def prepare_project(window: str, backend: str) -> tuple[Path, str]:
    project = json.loads(BASE_PROJECT.read_text(encoding="utf-8"))
    template_project = json.loads(TIMEGRID_PROJECT.read_text(encoding="utf-8"))
    template = copy.deepcopy(
        template_project["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]
    )
    name = case_name(window, backend)
    out = OUT_ROOT / window / ("mortar" if APPLY_MORTAR_BCS else backend)
    out.mkdir(parents=True, exist_ok=True)
    state_path = ROOT / "work" / "meshes" / MESH / f"{name}.state"
    shutil.copy2(REFERENCE_STATE, state_path)

    candidate = {
        **template,
        "mesh": MESH,
        "apply_mortar_bcs": APPLY_MORTAR_BCS,
        "initial_temperature": "T_0",
        "restart_from": None,
        "restart_file_base": REFERENCE_CASE,
        "restart_file_path": f"../work/meshes/{MESH}/{REFERENCE_CASE}.result",
        "preexisting_restart": True,
        "restart_time": 0.020,
        "timesteps": trim_schedule(hybrid_schedule(), window),
        "output_intervals": [],
        "series_file": f"{name}_series.csv",
        "iteration_series_file": f"{name}_iterations.csv",
        "state_file": f"work/meshes/{MESH}/{name}.state",
        "output_result": True,
        "output_file_path": f"../work/meshes/{MESH}/{name}.result",
        "vtu": False,
        "inner_circuit_step_commit": True,
        "solver_comment": (
            f"Phase24 Gate {4 if window != '1ms' else 5}: "
            f"{'mortar' if APPLY_MORTAR_BCS else 'no-mortar'} transient; "
            f"backend={backend}; window={window}"
        ),
        "comparison_time_grid": {
            "mode": (
                f"Phase24 {'mortar' if APPLY_MORTAR_BCS else 'no-mortar'} "
                f"Gate {4 if window != '1ms' else 5}"
            ),
            "event_s": PULSE_EVENT_S,
            "window": window,
            "mesh": MESH,
            "mortar": APPLY_MORTAR_BCS,
        },
        "phase24_smoke": {
            "purpose": "transient parity and stability",
            "independent_initial_temperature": False,
            "restart_reference": REFERENCE_CASE,
            "no_mortar": not APPLY_MORTAR_BCS,
            "memory_trace": True,
        },
    }
    candidate["output_intervals"] = [1] * len(candidate["timesteps"])
    solver = {
        **candidate.get("solver", {}),
        "linear_system": "mumps" if backend == "mumps" else HYPRE_SYSTEM,
        "nonlinear_max_iterations": NONLINEAR_MAX_ITERATIONS,
        "nonlinear_convergence_tolerance": NONLINEAR_TOLERANCE,
        "nonlinear_relaxation_factor": 1.0,
    }
    if backend == "hypre":
        solver.update(
            {
                "linear_system_max_iterations": HYPRE_MAX_ITERATIONS,
                "linear_system_convergence_tolerance": HYPRE_TOLERANCE,
                "hypre_gmres_dimension": 100,
                "boomer_amg_strong_threshold": HYPRE_STRONG_THRESHOLD,
            }
        )
    candidate["solver"] = solver
    project["cases"] = {name: candidate}
    project_path = out / "project.json"
    project_path.write_text(
        json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return project_path, name


def run_case(project_path: Path, name: str, solver: Path, runtime_bin: Path, toolchain_bin: Path) -> dict[str, Any]:
    out = project_path.parent
    sync = subprocess.run(
        [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project_path)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    (out / "sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
    command = [
        sys.executable,
        str(ROOT / "run.py"),
        name,
        "--project",
        str(project_path),
        "--mpi-procs",
        "1",
        "--elmer-solver",
        str(solver),
        "--runtime-bin",
        str(runtime_bin),
        "--toolchain-bin",
        str(toolchain_bin),
    ]
    log = out / "launcher.log"
    env = os.environ.copy()
    env["PHASE24_MEMORY_TRACE"] = "1"
    started = datetime.now(timezone.utc).isoformat()
    start = time.monotonic()
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, env=env)
    return {
        "command": command,
        "started": started,
        "finished": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - start,
        "exit_code": result.returncode,
        "launcher_log": str(log),
    }


def load_series(path: Path) -> tuple[list[float], list[float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    pairs = sorted(
        (float(row["time_s"]), float(row["tes_current_A"]) * 1.0e6)
        for row in rows
    )
    return [t for t, _ in pairs], [i for _, i in pairs]


def memory_trace(log_path: Path) -> dict[str, Any]:
    text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
    records: list[dict[str, int]] = []
    pattern = re.compile(
        r"PHASE24_MEMORY\s+[^\r\n]*?private_bytes=(\d+)\s+working_set=(\d+)\s+peak_working_set=(\d+)"
    )
    for match in pattern.finditer(text):
        records.append(
            {
                "private_bytes": int(match.group(1)),
                "working_set": int(match.group(2)),
                "peak_working_set": int(match.group(3)),
            }
        )
    private = [r["private_bytes"] for r in records]
    working = [r["working_set"] for r in records]
    return {
        "records": len(records),
        "private_bytes_initial": private[0] if private else None,
        "private_bytes_final": private[-1] if private else None,
        "private_bytes_peak": max(private) if private else None,
        "working_set_initial": working[0] if working else None,
        "working_set_final": working[-1] if working else None,
        "working_set_peak": max(working) if working else None,
        "private_bytes_monotonic_non_decreasing": bool(private) and all(a <= b for a, b in zip(private, private[1:])),
        "raw_records": records,
    }


def audit_solver(log_path: Path, manifest_path: Path) -> dict[str, Any]:
    text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
    residuals = [
        float(value.replace("D", "E").replace("d", "e"))
        for value in re.findall(r"SolveHypre: Required iterations \d+ \(method \d+\) to norm ([0-9.EeDd+-]+)", text)
    ]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    return {
        "all_done": "MAIN: *** Elmer Solver: ALL DONE ***" in text,
        "stop_1": text.count("STOP 1"),
        "mpi_abort": "MPI_ABORT" in text,
        "native_crash": any(marker in text for marker in ("Access violation", "Could not allocate memory", "Fortran runtime error")),
        "linear_residuals": residuals,
        "linear_residual_final": residuals[-1] if residuals else None,
        "manifest_exit_code": manifest.get("exit_code"),
        "solver_completed": manifest.get("solver_completed"),
    }


def compare_waveforms(hypre_series: Path, mumps_series: Path, out: Path, end_us: float) -> dict[str, Any]:
    sys.path.insert(0, str(ROOT))
    from scripts.analysis import compare_singlepixel_amgx_comsol as compare

    comsol = compare.read_comsol(COMSOL)
    hypre = compare.read_elmer(hypre_series)
    mumps = compare.read_elmer(mumps_series)
    hypre_metrics, hypre_aligned = compare.compare(comsol, hypre, end_us)
    mumps_metrics, mumps_aligned = compare.compare(comsol, mumps, end_us)
    compare.write_outputs(out / "hypre_vs_comsol", comsol, hypre, hypre_metrics, hypre_aligned, "HYPRE")
    compare.write_outputs(
        out / "mumps_vs_comsol",
        comsol,
        mumps,
        mumps_metrics,
        mumps_aligned,
        f"{'mortar' if APPLY_MORTAR_BCS else 'no-mortar'} MUMPS",
    )
    return {
        "hypre_vs_comsol": hypre_metrics,
        "mumps_vs_comsol": mumps_metrics,
        "hypre_baseline_minus_mumps_uA": hypre.baseline_uA - mumps.baseline_uA,
        "hypre_baseline_delta_percent_mumps": 100.0 * abs(hypre.baseline_uA - mumps.baseline_uA) / MUMPS_CURRENT_UA,
    }


def evaluate_gate4(metrics: dict[str, Any], audits: dict[str, Any], window: str = "100us") -> dict[str, Any]:
    hypre = metrics["hypre_vs_comsol"]
    crossings = hypre["crossing_us_at_comsol_peak_fraction"]
    t10 = crossings["0.1"].get("AMGX_minus_COMSOL")
    t50 = crossings["0.5"].get("AMGX_minus_COMSOL")
    time_ok = (
        True
        if window == "short"
        else all(value is not None and abs(value) <= TIME_CONSTANT_LIMIT_US for value in (t10, t50))
    )
    criteria = {
        "hypre_normal_exit": audits["hypre"]["manifest_exit_code"] == 0 and audits["hypre"]["all_done"],
        "mumps_normal_exit": audits["mumps"]["manifest_exit_code"] == 0 and audits["mumps"]["all_done"],
        "mumps_reference_recorded": math.isfinite(metrics["hypre_baseline_delta_percent_mumps"]),
        "comsol_waveform_max_within_0p05_uA": hypre["max_abs_difference_uA"] <= WAVEFORM_MAX_LIMIT_UA,
        "comsol_waveform_rmse_within_0p03_uA": hypre["rmse_uA"] <= WAVEFORM_RMSE_LIMIT_UA,
        "time_constants_within_1us_or_not_applicable": time_ok,
    }
    return {"status": "PASS" if all(criteria.values()) else "FAIL", "criteria": criteria, "metrics": metrics}


def evaluate_gate5(records: dict[str, Any]) -> dict[str, Any]:
    """Evaluate long-window stability independently for each executed backend."""
    backend_results: dict[str, Any] = {}
    for backend, record in records.items():
        audit = record["audit"]
        memory = record["memory"]
        initial = memory.get("working_set_initial")
        peak = memory.get("working_set_peak")
        criteria = {
            "normal_exit": audit.get("manifest_exit_code") == 0 and audit.get("all_done"),
            "no_stop_or_abort": not audit.get("stop_1") and not audit.get("mpi_abort") and not audit.get("native_crash"),
            "memory_trace_recorded": memory.get("records", 0) > 1,
            "private_usage_not_monotonic_increase": not memory.get("private_bytes_monotonic_non_decreasing", True),
            "working_set_peak_within_1p25x_initial": (
                initial is not None and peak is not None and peak <= 1.25 * initial
            ),
        }
        backend_results[backend] = {
            "status": "PASS" if all(criteria.values()) else "FAIL",
            "criteria": criteria,
            "memory": memory,
            "audit": audit,
        }
    status = "PASS" if backend_results and all(item["status"] == "PASS" for item in backend_results.values()) else "FAIL"
    return {"status": status, "backends": backend_results}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window", choices=("short", "40us", "100us", "1ms"), default="short")
    parser.add_argument("--backend", choices=("mumps", "hypre", "both"), default="both")
    parser.add_argument("--solver", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin\ElmerSolver.exe"))
    parser.add_argument("--runtime-bin", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin"))
    parser.add_argument("--toolchain-bin", type=Path, default=Path(r"C:\msys64\ucrt64\bin"))
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument(
        "--mortar",
        action="store_true",
        help="diagnostic only: enable Elmer Mortar BCs on the same mesh",
    )
    parser.add_argument("--linear-tolerance", type=float, default=None, help="override HYPRE linear convergence tolerance")
    parser.add_argument("--max-iterations", type=int, default=None, help="override HYPRE linear iteration limit")
    parser.add_argument("--tag", default="", help="suffix for an isolated transient diagnostic")
    parser.add_argument("--nonlinear-tolerance", type=float, default=None, help="override transient nonlinear convergence tolerance")
    parser.add_argument("--nonlinear-max-iterations", type=int, default=None, help="override transient nonlinear iteration limit")
    parser.add_argument("--base-project", type=Path, default=None, help="project containing the selected mesh registry")
    parser.add_argument("--mesh", default=None, help="mesh registry key")
    parser.add_argument("--mesh-tag", default=None, help="short case/artifact tag for the mesh")
    parser.add_argument("--reference-case", default=None, help="converged steady restart case")
    parser.add_argument("--reference-current-uA", type=float, default=None, help="same-mesh MUMPS current for diagnostics")
    parser.add_argument(
        "--refine-reference",
        action="store_true",
        help=(
            "run one full constrained MUMPS steady solve from the selected "
            "restart result/state, then use that refined checkpoint for all "
            "transient cases"
        ),
    )
    args = parser.parse_args()
    global APPLY_MORTAR_BCS, HYPRE_TOLERANCE, HYPRE_MAX_ITERATIONS, RUN_TAG
    global NONLINEAR_TOLERANCE, NONLINEAR_MAX_ITERATIONS
    global REFERENCE_CASE, REFERENCE_RESULT, REFERENCE_STATE
    APPLY_MORTAR_BCS = bool(args.mortar)
    if args.linear_tolerance is not None:
        HYPRE_TOLERANCE = args.linear_tolerance
    if args.max_iterations is not None:
        HYPRE_MAX_ITERATIONS = args.max_iterations
    RUN_TAG = re.sub(r"[^A-Za-z0-9_.-]+", "_", args.tag).strip("_")
    if args.nonlinear_tolerance is not None:
        NONLINEAR_TOLERANCE = args.nonlinear_tolerance
    if args.nonlinear_max_iterations is not None:
        NONLINEAR_MAX_ITERATIONS = args.nonlinear_max_iterations
    configure_mesh(
        base_project=args.base_project,
        mesh=args.mesh,
        mesh_tag=args.mesh_tag,
        reference_case=args.reference_case,
        reference_current_uA=args.reference_current_uA,
    )

    source_reference_case = REFERENCE_CASE
    refinement_record: dict[str, Any] | None = None
    refinement_paths: dict[str, str] | None = None
    if args.refine_reference:
        refinement_project, refined_case = prepare_refinement_project()
        refinement_paths = {
            "project": str(refinement_project),
            "name": refined_case,
            "result": str(
                ROOT / "work" / "meshes" / MESH / f"{refined_case}.result"
            ),
            "state": str(
                ROOT / "work" / "meshes" / MESH / f"{refined_case}.state"
            ),
            "log": str(ROOT / "results" / refined_case / "solver.log"),
            "manifest": str(ROOT / "results" / refined_case / "manifest.json"),
        }
        if not args.audit_only:
            refinement_record = run_case(
                refinement_project,
                refined_case,
                args.solver,
                args.runtime_bin,
                args.toolchain_bin,
            )
        else:
            refinement_record = {}
        refinement_record["audit"] = audit_solver(
            ROOT / "results" / refined_case / "solver.log",
            ROOT / "results" / refined_case / "manifest.json",
        )
        if refinement_record["audit"].get("manifest_exit_code") != 0 or not refinement_record["audit"].get("solver_completed"):
            print(
                f"restart refinement failed; transient cases were not started: "
                f"{refined_case}"
            )
            return 1
        REFERENCE_CASE = refined_case
        REFERENCE_RESULT = ROOT / "work" / "meshes" / MESH / f"{refined_case}.result"
        REFERENCE_STATE = ROOT / "work" / "meshes" / MESH / f"{refined_case}.state"

    backends = ("mumps", "hypre") if args.backend == "both" else (args.backend,)
    records: dict[str, Any] = {}
    paths: dict[str, dict[str, Path | str]] = {}
    for backend in backends:
        project_path, name = prepare_project(args.window, backend)
        out = project_path.parent
        paths[backend] = {
            "project": project_path,
            "name": name,
            "series": ROOT / "results" / name / f"{name}_series.csv",
            "log": ROOT / "results" / name / "solver.log",
            "manifest": ROOT / "results" / name / "manifest.json",
        }
        if not args.audit_only:
            records[backend] = run_case(project_path, name, args.solver, args.runtime_bin, args.toolchain_bin)
        records.setdefault(backend, {})
        records[backend]["audit"] = audit_solver(paths[backend]["log"], paths[backend]["manifest"])
        records[backend]["memory"] = memory_trace(paths[backend]["log"])
    out = OUT_ROOT / args.window / ("mortar" if APPLY_MORTAR_BCS else "nomortar")
    out.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "schema": "tes.phase24.gate4_5.nomortar.v1",
        "window": args.window,
        "policy": {
            "mesh": MESH,
            "mortar": APPLY_MORTAR_BCS,
            "initial_reference_case": source_reference_case,
            "reference_case": REFERENCE_CASE,
            "restart_refinement": bool(args.refine_reference),
            "hypre_system": HYPRE_SYSTEM,
            "hypre_max_iterations": HYPRE_MAX_ITERATIONS,
            "hypre_tolerance": HYPRE_TOLERANCE,
            "boomer_amg_strong_threshold": HYPRE_STRONG_THRESHOLD,
        },
        "paths": {backend: {key: str(value) for key, value in data.items()} for backend, data in paths.items()},
        "runs": records,
    }
    if refinement_record is not None:
        payload["restart_refinement"] = refinement_record
        payload["restart_refinement_paths"] = refinement_paths
    if {"mumps", "hypre"}.issubset(backends):
        end_us = {"short": 0.9, "40us": 40.0, "100us": 100.0, "1ms": 1000.0}[args.window]
        series_ready = all(Path(paths[backend]["series"]).is_file() for backend in ("mumps", "hypre"))
        if series_ready:
            metrics = compare_waveforms(paths["hypre"]["series"], paths["mumps"]["series"], out / "comparison", end_us)
            payload["gate4"] = evaluate_gate4(
                metrics,
                {k: records[k]["audit"] for k in ("mumps", "hypre")},
                args.window,
            )
        else:
            payload["gate4"] = {
                "status": "FAIL",
                "criteria": {"waveform_series_available": False},
                "reason": "mumps or hypre waveform series is missing; inspect per-backend manifest/log",
            }
    if args.window in ("100us", "1ms"):
        payload["gate5"] = evaluate_gate5(records)
    (out / "summary.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = [
        f"# Phase24 Gate 4/5 conformal {'mortar' if APPLY_MORTAR_BCS else 'no-mortar'}",
        "",
        f"- Window: `{args.window}`",
        f"- Mesh: `{MESH}`; mortar: `{APPLY_MORTAR_BCS}`",
        f"- Backend(s): `{', '.join(backends)}`",
        f"- Restart refinement: `{bool(args.refine_reference)}`",
        f"- Initial reference: `{source_reference_case}`",
        f"- Transient reference: `{REFERENCE_CASE}`",
    ]
    if "gate4" in payload:
        lines.extend([f"- Gate4: **{payload['gate4'].get('status', 'NOT RUN')}**"])
    if "gate5" in payload:
        gate5 = payload["gate5"]
        lines.extend([f"- Gate5: **{gate5.get('status', 'NOT RUN')}**", ""])
        for backend, result in gate5.get("backends", {}).items():
            memory = result.get("memory", {})
            lines.extend(
                [
                    f"## {backend}",
                    "",
                    f"- Status: **{result.get('status', 'NOT RUN')}**",
                    f"- Memory records: `{memory.get('records')}`",
                    f"- Working set initial/peak/final [bytes]: `{memory.get('working_set_initial')}` / `{memory.get('working_set_peak')}` / `{memory.get('working_set_final')}`",
                    f"- Private usage initial/peak/final [bytes]: `{memory.get('private_bytes_initial')}` / `{memory.get('private_bytes_peak')}` / `{memory.get('private_bytes_final')}`",
                    f"- Private usage monotonic non-decreasing: `{memory.get('private_bytes_monotonic_non_decreasing')}`",
                    "",
                ]
            )
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"window": args.window, "gate4": payload.get("gate4"), "gate5": payload.get("gate5"), "summary": str(out / "summary.json")}, ensure_ascii=False))
    statuses = [payload.get(key, {}).get("status") for key in ("gate4", "gate5") if key in payload]
    return 0 if statuses and all(status == "PASS" for status in statuses) else 1


if __name__ == "__main__":
    raise SystemExit(main())

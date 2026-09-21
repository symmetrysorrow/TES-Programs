"""Run and evaluate Phase24 Gate 6 on the conformal no-mortar route.

The CPU reference is the already qualified Phase24 refine20/no-mortar HYPRE
route.  The target is the native HIP HYPRE build under WSL.  The runner keeps
the physical project and mesh registry unchanged, changing only the HYPRE
device selector and the solver executable.  It records both an independent
steady solve and the short Gate 4 transient window.
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
BASE_PROJECT = ROOT / "artifacts" / "comparison" / "nomortar_refinement_probe" / "project.json"
TIMEGRID_PROJECT = ROOT / "elmer_project_comsol_timegrid.json"
MESH = "mesh_singlepixel_conformal_gpu_refine20"
MESH_TAG = "conformal_refine20"
REFERENCE_CASE = "case_tes_steady_singlepixel_conformal_gpu_refine20"
REFERENCE_CURRENT_UA = 144.71078126230953
CPU_STEADY_CURRENT_UA = 143.52217883888775
COMSOL_CURRENT_UA = 143.055049
PULSE_EVENT_S = 0.020020
LINEAR_SYSTEM = "iterative_hypre_flexgmres_boomeramg_gpu"
CPU_LINEAR_SYSTEM = "iterative_hypre_flexgmres_boomeramg"
LINEAR_MAX_ITERATIONS = 4000
LINEAR_TOLERANCE = 2.0e-8
GMRES_DIMENSION = 100
STRONG_THRESHOLD = 0.5
WAVEFORM_MAX_LIMIT_UA = 0.05
WAVEFORM_RMSE_LIMIT_UA = 0.03
STEADY_LIMIT_PERCENT = 0.05
SHORT_STAGES = [
    ["18[us]", 1],
    ["1[us]", 2],
    ["1[ns]", 1],
    ["10[ns]", 10],
    ["100[ns]", 9],
]
OUT_ROOT = ROOT / "artifacts" / "phase24_gate6_gpu_parity"


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def prepare_steady() -> tuple[Path, str]:
    project = json.loads(BASE_PROJECT.read_text(encoding="utf-8"))
    source_name = "case_tes_steady_singlepixel_conformal_gpu_fine"
    case_name = "case_phase24_g6_gpu_steady_conformal_refine20"
    case = copy.deepcopy(project["cases"][source_name])
    case.update(
        {
            "template": "steady",
            "mesh": MESH,
            "heat_source": "circuit_inner",
            "initial_temperature": "T_0",
            "apply_mortar_bcs": False,
            "restart_from": None,
            "restart_file_path": None,
            "preexisting_restart": False,
            "restart_file_base": None,
            "state_file": f"work/meshes/{MESH}/{case_name}.state",
            "series_file": f"{case_name}_series.csv",
            "iteration_series_file": f"{case_name}_iterations.csv",
            "output_result": True,
            "output_file_path": f"../work/meshes/{MESH}/{case_name}.result",
            "output_intervals": 1,
            "steady_state_max_iterations": 1,
            "vtu": False,
            "solver_comment": "Phase24 Gate 6: independent initial-T0 HIP HYPRE steady solve; conformal no-mortar",
            "comparison_time_grid": {
                "mode": "Phase24 Gate 6 GPU steady parity",
                "mesh": MESH,
                "mortar": False,
            },
            "phase24_smoke": {
                "purpose": "GPU parity steady qualification",
                "independent_initial_temperature": True,
                "no_restart": True,
                "no_mortar": True,
                "backend_selection": "HIP",
            },
            "phase24_hypre_backend": "device",
        }
    )
    case["solver"] = {
        **case.get("solver", {}),
        "linear_system": LINEAR_SYSTEM,
        "linear_system_max_iterations": LINEAR_MAX_ITERATIONS,
        "linear_system_convergence_tolerance": LINEAR_TOLERANCE,
        "hypre_gmres_dimension": GMRES_DIMENSION,
        "boomer_amg_strong_threshold": STRONG_THRESHOLD,
        "nonlinear_max_iterations": 120,
        "nonlinear_convergence_tolerance": 1.0e-8,
        "nonlinear_relaxation_factor": 1.0,
        "steady_state_convergence_tolerance": 1.0e-8,
    }
    project["cases"] = {case_name: case}
    out = OUT_ROOT / "steady" / "project.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return out, case_name


def hybrid_schedule() -> list[list[object]]:
    sys.path.insert(0, str(ROOT))
    from scripts.prep.run_singlepixel_prod_v2_original_timegrid import hybrid_timesteps

    source = json.loads(TIMEGRID_PROJECT.read_text(encoding="utf-8"))
    original = source["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]["timesteps"]
    return hybrid_timesteps(original)


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


def trim_schedule(window: str) -> list[list[object]]:
    if window == "short":
        return copy.deepcopy(SHORT_STAGES)
    target_s = PULSE_EVENT_S + 100.0e-6
    elapsed = 0.0
    result: list[list[object]] = []
    for token_obj, count_obj in hybrid_schedule():
        token = str(token_obj)
        count = int(count_obj)
        dt = seconds(token)
        available = target_s - (0.020000 + elapsed)
        if available <= 1.0e-15:
            break
        take = min(count, max(0, int(math.floor(available / dt + 1.0e-9))))
        if take:
            result.append([token, take])
            elapsed += dt * take
        if take < count:
            break
    if not result or elapsed < 100.0e-6 - 1.0e-9:
        raise ValueError(f"time grid does not cover 100us: elapsed={elapsed}")
    return result


def prepare_transient(window: str) -> tuple[Path, str]:
    project = json.loads(BASE_PROJECT.read_text(encoding="utf-8"))
    template_project = json.loads(TIMEGRID_PROJECT.read_text(encoding="utf-8"))
    template = copy.deepcopy(
        template_project["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]
    )
    case_name = f"case_phase24_g6_gpu_{window}_conformal_refine20"
    state_path = ROOT / "work" / "meshes" / MESH / f"{case_name}.state"
    reference_state = ROOT / "work" / "meshes" / MESH / f"{REFERENCE_CASE}.state"
    if not reference_state.is_file():
        raise FileNotFoundError(f"missing steady restart state: {reference_state}")
    shutil.copy2(reference_state, state_path)
    out_dir = OUT_ROOT / window
    candidate = {
        **template,
        "mesh": MESH,
        "apply_mortar_bcs": False,
        "initial_temperature": "T_0",
        "restart_from": None,
        "restart_file_base": REFERENCE_CASE,
        "restart_file_path": f"../work/meshes/{MESH}/{REFERENCE_CASE}.result",
        "preexisting_restart": True,
        "restart_time": 0.020,
        "timesteps": trim_schedule(window),
        "output_intervals": [1] * len(trim_schedule(window)),
        "series_file": f"{case_name}_series.csv",
        "iteration_series_file": f"{case_name}_iterations.csv",
        "state_file": f"work/meshes/{MESH}/{case_name}.state",
        "output_result": True,
        "output_file_path": f"../work/meshes/{MESH}/{case_name}.result",
        "vtu": False,
        "inner_circuit_step_commit": True,
        "solver_comment": f"Phase24 Gate 6: HIP HYPRE {window} transient parity; conformal no-mortar",
        "comparison_time_grid": {
            "mode": "Phase24 Gate 6 GPU transient parity",
            "event_s": PULSE_EVENT_S,
            "window": window,
            "mesh": MESH,
            "mortar": False,
        },
        "phase24_smoke": {
            "purpose": "GPU parity short transient",
            "independent_initial_temperature": False,
            "restart_reference": REFERENCE_CASE,
            "no_mortar": True,
            "backend_selection": "HIP",
        },
        "phase24_hypre_backend": "device",
    }
    candidate["solver"] = {
        **candidate.get("solver", {}),
        "linear_system": LINEAR_SYSTEM,
        "linear_system_max_iterations": LINEAR_MAX_ITERATIONS,
        "linear_system_convergence_tolerance": LINEAR_TOLERANCE,
        "hypre_gmres_dimension": GMRES_DIMENSION,
        "boomer_amg_strong_threshold": STRONG_THRESHOLD,
        "nonlinear_max_iterations": 120,
        "nonlinear_convergence_tolerance": 1.0e-8,
        "nonlinear_relaxation_factor": 1.0,
    }
    project["cases"] = {case_name: candidate}
    out = out_dir / "project.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return out, case_name


def prepare_completion() -> tuple[Path, str]:
    """Create the one-step tail needed to land exactly on 100 us.

    Elmer records 174 series rows for the 175-step 100-us schedule because
    the final circuit flush is written to the restart result, not the series
    CSV.  The existing result ends at 99.376 us after the pulse; one 0.625-us
    continuation supplies the final comparison sample at 100.001 us.
    """
    project = json.loads(BASE_PROJECT.read_text(encoding="utf-8"))
    template_project = json.loads(TIMEGRID_PROJECT.read_text(encoding="utf-8"))
    template = copy.deepcopy(
        template_project["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]
    )
    source_case = "case_phase24_g6_gpu_100us_conformal_refine20"
    case_name = "case_phase24_g6_gpu_100us_tail2_conformal_refine20"
    source_result = ROOT / "work" / "meshes" / MESH / f"{source_case}.result"
    iteration_path = ROOT / "results" / source_case / f"{source_case}_iterations.csv"
    if not source_result.is_file() or not iteration_path.is_file():
        raise FileNotFoundError("100us GPU result/iterations are required before tail completion")
    state_path = ROOT / "work" / "meshes" / MESH / f"{case_name}.state"
    with iteration_path.open(newline="", encoding="utf-8") as handle:
        iteration_rows = list(csv.DictReader(handle))
    if not iteration_rows:
        raise ValueError("100us GPU iteration CSV is empty")
    last = iteration_rows[-1]
    state_path.write_text(
        "  {tes_temperature_K} {raw_current_A} {tes_resistance_ohm} {raw_power_W} {raw_current_A}\n".format(**last),
        encoding="utf-8",
    )
    candidate = {
        **template,
        "mesh": MESH,
        "apply_mortar_bcs": False,
        "initial_temperature": "T_0",
        "restart_from": None,
        "restart_file_base": source_case,
        "restart_file_path": f"../work/meshes/{MESH}/{source_case}.result",
        "preexisting_restart": True,
        "restart_time": 0.020119376,
        "timesteps": [["0.625[us]", 1]],
        "output_intervals": [1],
        "series_file": f"{case_name}_series.csv",
        "iteration_series_file": f"{case_name}_iterations.csv",
        "state_file": f"work/meshes/{MESH}/{case_name}.state",
        "output_result": True,
        "output_file_path": f"../work/meshes/{MESH}/{case_name}.result",
        "vtu": False,
        "inner_circuit_step_commit": True,
        "solver_comment": "Phase24 Gate 6: exact-100us HIP HYPRE tail completion; conformal no-mortar",
        "comparison_time_grid": {
            "mode": "Phase24 Gate 6 GPU 100us tail completion",
            "event_s": PULSE_EVENT_S,
            "window": "100us",
            "mesh": MESH,
            "mortar": False,
        },
        "phase24_smoke": {
            "purpose": "complete final 100us comparison sample",
            "restart_reference": source_case,
            "no_mortar": True,
            "backend_selection": "HIP",
        },
        "phase24_hypre_backend": "device",
    }
    candidate["solver"] = {
        **candidate.get("solver", {}),
        "linear_system": LINEAR_SYSTEM,
        "linear_system_max_iterations": LINEAR_MAX_ITERATIONS,
        "linear_system_convergence_tolerance": LINEAR_TOLERANCE,
        "hypre_gmres_dimension": GMRES_DIMENSION,
        "boomer_amg_strong_threshold": STRONG_THRESHOLD,
        "nonlinear_max_iterations": 120,
        "nonlinear_convergence_tolerance": 1.0e-8,
        "nonlinear_relaxation_factor": 1.0,
    }
    project["cases"] = {case_name: candidate}
    out = OUT_ROOT / "100us" / "completion_project.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return out, case_name


def sync_project(project: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    (project.parent / "sync.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"sync_elmer_parameters.py failed ({result.returncode})")


def wsl_path(path: Path) -> str:
    resolved = absolute(path)
    drive = resolved.drive[0].lower()
    return f"/mnt/{drive}{resolved.as_posix()[2:]}"


def run_gpu(project: Path, case_name: str, solver: Path, out: Path) -> dict[str, Any]:
    # Linux Elmer resolves legacy ``<mesh>/../work/...`` output/restart paths
    # textually before normalizing ``..``.  The WSL GPU path therefore needs
    # the same empty root-level mesh anchor that the AMGX runner creates.
    (ROOT / MESH).mkdir(parents=True, exist_ok=True)
    solver_wsl = wsl_path(solver)
    project_wsl = wsl_path(project)
    repo_wsl = wsl_path(ROOT)
    prefix_wsl = wsl_path(solver.parent.parent)
    tools_wsl = wsl_path(ROOT.parent / "tools")
    hypre_wsl = f"{tools_wsl}/hypre-hip-v3-1-0-install"
    modules_wsl = f"{prefix_wsl}/share/elmersolver/include"
    lib_wsl = f"{prefix_wsl}/lib/elmersolver"
    udf_circuit = f"{repo_wsl}/tes_parallel_circuit.so"
    udf_pulse = f"{repo_wsl}/tes_transient_heat_source_t0.so"
    bash = f"""set -euo pipefail
export HIP_VISIBLE_DEVICES=0
export HIP_PATH=/opt/rocm/core-7.14
export ROCM_PATH=/opt/rocm/core-7.14
export PATH=/opt/rocm/core-7.14/bin:/opt/rocm/core-7.14/lib/llvm/bin:/opt/rocm-wsl/bin:/usr/lib/wsl/lib:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export ELMER_HOME='{prefix_wsl}'
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:{hypre_wsl}/lib:{lib_wsl}:{repo_wsl}:/opt/rocm-wsl/lib:/opt/rocm/core-7.14/lib
export PHASE24_MEMORY_TRACE=1
cd '{repo_wsl}'
gfortran -O2 -fPIC -shared -I'{modules_wsl}' tes_parallel_circuit.f90 -L'{lib_wsl}' -Wl,-rpath,'{lib_wsl}' -lelmersolver -o '{udf_circuit}'
gfortran -O2 -fPIC -shared -I'{modules_wsl}' tes_transient_heat_source.f90 -L'{lib_wsl}' -Wl,-rpath,'{lib_wsl}' -lelmersolver -o '{udf_pulse}'
python3 run.py '{case_name}' --project '{project_wsl}' --skip-sync --mpi-procs 1 --elmer-solver '{solver_wsl}' --runtime-bin ''
"""
    log = out / "launcher.log"
    started = datetime.now(timezone.utc).isoformat()
    start = time.monotonic()
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", bash],
            cwd=ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    return {
        "command": ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", bash],
        "started": started,
        "finished": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - start,
        "exit_code": result.returncode,
        "launcher_log": str(log),
    }


def audit(case_name: str) -> dict[str, Any]:
    result_dir = ROOT / "results" / case_name
    log_path = result_dir / "solver.log"
    manifest_path = result_dir / "manifest.json"
    text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    counters = re.findall(
        r"Phase24 GPU counters: backend matrix_H2D rhs_H2D solution_H2D solution_D2H migrations allocations=([^\r\n]+)",
        text,
    )
    backend_hits = re.findall(r"Phase24 HYPRE backend=([A-Za-z0-9_]+), persistent device execution enabled", text)
    device_lines = re.findall(r"Running on \"([^\r\n]+)\", gfx[0-9]+, Total VRAM:[^\r\n]+", text)
    return {
        "manifest_exit_code": manifest.get("exit_code"),
        "solver_completed": manifest.get("solver_completed"),
        "all_done": "MAIN: *** Elmer Solver: ALL DONE ***" in text,
        "stop_1": text.count("STOP 1"),
        "mpi_abort": "MPI_ABORT" in text,
        "native_crash": any(x in text for x in ("Access violation", "Could not allocate memory", "Fortran runtime error")),
        "hypre_backend_hits": backend_hits,
        "gpu_counter_records": counters,
        # The compact conformal SIF does not enable the optional Phase24
        # lifecycle counter block, but the HIP device line plus the native
        # backend marker prove actual device selection.  Counters remain
        # recorded when the selected build emits them.
        "gpu_execution_marker": bool(backend_hits and device_lines),
        "gpu_counter_records_present": bool(counters),
        "device_lines": device_lines,
        "cpu_fallback_marker": bool(re.search(r"HYPRE backend=CPU|falling back to CPU|CPU fallback", text, re.I)),
        "device_error_marker": bool(re.search(r"HIP (?:ERROR|error)|hipError|device resource leak", text)),
        "log_path": str(log_path),
        "manifest_path": str(manifest_path),
    }


def read_series(path: Path) -> tuple[list[float], list[float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    values = sorted((float(row["time_s"]), float(row["tes_current_A"]) * 1.0e6) for row in rows)
    return [x for x, _ in values], [y for _, y in values]


def merge_series(paths: list[Path], out: Path) -> Path:
    rows_by_time: dict[float, dict[str, str]] = {}
    fieldnames: list[str] | None = None
    for path in paths:
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if fieldnames is None:
                fieldnames = list(reader.fieldnames or [])
            for row in reader:
                rows_by_time[float(row["time_s"])] = row
    if not rows_by_time or not fieldnames:
        raise FileNotFoundError("no GPU series rows available for merge")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_by_time[key] for key in sorted(rows_by_time))
    return out


def series_row_from_last_iteration(path: Path, out: Path) -> Path:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"iteration CSV is empty: {path}")
    row = rows[-1]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "time_s",
                "tes_temperature_K",
                "tes_current_A",
                "tes_resistance_ohm",
                "tes_power_W",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "time_s": row["time_s"],
                "tes_temperature_K": row["tes_temperature_K"],
                "tes_current_A": row["raw_current_A"],
                "tes_resistance_ohm": row["tes_resistance_ohm"],
                "tes_power_W": row["raw_power_W"],
            }
        )
    return out


def direct_metrics(cpu_path: Path, gpu_path: Path) -> dict[str, float]:
    cpu_t, cpu_i = read_series(cpu_path)
    gpu_t, gpu_i = read_series(gpu_path)
    start = max(min(cpu_t), min(gpu_t))
    end = min(max(cpu_t), max(gpu_t))
    if end <= start:
        raise ValueError("CPU/GPU series do not overlap")
    step = min(0.05e-6, (end - start) / 100.0)
    count = max(2, int(math.ceil((end - start) / step)) + 1)
    grid = [start + (end - start) * i / (count - 1) for i in range(count)]
    import numpy as np

    diff = np.interp(grid, cpu_t, cpu_i) - np.interp(grid, gpu_t, gpu_i)
    return {
        "comparison_start_s": start,
        "comparison_end_s": end,
        "max_abs_difference_uA": float(np.max(np.abs(diff))),
        "rmse_uA": float(np.sqrt(np.mean(diff * diff))),
        "cpu_steady_last_uA": cpu_i[-1],
        "gpu_steady_last_uA": gpu_i[-1],
    }


def verify_gpu_sif(case_name: str) -> bool:
    path = ROOT / "generated" / "cases" / f"{case_name}.sif"
    text = path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""
    return all(
        token in text
        for token in (
            f'"{MESH}"',
            "Apply Mortar BCs = False",
            "HYPRE GPU = Logical True",
            "Linear System Max Iterations = 4000",
            "Linear System Convergence Tolerance = 2e-08",
        )
    )


def evaluate(steady_audit: dict[str, Any], transient_audit: dict[str, Any], steady_current: float | None, waveform: dict[str, Any] | None, cpu_gpu: dict[str, float] | None, sif_policy_ok: bool, completion_audit: dict[str, Any] | None = None) -> dict[str, Any]:
    steady_delta = None if steady_current is None else 100.0 * abs(steady_current - CPU_STEADY_CURRENT_UA) / CPU_STEADY_CURRENT_UA
    criteria = {
        "steady_gpu_normal_exit": steady_audit["manifest_exit_code"] == 0 and steady_audit["all_done"],
        "transient_gpu_normal_exit": transient_audit["manifest_exit_code"] == 0 and transient_audit["all_done"],
        "100us_tail_completion_normal_exit": (
            completion_audit is None
            or (completion_audit["manifest_exit_code"] == 0 and completion_audit["all_done"])
        ),
        "steady_current_within_cpu_hypre_0p05_percent": steady_delta is not None and steady_delta <= STEADY_LIMIT_PERCENT,
        "transient_waveform_max_within_gate4_0p05_uA": waveform is not None and waveform["max_abs_difference_uA"] <= WAVEFORM_MAX_LIMIT_UA,
        "transient_waveform_rmse_within_gate4_0p03_uA": waveform is not None and waveform["rmse_uA"] <= WAVEFORM_RMSE_LIMIT_UA,
        "gpu_execution_marker_present": steady_audit["gpu_execution_marker"] and transient_audit["gpu_execution_marker"],
        "no_cpu_fallback": not steady_audit["cpu_fallback_marker"] and not transient_audit["cpu_fallback_marker"],
        "no_device_error_or_native_crash": not any(
            item["device_error_marker"] or item["native_crash"] or item["mpi_abort"]
            for item in (steady_audit, transient_audit, *( [completion_audit] if completion_audit else [] ))
        ),
        "same_mesh_no_mortar_and_gpu_sif": sif_policy_ok,
    }
    return {
        "status": "PASS" if all(criteria.values()) else "FAIL",
        "criteria": criteria,
        "steady_gpu_current_uA": steady_current,
        "steady_cpu_hypre_current_uA": CPU_STEADY_CURRENT_UA,
        "steady_delta_percent": steady_delta,
        "cpu_gpu_transient": cpu_gpu,
        "gpu_vs_comsol": waveform,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("hip",), default="hip")
    parser.add_argument("--solver", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre-hip-v3-1-0-wsl\bin\ElmerSolver_mpi"))
    parser.add_argument("--window", choices=("short", "100us"), default="short")
    parser.add_argument("--rerun-steady", action="store_true", help="also rerun the independent GPU steady case")
    parser.add_argument("--rerun-transient", action="store_true", help="rerun the selected GPU transient case")
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    solver = absolute(args.solver)
    if not solver.is_file():
        raise SystemExit(f"GPU solver not found: {solver}")

    steady_project, steady_case = prepare_steady()
    transient_project, transient_case = prepare_transient(args.window)
    out = OUT_ROOT if args.window == "short" else OUT_ROOT / args.window
    existing_transient_audit = audit(transient_case)
    if not args.audit_only:
        if args.window == "short" or args.rerun_steady:
            sync_project(steady_project)
            steady_run = run_gpu(steady_project, steady_case, solver, steady_project.parent)
        else:
            steady_run = {"reused_existing_run": True}
        if (
            args.window == "short"
            or args.rerun_transient
            or not (existing_transient_audit["manifest_exit_code"] == 0 and existing_transient_audit["all_done"])
        ):
            sync_project(transient_project)
            transient_run = run_gpu(transient_project, transient_case, solver, transient_project.parent)
        else:
            transient_run = {"reused_existing_run": True}
    else:
        steady_run = {}
        transient_run = {}

    steady_audit = audit(steady_case)
    transient_audit = audit(transient_case)
    completion_project: Path | None = None
    completion_case: str | None = None
    completion_run: dict[str, Any] = {}
    completion_audit: dict[str, Any] | None = None
    if args.window == "100us":
        completion_project, completion_case = prepare_completion()
        completion_audit = audit(completion_case)
        if not args.audit_only and not (
            completion_audit["manifest_exit_code"] == 0 and completion_audit["all_done"]
        ):
            sync_project(completion_project)
            completion_run = run_gpu(completion_project, completion_case, solver, completion_project.parent)
            completion_audit = audit(completion_case)
    steady_series = ROOT / "results" / steady_case / f"{steady_case}_series.csv"
    transient_series = ROOT / "results" / transient_case / f"{transient_case}_series.csv"
    steady_current = None
    if steady_series.is_file():
        _, values = read_series(steady_series)
        if values:
            steady_current = values[-1]
    if steady_current is None:
        iteration_path = ROOT / "results" / steady_case / f"{steady_case}_iterations.csv"
        if iteration_path.is_file():
            with iteration_path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            if rows:
                raw = rows[-1].get("raw_current_A")
                if raw:
                    steady_current = float(raw) * 1.0e6

    waveform = None
    cpu_gpu = None
    comparison_series = transient_series
    if args.window == "100us" and completion_case:
        completion_series = ROOT / "results" / completion_case / f"{completion_case}_series.csv"
        if not completion_series.is_file():
            completion_iterations = ROOT / "results" / completion_case / f"{completion_case}_iterations.csv"
            if completion_iterations.is_file():
                completion_series = series_row_from_last_iteration(
                    completion_iterations,
                    out / "gpu_100us_tail_series.csv",
                )
        if completion_series.is_file():
            comparison_series = merge_series(
                [transient_series, completion_series],
                out / "gpu_100us_merged_series.csv",
            )
    if comparison_series.is_file():
        sys.path.insert(0, str(ROOT))
        from scripts.analysis import compare_singlepixel_amgx_comsol as compare

        requested_window_us = 0.9 if args.window == "short" else 100.0
        comsol_series = compare.read_comsol(ROOT / "docs" / "Single-Pixel.txt")
        gpu_series = compare.read_elmer(comparison_series)
        waveform, aligned = compare.compare(
            comsol_series,
            gpu_series,
            requested_window_us,
        )
        compare.write_outputs(
            out / "gpu_vs_comsol",
            comsol_series,
            gpu_series,
            waveform,
            aligned,
            "HIP HYPRE GPU",
        )
        cpu_name = "short" if args.window == "short" else "100us"
        cpu_path = ROOT / "results" / f"case_phase24_g45_conformal_refine20_{cpu_name}_hypre" / f"case_phase24_g45_conformal_refine20_{cpu_name}_hypre_series.csv"
        if cpu_path.is_file():
            cpu_gpu = direct_metrics(cpu_path, comparison_series)

    payload = {
        "schema": "tes.phase24.gate6.gpu.parity.v1",
        "window": args.window,
        "backend": "HIP",
        "policy": {
            "mesh": MESH,
            "mortar": False,
            "linear_system_cpu": CPU_LINEAR_SYSTEM,
            "linear_system_gpu": LINEAR_SYSTEM,
            "max_iterations": LINEAR_MAX_ITERATIONS,
            "tolerance": LINEAR_TOLERANCE,
            "gmres_dimension": GMRES_DIMENSION,
            "boomer_amg_strong_threshold": STRONG_THRESHOLD,
            "gpu_solver": str(solver),
        },
        "runs": {
            "steady": steady_run,
            f"transient_{args.window}": transient_run,
            "100us_tail_completion": completion_run,
        },
        "audits": {
            "steady": steady_audit,
            f"transient_{args.window}": transient_audit,
            "100us_tail_completion": completion_audit,
        },
        "provenance": {
            "steady_project_sha256": sha256(steady_project),
            "transient_project_sha256": sha256(transient_project),
            "steady_sif_sha256": sha256(ROOT / "generated" / "cases" / f"{steady_case}.sif"),
            "transient_sif_sha256": sha256(ROOT / "generated" / "cases" / f"{transient_case}.sif"),
            "mesh_header_sha256": sha256(ROOT / "work" / "meshes" / MESH / "mesh.header"),
            "solver_sha256": sha256(solver),
            "hip_device_validation": str(ROOT / "artifacts" / "hip_device_validation.json"),
        },
    }
    sif_policy_ok = verify_gpu_sif(steady_case) and verify_gpu_sif(transient_case)
    if completion_case:
        sif_policy_ok = sif_policy_ok and verify_gpu_sif(completion_case)
    payload["gate6"] = evaluate(
        steady_audit,
        transient_audit,
        steady_current,
        waveform,
        cpu_gpu,
        sif_policy_ok,
        completion_audit,
    )
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    gate = payload["gate6"]
    lines = [
        "# Phase24 Gate 6: HIP HYPRE GPU parity",
        "",
        f"- Status: **{gate['status']}**",
        f"- Mesh: `{MESH}`; mortar: `False`",
        "- Backend: `HIP HYPRE device`",
        f"- Steady GPU current [µA]: `{steady_current}`",
        f"- CPU HYPRE steady current [µA]: `{CPU_STEADY_CURRENT_UA}`",
        f"- Steady delta [%]: `{gate['steady_delta_percent']}`",
        f"- Requested transient window [µs]: `{100.0 if args.window == '100us' else 0.9}`",
        f"- Evaluated transient window [µs]: `{(waveform or {}).get('comparison_window_us', [None, None])[1]}`",
        f"- GPU vs COMSOL max difference [µA]: `{(waveform or {}).get('max_abs_difference_uA')}`",
        f"- GPU vs COMSOL RMSE [µA]: `{(waveform or {}).get('rmse_uA')}`",
        f"- CPU-HYPRE vs GPU max difference [µA]: `{(cpu_gpu or {}).get('max_abs_difference_uA')}`",
        f"- CPU-HYPRE vs GPU RMSE [µA]: `{(cpu_gpu or {}).get('rmse_uA')}`",
        "",
        "| criterion | result |",
        "|---|---|",
    ]
    lines.extend(f"| `{key}` | `{value}` |" for key, value in gate["criteria"].items())
    lines += [
        "",
        f"- GPU solver log (steady): `{steady_audit['log_path']}`",
        f"- GPU solver log (transient {args.window}): `{transient_audit['log_path']}`",
        f"- Machine-readable artifact: `{out / 'summary.json'}`",
        "",
    ]
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"status": gate["status"], "steady_current_uA": steady_current, "steady_delta_percent": gate["steady_delta_percent"], "gpu_execution": gate["criteria"]["gpu_execution_marker_present"]}, ensure_ascii=False))
    return 0 if gate["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

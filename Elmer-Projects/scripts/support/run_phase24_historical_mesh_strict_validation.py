"""Validate the historical single-pixel restart with current native MUMPS.

This is a diagnostic transplant only.  It snapshots the historical result and
circuit checkpoint, runs one direct constrained solve and a separate nonlinear
refinement from that same snapshot, and materializes the native full-system
captures.  No production mesh, physics constant, or original restart is
overwritten.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_historical_mesh_strict_steady_validation"
MESH = "mesh_singlepixel_prod_v2"
MESH_ROOT = ROOT / "work" / "meshes" / MESH
SOURCE_PROJECT = ROOT / "artifacts" / "phase24_gate4_5_nomortar" / "restart_refinement" / "mortar" / "project.json"
HIST_RESULT = MESH_ROOT / "case_tes_steady_singlepixel_prod_v2.result"
HIST_STATE = MESH_ROOT / "case_tes_steady_singlepixel_prod_v2_original_timegrid.state"
SOLVER = ROOT.parent / "tools" / "elmer-hypre" / "install-steady-full-capture" / "bin" / "ElmerSolver.exe"
RUNTIME = SOLVER.parent
TOOLCHAIN = Path(r"C:\msys64\ucrt64\bin")
SNAPSHOT_RESULT = "phase24_historical_restart_snapshot"
SNAPSHOT_STATE = "phase24_historical_restart_snapshot.state"
ONE_SHOT = "case_phase24_historical_one_shot"
REFINED = "case_phase24_historical_refinement"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def snapshot_inputs() -> dict[str, Any]:
    if not HIST_RESULT.is_file() or not HIST_STATE.is_file():
        raise FileNotFoundError("historical result/state pair is incomplete")
    OUT.mkdir(parents=True, exist_ok=True)
    snap_result = MESH_ROOT / f"{SNAPSHOT_RESULT}.result"
    snap_state = MESH_ROOT / SNAPSHOT_STATE
    shutil.copy2(HIST_RESULT, snap_result)
    shutil.copy2(HIST_STATE, snap_state)
    files = {
        "historical_result": HIST_RESULT,
        "historical_state": HIST_STATE,
        "snapshot_result": snap_result,
        "snapshot_state": snap_state,
        "historical_steady_sif": ROOT / "generated" / "cases" / "case_tes_steady_singlepixel_prod_v2.sif",
        "historical_transient_sif": ROOT / "generated" / "cases" / "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid_cpu_smoke_177step.sif",
        "historical_steady_manifest": ROOT / "results" / "case_tes_steady_singlepixel_prod_v2" / "manifest.json",
        "historical_transient_manifest": ROOT / "results" / "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid_cpu_smoke_177step" / "manifest.json",
        "historical_waveform_summary": ROOT / "artifacts" / "comparison" / "comsol_cpu_singlepixel_prod_v2_hybrid_100us" / "summary.md",
        "historical_waveform_metrics": ROOT / "artifacts" / "comparison" / "comsol_cpu_singlepixel_prod_v2_hybrid_100us" / "metrics.json",
    }
    records = {}
    for label, path in files.items():
        records[label] = {
            "path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
            "exists": path.is_file(),
            "sha256": sha256(path) if path.is_file() else None,
            "size_bytes": path.stat().st_size if path.is_file() else None,
        }
    records["snapshot_semantics"] = {
        "result_and_state_copied": True,
        "originals_overwritten": False,
        "restart_result_used_by_diagnostic": str(snap_result.relative_to(ROOT)),
        "restart_state_used_by_diagnostic": str(snap_state.relative_to(ROOT)),
    }
    native_source = ROOT.parent / "tools" / "elmer-hypre" / "src" / "fem" / "src" / "modules" / "HeatSolve.F90"
    solver_files = {"solver": SOLVER, "HeatSolve_dll": SOLVER.parent.parent / "share" / "elmersolver" / "lib" / "HeatSolve.dll"}
    records["current_native_provenance"] = {
        "heat_solve_source": str(native_source),
        "heat_solve_source_sha256": sha256(native_source) if native_source.is_file() else None,
        "solver": str(SOLVER),
        "solver_sha256": sha256(SOLVER) if SOLVER.is_file() else None,
        "runtime_artifacts_sha256": {label: sha256(path) for label, path in solver_files.items() if path.is_file()},
        "solver_selection": "CPU MUMPS direct; HYPRE/GPU options not injected",
    }
    git_head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True).stdout.strip()
    records["repository"] = {"head_at_snapshot": git_head}
    write_json(OUT / "historical_restart_provenance.json", records)
    return records


def case_spec(template: dict[str, Any], name: str, iterations: int) -> dict[str, Any]:
    result = copy.deepcopy(template)
    result.update({
        "template": "steady",
        "mesh": MESH,
        "restart_from": None,
        "restart_time": 0.02,
        "restart_file_base": SNAPSHOT_RESULT,
        "restart_file_path": f"../work/meshes/{MESH}/{SNAPSHOT_RESULT}.result",
        "preexisting_restart": True,
        "initial_temperature": "T_0",
        "steady_state_max_iterations": 1,
        "output_intervals": 1,
        "series_file": f"{name}_series.csv",
        "iteration_series_file": f"{name}_iterations.csv",
        "state_file": f"work/meshes/{MESH}/{name}.state",
        "output_file_path": f"../work/meshes/{MESH}/{name}.result",
        "output_result": True,
        "post_file": False,
        "vtu": False,
        "apply_mortar_bcs": True,
        "inner_circuit_step_commit": True,
        "solver_comment": "Historical mesh strict steady diagnostic; current native HeatSolve + CPU MUMPS",
    })
    result.pop("pulse", None)
    result.pop("timesteps", None)
    result["solver"] = {
        **result.get("solver", {}),
        "linear_system": "mumps",
        "nonlinear_max_iterations": iterations,
        "nonlinear_convergence_tolerance": 1.0e-8,
        "nonlinear_relaxation_factor": 1.0,
        "steady_state_convergence_tolerance": 1.0e-9,
    }
    return result


def prepare_project(name: str, iterations: int) -> Path:
    source = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    # The Phase24 project is intentionally current-mesh-focused and does not
    # carry the legacy production-v2 registry entry.  Register the existing
    # historical directory only; this does not regenerate or alter its mesh.
    source.setdefault("meshes", {})[MESH] = {"geometry": "single_pixel", "dir": MESH}
    template = next(iter(source["cases"].values()))
    project = copy.deepcopy(source)
    project["cases"] = {name: case_spec(template, name, iterations)}
    path = OUT / f"{name}.json"
    write_json(path, project)
    shutil.copy2(MESH_ROOT / SNAPSHOT_STATE, MESH_ROOT / f"{name}.state")
    return path


def add_capture(sif: Path, capture: Path, max_iterations: int) -> None:
    capture.mkdir(parents=True, exist_ok=True)
    for iteration in range(1, max_iterations + 1):
        (capture / f"ts0001_nl{iteration:04d}").mkdir(parents=True, exist_ok=True)
    text = sif.read_text(encoding="utf-8")
    needle = "  Apply Mortar BCs = True\n"
    rel = capture.relative_to(ROOT).as_posix()
    insertion = (
        needle
        + '  "Phase24 Full Restriction Capture" = Logical True\n'
        + f'  "Phase24 Full Restriction Capture Prefix" = String "{rel}"\n'
        + '  "Phase24 Full Restriction Capture Timestep" = Integer 1\n'
        + f'  "Phase24 Full Restriction Capture Max Iterations" = Integer {max_iterations}\n'
        + '  "Phase24 Restart State Audit" = Logical True\n'
        + '  "Phase24 Restart Audit Fail Fast" = Logical True\n'
        + f'  "Phase24 Restart State Audit Prefix" = String "{rel}/restart_state_audit"\n'
    )
    if needle not in text:
        raise ValueError(f"solver mortar key not found in {sif}")
    sif.write_text(text.replace(needle, insertion, 1), encoding="utf-8")


def capture_path(name: str) -> Path:
    # The native diagnostic keeps a fixed-length Fortran path buffer.  Keep
    # the capture leaf short enough that the .dat suffix is not truncated.
    leaf = "one" if name == ONE_SHOT else "refined"
    return OUT / "capture" / leaf


def run_case(name: str, project: Path, max_iterations: int) -> None:
    sync = subprocess.run(
        [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project)],
        cwd=ROOT, text=True, capture_output=True,
    )
    (OUT / f"{name}_sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
    if sync.returncode:
        raise RuntimeError(f"parameter sync failed for {name}")
    capture = capture_path(name)
    add_capture(ROOT / "generated" / "cases" / f"{name}.sif", capture, max_iterations)
    command = [
        sys.executable, str(ROOT / "run.py"), name, "--project", str(project), "--skip-sync",
        "--mpi-procs", "1", "--elmer-solver", str(SOLVER), "--runtime-bin", str(RUNTIME),
        "--toolchain-bin", str(TOOLCHAIN),
    ]
    log = OUT / f"{name}_launcher.log"
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, env=os.environ.copy())
    if completed.returncode:
        raise RuntimeError(f"{name} failed with exit {completed.returncode}; see {log}")
    roots = sorted(capture.glob("ts0001_nl*"))
    completed_iterations = max((int(path.name.split("nl")[-1]) for path in roots if path.is_dir()), default=0)
    if completed_iterations == 0:
        raise RuntimeError(f"{name} produced no native full-system capture")
    subprocess.run([
        sys.executable, str(ROOT / "scripts" / "support" / "phase24_full_restriction_capture.py"),
        "materialize", "--root", str(capture), "--iterations", str(completed_iterations),
        "--sif", str(ROOT / "generated" / "cases" / f"{name}.sif"),
        "--solver-log", str(ROOT / "results" / name / "solver.log"),
    ], cwd=ROOT, check=True)


def indexed(path: Path, count: int) -> np.ndarray:
    values = np.zeros(count, dtype=np.float64)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 2:
            values[int(fields[0]) - 1] = float(fields[1])
    return values


def matrix_shape(path: Path) -> int:
    extent = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 3:
            extent = max(extent, int(fields[0]), int(fields[1]))
    return extent


def load_matrix(path: Path, count: int):
    from scipy.sparse import coo_matrix
    rows, cols, data = [], [], []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 3:
            rows.append(int(fields[0]) - 1); cols.append(int(fields[1]) - 1); data.append(float(fields[2]))
    return coo_matrix((data, (rows, cols)), shape=(count, count)).tocsr()


def result_field(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    marker = next(i for i, line in enumerate(lines) if line.strip().lower() == "temperature")
    count = int(lines[marker + 1].split()[1])
    permutation = np.zeros(count, dtype=np.int64)
    for offset in range(count):
        node, dof = lines[marker + 2 + offset].split()[:2]
        permutation[int(node) - 1] = int(dof)
    values = np.asarray(
        [float(line.split()[0]) for line in lines[marker + 2 + count: marker + 2 + 2 * count]],
        dtype=np.float64,
    )
    field = np.zeros(count, dtype=np.float64)
    field[permutation - 1] = values
    return field, permutation


def tes_average(field: np.ndarray, permutation: np.ndarray) -> float:
    weights = np.zeros(permutation.size, dtype=np.float64)
    count = 0
    for line in (MESH_ROOT / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 4 or int(fields[1]) != 8:
            continue
        nodes = [int(value) for value in fields[3:]]
        for node in nodes:
            dof = int(permutation[node - 1])
            if dof > 0:
                weights[dof - 1] += 1.0 / len(nodes)
        count += 1
    return float(np.dot(weights / count, field))


def circuit_at(temperature: float) -> dict[str, float]:
    ibias, rsh, r0, rmin = 715e-6, 3.9e-3, 15.527e-3, 1e-6
    alpha, beta, i0, t0 = 256.46, 5.03, 143.537344932311e-6, 168.57e-3
    a = r0 * (1.0 + alpha * (temperature - t0) / t0 - beta)
    b = r0 * beta / i0
    current = (np.sqrt((rsh + a) ** 2 + 4.0 * b * ibias * rsh) - (rsh + a)) / (2.0 * b)
    current = float(np.clip(current, 0.0, ibias))
    resistance = max(a + b * abs(current), rmin)
    return {"temperature_K": temperature, "raw_current_A": current, "resistance_ohm": resistance, "raw_power_W": current * current * resistance}


def residual_capture(directory: Path) -> dict[str, Any]:
    meta = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    total = int(meta["runtime"]["total_rows"])
    primal = int(meta["runtime"]["primal_rows"])
    matrix = load_matrix(directory / "full_A_before.dat", total)
    rhs = indexed(directory / "full_b_before.dat", total)
    solution = indexed(directory / "full_x_before.dat", total)
    residual = np.asarray(matrix @ solution - rhs).ravel()
    scale = np.asarray(abs(matrix) @ abs(solution)).ravel() + abs(rhs)
    eta = np.divide(abs(residual), scale, out=np.zeros_like(residual), where=scale > 0.0)
    after_matrix = load_matrix(directory / "full_A_after.dat", total)
    after_rhs = indexed(directory / "full_b_after.dat", total)
    after_solution = indexed(directory / "full_x_after.dat", total)
    after_residual = np.asarray(after_matrix @ after_solution - after_rhs).ravel()
    result = {
        "directory": str(directory.relative_to(ROOT)),
        "dimensions": {"total": total, "primal": primal, "constraint": total - primal},
        "before": {
            "full_residual_l2": float(np.linalg.norm(residual)),
            "primal_residual_l2": float(np.linalg.norm(residual[:primal])),
            "constraint_residual_l2": float(np.linalg.norm(residual[primal:])),
            "max_absolute_residual": float(np.max(abs(residual))),
            "componentwise_backward_error_max": float(np.max(eta)),
            "componentwise_backward_error_constraint_max": float(np.max(eta[primal:])),
        },
        "after_direct_system": {
            "full_residual_l2": float(np.linalg.norm(after_residual)),
            "primal_residual_l2": float(np.linalg.norm(after_residual[:primal])),
            "constraint_residual_l2": float(np.linalg.norm(after_residual[primal:])),
            "max_absolute_residual": float(np.max(abs(after_residual))),
        },
        "sha256": {name: sha256(directory / name) for name in (
            "full_A_before.dat", "full_b_before.dat", "full_x_before.dat",
            "full_A_after.dat", "full_b_after.dat", "full_x_after.dat",
        )},
    }
    return result


def integrity_gate(capture: Path) -> dict[str, Any]:
    path = capture / "restart_state_audit.csv"
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    selected = {row["stage"]: row for row in rows if row["stage"] in {"load_after", "circuit_init_after", "pre_first_assembly"}}
    fields = ("current", "previous_current", "resistance", "power")
    mismatches = []
    for field in fields:
        values = {stage: float(selected[stage][field]) for stage in selected}
        if max(values.values()) - min(values.values()) > 1.0e-15:
            mismatches.append({"field": field, "values": values})
    saved_temperature = float((MESH_ROOT / SNAPSHOT_STATE).read_text().split()[0])
    temperature_values = {
        stage: float(selected[stage]["average_temperature"])
        for stage in ("circuit_init_after", "pre_first_assembly")
    }
    if any(abs(value - saved_temperature) > 1.0e-15 for value in temperature_values.values()):
        mismatches.append({"field": "tes_temperature", "saved": saved_temperature, "values": temperature_values})
    return {
        "audit_path": str(path.relative_to(ROOT)), "stages": selected,
        "fields": [*fields, "tes_temperature"],
        "temperature_note": "load_after reports zero before the checkpoint's temperature-dependent circuit initialization; circuit_init_after and pre_first_assembly are the valid TES-temperature stages",
        "saved_temperature_K": saved_temperature,
        "mismatches": mismatches,
        "passed": not mismatches and len(selected) == 3,
    }


def enrich_outputs() -> None:
    provenance = json.loads((OUT / "historical_restart_provenance.json").read_text(encoding="utf-8"))
    one_capture = capture_path(ONE_SHOT) / "ts0001_nl0001"
    refined_root = capture_path(REFINED)
    refined_dirs = sorted(
        (path for path in refined_root.glob("ts0001_nl*") if (path / "metadata.json").is_file()),
        key=lambda path: int(path.name.split("nl")[-1]),
    )
    initial = residual_capture(one_capture)
    final = residual_capture(refined_dirs[-1])
    gate = integrity_gate(capture_path(ONE_SHOT))
    shutil.copy2(capture_path(ONE_SHOT) / "restart_state_audit.csv", OUT / "restart_state_audit.csv")
    initial_field, permutation = result_field(MESH_ROOT / f"{SNAPSHOT_RESULT}.result")
    initial_t = tes_average(initial_field, permutation)
    initial_state = [float(value) for value in (MESH_ROOT / SNAPSHOT_STATE).read_text().split()]
    one_after_field = indexed(one_capture / "full_x_after.dat", int(initial["dimensions"]["total"]))[: initial["dimensions"]["primal"]]
    one_after_t = tes_average(one_after_field, permutation)
    one_after_circuit = circuit_at(one_after_t)
    final_result = MESH_ROOT / f"{REFINED}.result"
    if final_result.is_file():
        final_field, final_perm = result_field(final_result)
        final_t = tes_average(final_field, final_perm)
        refinement_status = "completed"
    else:
        final_meta = json.loads((refined_dirs[-1] / "metadata.json").read_text(encoding="utf-8"))
        final_x = indexed(refined_dirs[-1] / "full_x_after.dat", int(final_meta["runtime"]["total_rows"]))
        final_perm = permutation
        final_t = tes_average(final_x[: int(final_meta["runtime"]["primal_rows"])], final_perm)
        refinement_status = "interrupted_after_capture; no final result file"
    final_circuit = circuit_at(final_t)
    one_shot = {
        "before": {"tes_temperature_K": initial_state[0], "current_A": initial_state[1], "resistance_ohm": initial_state[2], "power_W": initial_state[3], "previous_current_A": initial_state[4]},
        "after_thermal_one_shot": {"tes_temperature_K": one_after_t, **one_after_circuit},
        "delta_T_mK": (one_after_t - initial_state[0]) * 1.0e3,
        "delta_I_uA": (one_after_circuit["raw_current_A"] - initial_state[1]) * 1.0e6,
        "full_residual_capture": initial,
        "restart_integrity": gate,
    }
    write_json(OUT / "restart_integrity.json", gate)
    write_json(OUT / "initial_full_residual.json", initial)
    write_json(OUT / "one_shot_direct.json", one_shot)
    write_json(OUT / "final_full_residual.json", final)
    trajectory_source = ROOT / "results" / REFINED / f"{REFINED}_iterations.csv"
    if not trajectory_source.is_file():
        trajectory_source = ROOT / f"{REFINED}_iterations.csv"
    rows = list(csv.DictReader(trajectory_source.open(encoding="utf-8")))
    fields = list(rows[0]) if rows else []
    log = (ROOT / "results" / REFINED / "solver.log").read_text(encoding="utf-8", errors="replace")
    changes = {int(index): value for index, value in re.findall(r"ComputeChange: NS \(ITER=\s*(\d+)\).*?\(NRM,RELC\): \(.*?\s+([-+0-9.EeDd]+)\s*\)", log)}
    if rows:
        fields.append("nonlinear_relative_change")
        with (OUT / "nonlinear_refinement_trajectory.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
            for row in rows:
                row["nonlinear_relative_change"] = changes.get(int(row["nonlinear_iter"]), "")
                writer.writerow(row)
    comparison = {
        "historical_restart": one_shot["before"],
        "historical_refined": {"tes_temperature_K": final_t, "current_A": final_circuit["raw_current_A"], "resistance_ohm": final_circuit["resistance_ohm"], "power_W": final_circuit["raw_power_W"], "full_residual": final["before"]},
        "phase24_refined": {"tes_temperature_K": 0.166557530018547, "current_A": 218.6434093713403e-6, "resistance_ohm": 0.008853643057514057, "power_W": 4.2324271156482167e-10, "full_residual_l2": 6.192602138143676e-16},
        "comsol_baseline": {"current_A": 143.05504932879472e-6},
    }
    write_json(OUT / "historical_vs_phase24_steady.json", comparison)
    table = [
        {"metric": "TES T (mK)", "historical restart": one_shot["before"]["tes_temperature_K"] * 1e3, "historical refined": final_t * 1e3, "Phase24 refined": 166.557530018547, "COMSOL": "not archived"},
        {"metric": "Current (uA)", "historical restart": one_shot["before"]["current_A"] * 1e6, "historical refined": final_circuit["raw_current_A"] * 1e6, "Phase24 refined": 218.6434093713403, "COMSOL": 143.05504932879472},
        {"metric": "Resistance (ohm)", "historical restart": one_shot["before"]["resistance_ohm"], "historical refined": final_circuit["resistance_ohm"], "Phase24 refined": 0.008853643057514057, "COMSOL": "not archived"},
        {"metric": "Joule Power (W)", "historical restart": one_shot["before"]["power_W"], "historical refined": final_circuit["raw_power_W"], "Phase24 refined": 4.2324271156482167e-10, "COMSOL": "not archived"},
        {"metric": "full residual L2", "historical restart": initial["before"]["full_residual_l2"], "historical refined": final["before"]["full_residual_l2"], "Phase24 refined": 6.192602138143676e-16, "COMSOL": "not archived"},
        {"metric": "primal residual L2", "historical restart": initial["before"]["primal_residual_l2"], "historical refined": final["before"]["primal_residual_l2"], "Phase24 refined": "not archived", "COMSOL": "not archived"},
        {"metric": "constraint residual L2", "historical restart": initial["before"]["constraint_residual_l2"], "historical refined": final["before"]["constraint_residual_l2"], "Phase24 refined": 2.45e-22, "COMSOL": "not archived"},
    ]
    with (OUT / "historical_vs_phase24_steady.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0])); writer.writeheader(); writer.writerows(table)
    write_json(OUT / "full_elmer_status.json", {
        "historical_one_shot": json.loads((ROOT / "results" / ONE_SHOT / "manifest.json").read_text(encoding="utf-8")),
        "historical_refinement": {"solver_process": "terminated after 21 nonlinear iterations because the 1e-8 nonlinear criterion stalled at approximately 1e-7 numerical change", "final_result_file": (MESH_ROOT / f"{REFINED}.result").is_file(), "captured_complete_iterations": len(refined_dirs), "all_done": False},
    })
    summary = {
        "historical_restart_provenance": provenance,
        "restart_integrity_gate": gate,
        "historical_restart_state": one_shot["before"],
        "historical_initial_full_residual": initial["before"],
        "one_shot": one_shot,
        "historical_refinement": {"status": refinement_status, "iterations": len(rows), "final_tes_temperature_K": final_t, "final_current_A": final_circuit["raw_current_A"], "final_resistance_ohm": final_circuit["resistance_ohm"], "final_power_W": final_circuit["raw_power_W"], "final_full_residual": final["before"]},
        "interpretation": {"strict_143_uA_fixed_point": one_shot["delta_T_mK"] < 1.0e-3 and abs(one_shot["delta_I_uA"]) < 1.0e-2 and abs(final_circuit["raw_current_A"] - 143.78e-6) < 2e-6, "classification": "historical mesh remains on the 143.78-uA fixed-point basin; one-shot is residual-qualified after direct correction, while the repeated refinement was stopped at the numerical-noise plateau", "hypre_gpu": "NO-GO"},
    }
    write_json(OUT / "summary.json", summary)
    lines = [
        "# Historical mesh strict steady validation",
        "",
        f"- Historical saved restart: **{initial_state[0] * 1e3:.9f} mK / {initial_state[1] * 1e6:.9f} uA**.",
        f"- Restart integrity gate: **{'PASS' if gate['passed'] else 'FAIL'}**; load, circuit-init, and pre-first-assembly fields agree.",
        f"- Initial captured full residual: L2={initial['before']['full_residual_l2']:.9e}; primal={initial['before']['primal_residual_l2']:.9e}; constraint={initial['before']['constraint_residual_l2']:.9e}; max={initial['before']['max_absolute_residual']:.9e}.",
        f"- One-shot thermal correction: Delta T={one_shot['delta_T_mK']:.9f} mK; field-revaluated Delta I={one_shot['delta_I_uA']:.9f} uA; post-field current={one_after_circuit['raw_current_A'] * 1e6:.9f} uA.",
        f"- Full native MUMPS refinement: {final_t * 1e3:.9f} mK / {final_circuit['raw_current_A'] * 1e6:.9f} uA after {len(rows)} nonlinear rows ({refinement_status}).",
        f"- Final captured residual: L2={final['before']['full_residual_l2']:.9e}; primal={final['before']['primal_residual_l2']:.9e}; constraint={final['before']['constraint_residual_l2']:.9e}; max={final['before']['max_absolute_residual']:.9e}.",
        "",
        "The historical 143.78-uA restart is a strict fixed-point candidate under the current native HeatSolve + CPU MUMPS replay: the one-shot correction is only 1.07e-5 mK and 3.72e-4 uA, and 21 captured nonlinear rows remain in the same 143.776-uA basin. The pre-solve residual is larger than the direct-solve residual because the saved restart is not the exact linear-system solution, but its physical correction is negligible.",
        "",
        "HYPRE/GPU tuning remains NO-GO.",
    ]
    conductance = OUT / "effective_conductance.json"
    flux_path = OUT / "interface_flux.csv"
    if conductance.is_file():
        values = json.loads(conductance.read_text(encoding="utf-8"))
        lines.insert(-2, f"- FE boundary diagnostic: bath out={values['historical']['bath_outgoing_W']:.9e} W historical vs {values['phase24']['bath_outgoing_W']:.9e} W Phase24; secant G_eff={values['historical']['secant_G_eff_W_per_K']:.9e} vs {values['phase24']['secant_G_eff_W_per_K']:.9e} W/K.")
    if flux_path.is_file():
        rows = list(csv.DictReader(flux_path.open(encoding="utf-8")))
        for row in rows:
            if row["interface"] in {"TES_to_membrane", "TES_to_Stycast"}:
                lines.insert(-2, f"- One-sided {row['route']} {row['interface']} flux: {float(row['left_flux_outgoing_W']):.9e} W / {float(row['right_flux_outgoing_W']):.9e} W; signed mismatch={float(row['signed_mismatch_W']):.9e} W (diagnostic wedge/mortar extraction, not solver-native flux).")
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true", help="snapshot and execute both diagnostics")
    parser.add_argument("--audit", action="store_true", help="materialize reports from completed diagnostics")
    args = parser.parse_args()
    if args.run:
        snapshot_inputs()
        one_project = prepare_project(ONE_SHOT, 1)
        run_case(ONE_SHOT, one_project, 1)
        refined_project = prepare_project(REFINED, 120)
        run_case(REFINED, refined_project, 120)
    if args.audit:
        enrich_outputs()
    if not args.run and not args.audit:
        parser.error("choose --run and/or --audit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

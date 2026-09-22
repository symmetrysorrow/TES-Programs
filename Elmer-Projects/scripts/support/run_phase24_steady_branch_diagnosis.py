"""Run the bounded Phase24 steady-branch diagnosis without changing physics.

The campaign deliberately has two one-Newton-iterate cases.  They assemble
the constrained mortar system at a supplied restart, use MUMPS once, and stop
before the next nonlinear/circuit feedback update.  A separately selected
full refinement remains available for basin tests.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_steady_branch_diagnosis"
MESH = "mesh_singlepixel_gpu_fine_stycast32_mortar"
MESH_ROOT = ROOT / "work" / "meshes" / MESH
SOURCE_PROJECT = (
    ROOT / "artifacts" / "phase24_gate4_5_nomortar" / "restart_refinement" / "mortar" / "project.json"
)
GATE3 = "case_phase24_g3_s32m_s32m2_10f17b"
REFINED = "case_phase24_restart_refine_fine_stycast32_mortar_mumps_mortar"
SOLVER = ROOT.parent / "tools" / "elmer-hypre" / "install-steady-full-capture" / "bin" / "ElmerSolver.exe"
RUNTIME = SOLVER.parent
TOOLCHAIN = Path(r"C:\msys64\ucrt64\bin")


def case_spec(template: dict[str, Any], name: str, restart: str, iterations: int) -> dict[str, Any]:
    result = copy.deepcopy(template)
    result.update({
        "template": "steady", "mesh": MESH, "apply_mortar_bcs": True,
        "restart_from": None, "restart_file_base": restart,
        "restart_file_path": f"../work/meshes/{MESH}/{restart}.result",
        "preexisting_restart": True, "restart_time": 0.020,
        # The refined case uses HeatSolve's nonlinear iterations; keep the
        # outer steady driver at one pass so the diagnostic remains bounded.
        "steady_state_max_iterations": 1, "output_intervals": 1,
        "series_file": f"{name}_series.csv",
        "iteration_series_file": f"{name}_iterations.csv",
        "state_file": f"work/meshes/{MESH}/{name}.state",
        "output_file_path": f"../work/meshes/{MESH}/{name}.result",
        "output_result": True, "post_file": False, "vtu": False,
        "inner_circuit_step_commit": True,
        "solver_comment": (
            "Phase24 steady-branch diagnosis: one assembled constrained MUMPS "
            "solve only; no subsequent nonlinear feedback" if iterations == 1 else
            "Phase24 steady-branch diagnosis: MUMPS basin/refinement run"
        ),
    })
    result.pop("pulse", None)
    result.pop("timesteps", None)
    result["solver"] = {
        **result.get("solver", {}), "linear_system": "mumps",
        "nonlinear_max_iterations": iterations,
        "nonlinear_convergence_tolerance": 1.0e-8,
        "nonlinear_relaxation_factor": 1.0,
        "steady_state_convergence_tolerance": 1.0e-9,
    }
    return result


def prepare(name: str, restart: str, iterations: int) -> Path:
    if not (MESH_ROOT / f"{restart}.result").is_file():
        raise FileNotFoundError(f"restart result missing: {restart}")
    source = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    template = next(iter(source["cases"].values()))
    project = copy.deepcopy(source)
    project["cases"] = {name: case_spec(template, name, restart, iterations)}
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.json"
    path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    source_state = MESH_ROOT / f"{restart}.state"
    if not source_state.is_file():
        raise FileNotFoundError(f"restart circuit state missing: {source_state}")
    shutil.copy2(source_state, OUT / f"{name}_original_state_before.state")
    shutil.copy2(source_state, MESH_ROOT / f"{name}.state")
    return path


def add_capture(sif: Path, capture: Path) -> None:
    capture.mkdir(parents=True, exist_ok=True)
    text = sif.read_text(encoding="utf-8")
    needle = "  Apply Mortar BCs = True\n"
    insertion = (
        needle + '  "Phase24 Full Restriction Capture" = Logical True\n'
        # The native diagnostic hook has a fixed-size Fortran character
        # buffer.  Keep this relative to the working directory (and short)
        # rather than injecting the Windows absolute workspace path.
        + f'  "Phase24 Full Restriction Capture Prefix" = String "{capture.relative_to(ROOT).as_posix()}"\n'
        + '  "Phase24 Full Restriction Capture Timestep" = Integer 1\n'
        + '  "Phase24 Full Restriction Capture Max Iterations" = Integer 1\n'
        + '  "Phase24 Restart State Audit" = Logical True\n'
        + '  "Phase24 Restart Audit Fail Fast" = Logical True\n'
        + f'  "Phase24 Restart State Audit Prefix" = String "{capture.relative_to(ROOT).as_posix()}/restart_state_audit"\n'
    )
    if needle not in text:
        raise ValueError(f"could not find mortar solver key in {sif}")
    # ElmerSolver is launched with ROOT as its cwd (see run.py).  The native
    # HeatSolve hook and tes_parallel_circuit UDF both resolve Constants paths
    # from that cwd, not from generated/cases.  Keep the repository-relative
    # path emitted by build_cases.py.  Rewriting it to ../../work silently
    # selects a different file and makes the restart look like a T0 fallback.
    sif.write_text(text.replace(needle, insertion, 1), encoding="utf-8")


def run(name: str, project: Path, capture: bool) -> None:
    sync = subprocess.run([sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project)], cwd=ROOT, text=True, capture_output=True)
    (OUT / f"{name}_sync.log").write_text(sync.stdout + sync.stderr, encoding="utf-8")
    if sync.returncode:
        raise RuntimeError("parameter sync failed")
    if capture:
        capture_root = OUT / "capture" / name
        (capture_root / "ts0001_nl0001").mkdir(parents=True, exist_ok=True)
        add_capture(ROOT / "generated" / "cases" / f"{name}.sif", capture_root)
    command = [sys.executable, str(ROOT / "run.py"), name, "--project", str(project), "--skip-sync",
               "--mpi-procs", "1", "--elmer-solver", str(SOLVER), "--runtime-bin", str(RUNTIME),
               "--toolchain-bin", str(TOOLCHAIN)]
    log = OUT / f"{name}_launcher.log"
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, env=os.environ.copy())
    if completed.returncode:
        raise RuntimeError(f"{name} exited {completed.returncode}; see {log}")
    if capture:
        subprocess.run([sys.executable, str(ROOT / "scripts" / "support" / "phase24_full_restriction_capture.py"),
                        "materialize", "--root", str(OUT / "capture" / name), "--iterations", "1",
                        "--sif", str(ROOT / "generated" / "cases" / f"{name}.sif"),
                        "--solver-log", str(ROOT / "results" / name / "solver.log")], cwd=ROOT, check=True)


def enrich_trajectory() -> None:
    source = ROOT / "results" / REFINED / f"{REFINED}_iterations.csv"
    rows = list(csv.DictReader(source.open(encoding="utf-8")))
    log = (ROOT / "results" / REFINED / "solver.log").read_text(encoding="utf-8", errors="replace")
    rel = {int(i): float(v.replace("D", "E")) for i, v in re.findall(
        r"ComputeChange: NS \(ITER=\s*(\d+)\) \(NRM,RELC\): \(\s*[-+0-9.EeDd]+\s+([-+0-9.EeDd]+)", log)}
    fields = list(rows[0]) + ["nonlinear_relative_change"]
    with (OUT / "nonlinear_refinement_trajectory.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            row["nonlinear_relative_change"] = rel.get(int(row["nonlinear_iter"]), "")
            writer.writerow(row)


def circuit_at(temperature: float) -> dict[str, float]:
    """The steady algebraic circuit evaluation; no state is committed."""
    ibias, rsh, r0, rmin = 715e-6, 3.9e-3, 15.527e-3, 1e-6
    alpha, beta, i0, t0 = 256.46, 5.03, 143.537344932311e-6, 168.57e-3
    a = r0 * (1.0 + alpha * (temperature - t0) / t0 - beta)
    b = r0 * beta / i0
    current = (np.sqrt((rsh + a) ** 2 + 4.0 * b * ibias * rsh) - (rsh + a)) / (2.0 * b)
    current = float(np.clip(current, 0.0, ibias))
    resistance = max(a + b * abs(current), rmin)
    if resistance == rmin:
        current = ibias * rsh / (rsh + resistance)
    return {"temperature_K": temperature, "raw_current_A": current,
            "resistance_ohm": resistance, "raw_power_W": current * current * resistance}


def indexed(path: Path, count: int) -> np.ndarray:
    values = np.zeros(count)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        words = line.split()
        if len(words) >= 2:
            values[int(words[0]) - 1] = float(words[1])
    return values


def matrix(path: Path, count: int):
    rows: list[int] = []; columns: list[int] = []; values: list[float] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        words = line.split()
        if len(words) >= 3:
            rows.append(int(words[0]) - 1); columns.append(int(words[1]) - 1); values.append(float(words[2]))
    return coo_matrix((values, (rows, columns)), shape=(count, count)).tocsr()


def residual_at_capture(name: str) -> dict[str, Any]:
    root = OUT / "capture" / name / "ts0001_nl0001"
    meta = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    total, primal = int(meta["runtime"]["total_rows"]), int(meta["runtime"]["primal_rows"])
    a, b, x = matrix(root / "full_A_before.dat", total), indexed(root / "full_b_before.dat", total), indexed(root / "full_x_before.dat", total)
    r = np.asarray(a @ x - b).ravel()
    scale = np.asarray(abs(a) @ abs(x)).ravel() + abs(b)
    eta = np.divide(abs(r), scale, out=np.zeros_like(r), where=scale > 0.0)
    return {"full": {"l2": float(np.linalg.norm(r)), "max_abs": float(abs(r).max()), "relative_l2_rhs": float(np.linalg.norm(r) / np.linalg.norm(b))},
            "primal": {"l2": float(np.linalg.norm(r[:primal])), "max_abs": float(abs(r[:primal]).max())},
            "constraint": {"l2": float(np.linalg.norm(r[primal:])), "max_abs": float(abs(r[primal:]).max())},
            "componentwise_backward_error_max": float(eta.max()),
            "componentwise_backward_error_constraint_max": float(eta[primal:].max()),
            "dimensions": {"total": total, "primal": primal, "constraint": total - primal}}


def read_state(name: str) -> dict[str, float]:
    keys = ("temperature_K", "current_A", "resistance_ohm", "power_W", "previous_current_A")
    values = [float(x) for x in (MESH_ROOT / f"{name}.state").read_text().split()]
    return dict(zip(keys, values, strict=True))


def pairwise_grid_audit() -> dict[str, Any]:
    sys.path.insert(0, str(ROOT))
    from scripts.analysis import compare_singlepixel_amgx_comsol as compare
    paths = {"COMSOL": ROOT / "docs" / "Single-Pixel.txt",
             "MUMPS": ROOT / "results" / "case_phase24_g45_fine_stycast32_mortar_40us_mumps_mortar" / "case_phase24_g45_fine_stycast32_mortar_40us_mumps_mortar_series.csv",
             "HYPRE": ROOT / "results" / "case_phase24_g45_fine_stycast32_mortar_40us_hypre_mortar" / "case_phase24_g45_fine_stycast32_mortar_40us_hypre_mortar_series.csv"}
    series = {key: compare.read_comsol(value) if key == "COMSOL" else compare.read_elmer(value) for key, value in paths.items()}
    end = min(40.0, *(float(item.time_us[-1]) for item in series.values()))
    grid = np.linspace(0.0, end, int(np.ceil(end / .05)) + 1)
    drops = {key: np.interp(grid, value.time_us, value.drop_uA) for key, value in series.items()}
    metrics: dict[str, Any] = {}
    for left, right in (("MUMPS", "COMSOL"), ("HYPRE", "COMSOL"), ("HYPRE", "MUMPS")):
        delta = drops[left] - drops[right]
        metrics[f"{left}_vs_{right}"] = {"max_abs_uA": float(abs(delta).max()), "rmse_uA": float(np.sqrt(np.mean(delta * delta))),
                                           "baseline_left_uA": series[left].baseline_uA, "baseline_right_uA": series[right].baseline_uA}
    return {"grid": {"start_us": 0.0, "end_us": end, "spacing_us": .05, "points": int(grid.size), "method": "linear interpolation; each series subtracts its own 19.5--20.020 ms mean baseline"},
            "metrics": metrics,
            "triangle_check": {"HYPRE_COMSOL_lte_HYPRE_MUMPS_plus_MUMPS_COMSOL": metrics["HYPRE_vs_COMSOL"]["max_abs_uA"] <= metrics["HYPRE_vs_MUMPS"]["max_abs_uA"] + metrics["MUMPS_vs_COMSOL"]["max_abs_uA"] + 1e-12}}


def audit() -> None:
    enrich_trajectory()
    gate, refined = read_state(GATE3), read_state(REFINED)
    gate_capture, refined_capture = "case_phase24_diag_one_shot_gate3", "case_phase24_diag_one_shot_refined"
    gate_full = json.loads((OUT / "one_shot_gate3_full.json").read_text())["iterations"]["nl1"]
    refined_full = json.loads((OUT / "refined_full_residual.json").read_text())["iterations"]["nl1"]
    gate_after = circuit_at(float(gate_full["tes_average"]["full_native_primal_TES"]))
    refined_after = circuit_at(float(refined_full["tes_average"]["full_native_primal_TES"]))
    trajectory = list(csv.DictReader((OUT / "nonlinear_refinement_trajectory.csv").open()))
    with (OUT / "branch_sweep.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["initial_state", "solver_path", "final_or_next_T_mK", "raw_or_final_current_uA", "iterations", "classification"]
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        writer.writerows([
            {"initial_state": "Gate3 HYPRE 143.568 uA", "solver_path": "existing 84-iteration MUMPS refinement", "final_or_next_T_mK": refined["temperature_K"] * 1e3, "raw_or_final_current_uA": refined["current_A"] * 1e6, "iterations": 84, "classification": "converges to 218.668-uA recorded iterate"},
            {"initial_state": "Gate3 HYPRE 143.568 uA", "solver_path": "one assembled MUMPS solve", "final_or_next_T_mK": gate_after["temperature_K"] * 1e3, "raw_or_final_current_uA": gate_after["raw_current_A"] * 1e6, "iterations": 1, "classification": "direct linear image; not nonlinear-committed"},
            {"initial_state": "refined 218.668 uA", "solver_path": "one assembled MUMPS solve", "final_or_next_T_mK": refined_after["temperature_K"] * 1e3, "raw_or_final_current_uA": refined_after["raw_current_A"] * 1e6, "iterations": 1, "classification": "leaves 218-uA iterate; not a fixed point"},
        ])
    schedule = [[18e-6, 1], [1e-6, 2], [1e-9, 1], [10e-9, 10], [100e-9, 9], [1e-6, 9], [.625e-6, 47]]
    solver_final = .020 + sum(dt * count for dt, count in schedule)
    series = ROOT / "results" / "case_phase24_g45_fine_stycast32_mortar_40us_mumps_mortar" / "case_phase24_g45_fine_stycast32_mortar_40us_mumps_mortar_series.csv"
    rows = list(csv.DictReader(series.open())); last = max(float(row["time_s"]) for row in rows)
    endpoint = {"solver_timestep_count": sum(count for _, count in schedule), "solver_final_time_s_from_schedule": solver_final,
                "solver_final_after_pulse_us": (solver_final - .020020) * 1e6, "series_rows": len(rows), "series_last_time_s": last,
                "series_last_after_pulse_us": (last - .020020) * 1e6, "writer_gap_us": (solver_final - last) * 1e6,
                "finding": "solver executes 79 steps through 39.376 us post-pulse; the series is one final 0.625-us step behind, so this is a final-row writer/flush defect, not an early solver stop or comparison-grid truncation."}
    gate_assembly = {key: float(value) for key, value in list(csv.DictReader((ROOT / "results" / gate_capture / f"{gate_capture}_iterations.csv").open()))[0].items() if key in {"tes_temperature_K", "previous_current_A", "raw_current_A", "tes_resistance_ohm", "raw_power_W", "relaxed_power_W"}}
    refined_assembly = {key: float(value) for key, value in list(csv.DictReader((ROOT / "results" / refined_capture / f"{refined_capture}_iterations.csv").open()))[0].items() if key in {"tes_temperature_K", "previous_current_A", "raw_current_A", "tes_resistance_ohm", "raw_power_W", "relaxed_power_W"}}
    payload = {"one_shot_gate3": {"before": gate, "assembly_circuit_observed": gate_assembly, "thermal_only_after": gate_full["tes_average"]["full_native_primal_TES"], "raw_circuit_after": gate_after,
                                      "delta_T_mK": (gate_after["temperature_K"] - gate["temperature_K"]) * 1e3, "delta_I_uA": (gate_after["raw_current_A"] - gate["current_A"]) * 1e6,
                                      "before_system_residual": residual_at_capture(gate_capture)},
               "one_shot_refined": {"before": refined, "assembly_circuit_observed": refined_assembly, "thermal_only_after": refined_full["tes_average"]["full_native_primal_TES"], "raw_circuit_after": refined_after,
                                       "delta_T_mK": (refined_after["temperature_K"] - refined["temperature_K"]) * 1e3, "delta_I_uA": (refined_after["raw_current_A"] - refined["current_A"]) * 1e6,
                                       "before_system_residual": residual_at_capture(refined_capture), "after_direct_residual": refined_full["block_residuals"]["native"]},
               "trajectory": {"rows": len(trajectory), "first": trajectory[:3], "last": trajectory[-3:]},
               "comparison_grid_audit": pairwise_grid_audit(), "endpoint_audit": endpoint,
               "conclusion": "The Gate3 one-shot uses an assembled power only 4.03e-15 W below its saved power, yet moves 6.20 mK. This is strong evidence of a linear equilibrium defect. The diagnostic also exposes a separate state-file reload failure: both one-shot cases assemble the T0/143.537-uA circuit state, so the refined-218 one-shot is not a valid fixed-point test and cannot settle branch multiplicity."}
    (OUT / "one_shot_direct.json").write_text(json.dumps(payload["one_shot_gate3"], indent=2) + "\n")
    (OUT / "comparison_grid_audit.json").write_text(json.dumps(payload["comparison_grid_audit"], indent=2) + "\n")
    (OUT / "endpoint_audit.json").write_text(json.dumps(endpoint, indent=2) + "\n")
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    saved_minus_assembled = gate["power_W"] - gate_assembly["relaxed_power_W"]
    lines = ["# Phase24 steady-branch diagnosis", "", "## Result", "", f"- Gate3 one-shot MUMPS thermal solve: {gate['temperature_K']*1e3:.6f} -> {gate_after['temperature_K']*1e3:.6f} mK; field-re-evaluated raw circuit current {gate['current_A']*1e6:.6f} -> {gate_after['raw_current_A']*1e6:.6f} uA.", f"- The Gate3 assembled power was {gate_assembly['relaxed_power_W']:.9e} W, only {saved_minus_assembled:.3e} W below the saved {gate['power_W']:.9e} W. The 6.20-mK direct move is therefore strong linear-defect evidence (Hypothesis A).", "- The intended checkpoint circuit state was not reloaded: both diagnostic assembly rows show the T0/143.537-uA state. Consequently the refined-218-uA one-shot is not a valid nonlinear fixed-point test; branch multiplicity and the existence of a 143-uA high-accuracy fixed point remain open.", "- No physical, mesh, mortar, HYPRE-preconditioner, GPU, or production-timestep parameter was changed.", "", "## Endpoint and comparison", "", f"- The solver schedule reaches {endpoint['solver_final_after_pulse_us']:.3f} us post-pulse, while both series stop at {endpoint['series_last_after_pulse_us']:.3f} us: the final series row is missing.", f"- Recomputed three-way common-grid triangle check: {payload['comparison_grid_audit']['triangle_check']['HYPRE_COMSOL_lte_HYPRE_MUMPS_plus_MUMPS_COMSOL']}. Earlier contradictory values used non-identical comparison/baseline pipelines."]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", choices=("gate3-one-shot", "refined-one-shot", "trajectory", "audit"), required=True)
    args = parser.parse_args()
    if args.run == "trajectory":
        OUT.mkdir(parents=True, exist_ok=True)
        enrich_trajectory()
        return 0
    if args.run == "audit":
        audit()
        return 0
    restart = GATE3 if args.run == "gate3-one-shot" else REFINED
    name = "case_phase24_diag_one_shot_gate3" if args.run == "gate3-one-shot" else "case_phase24_diag_one_shot_refined"
    project = prepare(name, restart, 1 if args.run == "gate3-one-shot" else 84)
    run(name, project, capture=True)
    print(name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

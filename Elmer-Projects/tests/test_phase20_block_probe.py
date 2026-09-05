from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.analysis.evaluate_solver_acceptance import evaluate
from scripts.analysis.summarize_block_schur_probe import summarize, summarize_trace, validate_same_workload
from scripts.prep.prepare_phase20_schur_sweep import build
from scripts.support.build_cases import solver1_block


ROOT = Path(__file__).resolve().parents[1]


def _peer_metrics() -> dict[str, object]:
    return {
        "full_absolute_residual": 1.0e-15,
        "full_relative_residual": 3.0e-11,
        "absolute_constraint_residual": 1.0e-15,
        "backward_error": 1.0e-15,
        "relative_primal_agreement_with_mumps": 2.2e-6,
        "tes_temperature_difference": 1.0e-6,
        "tes_current_difference": 1.0e-6,
        "constrained": True,
    }


def test_phase20_probe_cases_are_peer_cpu_gpu_lower_full_cases() -> None:
    project = json.loads((ROOT / "elmer_project_hypre_gpu_phase19.json").read_text())
    cases = project["cases"]
    names = [
        "case_p20_hypre_block_lower_cpu_probe",
        "case_p20_hypre_block_lower_gpu_probe",
        "case_p20_hypre_block_full_cpu_probe",
        "case_p20_hypre_block_full_gpu_probe",
    ]
    assert all(name in cases for name in names)
    assert len({cases[name]["mesh"] for name in names}) == 1
    assert len({tuple(cases[name]["timesteps"][0]) for name in names}) == 1
    assert len({cases[name]["restart_file_path"] for name in names}) == 1
    assert cases["case_p20_hypre_block_lower_cpu_probe"]["solver"]["block_schur_probe_workload_id"] == cases["case_p20_hypre_block_lower_gpu_probe"]["solver"]["block_schur_probe_workload_id"]
    assert cases["case_p20_hypre_block_full_cpu_probe"]["solver"]["block_schur_probe_workload_id"] == cases["case_p20_hypre_block_full_gpu_probe"]["solver"]["block_schur_probe_workload_id"]
    assert cases["case_p20_hypre_block_lower_cpu_probe"]["solver"]["block_schur_probe_workload_id"] != cases["case_p20_hypre_block_full_cpu_probe"]["solver"]["block_schur_probe_workload_id"]
    for name in names:
        solver = cases[name]["solver"]
        assert solver["linear_system_max_iterations"] == 15
        assert solver["linear_system_abort_not_converged"] is False
        assert solver["block_schur_probe"] is True
        assert solver["matrix_dump"] is False
        assert solver["block_schur_probe_lifecycle"] == "linear solve"
        sif = "\n".join(solver1_block(solver))
        assert "Linear System Save = True" not in sif
        assert "Block Schur Probe = Logical True" in sif
        assert "Block Schur Probe Workload ID" in sif
        assert "HYPRE GPU = " + ("True" if "gpu" in name else "False") in sif
    assert "Block Schur Probe" not in "\n".join(
        solver1_block(cases["case_p19_hypre_block_lower_cpu_time5us"]["solver"])
    )


def test_acceptance_profiles_keep_diagnostic_and_physical_gates_separate() -> None:
    diagnostic = evaluate(_peer_metrics())
    assert diagnostic["status"] == "PASS"
    assert diagnostic["production_ready"] is False
    assert diagnostic["physical_acceptance"]["status"] == "PASS"
    assert diagnostic["relative_residual_is_numerical_floor_warning"] is True

    production = evaluate(_peer_metrics(), profile="production")
    assert production["status"] == "PASS"
    assert production["production_ready"] is True

    missing_physics = evaluate({**_peer_metrics(), "tes_current_difference": None}, profile="production")
    assert missing_physics["status"] == "INCOMPLETE"
    assert "tes_current_difference" in missing_physics["physical_acceptance"]["missing_metrics"]


def test_acceptance_no_mortar_can_exempt_constraint_but_not_physics() -> None:
    metrics = {**_peer_metrics(), "constrained": False, "no_mortar": True}
    result = evaluate(metrics, profile="production")
    assert result["status"] == "PASS"
    assert result["no_mortar_constraint_exemption"] is True
    assert result["constraint_metric_required"] is False


def test_acceptance_nonfinite_or_breakdown_is_hard_failure() -> None:
    result = evaluate({"full_absolute_residual": float("nan"), "breakdown": True})
    assert result["status"] == "FAIL"
    assert "breakdown" in result["implementation_correctness"]["hard_fail_reasons"]


def test_acceptance_rejects_threshold_failures_and_invalid_timing() -> None:
    bad_numerics = evaluate({**_peer_metrics(), "full_relative_residual": 1.0e-3})
    assert bad_numerics["status"] == "FAIL"
    assert bad_numerics["numerical_convergence"]["status"] == "FAIL"

    bad_production = evaluate({**_peer_metrics(), "tes_current_difference": None,
                               "full_absolute_residual": 1.0e-3}, profile="production")
    assert bad_production["status"] == "FAIL"

    bad_timing = evaluate({**_peer_metrics(), "elapsed_wall_seconds": -1.0})
    assert bad_timing["status"] == "FAIL"
    assert bad_timing["performance_readiness"]["status"] == "FAIL"
    assert "elapsed_wall_seconds" in bad_timing["performance_readiness"]["invalid_metrics"]


def _write_probe_csvs(tmp_path: Path, *, missing_outer_residual: bool = False) -> tuple[Path, Path]:
    outer = tmp_path / "outer.csv"
    schur = tmp_path / "schur.csv"
    outer_fields = [
        "probe_version", "workload_id", "outer_iteration", "preconditioner_application",
        "solver_reported_iteration", "initial_residual", "current_residual",
        "solver_reported_residual", "schur_solve_count", "k_actions_total",
        "k_actions_primal_block_solve", "k_actions_matrix_free_schur",
        "k_actions_full_upper_correction", "k_actions_setup_rebuild",
        "elapsed_wall_seconds_cumulative", "elapsed_wall_seconds_per_call",
        "k_apply_seconds_cumulative", "schur_action_seconds_cumulative",
        "reached_tolerance", "hit_maxiter", "breakdown", "nonfinite",
    ]
    with outer.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=outer_fields)
        writer.writeheader()
        writer.writerows([
            {"probe_version": "phase20-v2", "workload_id": "peer-A", "outer_iteration": "1",
             "preconditioner_application": "1", "solver_reported_iteration": "1",
             "initial_residual": "1.0", "current_residual": "1.0", "schur_solve_count": "1",
             "k_actions_total": "6", "k_actions_primal_block_solve": "1",
             "k_actions_matrix_free_schur": "5", "k_actions_full_upper_correction": "0",
             "k_actions_setup_rebuild": "0", "elapsed_wall_seconds_cumulative": "2.0",
             "elapsed_wall_seconds_per_call": "2.0"},
            {"probe_version": "phase20-v2", "workload_id": "peer-A", "outer_iteration": "2",
             "preconditioner_application": "2", "solver_reported_iteration": "2",
             "initial_residual": "0.5", "current_residual": "0.25", "schur_solve_count": "1",
             "k_actions_total": "4", "k_actions_primal_block_solve": "1",
             "k_actions_matrix_free_schur": "2", "k_actions_full_upper_correction": "1",
             "k_actions_setup_rebuild": "0", "elapsed_wall_seconds_cumulative": "5.0",
             "elapsed_wall_seconds_per_call": "3.0"},
        ])
    if missing_outer_residual:
        text = outer.read_text(encoding="utf-8")
        text = text.replace(",1.0,1.0,", ",,,")
        text = text.replace(",0.5,0.25,", ",,,")
        outer.write_text(text, encoding="utf-8")

    schur_fields = [
        "probe_version", "workload_id", "outer_iteration", "schur_solve", "iterations",
        "initial_residual", "final_residual", "reached_tolerance", "hit_maxiter", "breakdown",
        "nonfinite", "k_actions_total", "k_actions_primal_block_solve",
        "k_actions_matrix_free_schur", "k_actions_full_upper_correction", "k_actions_setup_rebuild",
        "elapsed_wall_seconds_per_call", "elapsed_wall_seconds_cumulative",
        "k_apply_seconds_per_call", "k_apply_seconds_cumulative",
        "schur_action_seconds_per_call", "schur_action_seconds_cumulative", "gpu_synchronized",
    ]
    with schur.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=schur_fields)
        writer.writeheader()
        writer.writerows([
            {"probe_version": "phase20-v2", "workload_id": "peer-A", "outer_iteration": "1",
             "schur_solve": "1", "iterations": "5", "initial_residual": "1.0",
             "final_residual": "1e-3", "reached_tolerance": "F", "hit_maxiter": "T",
             "breakdown": "F", "nonfinite": "F", "k_actions_total": "6",
             "k_actions_primal_block_solve": "1", "k_actions_matrix_free_schur": "5",
             "elapsed_wall_seconds_per_call": "1.5", "elapsed_wall_seconds_cumulative": "1.5",
             "k_apply_seconds_per_call": "1.0", "k_apply_seconds_cumulative": "1.0",
             "schur_action_seconds_per_call": "1.2", "schur_action_seconds_cumulative": "1.2"},
            {"probe_version": "phase20-v2", "workload_id": "peer-A", "outer_iteration": "2",
             "schur_solve": "2", "iterations": "2", "initial_residual": "0.5",
             "final_residual": "1e-5", "reached_tolerance": "T", "hit_maxiter": "F",
             "breakdown": "F", "nonfinite": "F", "k_actions_total": "4",
             "k_actions_primal_block_solve": "1", "k_actions_matrix_free_schur": "2",
             "k_actions_full_upper_correction": "1", "elapsed_wall_seconds_per_call": "1.0",
             "elapsed_wall_seconds_cumulative": "2.5", "k_apply_seconds_per_call": "0.8",
             "k_apply_seconds_cumulative": "1.8", "schur_action_seconds_per_call": "0.9",
             "schur_action_seconds_cumulative": "2.1"},
        ])
    return outer, schur


def test_probe_summary_reports_wall_cost_stages_and_log_reductions(tmp_path: Path) -> None:
    outer, schur = _write_probe_csvs(tmp_path)
    result = summarize(outer, schur)
    assert result["status"] == "PASS"
    assert result["schur_solves"] == 2
    assert result["k_actions_total"] == 10.0
    assert result["k_actions_by_stage"]["full_factorization_upper_correction"] == 1.0
    assert result["outer_residual_reduction"] == 0.25
    assert result["actual_outer_iterations"] == 2
    assert result["elapsed_wall_seconds"] == 5.0
    assert result["elapsed_wall_seconds_per_call_sum"] == 5.0
    assert result["log10_reduction"] == math.log10(4.0)
    assert result["log10_reduction_per_k_action"] == math.log10(4.0) / 10.0


def test_probe_summary_rejects_missing_or_invalid_physical_rows(tmp_path: Path) -> None:
    outer, schur = _write_probe_csvs(tmp_path, missing_outer_residual=True)
    result = summarize(outer, schur)
    assert result["status"] == "INCOMPLETE"
    assert result["outer_last_residual"] is None
    assert "outer_residual_missing" in result["incomplete_reasons"]

    text = outer.read_text(encoding="utf-8").replace(",,,", ",-1.0,-1.0,")
    outer.write_text(text, encoding="utf-8")
    failed = summarize(outer, schur)
    assert failed["status"] == "FAIL"
    assert failed["invalid_residual_rows"]


def test_probe_summary_schema_and_header_only_are_incomplete(tmp_path: Path) -> None:
    outer, schur = _write_probe_csvs(tmp_path)
    header = outer.read_text(encoding="utf-8").splitlines()[0] + "\n"
    outer.write_text(header, encoding="utf-8")
    result = summarize(outer, schur)
    assert result["status"] == "INCOMPLETE"
    assert "outer_csv:HEADER_ONLY" in result["incomplete_reasons"]


def test_probe_summary_workload_matching_is_explicit(tmp_path: Path) -> None:
    outer, schur = _write_probe_csvs(tmp_path)
    first = summarize(outer, schur)
    second = summarize(outer, schur)
    assert validate_same_workload(first, second)["status"] == "PASS"
    altered = dict(second)
    altered["workload_signature"] = "different"
    assert validate_same_workload(first, altered)["status"] == "INCOMPLETE"


def test_probe_summary_rejects_missing_schur_and_bad_counter_contract(tmp_path: Path) -> None:
    outer, schur = _write_probe_csvs(tmp_path)
    missing_schur = summarize(outer, tmp_path / "does-not-exist.csv")
    assert missing_schur["status"] == "INCOMPLETE"
    assert "schur_csv:MISSING" in missing_schur["incomplete_reasons"]

    rows = list(csv.DictReader(outer.open(newline="", encoding="utf-8")))
    rows[0]["k_actions_total"] = "7"
    with outer.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    mismatch = summarize(outer, schur)
    assert mismatch["status"] == "FAIL"
    assert mismatch["k_stage_total_mismatch"] is True


def test_probe_summary_rejects_nonmonotonic_time_and_outer_sequence(tmp_path: Path) -> None:
    outer, schur = _write_probe_csvs(tmp_path)
    outer_rows = list(csv.DictReader(outer.open(newline="", encoding="utf-8")))
    outer_rows[1]["elapsed_wall_seconds_cumulative"] = "1.0"
    outer_rows[1]["solver_reported_iteration"] = "3"
    with outer.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(outer_rows[0]))
        writer.writeheader()
        writer.writerows(outer_rows)
    result = summarize(outer, schur)
    assert result["status"] == "FAIL"
    assert result["cumulative_timing_violations"]
    assert result["skipped_outer_iterations"] == [2]
    assert "outer_iteration_sequence_invalid" in result["incomplete_reasons"]


def test_inner_schur_trace_contract_and_stopping_reason(tmp_path: Path) -> None:
    trace = tmp_path / "schur_trace.csv"
    fields = [
        "probe_version", "workload_id", "preconditioner_application", "schur_solve",
        "inner_iteration", "residual_kind", "residual", "initial_residual",
        "residual_over_initial", "configured_tolerance", "rhs_norm",
        "convergence_threshold", "solution_norm", "update_norm", "k_actions_total",
        "stopping_reason",
    ]
    rows = [
        {"probe_version": "phase20-v3", "workload_id": "peer", "preconditioner_application": "7", "schur_solve": "7", "inner_iteration": "0", "residual_kind": "initial", "residual": "1.0", "initial_residual": "1.0", "residual_over_initial": "1.0", "configured_tolerance": "1e-4", "rhs_norm": "1.0", "convergence_threshold": "1e-4", "solution_norm": "0.0", "update_norm": "0.0", "k_actions_total": "1", "stopping_reason": "initial"},
        {"probe_version": "phase20-v3", "workload_id": "peer", "preconditioner_application": "7", "schur_solve": "7", "inner_iteration": "1", "residual_kind": "arnoldi_estimate", "residual": "1e-2", "initial_residual": "1.0", "residual_over_initial": "1e-2", "configured_tolerance": "1e-4", "rhs_norm": "1.0", "convergence_threshold": "1e-4", "solution_norm": "0.5", "update_norm": "0.5", "k_actions_total": "2", "stopping_reason": "iteration"},
        {"probe_version": "phase20-v3", "workload_id": "peer", "preconditioner_application": "7", "schur_solve": "7", "inner_iteration": "2", "residual_kind": "true_residual", "residual": "2e-3", "initial_residual": "1.0", "residual_over_initial": "2e-3", "configured_tolerance": "1e-4", "rhs_norm": "1.0", "convergence_threshold": "1e-4", "solution_norm": "0.6", "update_norm": "0.1", "k_actions_total": "3", "stopping_reason": "maxiter"},
    ]
    with trace.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    result = summarize_trace(trace)
    assert result["status"] == "PASS"
    assert result["stopping_reasons"]["maxiter"] == 1
    assert result["trajectories"][0]["final_residual_over_initial"] == 2e-3

    rows[-1]["convergence_threshold"] = "2e-4"
    with trace.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    assert summarize_trace(trace)["status"] == "FAIL"


def test_default_sweep_is_strictly_one_step_and_dump_free(tmp_path: Path) -> None:
    source = ROOT / "elmer_project_hypre_gpu_phase19.json"
    output = tmp_path / "sweep.json"
    project = build(source, output)
    assert len(project["cases"]) == 16
    assert project["phase20_sweep"]["bounded_outer_max_iterations"] == 15
    assert project["phase20_sweep"]["matrix_dump_default"] is False
    assert project["phase20_sweep"]["first_timestep_exactly_one_step"] is True
    assert all(case["solver"]["block_schur_probe"] for case in project["cases"].values())
    assert all(case["solver"]["matrix_dump"] is False for case in project["cases"].values())
    assert all(case["timesteps"][0][1] == 1 for case in project["cases"].values())
    assert all("Linear System Save = True" not in "\n".join(solver1_block(case["solver"])) for case in project["cases"].values())


def test_native_patch_has_real_timing_lifecycle_and_full_k_contract() -> None:
    patch = (ROOT / "docs/hypre_gpu_phase20_probe.patch").read_text()
    assert "CPU_TIME" not in patch
    assert "SYSTEM_CLOCK" in patch
    assert "Block Schur Probe Lifecycle" in patch
    assert "BlockSchurProbeActiveKey" in patch
    assert "k_actions_primal_block_solve" in patch
    assert "k_actions_matrix_free_schur" in patch
    assert "k_actions_full_upper_correction" in patch
    assert "k_actions_setup_rebuild" in patch
    assert "IOSTAT=io_status" in patch
    assert "Outer solver state is not owned by this hook" in patch
    assert "OuterNo,',0,'" not in patch
    timing_fix = (ROOT / "docs/hypre_phase20_timing_fix.patch").read_text()
    assert "BlockProbeSchurActionSeconds" in timing_fix
    assert "BlockProbeKSchurSeconds" in timing_fix
    contract = (ROOT / "docs/hypre_phase20_inner_contract.patch").read_text()
    assert "convergence_threshold" in contract
    assert "stopping_reason" in contract
    assert "_schur_trace.csv" in contract


def test_exit_code_artifacts_have_one_source_of_truth() -> None:
    artifact_dir = ROOT / "artifacts/hypre_phase20/gpu_correctness_20260905"
    parity_path = artifact_dir / "no_mortar_solution_parity.json"
    cpu_path = artifact_dir / "no_mortar_cpu.exitcode"
    gpu_path = artifact_dir / "no_mortar_gpu.exitcode"
    if not all(path.exists() for path in (parity_path, cpu_path, gpu_path)):
        return
    parity = json.loads(parity_path.read_text())
    assert parity["process_exit_code"]["cpu"] == int(cpu_path.read_text().strip())
    assert parity["process_exit_code"]["gpu"] == int(gpu_path.read_text().strip())


def test_native_patch_is_git_parseable_against_clean_baseline() -> None:
    baseline = Path(r"D:\Github\TES-Programs\tools\elmer-phase20-patch-base")
    if not (baseline / ".git").exists():
        return
    result = subprocess.run(
        ["git", "apply", "--numstat", "--check", str(ROOT / "docs/hypre_gpu_phase20_probe.patch")],
        cwd=baseline,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "BlockSolve.F90" in result.stdout

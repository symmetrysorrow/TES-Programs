from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.analysis.evaluate_solver_acceptance import evaluate
from scripts.analysis.summarize_block_schur_probe import summarize
from scripts.prep.prepare_phase20_schur_sweep import build
from scripts.support.build_cases import solver1_block


ROOT = Path(__file__).resolve().parents[1]


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
    for name in names:
        solver = cases[name]["solver"]
        assert solver["linear_system_max_iterations"] == 15
        assert solver["linear_system_abort_not_converged"] is False
        assert solver["block_schur_probe"] is True
        assert "HYPRE GPU = " + ("True" if "gpu" in name else "False") in "\n".join(
            solver1_block(solver)
        )
    assert "Block Schur Probe" not in "\n".join(solver1_block(cases["case_p19_hypre_block_lower_cpu_time5us"]["solver"]))


def test_acceptance_does_not_hard_fail_oracle_floor_relative_residual() -> None:
    result = evaluate({
        "full_absolute_residual": 1.1384928931822508e-16,
        "full_relative_residual": 3.180826027618793e-11,
        "absolute_constraint_residual": 1.1589146152036896e-24,
        "backward_error": 1.0e-16,
        "relative_primal_agreement_with_mumps": 2.2e-6,
    })
    assert result["status"] == "PASS"
    assert result["production_ready"] is True
    assert result["relative_residual_is_numerical_floor_warning"] is True


def test_acceptance_nonfinite_is_hard_failure() -> None:
    result = evaluate({
        "original_system_absolute_residual": float("nan"),
        "constraint_absolute_residual": 0.0,
        "backward_error": 0.0,
        "breakdown": True,
    })
    assert result["status"] == "FAIL"
    assert "breakdown" in result["hard_fail_reasons"]


def test_probe_summary_reports_cost_metrics(tmp_path: Path) -> None:
    outer = tmp_path / "outer.csv"
    schur = tmp_path / "schur.csv"
    with outer.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["outer_residual", "elapsed_wall_seconds"])
        writer.writeheader()
        writer.writerows([
            {"outer_residual": "1.0", "elapsed_wall_seconds": "2.0"},
            {"outer_residual": "0.25", "elapsed_wall_seconds": "5.0"},
        ])
    with schur.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["iterations", "final_residual", "reached_tolerance", "hit_maxiter", "k_actions"])
        writer.writeheader()
        writer.writerows([
            {"iterations": "5", "final_residual": "1e-3", "reached_tolerance": "F", "hit_maxiter": "T", "k_actions": "6"},
            {"iterations": "2", "final_residual": "1e-5", "reached_tolerance": "T", "hit_maxiter": "F", "k_actions": "3"},
        ])
    result = summarize(outer, schur)
    assert result["schur_solves"] == 2
    assert result["k_actions_total"] == 9
    assert result["outer_residual_reduction"] == 0.25
    assert result["residual_reduction_per_k_action"] == 0.75 / 9


def test_default_sweep_has_small_representative_product(tmp_path: Path) -> None:
    source = ROOT / "elmer_project_hypre_gpu_phase19.json"
    output = tmp_path / "sweep.json"
    project = build(source, output)
    assert len(project["cases"]) == 16
    assert project["phase20_sweep"]["bounded_outer_max_iterations"] == 15
    assert all(case["solver"]["block_schur_probe"] for case in project["cases"].values())

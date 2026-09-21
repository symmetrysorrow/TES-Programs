import json
from pathlib import Path

from scripts.support import run_phase24_hypre_max_iterations_diagnostic as diagnostic


def test_build_project_changes_only_max_iterations(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(diagnostic, "DIAGNOSTIC_DIR", tmp_path)
    monkeypatch.setattr(diagnostic, "PROJECT_PATH", tmp_path / "diagnostic.json")
    project = json.loads(diagnostic.build_project().read_text(encoding="utf-8"))
    base = project["cases"][diagnostic.BASE_CASE]
    case = project["cases"][diagnostic.CASE]

    assert case["timesteps"] == base["timesteps"][:5]
    assert case["output_intervals"] == base["output_intervals"][:5]
    assert case["bdf_order"] == base["bdf_order"] == 2
    assert case["solver"]["linear_system_convergence_tolerance"] == 1.0e-8
    assert case["solver"]["linear_system_max_iterations"] == 10000
    assert case["phase24_hypre_reuse"] is True
    assert case["phase24_preconditioner_lagging"] == "adaptive"


def test_flexgmres_sif_uses_configured_max_iterations(tmp_path: Path) -> None:
    from scripts.support.build_cases import solver1_block

    project_case = {
        "solver": {
            "nonlinear_max_iterations": 1,
            "nonlinear_convergence_tolerance": 1.0e-6,
            "nonlinear_relaxation_factor": 1.0,
            "steady_state_convergence_tolerance": 1.0e-9,
            "linear_system": "iterative_hypre_flexgmres_boomeramg",
            "linear_system_convergence_tolerance": 1.0e-8,
            "linear_system_max_iterations": 10000,
        }
    }
    rendered = "\n".join(solver1_block(project_case["solver"]))
    assert "Linear System Max Iterations = 10000" in rendered


def test_audit_log_detects_stagnating_10000_iteration_failure(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text(
        "PHASE24_HYPRE_SOLVE_FAILURE phase=backend_status method=901 "
        "solve_status=256 iterations=10000 final_relative_residual=1.5e-7 "
        "requested_tolerance=1.0e-8 max_iterations=10000 matrix_epoch=2 "
        "preconditioner_epoch=2 case=0 hypre_error=256 hypre_global_error=256\n",
        encoding="utf-8",
    )
    audit = diagnostic.audit_log(log)
    assert audit["telemetry_records"][0]["iterations"] == 10000
    assert audit["telemetry_records"][0]["final_relative_residual"] == 1.5e-7
    assert audit["telemetry_records"][0]["requested_tolerance"] == 1.0e-8

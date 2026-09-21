import json
from pathlib import Path

from scripts.support import run_phase24_gate3_hypre_steady as gate3


def test_build_project_is_independent_no_mortar_hypre_case(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(gate3, "DIAGNOSTIC_DIR", tmp_path)
    monkeypatch.setattr(gate3, "PROJECT_PATH", tmp_path / "gate3.json")
    project_path = gate3.build_project()
    project = json.loads(project_path.read_text(encoding="utf-8"))
    case = project["cases"][gate3.CASE_NAME]

    assert set(project["cases"]) == {gate3.CASE_NAME}
    assert case["mesh"] == "mesh_singlepixel_prod_v2"
    assert case["apply_mortar_bcs"] is False
    assert case["initial_temperature"] == "T_0"
    assert case["restart_from"] is None
    assert case["restart_file_path"] is None
    assert case["preexisting_restart"] is False
    assert case["solver"]["linear_system"] == "iterative_hypre_flexgmres_boomeramg"
    assert case["solver"]["linear_system_max_iterations"] == 4000
    assert case["solver"]["linear_system_convergence_tolerance"] == 1.0e-10
    assert case["solver"]["hypre_gmres_dimension"] == 100


def test_audit_solver_log_collects_linear_and_nonlinear_residuals(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text(
        "SolveHypre: Required iterations 12 to norm 9.0E-11\n"
        "ComputeChange: NS (ITER=1) (NRM,RELC): (  1.0E+00 2.0E-04 ) :: heat equation\n"
        "ComputeChange: SS (ITER=1) (NRM,RELC): (  2.0E-04 1.0E-09 ) :: heat equation\n"
        "MAIN: *** Elmer Solver: ALL DONE ***\n",
        encoding="utf-8",
    )

    audit = gate3.audit_solver_log(log)

    assert audit["all_done"] is True
    assert audit["hypre_markers"] == 1
    assert audit["linear_residual_final"] == 9.0e-11
    assert audit["nonlinear_residual_final"] == 2.0e-4


def test_gate3_requires_all_physical_and_solver_evidence() -> None:
    case = {
        "initial_temperature": "T_0",
        "apply_mortar_bcs": False,
        "restart_from": None,
        "restart_file_path": None,
        "preexisting_restart": False,
    }
    iteration = {
        "tes_temperature_K": 0.16857,
        "tes_resistance_ohm": 0.015527,
        "raw_power_W": 3.199e-10,
        "raw_current_A": 143.53734493231093e-6,
        "residual_W": 1.0e-15,
    }
    audit = {
        "all_done": True,
        "linear_residuals": [1.0e-10],
        "linear_residual_final": 1.0e-10,
        "nonlinear_residuals": [1.0e-9],
        "nonlinear_residual_final": 1.0e-9,
        "telemetry_records": [],
    }

    result = gate3.evaluate_gate3(
        run={"exit_code": 0},
        audit=audit,
        iteration=iteration,
        series_current_uA=None,
        project_case=case,
    )

    assert result["status"] == "PASS"
    assert result["criteria"]["independent_initial_temperature"] is True
    assert result["criteria"]["tes_circuit_residual_recorded"] is True


def test_gate3_rejects_restart_and_failed_process() -> None:
    result = gate3.evaluate_gate3(
        run={"exit_code": 1},
        audit={"all_done": False, "linear_residuals": [], "telemetry_records": []},
        iteration=None,
        series_current_uA=None,
        project_case={
            "initial_temperature": "T_0",
            "apply_mortar_bcs": False,
            "restart_from": "other_case",
            "restart_file_path": None,
            "preexisting_restart": False,
        },
    )

    assert result["status"] == "FAIL"
    assert result["criteria"]["normal_exit"] is False
    assert result["criteria"]["independent_initial_temperature"] is False

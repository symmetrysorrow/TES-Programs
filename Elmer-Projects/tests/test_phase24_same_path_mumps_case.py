import json
from pathlib import Path

from scripts.support import run_phase24_same_path_mumps_case as diagnostic


def test_build_project_changes_only_linear_backend(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(diagnostic, "DIAGNOSTIC_DIR", tmp_path)
    monkeypatch.setattr(diagnostic, "PROJECT_PATH", tmp_path / "diagnostic.json")
    project_path = diagnostic.build_project()
    project = json.loads(project_path.read_text(encoding="utf-8"))
    case = project["cases"][diagnostic.CASE_NAME]
    base = project["cases"][diagnostic.BASE_CASE]

    assert case["bdf_order"] == base["bdf_order"] == 2
    assert case["timesteps"] == base["timesteps"][:5]
    assert case["output_intervals"] == base["output_intervals"][:5]
    assert case["solver"]["linear_system"] == "mumps"
    assert case["phase24_hypre_reuse"] is True
    assert case["phase24_preconditioner_lagging"] == "adaptive"


def test_audit_solver_log_detects_completed_mumps_run(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text(
        "Linear System Direct Method = MUMPS\n"
        "Relative Change : 1.0E-4\n"
        "Relative Change : 0.0E+00\n"
        "MAIN: *** Elmer Solver: ALL DONE ***\n",
        encoding="utf-8",
    )
    audit = diagnostic.audit_solver_log(log)
    assert audit["all_done"] is True
    assert audit["mumps_markers"] == 1
    assert audit["relative_change_samples"] == 2
    assert audit["relative_change_zero"] == 1

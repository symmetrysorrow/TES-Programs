import json
from pathlib import Path

from scripts.support import run_phase24_lagging_isolation_case as diagnostic


def test_build_project_changes_only_lagging_policy(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(diagnostic, "DIAGNOSTIC_DIR", tmp_path)
    monkeypatch.setattr(diagnostic, "PROJECT_PATH", tmp_path / "diagnostic.json")
    project_path = diagnostic.build_project()
    project = json.loads(project_path.read_text(encoding="utf-8"))
    case = project["cases"][diagnostic.CASE_NAME]

    assert case["bdf_order"] == 2
    assert case["phase24_hypre_reuse"] is True
    assert case["phase24_preconditioner_lagging"] == "disabled"
    assert len(case["timesteps"]) == 5
    assert len(case["output_intervals"]) == 5
    assert case["timesteps"][-1][0] == "100[ns]"


def test_audit_solver_log_counts_zero_iteration_solves(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text(
        "Phase24 HYPRE lifecycle: full_setup=F\n"
        "SolveHypre: Required iterations 0, norm 1.0E-7\n"
        "Relative Change : 0.000000E+00\n"
        "ALL DONE\n",
        encoding="utf-8",
    )
    audit = diagnostic.audit_solver_log(log)
    assert audit["all_done"] is True
    assert audit["hypre_solve_calls"] == 1
    assert audit["hypre_required_iterations_zero"] == 1
    assert audit["relative_change_zero"] == 1

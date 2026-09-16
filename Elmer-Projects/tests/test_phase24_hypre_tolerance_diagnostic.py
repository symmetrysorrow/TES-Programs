import json
from pathlib import Path

from scripts.support import run_phase24_hypre_tolerance_diagnostic as diagnostic


def test_build_project_changes_only_linear_tolerance(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(diagnostic, "DIAGNOSTIC_DIR", tmp_path)
    monkeypatch.setattr(diagnostic, "PROJECT_PATH", tmp_path / "diagnostic.json")
    project_path = diagnostic.build_project()
    project = json.loads(project_path.read_text(encoding="utf-8"))
    base = project["cases"][diagnostic.BASE_CASE]

    for key, options in diagnostic.VARIANTS.items():
        case = project["cases"][options["case"]]
        assert case["bdf_order"] == base["bdf_order"] == 2
        assert case["timesteps"] == base["timesteps"][:5]
        assert case["output_intervals"] == base["output_intervals"][:5]
        assert case["solver"]["linear_system"] == base["solver"]["linear_system"]
        assert case["solver"]["linear_system_convergence_tolerance"] == options["tolerance"]
        assert case["phase24_hypre_reuse"] is True
        assert case["phase24_preconditioner_lagging"] == "adaptive"


def test_audit_log_counts_hypre_iterations(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text(
        "SolveHypre: Required iterations 12 to norm 1.0E-7\n"
        "SolveHypre: Required iterations 0 to norm 1.0E-7\n"
        "Relative Change : 1.0E-4\n"
        "MAIN: *** Elmer Solver: ALL DONE ***\n",
        encoding="utf-8",
    )
    audit = diagnostic.audit_log(log)
    assert audit["all_done"] is True
    assert audit["hypre_solve_calls"] == 2
    assert audit["required_iterations_zero"] == 1
    assert audit["required_iterations_sum"] == 12

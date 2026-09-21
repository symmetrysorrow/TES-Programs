import json
from pathlib import Path

from scripts.support import run_phase24_hypre_tol1e8_diagnostic as diagnostic


def test_build_project_uses_same_short_window_and_only_tol1e8(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(diagnostic, "DIAGNOSTIC_DIR", tmp_path)
    monkeypatch.setattr(diagnostic, "PROJECT_PATH", tmp_path / "diagnostic.json")
    project = json.loads(diagnostic.build_project().read_text(encoding="utf-8"))
    base = project["cases"][diagnostic.BASE_CASE]
    case = project["cases"][diagnostic.CASE]

    assert case["bdf_order"] == base["bdf_order"] == 2
    assert case["timesteps"] == base["timesteps"][:5]
    assert case["output_intervals"] == base["output_intervals"][:5]
    assert case["solver"]["linear_system"] == base["solver"]["linear_system"]
    assert case["solver"]["linear_system_convergence_tolerance"] == 1.0e-8
    assert case["phase24_hypre_reuse"] is True
    assert case["phase24_preconditioner_lagging"] == "adaptive"


def test_audit_log_records_tolerance_failure_and_selector_warning(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text(
        "SolveHypre: Required iterations 42 to norm 1.0E-8\n"
        "ERROR:: HypreSolver: HYPRE failed to satisfy the production linear tolerance\n"
        "NATIVE_XVEC_CAPTURE_CANDIDATE_SEEN_BUT_SELECTOR_MISMATCH\n",
        encoding="utf-8",
    )
    audit = diagnostic.audit_log(log)
    assert audit["hypre_solve_calls"] == 1
    assert audit["required_iterations_max"] == 42
    assert audit["tolerance_failure_marker"] is True
    assert audit["native_xvec_selector_warning"] is True

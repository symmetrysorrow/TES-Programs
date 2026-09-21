from pathlib import Path

from scripts.support import run_phase24_hypre_residual_trace as diagnostic


def test_parse_trace_and_failure_record(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text(
        "PHASE24_HYPRE_RESIDUAL_TRACE_ENABLED method=901 max_iterations=10000 tolerance=1e-8\n"
        "FlexGMRES Iteration = 10 residual = 2.0e-7\n"
        "PHASE24_HYPRE_SOLVE_FAILURE phase=backend_status method=901 iterations=10000 final_relative_residual=1.2e-7\n",
        encoding="utf-8",
    )
    audit = diagnostic.parse_trace(log)
    assert audit["trace_enabled_marker"] is True
    assert audit["trace_line_count"] >= 2
    assert audit["failure_records"][0]["iterations"] == 10000


def test_parse_trace_missing_marker(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text("STOP 1\n", encoding="utf-8")
    audit = diagnostic.parse_trace(log)
    assert audit["trace_missing"] is True
    assert audit["trace_line_count"] == 0

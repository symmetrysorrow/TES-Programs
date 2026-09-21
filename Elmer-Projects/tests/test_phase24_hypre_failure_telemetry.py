from pathlib import Path

from scripts.support import run_phase24_hypre_failure_telemetry as diagnostic


def test_parse_failure_telemetry(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text(
        "PHASE24_HYPRE_SOLVE_FAILURE phase=backend_status method=901 "
        "solve_status=1 iterations=2000 final_relative_residual=1.2e-8 "
        "requested_tolerance=1.0e-8 max_iterations=2000 matrix_epoch=2 "
        "preconditioner_epoch=2 case=1 hypre_error=0 hypre_global_error=0\n"
        "STOP 1\n",
        encoding="utf-8",
    )
    audit = diagnostic.parse_failure_telemetry(log)
    assert audit["telemetry_missing"] is False
    assert audit["stop_1"] == 1
    assert audit["telemetry_records"][0]["iterations"] == 2000
    assert audit["telemetry_records"][0]["final_relative_residual"] == 1.2e-8


def test_parse_missing_telemetry(tmp_path: Path) -> None:
    log = tmp_path / "solver.log"
    log.write_text("ERROR:: HYPRE failed to satisfy the production linear tolerance\n", encoding="utf-8")
    audit = diagnostic.parse_failure_telemetry(log)
    assert audit["telemetry_missing"] is True
    assert audit["telemetry_records"] == []

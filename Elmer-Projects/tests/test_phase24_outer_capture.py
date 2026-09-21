from __future__ import annotations

import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "support" / "phase24_outer_capture.py"
SPEC = importlib.util.spec_from_file_location("outer_capture", SCRIPT)
assert SPEC and SPEC.loader
outer_capture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(outer_capture)


def test_materialize_keeps_only_outer_capture_contract(tmp_path: Path) -> None:
    root = tmp_path / "capture"
    directory = root / "ts0001_nl0001"
    directory.mkdir(parents=True)
    (directory / "outer_before_a.dat").write_text("1 1 2.0\n2 1 -1.0\n", encoding="utf-8")
    (directory / "outer_before_b.dat").write_text("1 3.0\n2 0.0\n", encoding="utf-8")
    (directory / "outer_before_sol.dat").write_text("1 0.1\n", encoding="utf-8")
    (directory / "outer_after_sol.dat").write_text("1 0.2\n", encoding="utf-8")
    sif = tmp_path / "case.sif"
    sif.write_text('Apply Mortar BCs = True\n"Phase24 Vector Assembly" = Logical False\n', encoding="utf-8")
    log = tmp_path / "solver.log"
    log.write_text("PHASE24_OUTER_CAPTURE stage=before, timestep=1, nonlinear_iteration=1, time=0.02, dt=1e-6, bdf_order=1, matrix_epoch=7, rhs_epoch=8\n", encoding="utf-8")

    assert outer_capture.materialize(root, 1, 1, sif, log, None) == 0
    metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["capture_kind"] == "outer_heat_solve"
    assert metadata["matrix"] == {"total_rows": 2, "matrix_records": 2, "primal_rows": 1, "constraint_rows": 1}
    assert metadata["solver"]["matrix_epoch"] == 7
    assert (directory / "A.dat").read_text(encoding="utf-8") == "1 1 2.0\n2 1 -1.0\n"
    assert (directory / "x_after.dat").read_text(encoding="utf-8") == "1 0.2\n"

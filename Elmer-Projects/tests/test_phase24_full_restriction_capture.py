from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "support" / "phase24_full_restriction_capture.py"
SPEC = importlib.util.spec_from_file_location("full_capture", SCRIPT)
assert SPEC and SPEC.loader
full_capture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(full_capture)


def test_materialize_records_runtime_dimensions_hashes_and_provenance(tmp_path: Path) -> None:
    root = tmp_path / "capture"
    directory = root / "ts0001_nl0001"
    directory.mkdir(parents=True)
    (directory / "full_A_before.dat").write_text("1 1 2.0\n1 2 -1.0\n2 1 -1.0\n2 2 2.0\n3 1 1.0\n", encoding="utf-8")
    (directory / "full_A_after.dat").write_text("1 1 2.0\n1 2 -1.0\n2 1 -1.0\n2 2 2.0\n3 1 1.0\n", encoding="utf-8")
    (directory / "full_b_before.dat").write_text("1 3.0\n2 0.0\n3 0.0\n", encoding="utf-8")
    (directory / "full_b_after.dat").write_text("1 3.0\n2 0.0\n3 0.0\n", encoding="utf-8")
    (directory / "full_x_before.dat").write_text("1 0.1\n2 0.2\n3 0.0\n", encoding="utf-8")
    (directory / "full_x_after.dat").write_text("1 0.2\n2 0.1\n3 0.0\n", encoding="utf-8")
    (directory / "full_sizes_before.dat").write_text("3\n5\n2\n1\n", encoding="utf-8")
    (tmp_path / "case.sif").write_text(
        'Apply Mortar BCs = True\nLinear System Direct Method = MUMPS\n', encoding="utf-8"
    )
    (tmp_path / "solver.log").write_text(
        "PHASE24_FULL_RESTRICTION_CAPTURE stage=before, timestep=1, "
        "nonlinear_iteration=1, time=  2.0E-2, dt=  1.0E-6, bdf_order=1, "
        "solver=direct, direct_method=mumps, total_rows=3, primal_rows=2, "
        "constraint_rows=1, matrix_nnz=5, constraint_matrix_rows=1, "
        "eliminate=0, penalty=0, restriction_active=1\n",
        encoding="utf-8",
    )

    assert full_capture.materialize(root, 1, 1, tmp_path / "case.sif", tmp_path / "solver.log") == 0
    metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["capture_kind"] == "solve_with_linear_restriction_full_system"
    assert metadata["runtime"]["total_rows"] == 3
    assert metadata["runtime"]["primal_rows"] == 2
    assert metadata["runtime"]["constraint_rows"] == 1
    assert metadata["provenance"]["physical_time"] == 0.02
    assert metadata["provenance"]["eliminate_linear_constraints"] == 0
    assert metadata["sha256"]["matrix_before"] == metadata["sha256"]["matrix_after"]


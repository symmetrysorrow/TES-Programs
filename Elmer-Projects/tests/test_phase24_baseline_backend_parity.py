import json
from pathlib import Path

from scripts.support import run_phase24_baseline_backend_parity as diagnostic


def test_build_project_creates_one_solve_hypre_and_mumps_variants(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(diagnostic, "DIAGNOSTIC_DIR", tmp_path)
    monkeypatch.setattr(diagnostic, "PROJECT_PATH", tmp_path / "diagnostic.json")
    monkeypatch.setattr(diagnostic, "DUMP_DIR", tmp_path / "dumps")
    project_path = diagnostic.build_project()
    project = json.loads(project_path.read_text(encoding="utf-8"))

    assert set(diagnostic.VARIANTS) == {"hypre", "mumps"}
    for options in diagnostic.VARIANTS.values():
        case = project["cases"][options["case"]]
        assert case["timesteps"] == [["18[us]", 1]]
        assert case["output_intervals"] == [1]
        assert case["solver"]["nonlinear_max_iterations"] == 1
        assert case["solver"]["linear_system"] == options["linear_system"]
        assert case["solver"]["matrix_dump_solution"] is False
        assert case["phase24_hypre_reuse"] is True
        assert case["phase24_preconditioner_lagging"] == "adaptive"


def test_compare_dump_reports_numeric_difference(tmp_path: Path) -> None:
    left = tmp_path / "left.dat"
    right = tmp_path / "right.dat"
    left.write_text("1 1 2.0\n2 2 4.0\n", encoding="utf-8")
    right.write_text("1 1 2.5\n2 2 4.0\n", encoding="utf-8")
    report = diagnostic.compare_dump(left, right, 3)
    assert report["available"] is True
    assert report["union_records"] == 2
    assert report["nonzero_difference_records"] == 1
    assert report["max_absolute_difference"] == 0.5

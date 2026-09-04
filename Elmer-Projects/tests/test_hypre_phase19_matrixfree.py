import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_generated_factorized_cases_have_matrix_free_schur_configuration():
    project = json.loads((ROOT / "elmer_project_hypre_gpu_phase19.json").read_text())
    cases = project["cases"]
    for name in (
        "case_p19_hypre_block_lower_cpu_time5us",
        "case_p19_hypre_block_lower_gpu_time5us",
        "case_p19_hypre_block_full_cpu_time5us",
        "case_p19_hypre_block_full_gpu_time5us",
    ):
        assert cases[name]["solver"]["eliminate_linear_constraints"] is False
        assert cases[name]["solver"]["no_explicit_constrained_matrix"] is True
        assert "block_" in cases[name]["solver"]["linear_system"]


def test_matrix_free_source_has_explicit_sign_and_reuse_guard():
    source = (ROOT.parent / "tools" / "elmer-hypre" / "src" / "fem" / "src" / "BlockSolve.F90").read_text()
    assert "y = d_v - b_k_u" in source
    assert "BlockSchurMatrixFreeSolve" in source
    assert "No Precondition Recompute" in source
    assert "Block Matrix-free Schur" in source

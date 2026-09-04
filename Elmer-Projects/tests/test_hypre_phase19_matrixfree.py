import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.support.build_cases import solver1_block


ROOT = Path(__file__).resolve().parents[1]


FACTOR_CASES = {
    "lower_cpu": ("case_p19_hypre_block_lower_cpu_time5us", False, False),
    "lower_gpu": ("case_p19_hypre_block_lower_gpu_time5us", True, False),
    "full_cpu": ("case_p19_hypre_block_full_cpu_time5us", False, True),
    "full_gpu": ("case_p19_hypre_block_full_gpu_time5us", True, True),
}


def _render_factorized_solver(*, use_gpu: bool, full: bool) -> str:
    linear_system = f"iterative_hypre_block_{'full' if full else 'lower'}"
    if use_gpu:
        linear_system += "_gpu"
    solver = {
        "linear_system": linear_system,
        "nonlinear_max_iterations": 1,
        "nonlinear_convergence_tolerance": 1.0e-3,
        "nonlinear_relaxation_factor": 1.0,
        "steady_state_convergence_tolerance": 1.0e-8,
    }
    return "\n".join(solver1_block(solver, apply_mortar_bcs=True))


def test_generated_factorized_cases_have_matrix_free_schur_configuration():
    project = json.loads((ROOT / "elmer_project_hypre_gpu_phase19.json").read_text())
    cases = project["cases"]
    for _, (name, use_gpu, full) in FACTOR_CASES.items():
        assert cases[name]["solver"]["eliminate_linear_constraints"] is False
        assert cases[name]["solver"]["no_explicit_constrained_matrix"] is True
        assert cases[name]["solver"]["linear_system"] == (
            f"iterative_hypre_block_{'full' if full else 'lower'}"
            + ("_gpu" if use_gpu else "")
        )

        rendered = _render_factorized_solver(use_gpu=use_gpu, full=full)
        assert "Linear System Block Mode = True" in rendered
        assert "Linear System Use Hypre = True" in rendered
        assert "Linear System Iterative Method = FlexGMRES" in rendered
        assert "Block Preconditioner = True" in rendered
        assert "Block Matrix-free Schur = True" in rendered
        assert "Block Nested Primal AMG = True" in rendered
        assert "Block Lower Triangular = True" in rendered
        assert f"Block Full Factorization = {'True' if full else 'False'}" in rendered
        assert f"HYPRE GPU = {'True' if use_gpu else 'False'}" in rendered
        assert "Linear System Direct Method = Umfpack" not in rendered


def test_factorized_solver_rendering_is_not_the_default_umfpack_branch():
    for _, (_, use_gpu, full) in FACTOR_CASES.items():
        rendered = _render_factorized_solver(use_gpu=use_gpu, full=full)
        assert "Linear System Solver = Iterative" in rendered
        assert "Linear System Direct Method" not in rendered


def test_diagnostic_case_is_matrix_free_and_diagnostic_only():
    project = json.loads((ROOT / "elmer_project_hypre_gpu_phase19.json").read_text())
    solver = project["cases"]["case_p19_hypre_block_schur_diag_cpu_time5us"]["solver"]
    assert solver["linear_system"] == "iterative_hypre_block_lower"
    assert solver["block_schur_diagnostic"] is True
    assert solver["block_schur_diagnostic_direct"] is True
    assert solver["block_schur_diagnostic_only"] is True


def test_schur_diagonal_only_generation_path_is_available():
    solver = {
        "linear_system": "iterative_hypre_block_lower",
        "create_schur_matrix_approximation": False,
        "nonlinear_max_iterations": 1,
        "nonlinear_convergence_tolerance": 1.0e-3,
        "nonlinear_relaxation_factor": 1.0,
        "steady_state_convergence_tolerance": 1.0e-8,
    }
    rendered = "\n".join(solver1_block(solver))
    assert "Create Schur Matrix Approximation = False" in rendered


def test_matrix_free_source_has_explicit_sign_and_reuse_guard():
    source_path = ROOT.parent / "tools" / "elmer-hypre" / "src" / "fem" / "src" / "BlockSolve.F90"
    if not source_path.exists():
        pytest.skip("local Elmer source tree is optional integration-test input")
    source = source_path.read_text()
    assert "y = d_v - b_k_u" in source
    assert "BlockSchurMatrixFreeSolve" in source
    assert "No Precondition Recompute" in source
    assert "Block Matrix-free Schur" in source

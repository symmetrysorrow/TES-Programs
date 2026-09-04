"""Create CPU/GPU HYPRE Phase19 smoke cases from the validated MUMPS case."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hybrid_prism_phase19_time5us.json"
OUTPUT = ROOT / "elmer_project_hypre_gpu_phase19.json"
SOURCE_CASE = "case_p19_pulse_time5us"
MUMPS_DIAG_CASE = "case_p19_mumps_matrixdiag_time5us_smoke_7step"
CPU_CASE = "case_p19_hypre_flexgmres_boomeramg_cpu_time5us"
GPU_CASE = "case_p19_hypre_flexgmres_boomeramg_gpu_time5us"
MGR_CPU_CASE = "case_p19_hypre_flexgmres_mgr_cpu_time5us"
MGR_GPU_CASE = "case_p19_hypre_flexgmres_mgr_gpu_time5us"
CONTROL_CASE = "case_p19_hypre_flexgmres_boomeramg_cpu_nomortar_time5us"
GPU_CONTROL_CASE = "case_p19_hypre_flexgmres_boomeramg_gpu_nomortar_time5us"
BLOCK_CPU_CASE = "case_p19_hypre_block_diag_cpu_time5us"
BLOCK_GPU_CASE = "case_p19_hypre_block_diag_gpu_time5us"
LOWER_CPU_CASE = "case_p19_hypre_block_lower_cpu_time5us"
LOWER_GPU_CASE = "case_p19_hypre_block_lower_gpu_time5us"
FULL_CPU_CASE = "case_p19_hypre_block_full_cpu_time5us"
FULL_GPU_CASE = "case_p19_hypre_block_full_gpu_time5us"
SCHUR_DIAG_CASE = "case_p19_hypre_block_schur_diag_cpu_time5us"


def truncate(groups: list[list[object]], count: int) -> list[list[object]]:
    result: list[list[object]] = []
    remaining = count
    for step, interval in groups:
        taken = min(remaining, int(interval))
        if taken:
            result.append([step, taken])
            remaining -= taken
        if not remaining:
            return result
    raise ValueError("source time grid is shorter than the smoke run")


def hypre_case(source: dict, name: str, use_gpu: bool, use_mgr: bool = False) -> dict:
    result = copy.deepcopy(source)
    result.update(
        {
            "restart_from": None,
            "preexisting_restart": True,
            "restart_file_base": "case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight",
            "restart_file_path": "../work/meshes/mesh_hybrid_abs_tet_layers_prism_stack17_abs35r50_noextend/case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight.result",
            "state_file": "work/meshes/mesh_hybrid_abs_tet_layers_prism_stack17_abs35r50_noextend/phase19_steady.state",
            "series_file": f"{name}_series.csv",
            "iteration_series_file": f"{name}_iterations.csv",
            "output_file_path": f"../work/meshes/{source['mesh']}/{name}.result",
            "comparison_time_grid": {
                "mode": "Phase19 HYPRE CUDA/HIP benchmark",
                "purpose": "same mesh, restart, source, circuit, and 5-us rise grid as CPU/MUMPS",
                "reference_case": SOURCE_CASE,
                "linear_solver": "FlexGMRES + BoomerAMG (HYPRE CUDA/HIP)",
            },
        }
    )
    result["solver"] = dict(result["solver"])
    # Keep the mortar rows explicit until a validated MGR marker path (or an
    # exact reduced operator) is available.  Requesting elimination here is
    # misleading for this projector: Elmer reports zero eliminated rows.
    result["solver"]["eliminate_linear_constraints"] = False
    result["solver"]["linear_system"] = (
        "iterative_hypre_flexgmres_mgr_gpu"
        if use_mgr and use_gpu
        else "iterative_hypre_flexgmres_mgr"
        if use_mgr
        else "iterative_hypre_flexgmres_boomeramg_gpu"
        if use_gpu
        else "iterative_hypre_flexgmres_boomeramg"
    )
    result["solver"]["matrix_dump_prefix"] = name
    return result


def hypre_nomortar_case(source: dict, name: str, use_gpu: bool) -> dict:
    """Build an unconstrained SPD control using the same restart and physics.

    The mortar run is an explicit saddle-point system.  This control removes
    mortar assembly so BoomerAMG can be tested independently of saddle-point
    preconditioning; it must not be used as the constrained production case.
    """
    result = hypre_case(source, name, use_gpu=use_gpu)
    result["apply_mortar_bcs"] = False
    result["comparison_time_grid"] = {
        "mode": "Phase19 HYPRE CPU control",
        "purpose": "same restart and physics with mortar constraints disabled",
        "reference_case": SOURCE_CASE,
        "linear_solver": "FlexGMRES + BoomerAMG (HYPRE CPU)",
    }
    result["solver"] = dict(result["solver"])
    result["solver"]["eliminate_linear_constraints"] = False
    return result


def block_diag_case(source: dict, name: str, use_gpu: bool) -> dict:
    """Build the native Elmer block-driver Schur candidate.

    ``No Explicit Constrained Matrix`` is important here: it preserves the
    original ``ConstraintMatrix`` metadata so BlockSolve can construct the
    K/B/B^T blocks.  The ordinary collection-matrix path would flatten the
    saddle system before the block driver sees it.
    """
    result = hypre_case(source, name, use_gpu=use_gpu)
    result["comparison_time_grid"] = {
        "mode": "Phase19 native block Schur candidate",
        "purpose": "Elmer block driver with diagonal Schur approximation",
        "reference_case": SOURCE_CASE,
        "linear_solver": "native BlockSolve + HYPRE BoomerAMG diagonal blocks",
    }
    result["solver"] = dict(result["solver"])
    result["solver"]["linear_system"] = (
        "iterative_hypre_block_diag_gpu" if use_gpu
        else "iterative_hypre_block_diag"
    )
    result["solver"]["eliminate_linear_constraints"] = False
    result["solver"]["no_explicit_constrained_matrix"] = True
    result["solver"]["matrix_dump_prefix"] = name
    return result


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    for name, use_gpu, use_mgr in (
        (CPU_CASE, False, False), (GPU_CASE, True, False),
        (MGR_CPU_CASE, False, True), (MGR_GPU_CASE, True, True),
    ):
        full = hypre_case(project["cases"][SOURCE_CASE], name, use_gpu, use_mgr)
        project["cases"][name] = full

        smoke_name = f"{name}_smoke_7step"
        smoke = copy.deepcopy(full)
        smoke["timesteps"] = truncate(full["timesteps"], 7)
        smoke["output_intervals"] = [999999] * len(smoke["timesteps"])
        smoke["output_intervals"][-1] = 1
        smoke["series_file"] = f"{smoke_name}_series.csv"
        smoke["iteration_series_file"] = f"{smoke_name}_iterations.csv"
        smoke["output_file_path"] = f"../work/meshes/{full['mesh']}/{smoke_name}.result"
        smoke["solver"] = dict(smoke["solver"])
        smoke["solver"]["matrix_dump_prefix"] = smoke_name
        project["cases"][smoke_name] = smoke

        # One-step diagnostic: keep the assembled matrix/RHS at the restart
        # state while avoiding the many nonlinear iterations of the full
        # transient smoke run.  This is the fast CPU/GPU algebra comparison.
        diag_name = f"{name}_smoke_1step"
        diag = copy.deepcopy(full)
        diag["timesteps"] = truncate(full["timesteps"], 1)
        diag["output_intervals"] = [1]
        diag["series_file"] = f"{diag_name}_series.csv"
        diag["iteration_series_file"] = f"{diag_name}_iterations.csv"
        diag["output_file_path"] = f"../work/meshes/{full['mesh']}/{diag_name}.result"
        diag["solver"] = dict(diag["solver"])
        diag["solver"]["matrix_dump_prefix"] = diag_name
        diag["solver"]["nonlinear_max_iterations"] = 1
        diag["solver"]["nonlinear_convergence_tolerance"] = 1e-3
        diag["max_output_level"] = 10
        project["cases"][diag_name] = diag

    # Native block-driver candidate.  Keep this A/B pair one-step only until
    # the Schur approximation passes the original full-system residual gate.
    for block_name, use_gpu in ((BLOCK_CPU_CASE, False), (BLOCK_GPU_CASE, True)):
        block = block_diag_case(project["cases"][SOURCE_CASE], block_name, use_gpu)
        block["timesteps"] = truncate(block["timesteps"], 1)
        block["output_intervals"] = [1]
        block["series_file"] = f"{block_name}_series.csv"
        block["iteration_series_file"] = f"{block_name}_iterations.csv"
        block["output_file_path"] = f"../work/meshes/{block['mesh']}/{block_name}.result"
        block["solver"] = dict(block["solver"])
        block["solver"]["nonlinear_max_iterations"] = 1
        block["solver"]["nonlinear_convergence_tolerance"] = 1e-3
        block["solver"]["matrix_dump_prefix"] = block_name
        block["max_output_level"] = 10
        project["cases"][block_name] = block

    # Deterministic one-rank CPU action diagnostic.  The Elmer block driver
    # writes four matrix-free Schur action pairs; its K solve is switched to
    # SuperLU only for this oracle-quality operator check.
    schur_diag = block_diag_case(project["cases"][SOURCE_CASE], SCHUR_DIAG_CASE, False)
    schur_diag["comparison_time_grid"]["purpose"] = (
        "deterministic matrix-free Schur action versus independent SuperLU"
    )
    schur_diag["timesteps"] = truncate(schur_diag["timesteps"], 1)
    schur_diag["output_intervals"] = [1]
    schur_diag["series_file"] = f"{SCHUR_DIAG_CASE}_series.csv"
    schur_diag["iteration_series_file"] = f"{SCHUR_DIAG_CASE}_iterations.csv"
    schur_diag["output_file_path"] = f"../work/meshes/{schur_diag['mesh']}/{SCHUR_DIAG_CASE}.result"
    schur_diag["solver"] = dict(schur_diag["solver"])
    schur_diag["solver"]["block_schur_diagnostic"] = True
    schur_diag["solver"]["block_schur_diagnostic_direct"] = True
    schur_diag["solver"]["block_schur_diagnostic_prefix"] = SCHUR_DIAG_CASE
    schur_diag["solver"]["nonlinear_max_iterations"] = 1
    schur_diag["solver"]["nonlinear_convergence_tolerance"] = 1e-3
    schur_diag["max_output_level"] = 10
    project["cases"][SCHUR_DIAG_CASE] = schur_diag

    # Compare the diagonal candidate with the two factorized variants.  Keep
    # all three one-step cases on the same mesh/restart and tolerance.
    for block_name, use_gpu, variant in (
        (LOWER_CPU_CASE, False, "lower triangular"),
        (LOWER_GPU_CASE, True, "lower triangular"),
        (FULL_CPU_CASE, False, "full approximate factorization"),
        (FULL_GPU_CASE, True, "full approximate factorization"),
    ):
        block = block_diag_case(project["cases"][SOURCE_CASE], block_name, use_gpu)
        block["comparison_time_grid"]["purpose"] = (
            f"Elmer native {variant} with matrix-free Schur action"
        )
        block["solver"] = dict(block["solver"])
        block["solver"]["linear_system"] = (
            f"iterative_hypre_block_{'full' if variant.startswith('full') else 'lower'}"
            f"_{'gpu' if use_gpu else 'cpu'}"
        )
        block["timesteps"] = truncate(block["timesteps"], 1)
        block["output_intervals"] = [1]
        block["series_file"] = f"{block_name}_series.csv"
        block["iteration_series_file"] = f"{block_name}_iterations.csv"
        block["output_file_path"] = f"../work/meshes/{block['mesh']}/{block_name}.result"
        block["solver"]["nonlinear_max_iterations"] = 1
        block["solver"]["nonlinear_convergence_tolerance"] = 1e-3
        block["solver"]["matrix_dump_prefix"] = block_name
        block["max_output_level"] = 10
        project["cases"][block_name] = block

    # SPD control: validates the HYPRE/BoomerAMG path without explicit mortar
    # rows.  Keep this separate from the constrained CPU/GPU cases.
    for control_name, use_gpu in ((CONTROL_CASE, False), (GPU_CONTROL_CASE, True)):
        control = hypre_nomortar_case(project["cases"][SOURCE_CASE], control_name, use_gpu)
        control["timesteps"] = truncate(control["timesteps"], 1)
        control["output_intervals"] = [1]
        control["series_file"] = f"{control_name}_series.csv"
        control["iteration_series_file"] = f"{control_name}_iterations.csv"
        control["output_file_path"] = f"../work/meshes/{control['mesh']}/{control_name}.result"
        control["solver"] = dict(control["solver"])
        control["solver"]["nonlinear_max_iterations"] = 1
        control["solver"]["nonlinear_convergence_tolerance"] = 1e-3
        control["solver"]["matrix_dump_prefix"] = control_name
        control["max_output_level"] = 10
        project["cases"][control_name] = control

    mumps = hypre_case(project["cases"][SOURCE_CASE], MUMPS_DIAG_CASE, use_gpu=False)
    mumps["solver"]["linear_system"] = "mumps"
    mumps["timesteps"] = truncate(mumps["timesteps"], 7)
    mumps["output_intervals"] = [999999] * len(mumps["timesteps"])
    mumps["output_intervals"][-1] = 1
    project["cases"][MUMPS_DIAG_CASE] = mumps

    mdiag_name = "case_p19_mumps_matrixdiag_time5us_smoke_1step"
    mdiag = copy.deepcopy(mumps)
    mdiag["timesteps"] = truncate(mumps["timesteps"], 1)
    mdiag["output_intervals"] = [1]
    # Persist the direct MUMPS field so block-Schur diagnostics can compare
    # the reconstructed solution against the full saddle-point solve.
    mdiag["output_result"] = True
    mdiag["series_file"] = f"{mdiag_name}_series.csv"
    mdiag["iteration_series_file"] = f"{mdiag_name}_iterations.csv"
    mdiag["output_file_path"] = f"../work/meshes/{mumps['mesh']}/{mdiag_name}.result"
    mdiag["solver"] = dict(mdiag["solver"])
    mdiag["solver"]["matrix_dump_prefix"] = mdiag_name
    mdiag["solver"]["matrix_dump_solution"] = True
    mdiag["solver"]["nonlinear_max_iterations"] = 1
    mdiag["solver"]["nonlinear_convergence_tolerance"] = 1e-3
    mdiag["max_output_level"] = 10
    project["cases"][mdiag_name] = mdiag

    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()

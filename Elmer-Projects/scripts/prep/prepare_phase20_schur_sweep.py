"""Generate a small, reproducible Phase20 Schur parameter sweep manifest.

The default sweep is intentionally representative rather than the full
3x4x2 Cartesian product.  Use ``--all`` when a complete matrix is wanted.
Each generated case is a peer candidate and is bounded to the same first
step/outer iteration budget as the Phase20 probe cases.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "elmer_project_hypre_gpu_phase19.json"
DEFAULT_OUTPUT = ROOT / "elmer_project_hypre_gpu_phase20_sweep.json"
TOLERANCES = (1.0e-2, 1.0e-3, 1.0e-4)
MAX_ITERS = (5, 10, 20, 30)
AMG_CYCLES = (1, 2)


def tag(value: float) -> str:
    return f"{value:.0e}".replace("-0", "-").replace("+", "p")


def selected_points(all_points: bool) -> list[tuple[float, int, int]]:
    if all_points:
        return list(itertools.product(TOLERANCES, MAX_ITERS, AMG_CYCLES))
    return [(1.0e-2, 5, 1), (1.0e-3, 10, 1), (1.0e-4, 20, 1), (1.0e-4, 30, 2)]


def make_case(base: dict, name: str, variant: str, gpu: bool,
              tol: float, max_iters: int, cycles: int,
              *, matrix_dump: bool = False) -> dict:
    case = copy.deepcopy(base)
    case["comparison_time_grid"] = {
        "mode": "Phase20 bounded Schur sweep",
        "purpose": "measure outer convergence and total K cost, not minimum inner residual",
        "peer_candidate": True,
        "variant": variant,
        "gpu": gpu,
        "schur_tolerance": tol,
        "schur_max_iterations": max_iters,
        "k_amg_cycles": cycles,
    }
    case["solver"] = dict(case["solver"])
    case["solver"].update({
        "linear_system": f"iterative_hypre_block_{variant}" + ("_gpu" if gpu else ""),
        "linear_system_max_iterations": 15,
        "linear_system_abort_not_converged": False,
        "block_nested_primal_max_iterations": cycles,
        "block_schur_inner_tolerance": tol,
        "block_schur_max_iterations": max_iters,
        "block_schur_restart": min(30, max_iters),
        "block_schur_probe": True,
        "block_schur_probe_prefix": name,
        "block_schur_probe_workload_id": f"phase20-p19-hypre-block-{variant}-time5us",
        "block_schur_probe_lifecycle": "linear solve",
        "matrix_dump": matrix_dump,
        "matrix_dump_prefix": name if matrix_dump else None,
        "nonlinear_max_iterations": 1,
        "nonlinear_convergence_tolerance": 1.0e-3,
    })
    # Keep exactly one *step*, not merely the first stage.  The first source
    # stage can have a count greater than one in the parent project.
    first_size = case["timesteps"][0][0]
    case["timesteps"] = [[first_size, 1]]
    case["output_intervals"] = [1]
    case["series_file"] = f"{name}_series.csv"
    case["iteration_series_file"] = f"{name}_iterations.csv"
    case["output_file_path"] = f"../work/meshes/{case['mesh']}/{name}.result"
    case["max_output_level"] = 10
    return case


def build(input_path: Path, output_path: Path, all_points: bool = False,
          *, matrix_dump: bool = False) -> dict:
    project = json.loads(input_path.read_text(encoding="utf-8"))
    source = project["cases"]["case_p19_pulse_time5us"]
    cases: dict[str, dict] = {}
    for variant, gpu, point in itertools.product(("lower", "full"), (False, True), selected_points(all_points)):
        tol, max_iters, cycles = point
        name = f"case_p20_sweep_{variant}_{'gpu' if gpu else 'cpu'}_tol{tag(tol)}_m{max_iters}_k{cycles}"
        cases[name] = make_case(
            source, name, variant, gpu, tol, max_iters, cycles,
            matrix_dump=matrix_dump,
        )
    project["cases"] = cases
    project["phase20_sweep"] = {
        "same_mesh_restart_timestep": True,
        "bounded_outer_max_iterations": 15,
        "matrix_dump_default": False,
        "matrix_dump_enabled": matrix_dump,
        "first_timestep_exactly_one_step": True,
        "points": [
            {"schur_tolerance": t, "schur_max_iterations": m, "k_amg_cycles": c}
            for t, m, c in selected_points(all_points)
        ],
        "interpretation": "compare outer reduction, K actions, and wall time; do not optimize inner residual alone",
    }
    output_path.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    return project


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--all", action="store_true", help="generate the complete 3x4x2 sweep")
    parser.add_argument(
        "--matrix-dump", action="store_true",
        help="explicitly enable Linear System Save for diagnostic matrix capture",
    )
    args = parser.parse_args()
    project = build(args.input, args.output, args.all, matrix_dump=args.matrix_dump)
    print(f"{args.output}: {len(project['cases'])} cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

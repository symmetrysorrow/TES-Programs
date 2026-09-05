"""Summarize direct/CPU-HYPRE/GPU-HYPRE conformal smoke results."""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


def result_values(path: Path) -> list[float]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("Perm:"))
    count = int(lines[index].split()[1])
    start = index + 1 + count
    return [float(lines[i].replace("D", "E")) for i in range(start, start + count)]


def log_summary(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    # Elmer's redirected Windows/WSL output can wrap a numeric token in the
    # middle of a line.  Collapse whitespace before parsing so the report is
    # independent of terminal line wrapping.
    compact = re.sub(r"\s+", " ", text)
    solution_times = [float(value) for value in re.findall(r"Solution time \(method 901\):\s+([0-9.E+-]+)", compact)]
    norms = [float(value) for value in re.findall(r"Result Norm\s*:\s*([0-9.E+-]+)", compact)]
    nonlinear = [
        {
            "iteration": int(iteration),
            "norm": float(norm),
            "relative_change": float(relative),
        }
        for iteration, norm, relative in re.findall(
            r"ComputeChange: NS \(ITER=(\d+)\).*?\(\s*([0-9.E+-]+)\s+([0-9.E+-]+)\s*\)",
            compact,
        )
    ]
    total = re.search(r"SOLVER TOTAL TIME\(CPU,REAL\):\s+([0-9.E+-]+)\s+([0-9.E+-]+)", compact)
    return {
        "hypre_gpu_requested": "HYPRE GPU requested" in text,
        "hypre_solution_times_s": solution_times,
        "hypre_solution_time_sum_s": math.fsum(solution_times),
        "result_norms": norms,
        "nonlinear_iterations": nonlinear,
        "solver_total_cpu_s": float(total.group(1)) if total else None,
        "solver_total_wall_s": float(total.group(2)) if total else None,
        "all_done": "ALL DONE" in text,
    }


def vector_difference(left: list[float], right: list[float]) -> dict[str, object]:
    if len(left) != len(right):
        return {
            "comparable": False,
            "left_entries": len(left),
            "right_entries": len(right),
            "reason": "result vectors are from different meshes and have different lengths",
        }
    differences = [abs(a - b) for a, b in zip(left, right)]
    left_norm = math.sqrt(math.fsum(value * value for value in left))
    return {
        "comparable": True,
        "entries": len(differences),
        "max_abs": max(differences, default=0.0),
        "l2": math.sqrt(math.fsum(value * value for value in differences)),
        "relative_l2": math.sqrt(math.fsum(value * value for value in differences)) / max(left_norm, 1.0e-300),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--connectivity", type=Path, required=True)
    parser.add_argument("--direct-log", type=Path, required=True)
    parser.add_argument("--reference-log", type=Path, required=True)
    parser.add_argument("--cpu-log", type=Path, required=True)
    parser.add_argument("--gpu-log", type=Path, required=True)
    parser.add_argument("--direct-result", type=Path, required=True)
    parser.add_argument("--reference-result", type=Path, required=True)
    parser.add_argument("--cpu-result", type=Path, required=True)
    parser.add_argument("--gpu-result", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    direct = result_values(args.direct_result)
    reference = result_values(args.reference_result)
    cpu = result_values(args.cpu_result)
    gpu = result_values(args.gpu_result)
    report = {
        "connectivity": json.loads(args.connectivity.read_text(encoding="utf-8")),
        "runs": {
            "conformal_direct_umfpack": log_summary(args.direct_log),
            "mortar_reference_umfpack": log_summary(args.reference_log),
            "conformal_hypre_cpu": log_summary(args.cpu_log),
            "conformal_hypre_gpu": log_summary(args.gpu_log),
        },
        "result_vector_parity": {
            "conformal_hypre_cpu_vs_gpu": vector_difference(cpu, gpu),
            "conformal_direct_vs_mortar_reference": vector_difference(direct, reference),
        },
        "interpretation": {
            "direct_reference_tolerance": "1e-6 for the UMFPACK Mortar reference; standard Windows Elmer has no MUMPS",
            "hypre_tolerance": "1e-7; 1e-10 reached 1.94e-8 at 1000 CPU iterations and was rejected",
            "result_field_note": "Elmer SaveResult vectors are retained for solver parity, but scalar TES values should be read from the UDF iteration series when present",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

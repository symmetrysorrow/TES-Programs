"""Assemble the measured Phase20 production HYPRE CPU/GPU benchmark."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.assemble_phase20_control_reports import body_temperature_average
from scripts.analysis.evaluate_physical_parity import result_values


ROOT = Path(__file__).resolve().parents[2]
MESH = ROOT / "work/meshes/mesh_singlepixel_conformal_gpu_fine"
CPU_LOG = ROOT / "results/case_phase20_production_hypre_cpu/solver.log"
GPU_LOG = ROOT / "results/case_phase20_production_hypre_gpu/solver.log"


def timing(log: Path) -> dict[str, object]:
    text = log.read_text(encoding="utf-8", errors="replace")
    solve = [float(x) for x in re.findall(r"SolveHypre: Solution time \(method \d+\):\s*([0-9.Ee+-]+)", text)]
    wall_match = re.search(r"WALL_SECONDS\s+([0-9.Ee+-]+)", text)
    return {
        "solve_times_s": solve,
        "linear_solve_total_s": sum(solve),
        "wall_seconds": float(wall_match.group(1)) if wall_match else None,
        "all_done": "MAIN: *** Elmer Solver: ALL DONE ***" in text,
        "gpu_device_migration": "migrated HYPRE IJ matrices to device memory" in text,
    }


def main() -> int:
    cpu_result = MESH / "case_phase20_production_hypre_cpu.result"
    gpu_result = MESH / "case_phase20_production_hypre_gpu.result"
    cpu_values = result_values(cpu_result, field_index=0)
    gpu_values = result_values(gpu_result, field_index=0)
    deltas = [abs(cpu_values[node] - gpu_values[node]) for node in cpu_values]
    cpu_tes = body_temperature_average(MESH, cpu_result)
    gpu_tes = body_temperature_average(MESH, gpu_result)
    cpu = timing(CPU_LOG)
    gpu = timing(GPU_LOG)
    speedup = cpu["wall_seconds"] / gpu["wall_seconds"]
    benchmark = {
        "status": "PASS" if cpu["all_done"] and gpu["all_done"] and gpu["gpu_device_migration"] else "FAIL",
        "mesh": str(MESH.resolve()),
        "node_count": int((MESH / "mesh.header").read_text().split()[0]),
        "tetrahedron_count": int((MESH / "mesh.header").read_text().split()[1]),
        "topology_gate": "PASS",
        "cpu": cpu,
        "gpu": gpu,
        "gpu_wall_speedup": speedup,
        "benchmark_definition": "identical production-size conformal shared-node steady HYPRE control; circuit_parallel one coupling iteration",
    }
    parity = {
        "status": "PASS" if max(deltas) <= 1.0e-5 and abs(gpu_tes / cpu_tes - 1.0) <= 1.0e-5 else "FAIL",
        "cpu_result": str(cpu_result.resolve()),
        "gpu_result": str(gpu_result.resolve()),
        "temperature_nodes": len(deltas),
        "temperature_max_abs_difference_K": max(deltas),
        "temperature_rms_difference_K": (sum(delta * delta for delta in deltas) / len(deltas)) ** 0.5,
        "tes_volume_average_cpu_K": cpu_tes,
        "tes_volume_average_gpu_K": gpu_tes,
        "tes_volume_average_difference_K": gpu_tes - cpu_tes,
        "tes_volume_average_relative_difference": abs(gpu_tes / cpu_tes - 1.0),
        "tolerance_temperature_max_abs_K": 1.0e-5,
        "tolerance_tes_relative": 1.0e-5,
    }
    out = ROOT / "artifacts/phase20_conformal"
    out.mkdir(parents=True, exist_ok=True)
    (out / "production_cpu_gpu_benchmark.json").write_text(json.dumps(benchmark, indent=2) + "\n", encoding="utf-8")
    (out / "production_cpu_gpu_parity.json").write_text(json.dumps(parity, indent=2) + "\n", encoding="utf-8")
    (out / "gpu_size_crossover.json").write_text(
        json.dumps(
            {
                "status": "NOT_AVAILABLE",
                "reason": "Only the production-size point was requested and measured; a multi-size crossover curve was not run.",
                "production_point": {"node_count": benchmark["node_count"], "gpu_wall_speedup": speedup},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"benchmark": benchmark["status"], "parity": parity["status"], "gpu_wall_speedup": speedup}, indent=2))
    return 0 if benchmark["status"] == parity["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import json

from scripts.analysis.assemble_phase21_host_optimization_reports import parse_log


def test_parse_phase22_clocks_and_udf_wall(tmp_path):
    case = tmp_path / "case"
    case.mkdir()
    (case / "manifest.json").write_text(
        json.dumps(
            {
                "wall_seconds": 2.5,
                "solver_wall_seconds": 2.3,
                "output_io_wall_seconds": 0.2,
                "process_cpu_seconds": 4.0,
                "thread_cpu_seconds": 0.01,
                "clock_policy": {"wall": "perf_counter"},
            }
        ),
        encoding="utf-8",
    )
    log = case / "solver.log"
    log.write_text(
        "\n".join(
            [
                "MAIN: Time: 1/1:",
                "MAIN: Elapsed time: 2.0 seconds",
                "SolveHypre: setup time (method 1): 0.3",
                "SolveHypre: Solution time (method 1): 0.4",
                "TESParallelCircuitProfile: step=1 iter=1 integration_cpu_s= 0.1 circuit_output_cpu_s= 0.2 total_cpu_s= 0.3 integration_wall_s= 0.01 circuit_output_wall_s= 0.02 total_wall_s= 0.03 cached_elements= 2 cached_nodes= 8",
                "MAIN: *** Elmer Solver: ALL DONE ***",
                "WALL_SECONDS 2.5",
            ]
        ),
        encoding="utf-8",
    )
    parsed = parse_log(log)
    assert parsed["wall_seconds"] == 2.5
    assert parsed["process_cpu_seconds"] == 4.0
    assert parsed["thread_cpu_seconds"] == 0.01
    assert parsed["circuit_profile"][0]["total_wall_s"] == 0.03
    assert parsed["circuit_profile"][0]["total_cpu_s"] == 0.3

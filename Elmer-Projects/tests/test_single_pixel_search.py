from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "search" / "single_pixel_search.py"
SPEC = importlib.util.spec_from_file_location("single_pixel_search", MODULE_PATH)
assert SPEC and SPEC.loader
SEARCH = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SEARCH
SPEC.loader.exec_module(SEARCH)


def config() -> dict:
    return json.loads((ROOT / "single_pixel_search_config.json").read_text(encoding="utf-8"))


def test_candidate_count_and_ids() -> None:
    candidates = SEARCH.candidates_from_config(config())
    assert len(candidates) == 11
    assert candidates[0].candidate_id == "baseline"
    assert {candidate.candidate_id for candidate in candidates} >= {
        "pb_cp_m20",
        "pb_cp_p20",
        "stycast_k_m20",
        "pb_k_p20",
    }


def test_mutation_scales_only_selected_expression() -> None:
    cfg = config()
    base = json.loads((ROOT / cfg["base_project"]).read_text(encoding="utf-8"))
    candidate = SEARCH.Candidate("pb_cp_p20", {"pb_cp": 1.2}, False)
    project, metadata = SEARCH.mutate_project(base, cfg, candidate, 1.0, 0)
    assert project["materials"]["Pb"]["cp"]["expression"].endswith("*1.2")
    assert project["materials"]["TES"]["cp"]["expression"] == base["materials"]["TES"]["cp"]["expression"]
    assert set(project["cases"]) == {metadata["steady_case"], metadata["pulse_case"]}
    assert project["cases"][metadata["pulse_case"]]["restart_from"] == metadata["steady_case"]
    steady = project["cases"][metadata["steady_case"]]
    pulse = project["cases"][metadata["pulse_case"]]
    assert steady["output_file_path"].startswith("../work/meshes/")
    assert pulse["restart_file_path"] == steady["output_file_path"]


def test_multivariate_candidates_are_reproducible_and_bounded() -> None:
    cfg = config()
    first = SEARCH.multivariate_candidates_from_config(cfg, sample_count=8, seed=17)
    second = SEARCH.multivariate_candidates_from_config(cfg, sample_count=8, seed=17)
    assert first == second
    assert len(first) == 8
    assert len({candidate.candidate_id for candidate in first}) == 8
    for candidate in first:
        assert candidate.steady_sensitive
        assert set(candidate.factors) == set(cfg["multivariate_search"]["variables"])
        for variable_name, factor in candidate.factors.items():
            lower, upper = cfg["variables"][variable_name]["search_range"]
            assert lower < factor < upper


def test_multivariate_mutation_scales_all_selected_expressions() -> None:
    cfg = config()
    base = json.loads((ROOT / cfg["base_project"]).read_text(encoding="utf-8"))
    candidate = SEARCH.Candidate(
        "synthetic_multi", {"pb_cp": 1.1, "tes_cp": 0.9, "stycast_k": 1.3}, True
    )
    project, _ = SEARCH.mutate_project(base, cfg, candidate, 0.8, 3)
    assert project["materials"]["Pb"]["cp"]["expression"].endswith("*1.1")
    assert project["materials"]["TES"]["cp"]["expression"].endswith("*0.9")
    assert project["materials"]["Stycast"]["k"]["expression"].endswith("*1.3")
    assert project["parameter_expressions"]["G0"].endswith("*0.8")


def test_pulse_metrics_median_baseline_and_auto_polarity() -> None:
    time_ms = np.linspace(0.0, 2.0, 21)
    signal = np.full_like(time_ms, 10.0)
    signal[2] = 1000.0
    signal[time_ms >= 1.0] += np.linspace(0.0, 3.0, np.count_nonzero(time_ms >= 1.0))
    metrics = SEARCH.pulse_metrics(
        time_ms,
        signal,
        1.0,
        (0.0, 0.9),
        response_direction="auto",
        baseline_statistic="median",
    )
    assert metrics["baseline_signal"] == 10.0
    assert metrics["response_direction"] == "rise"
    assert metrics["peak_signal"] == 3.0


def test_score_recovers_known_time_shift(tmp_path: Path) -> None:
    cfg = config()
    cfg["reference"] = str(tmp_path / "reference.txt")
    cfg["reference_pulse_start_ms"] = 1.0
    cfg["baseline_window_ms"] = [0.0, 1.0]
    cfg["comparison_window_ms"] = [0.0, 0.5]
    cfg["score"]["max_shift_ms"] = 0.08
    cfg["score"]["shift_step_ms"] = 0.001
    cfg["score"]["grid_points"] = 600

    time_ms = np.linspace(0.0, 2.0, 2001)
    relative = time_ms - 1.0
    response = np.where(relative >= 0, (1 - np.exp(-relative / 0.04)) * np.exp(-relative / 0.5), 0.0)
    reference = np.column_stack(
        [time_ms, np.zeros_like(time_ms), np.zeros_like(time_ms), np.zeros_like(time_ms), 100.0 - response, np.zeros_like(time_ms)]
    )
    np.savetxt(tmp_path / "reference.txt", reference)

    shift_ms = 0.035
    sim_time_ms = time_ms
    shifted_relative = sim_time_ms - 1.0 - shift_ms
    sim_response = np.where(
        shifted_relative >= 0,
        (1 - np.exp(-shifted_relative / 0.04)) * np.exp(-shifted_relative / 0.5),
        0.0,
    )
    series = tmp_path / "series.csv"
    with series.open("w", encoding="utf-8") as file:
        file.write("time_s,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W\n")
        for time, value in zip(sim_time_ms, sim_response):
            file.write(f"{time/1e3},0.16,{(100.0-value)/1e6},0.01,1e-10\n")

    trace = tmp_path / "aligned.csv"
    score = SEARCH.score_series(
        cfg,
        series,
        SEARCH.Candidate("synthetic", {}, False),
        1.0,
        trace_path=trace,
    )
    # The scoring convention reports the shift applied to the simulated
    # trace.  A synthetically delayed simulation must therefore be moved
    # earlier by the same amount.
    assert abs(score["best_shift_ms"] + shift_ms) <= 0.002
    assert score["normalized_rmse"] < 0.01
    assert score["waveform_objective"] < 0.01
    assert set(score["region_rmse"]) >= {"rise", "peak"}
    assert trace.exists()
    assert "simulation_normalized" in trace.read_text(encoding="utf-8").splitlines()[0]


def test_analyze_sensitivity_writes_summary_and_correlation(tmp_path: Path) -> None:
    cfg = config()
    cfg["output_dir"] = str(tmp_path / "search")
    cfg["variables"] = {
        "thermal_x": {
            "path": ["materials", "Pb", "cp", "expression"],
            "factors": [0.8, 1.2],
            "search_range": [0.5, 2.0],
            "steady_sensitive": False,
        }
    }
    time_ms = np.linspace(0.0, 0.5, 101)
    reference = np.exp(-time_ms / 0.2)
    regions = np.where(time_ms < 0.1, "rise", np.where(time_ms < 0.2, "peak", "decay"))
    for factor, suffix, scale, objective in (
        (0.8, "m20", 0.9, 0.2),
        (1.2, "p20", 1.1, 0.1),
    ):
        candidate_id = f"thermal_x_{suffix}"
        directory = SEARCH.candidate_dir(cfg, candidate_id)
        trace_path = directory / "aligned_waveform.csv"
        SEARCH.write_aligned_trace(
            trace_path,
            time_ms,
            reference,
            reference * scale,
            regions,
        )
        SEARCH.write_json(
            directory / "score.json",
            {
                "candidate_id": candidate_id,
                "factors": {"thermal_x": factor},
                "aligned_trace": str(trace_path),
                "objective": objective,
            },
        )

    paths = SEARCH.analyze_sensitivity(cfg)
    assert paths["summary"].exists()
    assert paths["correlation"].exists()
    summary = paths["summary"].read_text(encoding="utf-8")
    assert "thermal_x" in summary
    assert "waveform_sensitivity_rms" in summary

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "search" / "post_multivariate_search.py"
SPEC = importlib.util.spec_from_file_location("post_multivariate_search", MODULE_PATH)
assert SPEC and SPEC.loader
SEARCH = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SEARCH
SPEC.loader.exec_module(SEARCH)


def config() -> dict:
    return json.loads(
        (ROOT / "post_multivariate_search_config.json").read_text(encoding="utf-8")
    )


def test_candidates_are_reproducible_and_bounded() -> None:
    cfg = config()
    first = SEARCH.candidates_from_config(cfg, sample_count=8, seed=17)
    second = SEARCH.candidates_from_config(cfg, sample_count=8, seed=17)
    assert first == second
    assert len(first) == 8
    assert len({candidate.candidate_id for candidate in first}) == 8
    for candidate in first:
        assert set(candidate.factors) == set(cfg["multivariate_search"]["variables"])
        for name, factor in candidate.factors.items():
            lower, upper = cfg["variables"][name]["search_range"]
            assert lower < factor < upper


def test_mutation_builds_dual_position_cases() -> None:
    cfg = config()
    base = json.loads((ROOT / cfg["base_project"]).read_text(encoding="utf-8"))
    candidate = SEARCH.Candidate(
        "synthetic", {"pb_cp": 1.1, "stycast_k": 0.8}
    )
    project, metadata = SEARCH.mutate_project(base, cfg, candidate, 1.2, 3)
    assert project["materials"]["Pb"]["cp"]["expression"].endswith("*1.1")
    assert project["materials"]["Stycast"]["k"]["expression"].endswith("*0.8")
    assert project["parameter_expressions"]["G0"].endswith("*1.2")
    assert set(project["cases"]) == {
        metadata["steady_case"],
        *(value["case"] for value in metadata["pulse_cases"].values()),
    }
    steady = project["cases"][metadata["steady_case"]]
    assert steady["mesh"] == "mesh_dual_20mm_localrefine"
    steady_result_path = (
        f"../work/meshes/{metadata['mesh_dir']}/{metadata['steady_case']}.result"
    )
    assert steady["output_file_path"] == steady_result_path
    for target in cfg["targets"]:
        pulse_name = metadata["pulse_cases"][target["name"]]["case"]
        pulse = project["cases"][pulse_name]
        assert pulse["restart_from"] == metadata["steady_case"]
        assert pulse["restart_file_path"] == steady_result_path
        assert pulse["vtu"] is False
        assert pulse["timesteps"] == cfg["cases"]["pulse_overrides"]["timesteps"]


def test_reference_key_selection_tracks_channel_fraction_quantiles() -> None:
    selected = list(range(1, 11))
    left = {key: float(key) for key in selected}
    right = {key: float(11 - key) for key in selected}
    targets = [
        {"name": "left", "selector": {"fraction_quantile": 0.9, "count": 2}},
        {"name": "center", "selector": {"fraction_quantile": 0.5, "count": 2}},
        {"name": "right", "selector": {"fraction_quantile": 0.1, "count": 2}},
    ]
    result = SEARCH.select_reference_keys(
        selected, left, right, targets, left_gain=5.5, right_gain=5.5
    )
    assert result["left"]["mean_left_fraction"] > result["center"]["mean_left_fraction"]
    assert result["center"]["mean_left_fraction"] > result["right"]["mean_left_fraction"]
    assert len(result["left"]["keys"]) == 2


def test_pair_score_recovers_shared_shift_and_peak_share() -> None:
    score = {
        "region_weights": {"rise": 2.0, "peak": 3.0, "decay": 2.0, "tail": 0.5}
    }
    grid = np.linspace(0.0, 20.0, 2001)
    filter_grid = np.linspace(-2.0, 22.0, 2401)

    def pulse(time: np.ndarray) -> np.ndarray:
        return np.where(
            time >= 0.0,
            (1.0 - np.exp(-time / 0.8)) * np.exp(-time / 7.0),
            0.0,
        )

    shift = 0.35
    reference_left = 0.7 * pulse(grid)
    reference_right = 0.3 * pulse(grid)
    simulation_left = 0.7 * pulse(filter_grid - shift)
    simulation_right = 0.3 * pulse(filter_grid - shift)
    envelope = reference_left + reference_right
    regions = SEARCH.waveform_regions(
        grid, envelope / np.max(envelope), float(grid[np.argmax(envelope)]),
        {"peak_half_width_ms": 0.5, "tail_threshold": 0.2},
    )
    data = {
        "comparison_grid": grid,
        "filter_grid": filter_grid,
        "reference_left": reference_left,
        "reference_right": reference_right,
        "simulation_left": simulation_left,
        "simulation_right": simulation_right,
        "regions": regions,
        "reference_left_fraction": 0.7,
        "simulation_left_fraction": 0.7,
    }
    metrics = SEARCH.score_target_shift(data, -shift, score)
    assert metrics["waveform_objective"] < 1e-3
    assert metrics["peak_share_error"] < 1e-12


def test_ch0_extrema_score_ignores_unselected_channel() -> None:
    score = {
        "region_weights": {"rise": 2.0, "peak": 3.0, "decay": 2.0, "tail": 0.5}
    }
    grid = np.linspace(0.0, 20.0, 2001)
    filter_grid = np.linspace(-2.0, 22.0, 2401)

    def pulse(time: np.ndarray) -> np.ndarray:
        return np.where(
            time >= 0.0,
            (1.0 - np.exp(-time / 0.8)) * np.exp(-time / 7.0),
            0.0,
        )

    shift = 0.35
    reference = pulse(grid)
    simulation = pulse(filter_grid - shift)
    regions = SEARCH.waveform_regions(
        grid,
        reference / np.max(reference),
        float(grid[np.argmax(reference)]),
        {"peak_half_width_ms": 0.5, "tail_threshold": 0.2},
    )
    data = {
        "comparison_mode": "ch0_extrema",
        "comparison_grid": grid,
        "filter_grid": filter_grid,
        "reference_selected": reference,
        "simulation_selected": simulation,
        "simulation_left": np.full_like(filter_grid, 999.0),
        "simulation_right": np.full_like(filter_grid, -999.0),
        "regions": regions,
    }
    metrics = SEARCH.score_target_shift(data, -shift, score)
    assert metrics["waveform_objective"] < 1e-3
    assert metrics["selected_rmse"] < 1e-3
    assert "peak_share_error" not in metrics


def test_ch0_extrema_peak_ratio_uses_high_and_low_only() -> None:
    metrics = SEARCH.relative_peak_ratio_metrics(
        [
            {
                "name": "high",
                "reference_selected_peak": 4.0,
                "simulation_selected_peak": 6.0,
            },
            {
                "name": "low",
                "reference_selected_peak": 2.0,
                "simulation_selected_peak": 3.0,
            },
        ],
        {"peak_ratio_targets": ["high", "low"]},
    )
    assert metrics["reference_peak_ratio"] == 2.0
    assert metrics["simulation_peak_ratio"] == 2.0
    assert metrics["peak_ratio_error"] == 0.0


def test_ch0_extrema_artifacts_do_not_replace_paired_score() -> None:
    cfg = config()
    assert SEARCH.score_filename(cfg) == "score_ch0_extrema.json"
    assert SEARCH.leaderboard_filename(cfg) == "leaderboard_ch0_extrema.csv"


def test_side_series_file_matches_dual_builder_convention() -> None:
    assert (
        SEARCH.side_series_file("postsearch_case_series.csv", "L")
        == "postsearch_case_L_series.csv"
    )

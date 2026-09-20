from __future__ import annotations

import importlib.util
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "support"
    / "run_phase24_restart_continuity_campaign.py"
)
SPEC = importlib.util.spec_from_file_location("phase24_restart_continuity_campaign", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
campaign = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(campaign)


def test_read_state_file_parses_five_value_checkpoint(tmp_path: Path) -> None:
    state = tmp_path / "steady.state"
    state.write_text(
        "1.6856910000000000E-01 1.4356758900000000E-04 "
        "1.0000000000000000E-02 2.0000000000000000E-10 "
        "1.4356758900000000E-04\n",
        encoding="utf-8",
    )

    parsed = campaign.read_state_file(state)

    assert parsed["exists"] is True
    assert parsed["temperature_K"] == 0.1685691
    assert abs(parsed["current_uA"] - 143.567589) < 1.0e-9
    assert parsed["sha256"]


def test_set_short_variant_is_pulse_off_and_non_destructive() -> None:
    base = {
        "template": "pulse",
        "state_file": "work/meshes/example/original.state",
        "timesteps": [["18[us]", 1], ["1[us]", 2]],
        "output_intervals": [1, 1],
        "pulse": {
            "energy": "1332[keV]",
            "start": "20.02[ms]",
            "duration": "1[ns]",
        },
        "solver": {
            "linear_system": "iterative_hypre_flexgmres_boomeramg",
            "nonlinear_max_iterations": 50,
        },
    }

    variant = campaign.set_short_variant(
        base,
        "case_restartdiag",
        backend="mumps",
        bdf_order=1,
        steps=5,
        state_file="artifacts/phase24_restart_continuity_campaign/state_snapshots/gate3.state",
        apply_mortar=True,
        dump_matrix=False,
    )

    assert base["pulse"]["energy"] == "1332[keV]"
    assert base["state_file"] == "work/meshes/example/original.state"
    assert variant["pulse"]["energy"] == 0.0
    assert variant["timesteps"] == [["18[us]", 5]]
    assert variant["output_intervals"] == [1]
    assert variant["bdf_order"] == 1
    assert variant["apply_mortar_bcs"] is True
    assert variant["solver"]["linear_system"] == "mumps"
    assert variant["state_file"].endswith("gate3.state")


def test_snapshot_state_uses_private_artifact_copy(tmp_path: Path, monkeypatch) -> None:
    source = tmp_path / "source.state"
    source.write_text("1 2 3 4 5\n", encoding="utf-8")
    snapshot_dir = tmp_path / "snapshots"
    monkeypatch.setattr(campaign, "STATE_SNAPSHOT_DIR", snapshot_dir)
    monkeypatch.setattr(campaign, "ROOT", tmp_path)

    stored = campaign.snapshot_state(source, "gate3")
    copied = tmp_path / stored

    assert copied == snapshot_dir / "gate3.state"
    assert copied.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")
    copied.write_text("changed\n", encoding="utf-8")
    assert source.read_text(encoding="utf-8") == "1 2 3 4 5\n"


def test_first_jump_location_prefers_pre_solver_state() -> None:
    metrics = {
        "state_current_uA": 165.0,
        "first_iteration_current_uA": 164.0,
        "first_accepted_current_uA": 163.0,
    }

    assert (
        campaign.first_jump_location(metrics, 143.567589)
        == "TES state file before solver"
    )


def test_first_jump_location_falls_through_to_first_accepted_step() -> None:
    reference = 143.567589
    metrics = {
        "state_current_uA": reference,
        "first_iteration_current_uA": reference,
        "first_accepted_current_uA": 145.0,
    }

    assert (
        campaign.first_jump_location(metrics, reference)
        == "first accepted timestep"
    )


def test_dict_diff_reports_nested_changes() -> None:
    left = {"solver": {"linear_system": "mumps"}, "bdf_order": 1}
    right = {"solver": {"linear_system": "hypre"}, "bdf_order": 2}

    diff = campaign.dict_diff(left, right)
    keys = {row["key"] for row in diff}

    assert "solver.linear_system" in keys
    assert "bdf_order" in keys


def test_source_keeps_observable_nomortar_hold_variant() -> None:
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert '"gate3_state_mumps_bdf1_hold5_nomortar"' in source
    assert '"steps": 5' in source


def test_diagnosis_uses_hold_case_for_mortar_comparison() -> None:
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert 'by_variant.get("gate3_state_mumps_bdf1_hold5")' in source
    assert 'by_variant.get("gate3_state_mumps_bdf1_hold5_nomortar")' in source
    assert "one-step variants may legitimately have row_count=0" in source


def test_restart_candidate_residual_splits_zero_diagonal_constraint_rows(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "firstsolve"
    # Saddle matrix [[2,0,1],[0,3,1],[1,1,0]], b=A*[1,2,0].
    (tmp_path / "firstsolve_a.dat").write_text(
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n",
        encoding="utf-8",
    )
    (tmp_path / "firstsolve_b.dat").write_text(
        "1 2\n2 6\n3 3\n",
        encoding="utf-8",
    )
    # SaveLinearSystem-style primal-only vector; multiplier is omitted.
    (tmp_path / "firstsolve_sol.dat").write_text(
        "1 1\n2 2\n",
        encoding="utf-8",
    )

    result = campaign.restart_candidate_residual(prefix)

    assert result["available"] is True
    assert result["rows"] == 3
    assert result["saved_solution_records"] == 2
    assert result["missing_solution_entries_zero_extended"] == 1
    assert result["primal_rows"] == 2
    assert result["constraint_rows"] == 1
    assert result["full"]["residual_l2"] == 0.0
    assert result["primal"]["residual_l2"] == 0.0
    assert result["constraint"]["residual_l2"] == 0.0


def test_restart_candidate_residual_detects_constraint_violation(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "firstsolve"
    (tmp_path / "firstsolve_a.dat").write_text(
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n",
        encoding="utf-8",
    )
    (tmp_path / "firstsolve_b.dat").write_text(
        "1 2\n2 6\n3 0\n",
        encoding="utf-8",
    )
    (tmp_path / "firstsolve_sol.dat").write_text(
        "1 1\n2 2\n",
        encoding="utf-8",
    )

    result = campaign.restart_candidate_residual(prefix)

    assert result["constraint"]["residual_l2"] == 3.0
    assert result["constraint"]["residual_max_abs"] == 3.0
    assert result["constraint"]["backward_error"] > 0.0


def test_matrix_block_decomposition_reports_saddle_blocks(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "blocks"
    (tmp_path / "blocks_a.dat").write_text(
        "1 1 2\n"
        "1 3 4\n"
        "2 2 3\n"
        "2 3 5\n"
        "3 1 4\n"
        "3 2 5\n",
        encoding="utf-8",
    )
    (tmp_path / "blocks_b.dat").write_text(
        "1 7\n2 8\n3 0\n",
        encoding="utf-8",
    )

    result = campaign.matrix_block_decomposition(prefix)

    assert result["available"] is True
    assert result["primal_rows"] == 2
    assert result["constraint_rows"] == 1
    assert result["blocks"]["K"]["records"] == 2
    assert result["blocks"]["Bt"]["records"] == 2
    assert result["blocks"]["B"]["records"] == 2
    assert result["blocks"]["D"]["records"] == 0
    assert result["B_minus_BtT_max_abs"] == 0.0


def test_compare_primal_systems_ignores_added_constraint_block(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    mortar = tmp_path / "mortar"
    nomortar = tmp_path / "nomortar"
    (tmp_path / "mortar_a.dat").write_text(
        "1 1 2\n"
        "1 3 4\n"
        "2 2 3\n"
        "2 3 5\n"
        "3 1 4\n"
        "3 2 5\n",
        encoding="utf-8",
    )
    (tmp_path / "mortar_b.dat").write_text(
        "1 7\n2 8\n3 0\n",
        encoding="utf-8",
    )
    (tmp_path / "nomortar_a.dat").write_text(
        "1 1 2\n"
        "2 2 3\n",
        encoding="utf-8",
    )
    (tmp_path / "nomortar_b.dat").write_text(
        "1 7\n2 8\n",
        encoding="utf-8",
    )

    result = campaign.compare_primal_systems(mortar, nomortar)

    assert result["available"] is True
    assert result["same_primal_row_set"] is True
    assert result["K_nonzero_difference_records"] == 0
    assert result["K_max_absolute_difference"] == 0.0
    assert result["rhs_nonzero_difference_records"] == 0
    assert result["rhs_max_absolute_difference"] == 0.0


def test_independent_direct_solve_is_best_effort_without_crashing(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "direct"
    (tmp_path / "direct_a.dat").write_text(
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n",
        encoding="utf-8",
    )
    (tmp_path / "direct_b.dat").write_text(
        "1 2\n2 6\n3 3\n",
        encoding="utf-8",
    )
    (tmp_path / "direct_sol.dat").write_text(
        "1 1\n2 2\n",
        encoding="utf-8",
    )

    result = campaign.independent_direct_solve(prefix)

    assert "available" in result
    if result["available"]:
        assert result["direct_relative_residual"] < 1.0e-10
        assert result["primal"]["delta_max_abs"] < 1.0e-10

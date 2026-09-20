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


def test_discover_linear_solve_dumps_orders_continuous_numbering(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "diag"
    for base in ("diag", "diag_2", "diag_10"):
        (tmp_path / f"{base}_a.dat").write_text("1 1 1\n", encoding="utf-8")
        (tmp_path / f"{base}_b.dat").write_text("1 1\n", encoding="utf-8")

    dumps = campaign.discover_linear_solve_dumps(prefix)

    assert [item["base"] for item in dumps] == ["diag", "diag_2", "diag_10"]
    assert [item["ordinal"] for item in dumps] == [1, 2, 3]


def test_compare_linear_solve_dump_pair_splits_rhs_and_matrix_blocks(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    left_a = tmp_path / "left_a.dat"
    left_b = tmp_path / "left_b.dat"
    right_a = tmp_path / "right_a.dat"
    right_b = tmp_path / "right_b.dat"

    left_a.write_text(
        "1 1 2\n"
        "1 3 4\n"
        "2 2 3\n"
        "2 3 5\n"
        "3 1 4\n"
        "3 2 5\n",
        encoding="utf-8",
    )
    right_a.write_text(
        "1 1 2.5\n"
        "1 3 4\n"
        "2 2 3\n"
        "2 3 5\n"
        "3 1 4\n"
        "3 2 5\n",
        encoding="utf-8",
    )
    left_b.write_text("1 7\n2 8\n3 0\n", encoding="utf-8")
    right_b.write_text("1 7.25\n2 8\n3 0.5\n", encoding="utf-8")

    result = campaign.compare_linear_solve_dump_pair(
        {
            "ordinal": 1,
            "A_path": left_a,
            "b_path": left_b,
            "x_path": None,
        },
        {
            "ordinal": 2,
            "A_path": right_a,
            "b_path": right_b,
            "x_path": None,
        },
    )

    assert result["available"] is True
    assert result["same_primal_constraint_partition"] is True
    assert result["A_blocks"]["K"]["difference_max_abs"] == 0.5
    assert result["A_blocks"]["B"]["difference_max_abs"] == 0.0
    assert result["A_blocks"]["Bt"]["difference_max_abs"] == 0.0
    assert result["b_primal"]["difference_max_abs"] == 0.25
    assert result["b_constraint"]["difference_max_abs"] == 0.5


def test_nonlinear_linear_system_sequence_compares_first_three_solves(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "diag"
    for ordinal, rhs in ((1, 1.0), (2, 1.0), (3, 2.0)):
        base = "diag" if ordinal == 1 else f"diag_{ordinal}"
        (tmp_path / f"{base}_a.dat").write_text("1 1 2\n", encoding="utf-8")
        (tmp_path / f"{base}_b.dat").write_text(f"1 {rhs}\n", encoding="utf-8")

    result = campaign.nonlinear_linear_system_sequence(prefix, limit=3)

    assert result["available"] is True
    assert result["dump_count_discovered"] == 3
    assert result["analyzed_count"] == 3
    assert "solve1__vs__solve2" in result["pairs"]
    assert "solve2__vs__solve3" in result["pairs"]
    assert result["pairs"]["solve1__vs__solve2"]["b_full"]["difference_l2"] == 0.0
    assert result["pairs"]["solve2__vs__solve3"]["b_full"]["difference_l2"] == 1.0


def test_direct_solve_transition_sensitivity_isolates_rhs_effect(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    a1 = tmp_path / "s1_a.dat"
    b1 = tmp_path / "s1_b.dat"
    a2 = tmp_path / "s2_a.dat"
    b2 = tmp_path / "s2_b.dat"
    matrix = (
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n"
    )
    a1.write_text(matrix, encoding="utf-8")
    a2.write_text(matrix, encoding="utf-8")
    b1.write_text("1 2\n2 6\n3 3\n", encoding="utf-8")
    b2.write_text("1 2.2\n2 6\n3 3\n", encoding="utf-8")

    result = campaign.direct_solve_transition_sensitivity(
        {"ordinal": 1, "A_path": a1, "b_path": b1},
        {"ordinal": 2, "A_path": a2, "b_path": b2},
    )

    assert "available" in result
    if result["available"]:
        primal = result["primal"]
        full = primal["full_A2b2_minus_A1b1"]
        rhs = primal["rhs_only_A1b2_minus_A1b1"]
        op = primal["operator_only_A2b1_minus_A1b1"]
        interaction = primal["interaction"]
        assert full["delta_l2"] > 0.0
        assert abs(full["delta_l2"] - rhs["delta_l2"]) < 1.0e-12
        assert op["delta_l2"] < 1.0e-12
        assert interaction["delta_l2"] < 1.0e-12
        assert result["primal_l2_ratios"]["rhs_only_to_full"] > 0.999999


def test_nonlinear_transition_direct_sensitivity_uses_adjacent_saved_solves(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "diag"
    matrix = (
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n"
    )
    for ordinal, rhs1 in ((1, 2.0), (2, 2.1), (3, 2.2)):
        base = "diag" if ordinal == 1 else f"diag_{ordinal}"
        (tmp_path / f"{base}_a.dat").write_text(matrix, encoding="utf-8")
        (tmp_path / f"{base}_b.dat").write_text(
            f"1 {rhs1}\n2 6\n3 3\n",
            encoding="utf-8",
        )

    result = campaign.nonlinear_transition_direct_sensitivity(prefix, limit=3)

    assert result["available"] is True
    assert result["dump_count_analyzed"] == 3
    assert "solve1__vs__solve2" in result["pairs"]
    assert "solve2__vs__solve3" in result["pairs"]


def test_matrix_dimension_comes_from_A_when_zero_rhs_tail_is_skipped(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "tail"
    (tmp_path / "tail_a.dat").write_text(
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n",
        encoding="utf-8",
    )
    # Row 3 is a zero-RHS multiplier row and is omitted by Save Skip Zeros.
    (tmp_path / "tail_b.dat").write_text(
        "1 2\n"
        "2 6\n",
        encoding="utf-8",
    )
    (tmp_path / "tail_sol.dat").write_text(
        "1 1\n"
        "2 2\n",
        encoding="utf-8",
    )

    dimension = campaign.matrix_dump_dimension(tmp_path / "tail_a.dat")
    residual = campaign.restart_candidate_residual(prefix)
    blocks = campaign.matrix_block_decomposition(prefix)

    assert dimension["rows"] == 3
    assert residual["rows"] == 3
    assert residual["rhs_saved_records"] == 2
    assert residual["rhs_implicit_zero_entries"] == 1
    assert residual["constraint_rows"] == 1
    assert blocks["rows"] == 3
    assert blocks["constraint_rows"] == 1
    assert blocks["rhs_implicit_zero_entries"] == 1


def test_direct_sensitivity_accepts_different_rhs_max_indices_with_same_A_dimension(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    matrix = (
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n"
    )
    a1 = tmp_path / "s1_a.dat"
    a2 = tmp_path / "s2_a.dat"
    b1 = tmp_path / "s1_b.dat"
    b2 = tmp_path / "s2_b.dat"
    a1.write_text(matrix, encoding="utf-8")
    a2.write_text(matrix, encoding="utf-8")
    # Same 3x3 physical systems, but solve1 happens to save a nonzero row-3
    # RHS while solve2 omits the now-zero trailing entry.
    b1.write_text("1 2\n2 6\n3 0.5\n", encoding="utf-8")
    b2.write_text("1 2.1\n2 6\n", encoding="utf-8")

    result = campaign.direct_solve_transition_sensitivity(
        {"ordinal": 1, "A_path": a1, "b_path": b1, "sizes_path": None},
        {"ordinal": 2, "A_path": a2, "b_path": b2, "sizes_path": None},
    )

    assert "available" in result
    if result["available"]:
        assert result["rows"] == 3
        assert result["constraint_rows"] == 1
        assert result["left_dimension"]["rows"] == 3
        assert result["right_dimension"]["rows"] == 3


def test_linear_dump_provenance_excludes_heterogeneous_numbered_dump(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "diag"

    (tmp_path / "diag_a.dat").write_text(
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n",
        encoding="utf-8",
    )
    (tmp_path / "diag_b.dat").write_text("1 2\n2 6\n3 0\n", encoding="utf-8")

    # Continuously-numbered second dump is structurally different: it must
    # not be interpreted as the next outer nonlinear iteration.
    (tmp_path / "diag_2_a.dat").write_text(
        "1 1 4\n"
        "2 2 5\n",
        encoding="utf-8",
    )
    (tmp_path / "diag_2_b.dat").write_text("1 1\n2 1\n", encoding="utf-8")

    log = tmp_path / "launcher.log"
    log.write_text(
        "HeatSolve: Assembly done\n"
        "SaveLinearSystem: Saving matrix to: diag_a.dat\n"
        "SaveLinearSystem: Saving matrix rhs to: diag_b.dat\n"
        "SaveLinearSystem: Saving matrix to: diag_2_a.dat\n"
        "SaveLinearSystem: Saving matrix rhs to: diag_2_b.dat\n"
        "ComputeChange: NS (ITER=1) (NRM,RELC): (1 1) :: heat equation\n",
        encoding="utf-8",
    )

    provenance = campaign.linear_dump_provenance(prefix, log)
    sequence = campaign.nonlinear_linear_system_sequence(
        prefix, provenance=provenance
    )
    sensitivity = campaign.nonlinear_transition_direct_sensitivity(
        prefix, provenance=provenance
    )

    assert provenance["dump_count"] == 2
    assert provenance["comparable_outer_candidate_ordinals"] == [1]
    assert provenance["heterogeneous_ordinals"] == [2]
    assert provenance["records"][0]["role"] == "reference_outer_candidate"
    assert (
        provenance["records"][1]["role"]
        == "heterogeneous_auxiliary_or_restricted"
    )
    assert provenance["records"][0]["log_event"]["next_computechange_iter"] == 1
    assert provenance["records"][1]["log_event"]["next_computechange_iter"] == 1
    assert sequence["available"] is False
    assert sequence["heterogeneous_count"] == 1
    assert sensitivity["available"] is False
    assert "fewer than two structurally homogeneous" in sensitivity["reason"]


def test_linear_dump_provenance_allows_same_shape_outer_candidates(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "diag"
    matrix = (
        "1 1 2\n"
        "1 3 1\n"
        "2 2 3\n"
        "2 3 1\n"
        "3 1 1\n"
        "3 2 1\n"
    )
    for ordinal, rhs in ((1, 2.0), (2, 2.1)):
        base = "diag" if ordinal == 1 else f"diag_{ordinal}"
        (tmp_path / f"{base}_a.dat").write_text(matrix, encoding="utf-8")
        (tmp_path / f"{base}_b.dat").write_text(
            f"1 {rhs}\n2 6\n3 0\n",
            encoding="utf-8",
        )

    provenance = campaign.linear_dump_provenance(prefix)
    sequence = campaign.nonlinear_linear_system_sequence(
        prefix, provenance=provenance
    )

    assert provenance["heterogeneous_count"] == 0
    assert provenance["comparable_outer_candidate_ordinals"] == [1, 2]
    assert sequence["available"] is True
    assert "solve1__vs__solve2" in sequence["pairs"]


def test_discover_linear_solve_dumps_excludes_textual_sibling_variant(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "case_x"
    for base in ("case_x", "case_x_1", "case_x_2", "case_x_nomortar", "case_x_hold5"):
        (tmp_path / f"{base}_a.dat").write_text("1 1 1\n", encoding="utf-8")
        (tmp_path / f"{base}_b.dat").write_text("1 1\n", encoding="utf-8")

    dumps = campaign.discover_linear_solve_dumps(prefix)

    assert [item["base"] for item in dumps] == ["case_x", "case_x_1", "case_x_2"]
    assert [item["continuous_number"] for item in dumps] == [None, 1, 2]
    assert dumps[0]["ignored_sibling_bases"] == ["case_x_hold5", "case_x_nomortar"]


def test_discover_linear_solve_dumps_does_not_treat_numeric_case_prefix_as_suffix(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "case_1"
    for base in ("case_1", "case_10", "case_1_2"):
        (tmp_path / f"{base}_a.dat").write_text("1 1 1\n", encoding="utf-8")
        (tmp_path / f"{base}_b.dat").write_text("1 1\n", encoding="utf-8")

    dumps = campaign.discover_linear_solve_dumps(prefix)

    assert [item["base"] for item in dumps] == ["case_1", "case_1_2"]
    assert "case_10" in dumps[0]["ignored_sibling_bases"]


def test_linear_dump_provenance_reads_solver_log_computechange_context(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    prefix = tmp_path / "diag"
    for base in ("diag", "diag_1"):
        (tmp_path / f"{base}_a.dat").write_text("1 1 2\n", encoding="utf-8")
        (tmp_path / f"{base}_b.dat").write_text("1 1\n", encoding="utf-8")

    log = tmp_path / "solver.log"
    log.write_text(
        "SaveLinearSystem: Saving matrix to: diag_a.dat\n"
        "SaveLinearSystem: Saving matrix rhs to: diag_b.dat\n"
        "ComputeChange: NS (ITER=1) (NRM,RELC): (1 1) :: heat equation\n"
        "SaveLinearSystem: Saving matrix to: diag_1_a.dat\n"
        "SaveLinearSystem: Saving matrix rhs to: diag_1_b.dat\n"
        "ComputeChange: NS (ITER=2) (NRM,RELC): (1 1) :: heat equation\n",
        encoding="utf-8",
    )

    provenance = campaign.linear_dump_provenance(prefix, log)

    assert provenance["log_save_event_count"] == 2
    assert provenance["records"][0]["log_event"]["next_computechange_iter"] == 1
    assert provenance["records"][1]["log_event"]["previous_computechange_iter"] == 1
    assert provenance["records"][1]["log_event"]["next_computechange_iter"] == 2

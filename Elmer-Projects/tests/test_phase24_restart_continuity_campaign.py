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

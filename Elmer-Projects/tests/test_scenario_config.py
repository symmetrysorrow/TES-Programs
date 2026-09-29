from __future__ import annotations

from pathlib import Path

import pytest

from scripts.support import scenario_config
from scripts.support.build_cases import body_force_blocks
from scripts.support.scenario_config import ROOT, compile_project, load_scenario

CASES = ROOT / "cases"
MODELS = ROOT / "models"


def changed_file(tmp_path: Path, source: Path, before: str, after: str) -> Path:
    text = source.read_text(encoding="utf-8")
    assert before in text
    target = tmp_path / source.name
    target.write_text(text.replace(before, after, 1), encoding="utf-8")
    return target


def test_group_expansion_and_mesh_reuse() -> None:
    nominal = load_scenario(CASES / "post_four_tes_nominal.toml")
    lower = load_scenario(CASES / "post_four_tes_low_bias.toml")
    assert [item["id"] for item in nominal["tes"]] == ["west.a", "west.b", "east.a", "east.b"]
    assert [item["x"] for item in nominal["tes"]] == pytest.approx([-0.007, -0.003, 0.003, 0.007])
    assert nominal["mesh_hash"] == lower["mesh_hash"]
    assert nominal["steady_hash"] != lower["steady_hash"]
    assert nominal["tes"][3]["sources"]["I_bias"] == "TES east.b"
    assert nominal["tes"][2]["sources"]["R_sh"] == "group east"
    project, mesh, steady, pulse = compile_project(nominal)
    assert project["cases"][steady]["mesh"] == mesh
    assert len(project["cases"][steady]["tes_circuits"]) == 4
    assert project["cases"][pulse]["restart_file_base"] == steady


def test_single_pixel_alpha_changes_steady_and_reuses_mesh() -> None:
    baseline = load_scenario(CASES / "single_pixel.toml")
    changed = load_scenario(CASES / "single_pixel_alpha240.toml")
    assert baseline["mesh_hash"] == changed["mesh_hash"]
    assert baseline["steady_hash"] != changed["steady_hash"]
    assert len(changed["tes"]) == 1
    project, mesh, steady, pulse = compile_project(changed)
    assert mesh == "mesh_hybrid_fullconf_h8"
    assert project["cases"][steady]["tes_circuits"] == [{"alpha": 240}]
    assert project["cases"][pulse]["bdf_order"] == 2
    assert len(project["cases"][pulse]["timesteps"]) == 26


def test_pulse_change_reuses_steady(tmp_path: Path) -> None:
    original = CASES / "single_pixel.toml"
    changed = changed_file(tmp_path, original, 'energy = "1332[keV]"', 'energy = "1000[keV]"')
    before, after = load_scenario(original), load_scenario(changed)
    assert before["mesh_hash"] == after["mesh_hash"]
    assert before["steady_hash"] == after["steady_hash"]
    assert before["pulse_hash"] != after["pulse_hash"]


def test_unknown_tes_override_is_rejected(tmp_path: Path) -> None:
    original = CASES / "post_four_tes_nominal.toml"
    changed = changed_file(tmp_path, original, '"east.b"', '"missing.tes"')
    with pytest.raises(ValueError, match="unknown TES override"):
        load_scenario(changed)


def test_geometry_change_invalidates_mesh(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    changed_file(tmp_path, MODELS / "post_four_tes.toml", 'absorber_length = "16[mm]"', 'absorber_length = "17[mm]"')
    before = load_scenario(CASES / "post_four_tes_nominal.toml")
    monkeypatch.setattr(scenario_config, "MODEL_DIR", tmp_path)
    after = load_scenario(CASES / "post_four_tes_nominal.toml")
    assert before["mesh_hash"] != after["mesh_hash"]


def test_mesh_size_change_invalidates_mesh(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    changed_file(tmp_path, MODELS / "post_four_tes.toml", 'tes_local_size = "8[um]"', 'tes_local_size = "7[um]"')
    before = load_scenario(CASES / "post_four_tes_nominal.toml")
    monkeypatch.setattr(scenario_config, "MODEL_DIR", tmp_path)
    after = load_scenario(CASES / "post_four_tes_nominal.toml")
    assert before["mesh_hash"] != after["mesh_hash"]


def test_prebuilt_single_pixel_rejects_geometry_edit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    changed_file(tmp_path, MODELS / "single_pixel.toml", 'absorber_length = "1[mm]"', 'absorber_length = "2[mm]"')
    monkeypatch.setattr(scenario_config, "MODEL_DIR", tmp_path)
    with pytest.raises(ValueError, match="prebuilt single-pixel mesh"):
        load_scenario(CASES / "single_pixel.toml")


def test_one_numbered_tes_uses_numbered_heat_source() -> None:
    blocks = body_force_blocks("circuit_inner", ["TES_T1"], False)
    assert any("TESParallelHeatSource1" in line for line in blocks)


def test_misspelled_case_key_is_rejected(tmp_path: Path) -> None:
    changed = changed_file(tmp_path, CASES / "single_pixel.toml", 'model = "single_pixel"', 'modle = "single_pixel"')
    with pytest.raises(ValueError, match="unknown keys"):
        load_scenario(changed)

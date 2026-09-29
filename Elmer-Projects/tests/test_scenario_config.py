from __future__ import annotations

from pathlib import Path

import pytest

from support.build_cases import body_force_blocks
from support.scenario_config import ROOT, compile_project, load_scenario

PROJECTS = ROOT / "projects"
MODELS = PROJECTS / "models"


def changed_tree(
    tmp_path: Path,
    project_name: str,
    model_name: str,
    *,
    project_replace: tuple[str, str] | None = None,
    model_replace: tuple[str, str] | None = None,
) -> Path:
    root = tmp_path / "projects"
    models = root / "models"
    models.mkdir(parents=True)

    project_text = (PROJECTS / f"{project_name}.toml").read_text(encoding="utf-8")
    if project_replace:
        before, after = project_replace
        assert before in project_text
        project_text = project_text.replace(before, after, 1)
    target_project = root / f"{project_name}.toml"
    target_project.write_text(project_text, encoding="utf-8")

    model_text = (MODELS / f"{model_name}.toml").read_text(encoding="utf-8")
    if model_replace:
        before, after = model_replace
        assert before in model_text
        model_text = model_text.replace(before, after, 1)
    (models / f"{model_name}.toml").write_text(model_text, encoding="utf-8")
    return target_project


def test_group_expansion_and_explicit_circuits_reuse_mesh() -> None:
    nominal = load_scenario(PROJECTS / "post_four_tes_nominal.toml")
    lower = load_scenario(PROJECTS / "post_four_tes_low_bias.toml")
    assert [item["id"] for item in nominal["tes"]] == ["west.a", "west.b", "east.a", "east.b"]
    assert [item["x"] for item in nominal["tes"]] == pytest.approx([-0.007, -0.003, 0.003, 0.007])
    assert nominal["mesh_hash"] == lower["mesh_hash"]
    assert nominal["steady_hash"] != lower["steady_hash"]
    assert nominal["tes"][3]["sources"]["I_bias"] == "project post_four_tes_nominal.toml"
    assert nominal["tes"][2]["circuit"]["R_sh"] == "3.8[mohm]"
    project, mesh, steady, pulse = compile_project(nominal)
    assert project["cases"][steady]["mesh"] == mesh
    assert len(project["cases"][steady]["tes_circuits"]) == 4
    assert project["cases"][pulse]["restart_file_base"] == steady


def test_single_pixel_alpha_changes_steady_and_reuses_mesh() -> None:
    baseline = load_scenario(PROJECTS / "single_pixel.toml")
    changed = load_scenario(PROJECTS / "single_pixel_alpha240.toml")
    assert baseline["mesh_hash"] == changed["mesh_hash"]
    assert baseline["steady_hash"] != changed["steady_hash"]
    assert len(changed["tes"]) == 1
    project, mesh, steady, pulse = compile_project(changed)
    assert mesh == "mesh_hybrid_fullconf_h8"
    circuit = project["cases"][steady]["tes_circuits"][0]
    assert circuit["alpha"] == 240
    assert circuit["I_bias"] == "715[uA]"
    assert "I_0" in circuit
    # The generated JSON is a flattened adapter for the legacy single-TES
    # builder; user-facing ownership still remains in the project TOML.
    assert project["parameter_expressions"]["alpha"] == "240"
    assert project["parameter_expressions"]["I_bias"] == "715[uA]"
    assert project["cases"][pulse]["bdf_order"] == 2
    assert len(project["cases"][pulse]["timesteps"]) == 26


def test_pulse_change_reuses_steady(tmp_path: Path) -> None:
    changed = changed_tree(
        tmp_path,
        "single_pixel",
        "single_pixel",
        project_replace=('energy = "1332[keV]"', 'energy = "1000[keV]"'),
    )
    before, after = load_scenario(PROJECTS / "single_pixel.toml"), load_scenario(changed)
    assert before["mesh_hash"] == after["mesh_hash"]
    assert before["steady_hash"] == after["steady_hash"]
    assert before["pulse_hash"] != after["pulse_hash"]


def test_unknown_tes_id_is_rejected_with_available_ids(tmp_path: Path) -> None:
    changed = changed_tree(
        tmp_path,
        "post_four_tes_nominal",
        "post_four_tes",
        project_replace=('"east.b"', '"missing.tes"'),
    )
    with pytest.raises(ValueError, match="unknown circuit TES IDs.*available IDs"):
        load_scenario(changed)


def test_missing_final_circuit_value_is_rejected(tmp_path: Path) -> None:
    changed = changed_tree(
        tmp_path,
        "single_pixel",
        "single_pixel",
        project_replace=("beta = 5.03\n", ""),
    )
    with pytest.raises(ValueError, match="missing final values"):
        load_scenario(changed)


def test_model_cannot_reintroduce_project_owned_values(tmp_path: Path) -> None:
    changed = changed_tree(
        tmp_path,
        "single_pixel",
        "single_pixel",
        model_replace=("[parameters]\n", '[parameters]\nI_bias = "715[uA]"\n'),
    )
    with pytest.raises(ValueError, match="project-owned circuit keys"):
        load_scenario(changed)


def test_geometry_change_invalidates_mesh(tmp_path: Path) -> None:
    changed = changed_tree(
        tmp_path,
        "post_four_tes_nominal",
        "post_four_tes",
        model_replace=('absorber_length = "16[mm]"', 'absorber_length = "17[mm]"'),
    )
    before = load_scenario(PROJECTS / "post_four_tes_nominal.toml")
    after = load_scenario(changed)
    assert before["mesh_hash"] != after["mesh_hash"]


def test_mesh_size_change_invalidates_mesh(tmp_path: Path) -> None:
    changed = changed_tree(
        tmp_path,
        "post_four_tes_nominal",
        "post_four_tes",
        model_replace=('tes_local_size = "8[um]"', 'tes_local_size = "7[um]"'),
    )
    before = load_scenario(PROJECTS / "post_four_tes_nominal.toml")
    after = load_scenario(changed)
    assert before["mesh_hash"] != after["mesh_hash"]


def test_prebuilt_single_pixel_rejects_geometry_edit(tmp_path: Path) -> None:
    changed = changed_tree(
        tmp_path,
        "single_pixel",
        "single_pixel",
        model_replace=('absorber_length = "1[mm]"', 'absorber_length = "2[mm]"'),
    )
    with pytest.raises(ValueError, match="prebuilt single-pixel mesh"):
        load_scenario(changed)


def test_one_numbered_tes_uses_numbered_heat_source() -> None:
    blocks = body_force_blocks("circuit_inner", ["TES_T1"], False)
    assert any("TESParallelHeatSource1" in line for line in blocks)


def test_misspelled_project_key_is_rejected(tmp_path: Path) -> None:
    changed = changed_tree(
        tmp_path,
        "single_pixel",
        "single_pixel",
        project_replace=('model = "models/single_pixel.toml"', 'modle = "models/single_pixel.toml"'),
    )
    with pytest.raises(ValueError, match="unknown keys"):
        load_scenario(changed)

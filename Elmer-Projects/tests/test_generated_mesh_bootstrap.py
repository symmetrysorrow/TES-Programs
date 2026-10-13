"""First-stage contract tests for opt-in generated meshes and CPU onboarding."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import main
from support import cpu_setup, scenario_cli
from support.scenario_config import PROJECT_DIR, compile_project, load_scenario


def test_single_pixel_generated_is_opt_in_and_preserves_fixed_reference():
    reference = load_scenario(PROJECT_DIR / "single_pixel_alpha240.toml")
    generated = load_scenario(PROJECT_DIR / "single_pixel_generated.toml")

    assert reference["source_mesh"] == "mesh_hybrid_fullconf_h8"
    assert generated["source_mesh"] is None
    assert generated["mesh"]["tes_local_size"] == pytest.approx(8e-6)
    assert generated["physics"]["tes_circuits"] == reference["physics"]["tes_circuits"]
    assert generated["mesh_hash"] != reference["mesh_hash"]
    _, mesh_name, _, _ = compile_project(generated)
    assert mesh_name != "mesh_hybrid_fullconf_h8"


def test_generated_mesh_is_created_from_config_and_requires_matching_provenance(tmp_path, monkeypatch):
    scenario = load_scenario(PROJECT_DIR / "single_pixel_generated.toml")
    _, name, _, _ = compile_project(scenario)
    monkeypatch.setattr(scenario_cli, "ROOT", tmp_path)
    monkeypatch.setattr(scenario_cli.shutil, "which", lambda name: "ElmerGrid" if name == "ElmerGrid" else None)
    commands = []

    def fake_run(command, *, cwd=None):
        commands.append(command)
        if "--output" in command:
            path = Path(command[command.index("--output") + 1])
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("Gmsh test fixture", encoding="utf-8")
        elif "-out" in command:
            out_dir = Path(command[command.index("-out") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            for filename in ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary", "mesh.names"):
                (out_dir / filename).write_text("test", encoding="utf-8")
        else:
            raise AssertionError(f"unexpected command: {command}")

    monkeypatch.setattr(scenario_cli, "run_command", fake_run)
    scenario_cli.mesh_command(scenario, tmp_path / "generated_project.json", name)
    mesh_dir = tmp_path / "workspace" / "work" / "meshes" / name
    provenance = json.loads((mesh_dir / "PROVENANCE.json").read_text(encoding="utf-8"))
    assert provenance["mesh_hash"] == scenario["mesh_hash"]
    assert len(commands) == 2
    assert "--tes-positions" in commands[0]
    assert "--conformal-tes-stack" in commands[0]
    assert "--conformal-abs" in commands[0]

    scenario_cli.mesh_command(scenario, tmp_path / "generated_project.json", name)
    assert len(commands) == 2

    provenance["mesh_hash"] = "wrong-cache-id"
    (mesh_dir / "PROVENANCE.json").write_text(json.dumps(provenance), encoding="utf-8")
    with pytest.raises(RuntimeError, match="does not match its TOML configuration"):
        scenario_cli.mesh_command(scenario, tmp_path / "generated_project.json", name)


def test_main_setup_dispatch_is_a_real_command(monkeypatch):
    calls = []
    monkeypatch.setattr(cpu_setup, "setup_cpu", lambda: calls.append("called"))
    assert main.main(["setup"]) == 0
    assert calls == ["called"]

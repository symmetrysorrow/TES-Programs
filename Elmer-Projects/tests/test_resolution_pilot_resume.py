import json
import csv

import pytest

from scripts.prep import run_singlepixel_resolution_pilot as pilot


def write_iterations(path, *, nonlinear_iter=3, temperature=0.1685, current=145e-6, resistance=0.015) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["time_step", "nonlinear_iter", "tes_temperature_K", "raw_current_A", "tes_resistance_ohm"])
        writer.writeheader()
        writer.writerow({"time_step": 1, "nonlinear_iter": nonlinear_iter, "tes_temperature_K": temperature, "raw_current_A": current, "tes_resistance_ohm": resistance})


def test_case_complete_requires_success_series_and_matching_project(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "ROOT", tmp_path)
    project = tmp_path / "pilot.json"
    project.write_text('{"end_us": 105}\n', encoding="utf-8")
    case = "case_pilot_tes_layers_1_pulse"
    result = tmp_path / "results" / case
    result.mkdir(parents=True)
    (result / f"{case}_series.csv").write_text("time_s,tes_current_A\n", encoding="utf-8")
    (result / "manifest.json").write_text(json.dumps({
        "exit_code": 0,
        "inputs_sha256": {project.name: pilot.sha256(project)},
        "restart_from": case.replace("_pulse", "_steady"),
    }), encoding="utf-8")
    write_iterations(result / f"{case}_iterations.csv")
    steady = tmp_path / "results" / case.replace("_pulse", "_steady")
    steady.mkdir()
    write_iterations(steady / f"{steady.name}_iterations.csv")
    assert pilot.case_complete(case, project)

    project.write_text('{"end_us": 225}\n', encoding="utf-8")
    assert not pilot.case_complete(case, project)


def test_case_complete_rejects_iteration_cap_and_wrong_steady_branch(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "ROOT", tmp_path)
    project = tmp_path / "pilot.json"
    project.write_text("{}\n", encoding="utf-8")
    case = "case_pilot_tes_layers_2_pulse"
    result = tmp_path / "results" / case
    result.mkdir(parents=True)
    (result / f"{case}_series.csv").write_text("time_s,tes_current_A\n", encoding="utf-8")
    (result / "manifest.json").write_text(json.dumps({
        "exit_code": 0,
        "inputs_sha256": {project.name: pilot.sha256(project)},
        "restart_from": case.replace("_pulse", "_steady"),
    }), encoding="utf-8")
    write_iterations(result / f"{case}_iterations.csv", nonlinear_iter=pilot.PILOT_PULSE_NONLINEAR_MAX_ITERATIONS)
    steady = tmp_path / "results" / case.replace("_pulse", "_steady")
    steady.mkdir()
    write_iterations(steady / f"{steady.name}_iterations.csv", temperature=0.15, current=715e-6, resistance=1e-6)
    reason = pilot.case_failure_reason(case, project)
    assert reason is not None
    assert "iteration cap" in reason


def test_steady_state_validation_rejects_low_temperature_branch(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "ROOT", tmp_path)
    case = "case_pilot_stack_steady"
    result = tmp_path / "results" / case
    result.mkdir(parents=True)
    write_iterations(result / f"{case}_iterations.csv", temperature=0.15, current=715e-6, resistance=1e-6)
    reason = pilot.steady_state_failure_reason(case)
    assert reason is not None
    assert "temperature" in reason


def test_build_project_uses_pilot_nonlinear_limits_and_tolerance(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "resolve_elmergrid", lambda: "ElmerGrid.exe")
    project, _ = pilot.build_project("tes_layers", [1.0, 2.0], 0.625, 105.0, None, "direct", tmp_path / "pilot.json")
    base = f"case_pilot_tes_layers_1_{pilot.PILOT_IMPLEMENTATION_SUFFIX}"
    assert project["cases"][f"{base}_steady"]["solver"]["nonlinear_max_iterations"] == 60
    assert project["cases"][f"{base}_pulse"]["solver"]["nonlinear_max_iterations"] == 25
    assert project["cases"][f"{base}_steady"]["solver"]["nonlinear_convergence_tolerance"] == 1e-7
    assert project["cases"][f"{base}_pulse"]["solver"]["nonlinear_convergence_tolerance"] == 1e-7


def test_steady_iteration_cap_does_not_reject_iteration_25(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "ROOT", tmp_path)
    case = "case_pilot_sio2_1_layers_1_alg5qg1_steady"
    result = tmp_path / "results" / case
    result.mkdir(parents=True)
    write_iterations(result / f"{case}_iterations.csv", nonlinear_iter=25)
    assert pilot.iteration_failure_reason(case) is None


def test_project_hash_changes_when_steady_cap_changes(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "resolve_elmergrid", lambda: "ElmerGrid.exe")
    output = tmp_path / "pilot.json"
    pilot.build_project("tes_layers", [1.0, 2.0], 0.625, 105.0, None, "direct", output)
    original_hash = pilot.sha256(output)
    monkeypatch.setattr(pilot, "PILOT_STEADY_NONLINEAR_MAX_ITERATIONS", 61)
    pilot.build_project("tes_layers", [1.0, 2.0], 0.625, 105.0, None, "direct", output)
    assert pilot.sha256(output) != original_hash


def test_project_hash_changes_when_pilot_tolerance_changes(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "resolve_elmergrid", lambda: "ElmerGrid.exe")
    output = tmp_path / "pilot.json"
    pilot.build_project("tes_layers", [1.0, 2.0], 0.625, 105.0, None, "direct", output)
    original_hash = pilot.sha256(output)
    monkeypatch.setattr(pilot, "PILOT_NONLINEAR_CONVERGENCE_TOLERANCE", 3e-7)
    pilot.build_project("tes_layers", [1.0, 2.0], 0.625, 105.0, None, "direct", output)
    assert pilot.sha256(output) != original_hash


def test_build_project_namespaces_all_pilot_artifacts(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "resolve_elmergrid", lambda: "ElmerGrid.exe")
    project, cases = pilot.build_project("tes_layers", [1.0, 2.0], 0.625, 105.0, None, "direct", tmp_path / "pilot.json")
    pulse, mesh = cases[0]
    suffix = pilot.PILOT_IMPLEMENTATION_SUFFIX
    assert pulse.endswith(f"_{suffix}_pulse")
    assert mesh.endswith(f"_{suffix}")
    commands = project["meshes"][mesh]["recipe"]["commands"]
    assert f"_{suffix}.msh" in commands[0]
    assert f"work/meshes/{mesh}" in commands[1]


def test_pilot_output_paths_are_versioned(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "OUT_DIR", tmp_path)
    project, manifest = pilot.pilot_output_paths("stack_local_size_um", "stack_refine")
    suffix = pilot.PILOT_IMPLEMENTATION_SUFFIX
    assert project.name == f"singlepixel_stack_refine_stack_local_size_um_{suffix}_pilot.json"
    assert manifest.name == f"singlepixel_stack_refine_stack_local_size_um_{suffix}_pilot_manifest.json"


def test_mesh_complete_requires_elmer_mesh_header(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "ROOT", tmp_path)
    assert not pilot.mesh_complete("mesh_pilot")
    header = tmp_path / "work" / "meshes" / "mesh_pilot" / "mesh.header"
    header.parent.mkdir(parents=True)
    header.write_text("1 1 1\n", encoding="utf-8")
    assert pilot.mesh_complete("mesh_pilot")


def test_mesh_complete_rejects_existing_bad_raw_mesh(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "ROOT", tmp_path)
    header = tmp_path / "work" / "meshes" / "mesh_pilot" / "mesh.header"
    header.parent.mkdir(parents=True)
    header.write_text("1 1 1\n", encoding="utf-8")
    raw = tmp_path / "gmsh" / "pilot.msh"
    monkeypatch.setattr(pilot, "inspect_mesh_file", lambda path: {"reasons": ["face ratio 21.3 exceeds 4"]})
    try:
        pilot.mesh_complete("mesh_pilot", raw)
    except RuntimeError as error:
        assert "mesh_pilot" in str(error)
        assert "face ratio 21.3" in str(error)
    else:
        raise AssertionError("bad existing mesh was reused")


def test_raw_mesh_from_recipe_resolves_output_under_root(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(pilot, "ROOT", tmp_path)
    path = pilot.raw_mesh_from_recipe([
        'python generate_hybrid_prism_geometry.py project.json --output "gmsh/pilot mesh.msh" --mesh-algorithm 5',
        "ElmerGrid 14 2 gmsh/pilot.msh",
    ])
    assert path == tmp_path / "gmsh" / "pilot mesh.msh"


def test_pilot_mesh_command_pins_algorithm_5(monkeypatch) -> None:
    monkeypatch.setattr(pilot, "resolve_elmergrid", lambda: "ElmerGrid.exe")
    generator, _ = pilot.mesh_command("pilot.msh", "mesh_pilot", [])
    index = generator.index("--mesh-algorithm")
    assert generator[index + 1] == "5"


def test_mesh_recipe_quotes_elmergrid_path_with_spaces(monkeypatch) -> None:
    elmergrid = r"C:\Program Files\Elmer 26.1-Release\bin\ElmerGrid.exe"
    monkeypatch.setattr(pilot.shutil, "which", lambda name: elmergrid if name == "ElmerGrid" else None)
    _, grid = pilot.mesh_command("pilot mesh.msh", "mesh pilot", [])
    command = pilot.recipe_command(grid)
    assert grid[0] == elmergrid
    assert command.startswith(f'"{elmergrid}" ')
    assert '"gmsh/pilot mesh.msh"' in command
    assert '"work/meshes/mesh pilot"' in command


def test_resolve_elmergrid_fails_clearly_when_absent(monkeypatch) -> None:
    monkeypatch.setattr(pilot.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="ElmerGrid executable was not found"):
        pilot.resolve_elmergrid()

import json

import run


def test_restart_result_is_not_reused_after_failed_run(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    restart = tmp_path / "work" / "meshes" / "mesh" / f"{case}.result"
    restart.parent.mkdir(parents=True)
    restart.write_text("partial restart", encoding="utf-8")
    manifest = tmp_path / "results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    project = tmp_path / "project.json"
    project.write_text("{}", encoding="utf-8")
    manifest.write_text(json.dumps({"exit_code": 3, "errors": ["MUMPS failed"]}), encoding="utf-8")

    assert not run.restart_result_is_reusable(restart, case, project)


def test_restart_result_reuses_success_or_untracked_legacy_result(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    restart = tmp_path / "work" / "meshes" / "mesh" / f"{case}.result"
    restart.parent.mkdir(parents=True)
    restart.write_text("complete restart", encoding="utf-8")
    project = tmp_path / "project.json"
    project.write_text("{}", encoding="utf-8")

    assert run.restart_result_is_reusable(restart, case, project)

    manifest = tmp_path / "results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({
        "exit_code": 0,
        "errors": [],
        "inputs_sha256": {project.name: run.sha256(project)},
    }), encoding="utf-8")
    assert run.restart_result_is_reusable(restart, case, project)


def test_restart_result_is_not_reused_after_project_change(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    restart = tmp_path / "work" / "meshes" / "mesh" / f"{case}.result"
    restart.parent.mkdir(parents=True)
    restart.write_text("complete restart", encoding="utf-8")
    project = tmp_path / "project.json"
    project.write_text('{"linear_system": "mumps"}', encoding="utf-8")
    manifest = tmp_path / "results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({
        "exit_code": 0,
        "errors": [],
        "inputs_sha256": {project.name: run.sha256(project)},
    }), encoding="utf-8")

    project.write_text('{"linear_system": "iterative_hypre_boomeramg"}', encoding="utf-8")
    assert not run.restart_result_is_reusable(restart, case, project)

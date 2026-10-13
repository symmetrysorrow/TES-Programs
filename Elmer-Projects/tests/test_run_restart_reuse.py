import json

from src import run


def manifest_inputs(tmp_path, case, project):
    sif = tmp_path / "workspace/generated/cases" / f"{case}.sif"
    sif.parent.mkdir(parents=True, exist_ok=True)
    sif.write_text("Solver = test", encoding="utf-8")
    return {project.name: run.sha256(project), sif.name: run.sha256(sif)}


def test_restart_result_is_not_reused_after_failed_run(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    restart = tmp_path / "workspace/work" / "meshes" / "mesh" / f"{case}.result"
    restart.parent.mkdir(parents=True)
    restart.write_text("partial restart", encoding="utf-8")
    manifest = tmp_path / "outputs/results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    project = tmp_path / "project.json"
    project.write_text("{}", encoding="utf-8")
    manifest.write_text(json.dumps({
        "exit_code": 3, "errors": ["MUMPS failed"],
        "inputs_sha256": manifest_inputs(tmp_path, case, project),
    }), encoding="utf-8")

    assert not run.restart_result_is_reusable(restart, case, project)


def test_restart_result_reuses_success_with_provenance_only(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    restart = tmp_path / "workspace/work" / "meshes" / "mesh" / f"{case}.result"
    restart.parent.mkdir(parents=True)
    restart.write_text("complete restart", encoding="utf-8")
    project = tmp_path / "project.json"
    project.write_text("{}", encoding="utf-8")
    inputs = manifest_inputs(tmp_path, case, project)

    assert not run.restart_result_is_reusable(restart, case, project)

    manifest = tmp_path / "outputs/results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({
        "exit_code": 0,
        "errors": [],
        "solver_completed": True,
        "runtime_sif_sha256": inputs[f"{case}.sif"],
        "mpi_procs": 1,
        "backend": "native",
        "solver": str((tmp_path / "solver").resolve()),
        "inputs_sha256": inputs,
    }), encoding="utf-8")
    assert run.restart_result_is_reusable(
        restart, case, project, backend="native", solver=str(tmp_path / "solver")
    )
    assert not run.restart_result_is_reusable(
        restart, case, project, mpi_procs=2, backend="native",
        solver=str(tmp_path / "solver"),
    )

    # A runtime SIF override is part of solver provenance, even when
    # the generated case SIF itself remains unchanged.
    runtime_sif = manifest.parent / "runtime.sif"
    runtime_sif.write_text("runtime override", encoding="utf-8")
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["runtime_sif"] = str(runtime_sif.relative_to(tmp_path))
    payload["runtime_sif_sha256"] = run.sha256(runtime_sif)
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    assert run.restart_result_is_reusable(restart, case, project, backend="native")
    runtime_sif.write_text("mutated override", encoding="utf-8")
    assert not run.restart_result_is_reusable(restart, case, project, backend="native")


def test_restart_result_is_not_reused_after_project_change(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    restart = tmp_path / "workspace/work" / "meshes" / "mesh" / f"{case}.result"
    restart.parent.mkdir(parents=True)
    restart.write_text("complete restart", encoding="utf-8")
    project = tmp_path / "project.json"
    project.write_text('{"linear_system": "mumps"}', encoding="utf-8")
    inputs = manifest_inputs(tmp_path, case, project)
    manifest = tmp_path / "outputs/results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({
        "exit_code": 0,
        "errors": [],
        "solver_completed": True,
        "runtime_sif_sha256": inputs[f"{case}.sif"],
        "mpi_procs": 1,
        "backend": "native",
        "solver": str((tmp_path / "solver").resolve()),
        "inputs_sha256": inputs,
    }), encoding="utf-8")

    project.write_text('{"linear_system": "iterative_hypre_boomeramg"}', encoding="utf-8")
    assert not run.restart_result_is_reusable(restart, case, project)


def test_restart_requires_matching_mpi_count_and_every_rank(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    mesh = tmp_path / "workspace/work/meshes/mesh"
    mesh.mkdir(parents=True)
    project = tmp_path / "project.json"
    project.write_text("{}", encoding="utf-8")
    inputs = manifest_inputs(tmp_path, case, project)
    for rank in range(4):
        (mesh / f"{case}.result.{rank}").write_bytes(b"restart")
    manifest = tmp_path / "outputs/results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({
        "exit_code": 0, "errors": [], "solver_completed": True,
        "runtime_sif_sha256": inputs[f"{case}.sif"], "mpi_procs": 4, "backend": "elmergpu",
        "solver": "/opt/elmer-gpu/elmer/bin/ElmerSolver_mpi",
        "inputs_sha256": inputs,
    }), encoding="utf-8")

    assert run.restart_result_is_reusable(
        mesh / f"{case}.result", case, project, mpi_procs=4,
        backend="elmergpu", solver="/opt/elmer-gpu/elmer/bin/ElmerSolver_mpi",
    )
    # The historical manifest name and the current logical backend identify
    # the same pinned WSL solver path; the alias itself must not invalidate it.
    assert run.restart_result_is_reusable(
        mesh / f"{case}.result", case, project, mpi_procs=4,
        backend="elmer_wsl", solver="/opt/elmer-gpu/elmer/bin/ElmerSolver_mpi",
    )
    assert run.normalize_backend_name("elmergpu") == "elmer_wsl"
    assert run.normalize_backend_name("native") == "native"
    assert not run.restart_result_is_reusable(
        mesh / f"{case}.result", case, project, mpi_procs=2,
        backend="elmergpu", solver="/opt/elmer-gpu/elmer/bin/ElmerSolver_mpi",
    )
    # A parallel restart set cannot satisfy a serial request, even if its
    # manifest is otherwise successful.
    assert not run.restart_result_is_reusable(
        mesh / f"{case}.result", case, project, mpi_procs=1,
        backend="elmergpu", solver="/opt/elmer-gpu/elmer/bin/ElmerSolver_mpi",
    )
    (mesh / f"{case}.result.2").unlink()
    assert not run.restart_result_is_reusable(
        mesh / f"{case}.result", case, project, mpi_procs=4,
        backend="elmergpu", solver="/opt/elmer-gpu/elmer/bin/ElmerSolver_mpi",
    )


def test_restart_rejects_backend_solver_and_state_mismatch(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    mesh = tmp_path / "workspace/work/meshes/mesh"
    mesh.mkdir(parents=True)
    result = mesh / f"{case}.result"
    result.write_bytes(b"restart")
    state = mesh / f"{case}_1.state"
    state.write_bytes(b"state")
    project = tmp_path / "project.json"
    project.write_text("{}", encoding="utf-8")
    inputs = manifest_inputs(tmp_path, case, project)
    manifest = tmp_path / "outputs/results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({
        "exit_code": 0, "errors": [], "solver_completed": True,
        "runtime_sif_sha256": inputs[f"{case}.sif"], "mpi_procs": 1, "backend": "native",
        "solver": str((tmp_path / "solver").resolve()),
        "inputs_sha256": inputs,
    }), encoding="utf-8")
    assert not run.restart_result_is_reusable(result, case, project, backend="elmergpu")
    assert not run.restart_result_is_reusable(result, case, project, solver="other-solver")
    state.unlink()
    assert not run.restart_result_is_reusable(
        result, case, project, required_paths=[state]
    )


def test_restart_rejects_rebuilt_accessible_native_solver(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run, "ROOT", tmp_path)
    case = "case_steady"
    mesh = tmp_path / "workspace/work/meshes/mesh"
    mesh.mkdir(parents=True)
    result = mesh / f"{case}.result"
    result.write_bytes(b"complete")
    project = tmp_path / "project.json"
    project.write_text("{}", encoding="utf-8")
    inputs = manifest_inputs(tmp_path, case, project)
    solver = tmp_path / "ElmerSolver"
    solver.write_bytes(b"first build")
    manifest = tmp_path / "outputs/results" / case / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({
        "exit_code": 0, "errors": [], "solver_completed": True,
        "runtime_sif_sha256": inputs[f"{case}.sif"],
        "mpi_procs": 1, "backend": "native",
        "solver": str(solver.resolve()),
        "runtime_artifacts_sha256": {"solver": run.sha256(solver)},
        "inputs_sha256": inputs,
    }), encoding="utf-8")
    assert run.restart_result_is_reusable(
        result, case, project, backend="native", solver=str(solver)
    )
    solver.write_bytes(b"updated build")
    assert not run.restart_result_is_reusable(
        result, case, project, backend="native", solver=str(solver)
    )

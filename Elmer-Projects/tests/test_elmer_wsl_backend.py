import json
from types import SimpleNamespace

import pytest

from support import elmer_wsl_backend
from support.build_cases import solver1_block


def test_fast_settings_apply_to_transient_and_are_idempotent(tmp_path):
    project_path = tmp_path / "project.json"
    project_path.write_text(json.dumps({"cases": {"pulse": {
        "template": "pulse", "pulse": {"start": "20[ms]", "duration": "1[ns]"},
        "solver": {}, "extra_solver_keywords": [],
    }}}), encoding="utf-8")

    elmer_wsl_backend.configure_fast_case(project_path, "pulse")
    first = json.loads(project_path.read_text(encoding="utf-8"))
    elmer_wsl_backend.configure_fast_case(project_path, "pulse")
    second = json.loads(project_path.read_text(encoding="utf-8"))

    keywords = second["cases"]["pulse"]["extra_solver_keywords"]
    assert first == second
    assert any(item.startswith("Phase24 BDF1 Step Start Times(2) =") for item in keywords)
    assert any(item.startswith("TES Inner Circuit Newton =") for item in keywords)
    assert len({item.partition("=")[0].strip().casefold() for item in keywords}) == len(keywords)


def test_partition_reuses_complete_partition_and_rejects_incomplete(tmp_path, monkeypatch):
    mesh = tmp_path / "mesh"
    partition = mesh / "partitioning.2"
    partition.mkdir(parents=True)
    for rank in (1, 2):
        for suffix in ("nodes", "elements", "boundary"):
            (partition / f"part.{rank}.{suffix}").write_text("ok", encoding="ascii")
    monkeypatch.setattr(elmer_wsl_backend.subprocess, "run", lambda *a, **kw: pytest.fail("unexpected partition"))

    elmer_wsl_backend.ensure_partition(mesh, 2)
    (partition / "part.2.boundary").unlink()
    with pytest.raises(RuntimeError, match="incomplete 2-rank mesh partition"):
        elmer_wsl_backend.ensure_partition(mesh, 2)


def test_partition_dry_run_does_not_create_or_run(tmp_path, monkeypatch):
    mesh = tmp_path / "mesh"
    monkeypatch.setattr(elmer_wsl_backend.subprocess, "run", lambda *a, **kw: pytest.fail("dry run launched ElmerGrid"))

    elmer_wsl_backend.ensure_partition(mesh, 4, dry_run=True)
    assert not mesh.exists()


@pytest.mark.parametrize(
    ("linear_system", "expected"),
    [
        ("iterative_hypre_pcg_boomeramg", "HYPRE GPU = Logical False"),
        ("iterative_hypre_pcg_boomeramg_gpu", "HYPRE GPU = Logical True"),
    ],
)
def test_cpu_and_gpu_hypre_are_written_to_sif(linear_system, expected):
    sif = "\n".join(solver1_block({
        "linear_system": linear_system,
        "nonlinear_max_iterations": 10,
        "nonlinear_convergence_tolerance": 1e-8,
        "nonlinear_relaxation_factor": 1,
        "steady_state_convergence_tolerance": 1e-8,
    }))

    assert expected in sif


def test_wsl_launcher_keeps_existing_distro_but_uses_canonical_backend(tmp_path, monkeypatch):
    command_calls = []

    def fake_run(command, **kwargs):
        command_calls.append(command)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(elmer_wsl_backend, "ROOT", tmp_path)
    monkeypatch.setattr(elmer_wsl_backend.subprocess, "run", fake_run)
    project = tmp_path / "steady.json"
    record = elmer_wsl_backend.run_wsl(project, "case_steady", "prepared_mesh")

    command = command_calls[0]
    assert command[:5] == ["wsl.exe", "-d", "ElmerGPU", "-u", "root"]
    assert "--backend elmer_wsl" in command[-1]
    assert "/opt/tes-elmer-cpu/elmer/bin/ElmerSolver" in command[-1]
    assert "/opt/tes-elmer-cpu/venv/bin/python" in command[-1]
    assert "ROCM_PATH=" not in command[-1]
    assert record["exit_code"] == 0
    assert record["launcher_log"].endswith(".wsl_launcher.log")

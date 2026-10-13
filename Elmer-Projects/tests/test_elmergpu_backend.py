import json

import pytest

from support import elmergpu_backend
from support.build_cases import solver1_block


def test_fast_settings_apply_to_transient_and_are_idempotent(tmp_path):
    project_path = tmp_path / "project.json"
    project_path.write_text(json.dumps({"cases": {"pulse": {
        "template": "pulse", "pulse": {"start": "20[ms]", "duration": "1[ns]"},
        "solver": {}, "extra_solver_keywords": [],
    }}}), encoding="utf-8")

    elmergpu_backend.configure_fast_case(project_path, "pulse")
    first = json.loads(project_path.read_text(encoding="utf-8"))
    elmergpu_backend.configure_fast_case(project_path, "pulse")
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
    monkeypatch.setattr(elmergpu_backend.subprocess, "run", lambda *a, **kw: pytest.fail("unexpected partition"))

    elmergpu_backend.ensure_partition(mesh, 2)
    (partition / "part.2.boundary").unlink()
    with pytest.raises(RuntimeError, match="incomplete 2-rank mesh partition"):
        elmergpu_backend.ensure_partition(mesh, 2)


def test_partition_dry_run_does_not_create_or_run(tmp_path, monkeypatch):
    mesh = tmp_path / "mesh"
    monkeypatch.setattr(elmergpu_backend.subprocess, "run", lambda *a, **kw: pytest.fail("dry run launched ElmerGrid"))

    elmergpu_backend.ensure_partition(mesh, 4, dry_run=True)
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

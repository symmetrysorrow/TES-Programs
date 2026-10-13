from types import SimpleNamespace

from support import run_phase24_gpu_elmergpu as legacy_runner


def test_legacy_sync_invokes_parameter_generator_and_writes_log(monkeypatch, tmp_path):
    project = tmp_path / "project.json"
    project.write_text("{}", encoding="utf-8")
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout="synced", stderr="")

    monkeypatch.setattr(legacy_runner.subprocess, "run", fake_run)
    legacy_runner.sync(project)

    assert calls[0][0][1].endswith("sync_elmer_parameters.py")
    assert calls[0][0][-1] == str(project)
    assert calls[0][1]["cwd"] == legacy_runner.ROOT
    assert (tmp_path / "sync.log").read_text(encoding="utf-8") == "synced"

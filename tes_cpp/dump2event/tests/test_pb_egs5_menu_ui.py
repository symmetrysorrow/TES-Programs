"""Small, independent tests for the interactive opt-in Pb physics menu."""

import json
from pathlib import Path
from types import SimpleNamespace

import h5py
import pytest

from PoST_Simulations.lib.pb_egs5_ui import _inspect_v2, run_pb_egs5_ui


class Answer:
    def __init__(self, value):
        self.value = value

    def ask(self):
        return self.value


class FakeQuestionary:
    def __init__(self, paths, *, confirmed=True):
        self.paths = iter(paths)
        self.confirmed = confirmed
        self.default_paths = []

    def path(self, message, *, default=""):
        self.default_paths.append(default)
        return Answer(next(self.paths))

    def confirm(self, message, *, default=False):
        assert default is False
        return Answer(self.confirmed)


def _inputs(tmp_path):
    raw = tmp_path / "input.dumpall"
    raw.write_text("not used by this mock")
    phits = tmp_path / "source.inp"
    phits.write_text("input")
    config = tmp_path / "physics.json"
    config.write_text(json.dumps({"phits_input": "source.inp"}))
    table = tmp_path / "stopping.csv"
    table.write_text("private-table-test")
    output = tmp_path / "event_pb_egs5_v2.h5"
    return raw, config, table, output


def _fake_hdf5(path):
    with h5py.File(path, "w") as handle:
        handle.attrs["format"] = "tes-dump2event"
        handle.attrs["format_version"] = 2
        events = handle.create_group("events")
        events.create_dataset("event_id", data=[1, 2])
        events.create_dataset("unresolved_count", data=[0, 0])
        events.create_dataset("complete", data=[True, True])
        particles = handle.create_group("particles")
        particles.create_dataset("particle_id", data=[1, 2])
        tracks = handle.create_group("tracks")
        tracks.create_dataset("ncol", data=[14, 14])
        deposits = handle.create_group("deposits")
        deposits.create_dataset("addition_id", data=[b"1", b"2"])
        prov = handle.create_group("provenance")
        prov.attrs["phits_version"] = "3.341"
        prov.attrs["algorithm_version"] = "test"


def test_pb_ui_routes_explicit_api_and_does_not_run_pulse(tmp_path, monkeypatch):
    import PoST_Simulations.lib.pb_egs5_ui as ui
    raw, config, table, output = _inputs(tmp_path)
    monkeypatch.setattr(ui, "_bundled_table_present", lambda: False)
    fake = FakeQuestionary([raw, table, output])
    calls = []

    def convert(*args, **kwargs):
        calls.append((args, kwargs))
        _fake_hdf5(output)
        return SimpleNamespace(output_path=output)

    result = run_pb_egs5_ui(fake, convert)
    assert result.output_path == output
    assert calls == [((raw.resolve(), output), {
        "physics": "pb-egs5",
        "stopping_table": table.resolve(),
    })]
    assert _inspect_v2(output)["deposits"] == 2
    assert "event_pb_egs5_v2.h5" in fake.default_paths[-1]


def test_pb_ui_rejects_existing_output_without_running_converter(tmp_path, monkeypatch):
    import PoST_Simulations.lib.pb_egs5_ui as ui
    raw, config, table, output = _inputs(tmp_path)
    monkeypatch.setattr(ui, "_bundled_table_present", lambda: False)
    output.write_text("do not replace")
    with pytest.raises(FileExistsError):
        run_pb_egs5_ui(FakeQuestionary([raw, table, output]),
                       lambda *args, **kwargs: pytest.fail("unexpected converter call"))
    assert output.read_text() == "do not replace"


def test_pb_ui_leaves_phits_pair_validation_to_fail_closed_converter(tmp_path, monkeypatch):
    import PoST_Simulations.lib.pb_egs5_ui as ui
    raw, _, table, output = _inputs(tmp_path)
    monkeypatch.setattr(ui, "_bundled_table_present", lambda: False)
    def reject_pair(*_args, **_kwargs):
        raise ValueError("no matching input/log pair")
    with pytest.raises(ValueError, match="no matching input/log pair"):
        run_pb_egs5_ui(FakeQuestionary([raw, table, output]), reject_pair)


def test_pb_ui_cancelled_before_conversion(tmp_path, monkeypatch):
    import PoST_Simulations.lib.pb_egs5_ui as ui
    raw, config, table, output = _inputs(tmp_path)
    monkeypatch.setattr(ui, "_bundled_table_present", lambda: False)
    assert run_pb_egs5_ui(
        FakeQuestionary([raw, table, output], confirmed=False),
        lambda *args, **kwargs: pytest.fail("unexpected converter call"),
    ) is None
    assert not output.exists()


def test_pb_ui_rejects_wrong_schema_and_unresolved_hdf5(tmp_path):
    output = tmp_path / "bad.h5"
    with h5py.File(output, "w") as handle:
        handle.attrs["format"] = "tes-dump2event"
        handle.attrs["format_version"] = 1
    with pytest.raises(ValueError, match="HDF5 v2"):
        _inspect_v2(output)
    output.unlink()
    _fake_hdf5(output)
    with h5py.File(output, "r+") as handle:
        handle["events"]["unresolved_count"][0] = 1
    with pytest.raises(ValueError, match="未解決"):
        _inspect_v2(output)



def test_pb_ui_uses_bundled_table_without_picker_when_installed(tmp_path, monkeypatch):
    import PoST_Simulations.lib.pb_egs5_ui as ui
    raw, config, table, output = _inputs(tmp_path)
    assert ui._bundled_table_present() is True
    dialog = FakeQuestionary([raw, output])
    kwargs_seen = []

    def convert(*args, **kwargs):
        kwargs_seen.append(kwargs)
        _fake_hdf5(output)
        return SimpleNamespace(output_path=output)

    run_pb_egs5_ui(dialog, convert)
    assert len(dialog.default_paths) == 2  # only Dumpall and output
    assert "stopping_table" not in kwargs_seen[0]
    assert "physics_config" not in kwargs_seen[0]


def test_pb_ui_detects_scikit_build_installed_table(monkeypatch, tmp_path):
    import PoST_Simulations.lib.pb_egs5_ui as ui
    prefix = tmp_path / "prefix"
    module_dir = prefix / "python" / "tes_cpp" / "lib"
    module_dir.mkdir(parents=True)
    monkeypatch.setattr(ui, "__file__", str(module_dir / "pb_egs5_ui.py"))
    table = prefix / "tes_cpp" / "bin" / "pb_egs5_phits3341.csv"
    table.parent.mkdir(parents=True)
    table.write_bytes(b"software test placeholder")
    class FakeDigest:
        def hexdigest(self):
            return "20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9"
    monkeypatch.setattr(ui.hashlib, "sha256", lambda _data: FakeDigest())
    assert ui._bundled_table_present() is True


def test_modified_private_table_environment_setting_fails_closed(tmp_path, monkeypatch):
    import PoST_Simulations.lib.pb_egs5_ui as ui
    _, _, table, _ = _inputs(tmp_path)
    fake_ui = tmp_path / "isolated" / "lib" / "pb_egs5_ui.py"
    fake_ui.parent.mkdir(parents=True)
    fake_ui.write_text("", encoding="utf-8")
    monkeypatch.setattr(ui, "__file__", str(fake_ui))
    table.write_bytes(b"private but not valid pinned bytes")
    monkeypatch.setenv("TES_CPP_PB_HATCH_TABLE", str(table))
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        ui._bundled_table_present()

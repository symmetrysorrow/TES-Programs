"""Numerical and end-to-end regressions for Pb/EGS5 HDF5 v2 -> Pulse_Ms."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tes_cpp" / "python"))
sys.path.insert(0, str(ROOT / "PoST_Simulations"))

from tes_cpp.event2pulse import (
    SEGMENT_ALLOCATION_MODEL,
    allocate_deposit,
    event_format_version,
    iter_v2_block_energies,
)
from tes_cpp.event_hdf5 import EventWriter
from tes_cpp.event_hdf5_v2 import EventV2Writer


PROVENANCE = {
    "phits_version": "3.341",
    "dumpall_sha256": "a" * 64,
    "physics_config_sha256": "b" * 64,
    "stopping_table_sha256": "c" * 64,
    "material": "Pb",
    "medium": "2",
    "cutoff_MeV": "0.001",
    "algorithm_version": "test",
    "position_semantics": "point score and segment transport chord are distinct",
}


def _point(identifier, energy, x):
    return {
        "addition_id": identifier, "particle_id": -1, "process": "point",
        "region_id": 1000, "energy_MeV": energy, "kind": "POINT",
        "score_xyz_cm": (x, 0, 0), "source_record_id": identifier,
    }


def _segment(identifier, energy, x0, x1, *, score=0.0):
    return {
        "addition_id": identifier, "particle_id": -1, "process": "egs5-continuous",
        "region_id": 1000, "energy_MeV": energy, "kind": "SEGMENT",
        "start_xyz_cm": (x0, 0, 0), "end_xyz_cm": (x1, 0, 0),
        "score_xyz_cm": (score, 0, 0), "source_record_id": identifier,
    }


def _write_v2(path, deposits, *, complete=True):
    with EventV2Writer(path, PROVENANCE) as writer:
        writer.append_event(42, [], deposits, complete=complete)
        writer.append_event(43, [], [], complete=True, all_absorption="undetermined")


def test_segment_crossing_blocks_conserves_energy_and_ignores_score_point():
    energies = np.zeros(4)
    allocate_deposit(energies, _segment("a", 1.0, -0.15, 0.05, score=0.19),
                     length_mm=4)
    np.testing.assert_allclose(energies, [0.25, 0.5, 0.25, 0], atol=1e-14)
    assert energies.sum() == pytest.approx(1)
    mirrored = np.zeros(4)
    allocate_deposit(mirrored, _segment("a", 1.0, -0.15, 0.05),
                     length_mm=4, mirror_x=True)
    np.testing.assert_allclose(mirrored, energies[::-1], atol=1e-14)


def test_point_boundary_and_parallel_segment_and_outside_geometry():
    energy = np.zeros(4)
    allocate_deposit(energy, _point("p", 0.1, 0.0), length_mm=4)
    allocate_deposit(energy, _segment("s", 0.2, 0.0, 0.0), length_mm=4)
    np.testing.assert_allclose(energy, [0, 0, 0.3, 0])
    with pytest.raises(ValueError, match="outside absorber"):
        allocate_deposit(energy, _segment("bad", 0.5, -0.3, -0.1), length_mm=4)
    with pytest.raises(ValueError, match="score position"):
        allocate_deposit(energy, dict(_point("bad2", 0.5, 0), score_xyz_cm=None), length_mm=4)
    with pytest.raises(ValueError, match="unresolved"):
        allocate_deposit(energy, dict(_point("bad3", 0.5, 0), kind="UNRESOLVED"), length_mm=4)


def test_v2_streaming_preserves_events_zero_rows_and_energy(tmp_path):
    path = tmp_path / "event.h5"
    _write_v2(path, [_segment("seg", 0.4, -0.2, 0.2),
                     _point("point", 0.6, -0.19)])
    assert event_format_version(path) == 2
    rows = list(iter_v2_block_energies(path, 4, 4))
    assert [event_id for event_id, _ in rows] == [42, 43]
    np.testing.assert_allclose(rows[0][1], [0.7, 0.1, 0.1, 0.1], atol=1e-14)
    np.testing.assert_array_equal(rows[1][1], np.zeros(4))


def test_incomplete_and_unresolved_are_not_coerced_to_zero(tmp_path):
    path = tmp_path / "bad.h5"
    unresolved = {
        "addition_id": "u", "particle_id": -1, "process": "unknown",
        "region_id": 1000, "energy_MeV": None, "kind": "UNRESOLVED",
        "source_record_id": "u",
    }
    _write_v2(path, [unresolved], complete=False)
    with pytest.raises(ValueError, match="incomplete or has unresolved"):
        list(iter_v2_block_energies(path, 4, 4))


def test_reader_rejects_wrong_deposit_event_identity(tmp_path):
    path = tmp_path / "altered.h5"
    _write_v2(path, [_point("p", 0.2, 0)])
    with h5py.File(path, "r+") as handle:
        handle["deposits/event_id"][0] = 99
    with pytest.raises(ValueError, match="event/deposit ID mismatch"):
        list(iter_v2_block_energies(path, 4, 4))


def _prepare_pulse_run(tmp_path, monkeypatch, event_path):
    import PoST_Simulation as sim

    para = {
        "position": [1], "n_abs": 4, "length": 4.0, "E": 1000.0,
        "UseMirroredEvents": False,
    }
    (tmp_path / "input.json").write_text(json.dumps(para), encoding="utf-8")
    monkeypatch.setattr(sim, "output", str(tmp_path))
    monkeypatch.setattr(sim, "event_root", str(event_path))
    monkeypatch.setattr(sim, "use_mirrored_events", None)
    reference = {
        i: {"time": np.array([0.0, 0.001]),
            "ch0": np.array([float(i), 0.0]),
            "ch1": np.array([0.0, float(i) * 2.0])}
        for i in range(1, 5)
    }
    monkeypatch.setattr(sim, "LoadPulses", lambda: reference)
    return sim


def test_pulse_ms_reads_direct_v2_and_records_approximation(tmp_path, monkeypatch):
    event_path = tmp_path / "event_pb_egs5_v2.h5"
    _write_v2(event_path, [_segment("seg", 0.4, -0.2, 0.2),
                           _point("point", 0.6, -0.19)])
    sim = _prepare_pulse_run(tmp_path, monkeypatch, event_path)
    sim.Pulse_Ms()
    with h5py.File(tmp_path / "1" / "pulse_MS.h5") as output:
        assert list(output["event_id"].asstr()) == ["42", "43"]
        np.testing.assert_allclose(output["ch0"][0], [1.6, 0.0], atol=1e-14)
        np.testing.assert_allclose(output["ch1"][0], [0.0, 3.2], atol=1e-14)
        np.testing.assert_array_equal(output["ch0"][1], np.zeros(2))
        assert output.attrs["source_event_format_version"] == 2
        assert output.attrs["energy_allocation_model"] == SEGMENT_ALLOCATION_MODEL


def test_pulse_ms_fails_atomically_for_unresolved_v2(tmp_path, monkeypatch):
    event_path = tmp_path / "bad.h5"
    _write_v2(event_path, [{"addition_id": "u", "particle_id": -1,
        "process": "unknown", "region_id": 1000,
        "energy_MeV": None, "kind": "UNRESOLVED", "source_record_id": "u"}],
        complete=False)
    sim = _prepare_pulse_run(tmp_path, monkeypatch, event_path)
    existing = tmp_path / "1" / "pulse_MS.h5"
    existing.parent.mkdir()
    existing.write_bytes(b"do not replace")
    with pytest.raises(ValueError, match="incomplete"):
        sim.Pulse_Ms()
    assert existing.read_bytes() == b"do not replace"
    assert not list(existing.parent.glob(".pulse_MS.*.h5"))


def test_old_flat_v1_events_remain_usable(tmp_path, monkeypatch):
    path = tmp_path / "old.h5"
    with EventWriter(path) as writer:
        writer.append(7, {"1": {"x_deposit": [-0.15], "y_deposit": [0],
                               "z_deposit": [0], "E_deposit": [1.0],
                               "x": [], "y": [], "z": [], "E": []}})
    sim = _prepare_pulse_run(tmp_path, monkeypatch, path)
    assert event_format_version(path) == 1
    sim.Pulse_Ms()
    with h5py.File(tmp_path / "1" / "pulse_MS.h5") as output:
        assert list(output["event_id"].asstr()) == ["7"]
        np.testing.assert_allclose(output["ch0"][0], [1.0, 0.0])
        assert output.attrs["energy_allocation_model"] == "legacy_point_block"


def test_interactive_action_chain_passes_new_v2_path(tmp_path, monkeypatch):
    from types import SimpleNamespace
    import PoST_Simulation as sim

    (tmp_path / "input.json").write_text(json.dumps({"position": [2]}), encoding="utf-8")
    monkeypatch.setattr(sim, "output", str(tmp_path))
    monkeypatch.setattr(sim, "event_root", None)
    selected = tmp_path / "event_v2.h5"
    called = []

    monkeypatch.setattr(sim, "Dump2Event", lambda: SimpleNamespace(output_path=selected))
    monkeypatch.setattr(sim, "Pulse_Ms", lambda: called.append(sim.event_root))
    sim._run_selected_actions(["dump_to_event", "pulse_ms"])
    assert called == [str(selected)]


def test_interactive_action_chain_rejects_ambiguous_multi_position(tmp_path, monkeypatch):
    import PoST_Simulation as sim

    (tmp_path / "input.json").write_text(json.dumps({"position": [1, 2]}), encoding="utf-8")
    monkeypatch.setattr(sim, "output", str(tmp_path))
    monkeypatch.setattr(sim, "Dump2Event", lambda: pytest.fail("must not convert multiple positions"))
    with pytest.raises(ValueError, match="1位置"):
        sim._run_selected_actions(["dump_to_event", "pulse_ms"])


def test_root_directory_accepts_new_default_converter_filename(tmp_path):
    import PoST_Simulation as sim

    directory = tmp_path / "events" / "1"
    directory.mkdir(parents=True)
    path = directory / "event_pb_egs5_v2.h5"
    _write_v2(path, [_point("p", 1, 0)])
    assert sim._resolve_ms_event_path(tmp_path / "events", 1, allow_single_file=False) == path

    (directory / "event.h5").write_bytes(b"ambiguous")
    with pytest.raises(ValueError, match="ambiguous"):
        sim._resolve_ms_event_path(tmp_path / "events", 1, allow_single_file=False)

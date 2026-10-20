from __future__ import annotations

import h5py
import pytest

from tes_cpp.event_hdf5 import iter_deposits, iter_events, require_legacy_point_events
from tes_cpp.event_hdf5_v2 import EventV2Writer, iter_deposits as iter_v2_deposits


PROVENANCE = {
    "phits_version": "3.341",
    "dumpall_sha256": "a" * 64,
    "physics_config_sha256": "b" * 64,
    "stopping_table_sha256": "c" * 64,
    "material": "Pb",
    "medium": "2",
    "cutoff_MeV": "0.001",
    "algorithm_version": "test-1",
    "position_semantics": "PHITS score point and transport chord stored separately",
}


def _point():
    return {"addition_id": "add-1", "particle_id": 7, "process": "test-point",
            "region_id": 1000, "energy_MeV": 0.25, "kind": "POINT",
            "score_xyz_cm": (1.0, 2.0, 3.0), "source_record_id": "line-12",
            "position_semantics": "validated PHITS score point"}


def _segment():
    return {"addition_id": "add-2", "particle_id": 7, "process": "continuous-loss",
            "region_id": 1000, "energy_MeV": 0.5, "kind": "SEGMENT",
            "start_xyz_cm": (0.0, 1.0, 2.0), "end_xyz_cm": (1.0, 2.0, 3.0),
            "score_xyz_cm": (0.5, 1.5, 2.5), "source_record_id": "line-13",
            "position_semantics": "transport segment; score point is representative"}


def test_v2_round_trip_preserves_kinds_offsets_and_metadata(tmp_path):
    path = tmp_path / "event.h5"
    particles = [{"particle_id": 7, "parent_particle_id": None, "ityp": 12,
                  "tracks": [{"source_record_id": "line-12", "ncol": 14,
                              "xyz_cm": (1.0, 2.0, 3.0), "energy_MeV": 0.75,
                              "time": 2.5, "region_id": 1000}]}]
    with EventV2Writer(path, PROVENANCE) as writer:
        result = writer.append_event(44, particles, [_point(), _segment()], complete=True,
                                     all_absorption="confirmed")
        assert result["deposit_total_MeV"] == 0.75
        writer.append_event(45, [{"particle_id": 2, "parent_particle_id": 7, "ityp": 12}], [],
                            complete=True, all_absorption="undetermined")
    assert list(iter_v2_deposits(path)) == [
        (44, {"addition_id": "add-1", "particle_id": 7, "process": "test-point",
              "region_id": 1000, "energy_MeV": 0.25, "kind": "POINT",
              "start_xyz_cm": (0.0, 0.0, 0.0), "end_xyz_cm": (0.0, 0.0, 0.0),
              "score_xyz_cm": (1.0, 2.0, 3.0), "score_position_valid": True,
              "source_record_id": "line-12", "position_semantics": "validated PHITS score point"}),
        (44, {"addition_id": "add-2", "particle_id": 7, "process": "continuous-loss",
              "region_id": 1000, "energy_MeV": 0.5, "kind": "SEGMENT",
              "start_xyz_cm": (0.0, 1.0, 2.0), "end_xyz_cm": (1.0, 2.0, 3.0),
              "score_xyz_cm": (0.5, 1.5, 2.5), "score_position_valid": True,
              "source_record_id": "line-13", "position_semantics": "transport segment; score point is representative"}),
    ]
    assert list(iter_events(path)) == [(44, {"particle_range": (0, 1), "deposit_range": (0, 2),
                                               "deposit_total_MeV": 0.75, "complete": True,
                                               "all_absorption": "confirmed", "unresolved_count": 0}),
                                        (45, {"particle_range": (1, 2), "deposit_range": (2, 2),
                                               "deposit_total_MeV": 0.0, "complete": True,
                                               "all_absorption": "undetermined", "unresolved_count": 0})]
    with h5py.File(path, "r") as file:
        assert file.attrs["format_version"] == 2
        assert list(file["events/particle_offsets"][:]) == [0, 1, 2]
        assert list(file["events/deposit_offsets"][:]) == [0, 2, 2]
        assert list(file["particles/track_offsets"][:]) == [0, 1, 1]
        assert file["tracks/source_record_id"][0].decode() == "line-12"
        assert file["provenance"].attrs["phits_version"] == "3.341"


def test_unresolved_energy_is_not_encoded_as_known_zero(tmp_path):
    path = tmp_path / "event.h5"
    unresolved = {"addition_id": "add-u", "particle_id": -1, "process": "unknown",
                  "region_id": 1000, "energy_MeV": None, "kind": "UNRESOLVED",
                  "source_record_id": "line-99", "position_semantics": "unknown"}
    with EventV2Writer(path, PROVENANCE) as writer:
        writer.append_event(1, [], [unresolved], complete=False)
    row = list(iter_v2_deposits(path))[0][1]
    assert row["energy_MeV"] is None
    assert row["kind"] == "UNRESOLVED"
    assert list(iter_events(path))[0][1]["unresolved_count"] == 1


def test_failed_write_does_not_replace_existing_output(tmp_path):
    path = tmp_path / "event.h5"
    path.write_bytes(b"keep-existing")
    with pytest.raises(ValueError, match="nonnegative"):
        with EventV2Writer(path, PROVENANCE) as writer:
            bad = dict(_point(), energy_MeV=-0.1)
            writer.append_event(1, [{"particle_id": 7, "ityp": 12}], [bad], complete=True)
    assert path.read_bytes() == b"keep-existing"
    assert not list(tmp_path.glob(".event.h5.*.tmp"))


def test_unknown_version_is_rejected_before_schema_guess(tmp_path):
    path = tmp_path / "future.h5"
    with h5py.File(path, "w") as file:
        file.attrs["format"] = "tes-dump2event"
        file.attrs["format_version"] = 77
        file.create_group("events")
    with pytest.raises(ValueError, match="format_version=77"):
        list(iter_events(path))
    with pytest.raises(ValueError, match="format_version=77"):
        list(iter_deposits(path))
    with pytest.raises(ValueError, match="format_version=77"):
        require_legacy_point_events(path)


def test_legacy_cpp_hierarchy_is_recognized_only_without_version(tmp_path):
    path = tmp_path / "legacy-cpp.h5"
    with h5py.File(path, "w") as file:
        file.attrs["format"] = "tes-dump2event"
        particle = file.create_group("events/11/2")
        particle.attrs["ityp"] = 12
        for name, values in {"x": [1.0], "y": [2.0], "z": [3.0], "E": [4.0],
                             "x_deposit": [1.5], "y_deposit": [2.5],
                             "z_deposit": [3.5], "E_deposit": [0.4]}.items():
            particle.create_dataset(name, data=values)
    event_id, particles = next(iter_events(path))
    assert event_id == "11"
    assert particles["2"]["ityp"] == 12
    assert particles["2"]["E_deposit"] == [0.4]


def test_v2_is_not_admitted_to_legacy_pulse_reader(tmp_path):
    path = tmp_path / "event.h5"
    with EventV2Writer(path, PROVENANCE):
        pass
    with pytest.raises(ValueError, match="SEGMENT allocation"):
        require_legacy_point_events(path)
    assert list(iter_events(path)) == []

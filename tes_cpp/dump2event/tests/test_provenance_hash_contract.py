import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import build_physical_hdf5_v2 as build
from tes_cpp.event_hdf5_v2 import EventV2Writer


SOURCE = """[Parameters]
bitrseed = 0101101011001010010111000011100101101001110000111010101010101010
icntl=0
maxbch=1
negs=1
e-mode=2
emin(12)=0.001
emin(13)=0.001
emin(14)=0.001
[Source]
totfact=1.0
wgt=1.0
s-type=1
proj=photon
e-type=8
ne=1
0.511 1
[Material]
mat[2] Pb 1
[Surface]
100 rpp -1 1 -0.05 0.05 -0.1 0
[Cell]
1000 2 -9.78 -100
[t-deposit]
mesh=reg
reg=1000
part=all
deposit=0
unit=3
output=deposit
e-type=2
ne=150
emin=0.0
emax=1.5
"""

SETTINGS = {
    "phits_version": "3.341", "phits_input": "input.inp", "region_id": 1000,
    "region_material": "Pb", "dumpall_material_id": 2, "egs5_medium": 2,
    "electron_cutoff_MeV": 0.001, "positron_cutoff_MeV": 0.001,
    "weight_mode": "unit", "energy_unit": "MeV",
    "deposit_mode": "t-deposit-output-deposit",
}


def _files(tmp_path, json_text=None):
    source = tmp_path / "input.inp"
    source.write_text(SOURCE, encoding="utf-8")
    config = tmp_path / "config.json"
    config.write_text(json_text or json.dumps(SETTINGS), encoding="utf-8")
    dumpall = tmp_path / "dumpall.dat"
    dumpall.write_bytes(b"raw dump")
    table = tmp_path / "pb.csv"
    table.write_bytes(b"stopping coefficients")
    return source, config, dumpall, table


def test_mono_gamma_source_guard_accepts_only_supported_energy_and_explicit_photon(tmp_path):
    path = tmp_path / "trace.inp"
    path.write_text(SOURCE, encoding="utf-8")
    build.validate_source_input(path)
    with pytest.raises(ValueError, match="proj=photon"):
        build.source_definition(["totfact=1", "e-type=8", "ne=1", "0.511 1"])
    with pytest.raises(ValueError, match="0.511/2.0 MeV"):
        build.source_definition(["totfact=1", "proj=photon", "e-type=8", "ne=1", "0.662 1"])


def test_raw_contract_hash_is_distinct_from_canonical_resolved_hash(tmp_path):
    source, config_a, dumpall, table = _files(tmp_path, '{"phits_input":"input.inp", "region_id":1000}\n')
    settings_a = dict(SETTINGS)
    meta_a = build.provenance_metadata(settings_a, config_a, source, dumpall, table)
    config_b = tmp_path / "config_b.json"
    config_b.write_text('{\n  "region_id": 1000,\n  "phits_input": "./input.inp"\n}\n', encoding="utf-8")
    meta_b = build.provenance_metadata(dict(reversed(list(SETTINGS.items()))), config_b, source, dumpall, table)
    assert meta_a["physics_contract_sha256"] != meta_b["physics_contract_sha256"]
    assert meta_a["resolved_physics_sha256"] == meta_b["resolved_physics_sha256"]
    assert meta_a["physics_config_sha256"] == meta_a["resolved_physics_sha256"]


def test_resolved_physics_and_file_hashes_change_at_the_right_layer(tmp_path):
    source, config, dumpall, table = _files(tmp_path)
    base = build.provenance_metadata(dict(SETTINGS), config, source, dumpall, table)
    changed_settings = dict(SETTINGS, electron_cutoff_MeV=0.002)
    changed = build.provenance_metadata(changed_settings, config, source, dumpall, table)
    assert changed["resolved_physics_sha256"] != base["resolved_physics_sha256"]
    assert changed["physics_contract_sha256"] == base["physics_contract_sha256"]
    changed_source = tmp_path / "changed.inp"
    changed_source.write_text(SOURCE + "# input byte change\n", encoding="utf-8")
    changed_input = build.provenance_metadata(dict(SETTINGS), config, changed_source, dumpall, table)
    assert changed_input["phits_input_sha256"] != base["phits_input_sha256"]
    assert changed_input["resolved_physics_sha256"] == base["resolved_physics_sha256"]
    changed_table = tmp_path / "other.csv"
    changed_table.write_bytes(b"different coefficients")
    changed_table_meta = build.provenance_metadata(dict(SETTINGS), config, source, dumpall, changed_table)
    assert changed_table_meta["stopping_table_sha256"] != base["stopping_table_sha256"]


def test_missing_input_or_invalid_v2_hash_fails_closed(tmp_path):
    source, config, dumpall, table = _files(tmp_path)
    with pytest.raises(FileNotFoundError):
        build.provenance_metadata(dict(SETTINGS), config, source.with_name("missing.inp"), dumpall, table)
    provenance = {
        "phits_version": "3.341", "dumpall_sha256": "a" * 64,
        "physics_config_sha256": "b" * 64, "resolved_physics_sha256": "c" * 64,
        "physics_contract_sha256": "d" * 63, "stopping_table_sha256": "e" * 64,
        "material": "Pb", "medium": "2", "cutoff_MeV": "0.001",
        "algorithm_version": "test", "position_semantics": "test",
        "provenance_schema_revision": 2,
    }
    with pytest.raises(ValueError, match="64-character SHA-256"):
        EventV2Writer(tmp_path / "bad.h5", provenance)


def test_revision_two_requires_matching_compatibility_alias_and_roundtrips(tmp_path):
    resolved = "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a"
    contract = "b" * 64
    provenance = {
        "phits_version": "3.341", "dumpall_sha256": "c" * 64,
        "physics_config_sha256": resolved, "resolved_physics_sha256": resolved,
        "physics_contract_sha256": contract, "phits_input_sha256": "d" * 64,
        "source_input_sha256": "d" * 64, "stopping_table_sha256": "e" * 64,
        "provenance_schema_revision": 2, "material": "Pb", "medium": "2",
        "cutoff_MeV": "0.001", "algorithm_version": "test",
        "position_semantics": "test", "resolved_physics_json": "{}",
    }
    path = tmp_path / "revision2.h5"
    with EventV2Writer(path, provenance):
        pass
    import h5py
    with h5py.File(path, "r") as h5:
        attrs = h5["provenance"].attrs
        assert attrs["physics_config_sha256"] == attrs["resolved_physics_sha256"]
        assert attrs["physics_contract_sha256"] == contract
        assert attrs["provenance_schema_revision"] == 2

    invalid = dict(provenance, physics_config_sha256="f" * 64)
    with pytest.raises(ValueError, match="compatibility alias"):
        EventV2Writer(tmp_path / "invalid.h5", invalid)


def test_executable_provenance_is_recorded_from_exact_files(tmp_path):
    source, config, dumpall, table = _files(tmp_path)
    converter = tmp_path / "dump2event.exe"
    reconstructor = tmp_path / "reconstructor.exe"
    exporter = tmp_path / "exporter.exe"
    converter.write_bytes(b"converter")
    reconstructor.write_bytes(b"reconstructor")
    exporter.write_bytes(b"exporter")
    metadata = build.provenance_metadata(
        dict(SETTINGS), config, source, dumpall, table, converter, reconstructor, exporter
    )
    assert metadata["converter_executable_sha256"] == build.sha256(converter)
    assert metadata["reconstructor_executable_sha256"] == build.sha256(reconstructor)
    assert metadata["transport_exporter_sha256"] == build.sha256(exporter)
    assert metadata["driver_sha256"] == build.sha256(Path(build.__file__))


def test_explicit_phits_executable_hash_is_validated_and_recorded(tmp_path):
    source, config, dumpall, table = _files(tmp_path)
    config.write_text(json.dumps(dict(SETTINGS, phits_executable_sha256="a" * 64)), encoding="utf-8")
    metadata = build.provenance_metadata(dict(SETTINGS), config, source, dumpall, table)
    assert metadata["phits_executable_sha256"] == "a" * 64
    config.write_text(json.dumps(dict(SETTINGS, phits_executable_sha256="not-a-hash")), encoding="utf-8")
    with pytest.raises(ValueError, match="phits_executable_sha256"):
        build.provenance_metadata(dict(SETTINGS), config, source, dumpall, table)

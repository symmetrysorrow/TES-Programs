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
    "electron_cutoff_MeV": 0.001, "positron_cutoff_MeV": 0.001, "esmin_MeV": 1.0e-3,
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


@pytest.mark.parametrize("energy", [0.0005, 0.511, 0.662, 1.173, 1.5, 1.999999, 2.0])
def test_mono_gamma_accepts_positive_energy_at_or_below_2mev(tmp_path, energy):
    source = SOURCE.replace("0.511 1", f"{energy} 1")
    path = tmp_path / "trace.inp"
    path.write_text(source, encoding="utf-8")
    build.validate_source_input(path)
    parsed = build.source_definition(
        ["totfact=1", "wgt=1", "proj=photon", "e-type=8", "ne=1", f"{energy} 1"]
    )
    assert parsed["source_particle"] == "photon"
    assert parsed["source_energies_MeV"] == [energy]


@pytest.mark.parametrize("bad_energy", [
    "0", "-0.1", "2.000001", "3", "nan", "inf", "-inf", "1e309", "not-a-number"
])
def test_mono_gamma_rejects_nonphysical_or_above_2mev_energy(bad_energy):
    source = ["totfact=1", "proj=photon", "e-type=8", "ne=1", f"{bad_energy} 1"]
    with pytest.raises(ValueError, match="energy"):
        build.source_definition(source)


def test_source_guards_energy_line_counts_and_unit_weights():
    valid = ["totfact=1", "proj=photon", "e-type=8", "ne=2", "0.662 1", "1.8 1"]
    parsed = build.source_definition(valid)
    assert parsed["source_mode"] == "discrete-gamma-lines"
    assert parsed["source_energies_MeV"] == [0.662, 1.8]
    with pytest.raises(ValueError, match="proj=photon"):
        build.source_definition(["totfact=1", "e-type=8", "ne=1", "0.662 1"])
    with pytest.raises(ValueError, match="ne"):
        build.source_definition(["totfact=1", "proj=photon", "e-type=8", "ne=1",
                                 "0.662 1", "1.8 1"])
    for bad_weight in ("0", "0.5", "nan", "inf"):
        with pytest.raises(ValueError, match="unit-intensity"):
            build.source_definition(["totfact=1", "proj=photon", "e-type=8",
                                     "ne=1", f"0.662 {bad_weight}"])
    for bad_totfact in ("0", "2", "nan", "inf"):
        with pytest.raises(ValueError, match="unit"):
            build.source_definition(["totfact=" + bad_totfact, "proj=photon",
                                     "e-type=8", "ne=1", "0.662 1"])
    with pytest.raises(ValueError, match="unit"):
        build.source_definition(["totfact=1", "wgt=0.5", "proj=photon",
                                 "e-type=8", "ne=1", "0.662 1"])
    with pytest.raises(ValueError, match="e-type=8"):
        build.source_definition(["totfact=1", "proj=photon", "e-type=2",
                                 "ne=1", "0.662 1"])


def test_legacy_photon_modes_keep_provenance_names():
    src = lambda rows: build.source_definition(rows)
    assert src(["totfact=1", "0.511 1", "proj=photon", "e-type=8",
                "ne=1"])["source_mode"] == "mono-gamma-511keV"
    assert src(["totfact=1", "2.0 1", "proj=photon", "e-type=8",
                "ne=1"])["source_mode"] == "mono-gamma-2MeV"
    # Historic Co-60 source syntax was allowed without proj/e-type/ne.
    assert src(["totfact=1", "1.173 1", "1.332 1"])["source_mode"] == "co60-lines"


def test_tally_histogram_bounds_are_separate_from_incident_energy(tmp_path):
    path = tmp_path / "trace.inp"
    # Preserve the previously verified 2-MeV, 1.5-MeV-histogram setup.
    path.write_text(SOURCE.replace("0.511 1", "2.0 1"), encoding="utf-8")
    build.validate_source_input(path)
    # A histogram spanning the full available photon energy is also valid.
    path.write_text(SOURCE.replace("0.511 1", "2.0 1").replace(
        "emax=1.5", "emax=2.0"), encoding="utf-8")
    build.validate_source_input(path)
    # Any histogram binning is optional for Dumpall reconstruction.
    path.write_text(SOURCE.replace("0.511 1", "2.0 1").replace(
        "emax=1.5", "emax=1.7"), encoding="utf-8")
    build.validate_source_input(path)
    path.write_text(SOURCE.replace("0.511 1", "0.662 1").replace(
        "emax=1.5", "emax=0.7"), encoding="utf-8")
    build.validate_source_input(path)


def test_mono_gamma_energy_changes_resolved_physics_hash(tmp_path):
    source, config, dumpall, table = _files(tmp_path)
    base = build.provenance_metadata(dict(SETTINGS), config, source, dumpall, table)
    source.write_text(SOURCE.replace("0.511 1", "0.662 1"), encoding="utf-8")
    changed = build.provenance_metadata(dict(SETTINGS), config, source, dumpall, table)
    assert changed["resolved_physics_sha256"] != base["resolved_physics_sha256"]
    assert changed["source_energy_MeV"] == "0.662"


def test_phits_341_esmin_default_and_explicit_score_gate_resolution(tmp_path):
    source = tmp_path / "default_esmin.inp"
    source.write_text(SOURCE, encoding="utf-8")
    assert build.resolved_esmin_MeV(source, "3.341") == 1.0e-3
    assert min(build.resolved_esmin_MeV(source, "3.341"), 1.0e-6) == 1.0e-6
    explicit = tmp_path / "explicit_esmin.inp"
    explicit.write_text(SOURCE.replace("icntl=0", "icntl=0\nesmin=5.0e-7"), encoding="utf-8")
    assert build.resolved_esmin_MeV(explicit, "3.341") == 5.0e-7
    assert min(build.resolved_esmin_MeV(explicit, "3.341"), 1.0e-6) == 5.0e-7
    with pytest.raises(ValueError, match="only for PHITS 3.341"):
        build.resolved_esmin_MeV(source, "3.342")


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


def _without_tally_or_seed(source_text: str) -> str:
    source_text = source_text.split("[t-deposit]", 1)[0]
    return "\n".join(
        line for line in source_text.splitlines()
        if not line.lstrip().lower().startswith("bitrseed")
    ) + "\n"


def test_dumpall_contract_does_not_require_t_deposit_seed_batch_count_or_icntl(tmp_path):
    source = tmp_path / "minimal.inp"
    minimal = _without_tally_or_seed(SOURCE)
    minimal = "\n".join(
        line for line in minimal.splitlines()
        if not line.lstrip().lower().startswith("icntl")
    ) + "\n"
    source.write_text(minimal, encoding="utf-8")
    build.validate_source_input(source)
    effective = build.physics_config(source)
    assert effective["source_seed_bits"] is None
    assert not build._input_section_present(source, "t-deposit")

    settings = dict(SETTINGS)
    dummy_config = tmp_path / "contract.json"
    dummy_config.write_text(json.dumps(settings), encoding="utf-8")
    dump = tmp_path / "dumpall.dat"
    dump.write_bytes(b"raw mock dump")
    table = tmp_path / "table.csv"
    table.write_bytes(b"mock stopping power bytes")
    metadata = build.provenance_metadata(settings, dummy_config, source, dump, table)
    assert metadata["t_deposit_status"] == "absent"
    assert metadata["source_seed_status"] == "not_explicit_in_input"
    assert metadata["phits_input_sha256"] == build.sha256(source)


def test_t_deposit_histogram_not_used_or_required_even_if_different(tmp_path):
    source = tmp_path / "other_tally.inp"
    source.write_text(SOURCE.replace("emax=1.5", "emax=8.0").replace(
        "unit=3", "unit=1").replace("reg=1000", "reg=9999"), encoding="utf-8")
    build.validate_source_input(source)
    source.write_text(SOURCE.split("[t-deposit]")[0], encoding="utf-8")
    build.validate_source_input(source)


@pytest.mark.parametrize("required", [
    "negs=1", "e-mode=2", "emin(12)=0.001",
    "emin(13)=0.001", "emin(14)=0.001",
])
def test_transport_mode_and_cutoff_settings_still_required(tmp_path, required):
    path = tmp_path / "wrong.inp"
    path.write_text(_without_tally_or_seed(SOURCE).replace(required, ""), encoding="utf-8")
    with pytest.raises(ValueError, match="source input must set"):
        build.validate_source_input(path)


def test_explicit_incompatible_mode_and_nonunit_source_still_fail(tmp_path):
    path = tmp_path / "wrong.inp"
    path.write_text(_without_tally_or_seed(SOURCE).replace("icntl=0", "icntl=1"),
                    encoding="utf-8")
    with pytest.raises(ValueError, match="icntl"):
        build.validate_source_input(path)
    path.write_text(_without_tally_or_seed(SOURCE).replace("wgt=1.0", "wgt=0.5"),
                    encoding="utf-8")
    with pytest.raises(ValueError, match="unit"):
        build.validate_source_input(path)


def test_source_ne_and_totfact_can_be_omitted_when_discrete_lines_are_present():
    lines = ["proj=photon", "e-type=8", "0.663 1"]
    result = build.source_definition(lines)
    assert result["source_energies_MeV"] == [0.663]
    with pytest.raises(ValueError, match="ne"):
        build.source_definition(lines + ["ne=2"])
    with pytest.raises(ValueError, match="unit"):
        build.source_definition(lines + ["totfact=2"])


def test_optional_seed_changes_byte_hash_not_semantic_physics_hash(tmp_path):
    source, config, dumpall, table = _files(tmp_path)
    old = build.provenance_metadata(dict(SETTINGS), config, source, dumpall, table)
    source.write_text(SOURCE.replace(
        "bitrseed = 0101101011001010010111000011100101101001110000111010101010101010\n",
        ""), encoding="utf-8")
    new = build.provenance_metadata(dict(SETTINGS), config, source, dumpall, table)
    assert new["source_seed_status"] == "not_explicit_in_input"
    assert new["source_input_sha256"] != old["source_input_sha256"]
    assert new["resolved_physics_sha256"] == old["resolved_physics_sha256"]

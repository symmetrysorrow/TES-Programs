"""Pb/Bi contract tests. Bi data are synthetic fixtures, NOT PHITS physics evidence."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import build_physical_hdf5_v2 as builder


SOURCE = """[Parameters]
icntl=0
maxbch=1
negs=1
e-mode=2
emin(12)=0.001
emin(13)=0.001
emin(14)=0.001
bitrseed=0101101011001010010111000011100101101001110000111011001010010111
[Source]
totfact=1.0
1.173 1
1.332 1
[Material]
mat[17] Bi 1
[Cell]
2007 17 -10.12345 -100
[t-deposit]
mesh=reg
reg=2007
part=all
deposit=0
unit=3
output=deposit
e-type=2
ne=150
emin=0.0
emax=1.5
"""


def fixture(tmp_path, *, source=SOURCE):
    path = tmp_path / "bi_synthetic.inp"
    path.write_text(source, encoding="utf-8")
    config = {
        "phits_version": "3.341",
        "phits_input": path.name,
        "region_material": "Bi",
        "region_id": 2007,
        "dumpall_material_id": 17,
        "egs5_medium": 5,
        # Intentionally synthetic test density, NOT a physical Bi value.
        "cell_density_g_cm3": 10.12345,
        "stopping_table_sha256": "a" * 64,
        "electron_cutoff_MeV": 0.001,
        "positron_cutoff_MeV": 0.001,
        "weight_mode": "unit",
        "energy_unit": "MeV",
        "deposit_mode": "t-deposit-output-deposit",
        "esmin_MeV": 0.001,
        "dead_secondary_score_threshold_MeV": 1.0e-6,
    }
    contract = tmp_path / "bi_contract.json"
    contract.write_text(json.dumps(config), encoding="utf-8")
    return contract, config, path


def test_bi_contract_parses_independent_material_medium_region(tmp_path):
    contract, data, path = fixture(tmp_path)
    parsed, resolved_path = builder.load_physics_config(contract, allow_pending_bi=True)
    assert resolved_path == path
    assert parsed == data
    settings = builder.material_settings(parsed)
    assert settings["region_material"] == "Bi"
    assert settings["region_id"] == 2007
    assert settings["dumpall_material_id"] == 17
    assert settings["egs5_medium"] == 5
    assert settings["validation_status"] == "BI_PENDING_PHITS_ORACLE"


def test_bi_normal_hdf5_output_is_not_approved(tmp_path):
    contract, _, _ = fixture(tmp_path)
    with pytest.raises(ValueError, match="not yet PHITS-oracle validated"):
        builder.load_physics_config(contract)


@pytest.mark.parametrize("change", [
    {"region_material": "Cu"},
    {"region_material": "Pb"},
    {"region_id": 0},
    {"dumpall_material_id": 0},
    {"egs5_medium": -1},
    {"cell_density_g_cm3": None},
    {"cell_density_g_cm3": -1.0},
    {"stopping_table_sha256": ""},
    {"stopping_table_sha256": "Pb HATCH"},
    {"weight_mode": "weighted"},
])
def test_bi_contract_rejects_bad_scope(tmp_path, change):
    contract, data, _ = fixture(tmp_path)
    data.update(change)
    contract.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises((ValueError, TypeError)):
        builder.load_physics_config(contract, allow_pending_bi=True)


@pytest.mark.parametrize("source", [
    SOURCE.replace("mat[17] Bi 1", "mat[17] Pb 1"),
    SOURCE.replace("2007 17 -10.12345", "2007 5 -10.12345"),
    SOURCE.replace("2007 17 -10.12345", "2007 17 -9.78"),
    SOURCE.replace("reg=2007", "reg=1000"),
    SOURCE.replace("emin(13)=0.001", "emin(13)=0.01"),
])
def test_bi_source_input_mismatch_rejected(tmp_path, source):
    contract, _, _ = fixture(tmp_path, source=source)
    with pytest.raises(ValueError):
        builder.load_physics_config(contract, allow_pending_bi=True)


def test_pb_contract_unchanged(tmp_path):
    contract, data, path = fixture(tmp_path)
    pb_path = path.with_name("pb.inp")
    pb_path.write_text(SOURCE.replace("mat[17] Bi 1", "mat[2] Pb 1")
                       .replace("2007 17 -10.12345", "1000 2 -9.78")
                       .replace("reg=2007", "reg=1000"), encoding="utf-8")
    data.update({"region_material": "Pb", "region_id": 1000,
                 "dumpall_material_id": 2, "egs5_medium": 2, "phits_input": "pb.inp"})
    data.pop("cell_density_g_cm3")
    data.pop("stopping_table_sha256")
    contract.write_text(json.dumps(data), encoding="utf-8")
    parsed, _ = builder.load_physics_config(contract)
    assert builder.material_settings(parsed)["validation_status"] == "PB_SAME_RUN_VALIDATED"


def test_bi_hatch_hash_must_match_exact_private_csv_bytes(tmp_path):
    import hashlib
    _, data, _ = fixture(tmp_path)
    table = tmp_path / "fake_bi_coefficients.csv"
    table.write_bytes(b"synthetic test-only coefficients")
    data["stopping_table_sha256"] = hashlib.sha256(table.read_bytes()).hexdigest()
    assert builder.validate_material_stopping_table(data, table) == data["stopping_table_sha256"]
    table.write_bytes(b"tampered synthetic bytes")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        builder.validate_material_stopping_table(data, table)

"""The optional in-repo Pb data MUST exactly match the previously validated bytes."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import build_physical_hdf5_v2 as pb


def test_explicit_modified_table_is_rejected(tmp_path):
    table = tmp_path / "private.csv"
    table.write_text("explicit data", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum mismatch"):
        pb.resolve_stopping_table(table)


def test_explicit_missing_table_reports_path(tmp_path):
    with pytest.raises(FileNotFoundError, match="explicit stopping table"):
        pb.resolve_stopping_table(tmp_path / "missing.csv")


def test_private_environment_table_is_checksum_pinned(monkeypatch, tmp_path):
    table = tmp_path / "private.csv"
    table.write_text("not the verified historical table", encoding="utf-8")
    monkeypatch.setenv("TES_CPP_PB_HATCH_TABLE", str(table))
    with pytest.raises(ValueError, match="checksum mismatch"):
        pb.resolve_stopping_table(None)


def test_missing_default_is_fail_closed(monkeypatch, tmp_path):
    driver = tmp_path / "validation" / "build_physical_hdf5_v2.py"
    driver.parent.mkdir(parents=True)
    driver.write_text("")
    monkeypatch.setattr(pb, "__file__", str(driver))
    monkeypatch.delenv("TES_CPP_PB_HATCH_TABLE", raising=False)
    with pytest.raises(FileNotFoundError, match="verified Pb HATCH"):
        pb.resolve_stopping_table(None)


def test_corrupt_bundled_table_is_fail_closed(monkeypatch, tmp_path):
    driver = tmp_path / "validation" / "build_physical_hdf5_v2.py"
    driver.parent.mkdir(parents=True)
    driver.write_text("")
    table_dir = tmp_path / "data"
    table_dir.mkdir()
    wrong = table_dir / "pb_egs5_phits3341.csv"
    wrong.write_text("synthetic data is NOT a HATCH table", encoding="utf-8")
    monkeypatch.setattr(pb, "__file__", str(driver))
    monkeypatch.delenv("TES_CPP_PB_HATCH_TABLE", raising=False)
    with pytest.raises(ValueError, match="checksum mismatch"):
        pb.resolve_stopping_table(None)


def test_exact_sha_pin_unchanged():
    assert pb.BUNDLED_PB_TABLE_SHA256 == (
        "20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9"
    )


def test_cmake_installed_bin_path_is_discovered(monkeypatch, tmp_path):
    prefix = tmp_path / "prefix"
    driver_dir = prefix / "share" / "tes_cpp_dump2event"
    driver_dir.mkdir(parents=True)
    driver = driver_dir / "build_physical_hdf5_v2.py"
    driver.write_text("", encoding="utf-8")
    table = prefix / "bin" / pb.BUNDLED_PB_TABLE_NAME
    table.parent.mkdir(parents=True)
    table.write_text("test-only path fixture; never physics data", encoding="utf-8")
    monkeypatch.setattr(pb, "__file__", str(driver))
    monkeypatch.delenv("TES_CPP_PB_HATCH_TABLE", raising=False)
    # Stub only the byte-digest function to test installation path discovery.
    # This synthetic fixture is never accepted or used as physical input.
    monkeypatch.setattr(pb, "sha256", lambda _: pb.BUNDLED_PB_TABLE_SHA256)
    assert pb.resolve_stopping_table(None) == table


def test_importer_refuses_fake_coefficients(tmp_path):
    import importlib.util
    importer_path = (Path(__file__).resolve().parents[1] /
                     "tools" / "import_verified_pb_table.py")
    spec = importlib.util.spec_from_file_location("_pb_coeff_importer", importer_path)
    importer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(importer)
    fake = tmp_path / "fake.csv"
    fake.write_text("medium,medname,lelke\n2,PB,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256"):
        importer.validate_original(fake)



def test_bundled_table_precedes_stale_environment_override(monkeypatch, tmp_path):
    driver = tmp_path / "validation" / "build_physical_hdf5_v2.py"
    driver.parent.mkdir(parents=True)
    driver.write_text("", encoding="utf-8")
    packaged = driver.parent / pb.BUNDLED_PB_TABLE_NAME
    packaged.write_text("synthetic package fixture", encoding="utf-8")
    stale = tmp_path / "stale.csv"
    stale.write_text("an invalid legacy private path", encoding="utf-8")
    monkeypatch.setattr(pb, "__file__", str(driver))
    monkeypatch.setenv("TES_CPP_PB_HATCH_TABLE", str(stale))
    # Stub digest only on the selected packaged fixture; never used in physics.
    original_sha = pb.sha256
    monkeypatch.setattr(
        pb, "sha256",
        lambda path: pb.BUNDLED_PB_TABLE_SHA256 if Path(path) == packaged else original_sha(path),
    )
    assert pb.resolve_stopping_table(None) == packaged


def test_private_importer_requires_confirmed_github_visibility(monkeypatch):
    import importlib.util
    from subprocess import CompletedProcess
    target = Path(__file__).resolve().parents[1] / "tools" / "import_verified_pb_table.py"
    spec = importlib.util.spec_from_file_location("pb_private_import", target)
    importer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(importer)
    monkeypatch.setattr(importer.subprocess, "run", lambda *args, **kwargs:
                        CompletedProcess(args[0], 0, stdout="false\\n"))
    with pytest.raises(RuntimeError, match="NOT private"):
        importer.require_private_github_repository()
    monkeypatch.setattr(importer.subprocess, "run", lambda *args, **kwargs:
                        CompletedProcess(args[0], 0, stdout="true\\n"))
    importer.require_private_github_repository()


def test_private_importer_rejects_missing_github_cli(monkeypatch):
    import importlib.util
    target = Path(__file__).resolve().parents[1] / "tools" / "import_verified_pb_table.py"
    spec = importlib.util.spec_from_file_location("pb_private_import_missing_gh", target)
    importer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(importer)

    def missing_gh(*args, **kwargs):
        raise FileNotFoundError("gh missing")

    monkeypatch.setattr(importer.subprocess, "run", missing_gh)
    with pytest.raises(RuntimeError, match="Cannot confirm GitHub privacy"):
        importer.require_private_github_repository()

"""The project CLI must keep new and legacy workflows separate."""

import main


def test_current_project_dry_run_shows_automatic_mesh_step(capsys):
    assert main.main(["run", "single_pixel_alpha240", "--dry-run"]) == 0
    output = capsys.readouterr().out
    assert "Mesh:" in output
    assert "projects\\models\\single_pixel.toml" in output or "projects/models/single_pixel.toml" in output
    assert "Solve: steady" in output


def test_current_project_show_uses_same_entry_point(capsys):
    assert main.main(["show", "single_pixel_alpha240"]) == 0
    output = capsys.readouterr().out
    assert "pixel.tes" in output
    assert "240.00" in output


def test_list_never_reads_legacy_json(monkeypatch, capsys):
    def fail(_):
        raise AssertionError("legacy JSON should not be read for the normal project list")

    monkeypatch.setattr(main, "read_legacy_project", fail)
    assert main.main(["list"]) == 0
    output = capsys.readouterr().out
    assert "single_pixel_alpha240" in output
    assert "Legacy JSON" not in output


def test_normal_parser_has_no_legacy_json_options():
    parser = main.build_parser()
    help_text = parser.format_help()
    assert "--project" not in help_text
    assert "--all" not in help_text
    assert "--record-only" not in help_text


def test_legacy_toml_is_not_accepted_by_normal_run(capsys):
    assert main.main(["run", "cases/single_pixel_h8.toml", "--dry-run"]) == 2
    error = capsys.readouterr().err
    assert "schema_version = 3" in error
    assert "main.py legacy" in error


def test_legacy_json_name_does_not_fallback_from_normal_show(capsys):
    assert main.main(["show", "case_tes_shunt_internal"]) == 2
    error = capsys.readouterr().err
    assert "unknown project" in error
    assert "main.py legacy" in error


def test_legacy_namespace_has_its_own_help(capsys):
    assert main.main(["legacy"]) == 0
    output = capsys.readouterr().out
    assert "python main.py legacy list" in output
    assert "python main.py legacy toml show" in output

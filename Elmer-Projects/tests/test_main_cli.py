"""The case-file CLI must not accidentally start a solver during inspection."""

import main


def test_current_toml_dry_run_shows_automatic_mesh_step(capsys):
    assert main.main(["run", "single_pixel_alpha240", "--dry-run"]) == 0
    output = capsys.readouterr().out
    assert "Mesh:" in output
    assert "models\\single_pixel.toml" in output or "models/single_pixel.toml" in output
    assert "Solve: steady" in output


def test_current_toml_show_uses_same_entry_point(capsys):
    assert main.main(["show", "single_pixel_alpha240"]) == 0
    output = capsys.readouterr().out
    assert "pixel.tes" in output
    assert "240.00" in output


def test_legacy_toml_dry_run_never_starts_solver(capsys):
    assert main.main(["run", "cases/single_pixel_h8.toml", "--dry-run"]) == 2
    assert "no calculation was started" in capsys.readouterr().err

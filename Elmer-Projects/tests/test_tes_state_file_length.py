import pytest

from scripts.support.build_cases import (
    TES_STATE_FILE_MAX_LEN,
    _side_state_file,
    validate_tes_state_file,
)


def test_tes_state_file_length_accepts_fortran_limit() -> None:
    validate_tes_state_file("x" * TES_STATE_FILE_MAX_LEN)


def test_tes_state_file_length_rejects_overflow() -> None:
    with pytest.raises(ValueError, match="TES State File exceeds 128 characters"):
        validate_tes_state_file("x" * (TES_STATE_FILE_MAX_LEN + 1))


def test_dual_state_file_is_root_relative_mesh_path() -> None:
    assert _side_state_file("mesh_dual", "steady_case", "L") == (
        "work/meshes/mesh_dual/steady_case_L.state"
    )

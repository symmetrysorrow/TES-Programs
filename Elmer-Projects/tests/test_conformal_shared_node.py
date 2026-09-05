from __future__ import annotations

import json
from pathlib import Path
import pytest

from scripts.analysis.check_conformal_interfaces import inspect
from scripts.support.build_cases import _MORTAR_PAIRS, bodies_and_bcs
from scripts.support.mesh_names import parse_mesh_names


ROOT = Path(__file__).resolve().parents[1]


def test_mortar_pair_contract_is_unchanged() -> None:
    assert _MORTAR_PAIRS == [
        ("TES", "zmin", "Membrane_SiNx", "zmax", True),
        ("Stycast", "zmin", "TES", "zmax", True),
        ("Stycast", "zmax", "abs", "zmin", False),
    ]


def test_mortar_bcs_are_removed_only_when_requested() -> None:
    mesh_names = ROOT / "work/meshes/mesh_refined_3x/mesh.names"
    if not mesh_names.exists():
        mesh_names = ROOT / "work/meshes/mesh_physical_parity_conformal/mesh.names"
    if not mesh_names.exists():
        pytest.skip("runtime mesh is not present in this clean source-only worktree")
    names = parse_mesh_names(mesh_names)
    with_mortar = "\n".join(bodies_and_bcs(names, False, apply_mortar_bcs=True))
    without_mortar = "\n".join(bodies_and_bcs(names, False, apply_mortar_bcs=False))
    assert "Mortar BC =" in with_mortar
    assert "Galerkin Projector" in with_mortar
    assert "Mortar BC =" not in without_mortar
    assert "Galerkin Projector" not in without_mortar


def test_conformal_project_selects_independent_shared_node_route() -> None:
    project = json.loads(
        (ROOT / "elmer_project_singlepixel_conformal_gpu.json").read_text(
            encoding="utf-8"
        )
    )
    overrides = project["elmer_overrides"]
    assert overrides["conformal_shared_node_interfaces"] is True
    assert overrides["conformal_mortar_interfaces"] is False
    assert overrides["fragment_mortar_interfaces"] is True
    assert project["cases"]["case_tes_steady_singlepixel_conformal_gpu"][
        "apply_mortar_bcs"
    ] is False


def test_connectivity_checker_requires_same_surface_partition(tmp_path: Path) -> None:
    mesh = tmp_path / "mesh"
    mesh.mkdir()
    (mesh / "mesh.names").write_text(
        "! ----- names for bodies -----\n"
        "$ TES = 101\n$ Stycast = 102\n$ abs = 100\n$ Membrane_SiNx = 103\n"
        "! ----- names for boundaries -----\n"
        "$ TES__zmax = 1105\n$ Stycast__zmin = 1204\n"
        "$ Stycast__zmax = 1205\n$ abs__zmin = 1004\n"
        "$ Membrane_SiNx__zmax = 1305\n$ TES__zmin = 1104\n",
        encoding="utf-8",
    )
    nodes = [
        "1 101 0 0 0", "2 101 1 0 0", "3 101 0 1 0",
        "4 101 0 0 1", "5 102 0 0 0", "6 102 1 0 0",
        "7 102 0 1 0", "8 102 0 0 -1", "9 100 0 0 -1",
        "10 100 1 0 -1", "11 100 0 1 -1", "12 100 0 0 -2",
        "13 103 0 0 1", "14 103 1 0 1", "15 103 0 1 1",
        "16 103 0 0 2",
    ]
    (mesh / "mesh.nodes").write_text("\n".join(nodes) + "\n", encoding="utf-8")
    (mesh / "mesh.boundary").write_text(
        "1 1104 101 504 303 1 2 3\n"
        "2 1305 103 504 303 1 2 3\n"
        "3 1105 101 504 303 5 6 7\n"
        "4 1204 102 504 303 5 6 7\n"
        "5 1205 102 504 303 9 10 11\n"
        "6 1004 100 504 303 9 10 11\n",
        encoding="utf-8",
    )
    (mesh / "mesh.elements").write_text(
        "1 101 504 1 2 3 4\n"
        "2 102 504 5 6 7 8\n"
        "3 100 504 9 10 11 12\n"
        "4 103 504 13 14 15 16\n",
        encoding="utf-8",
    )
    report = inspect(mesh)
    assert report["status"] == "PASS"
    assert report["interfaces"][0]["status"] == "PASS"
    assert report["interfaces"][1]["status"] == "PASS"
    assert report["interfaces"][2]["status"] == "PASS"
    (mesh / "mesh.boundary").write_text(
        (mesh / "mesh.boundary").read_text(encoding="utf-8").replace(
            "6 1004 100 504 303 9 10 11",
            "6 1004 100 504 303 9 10 12",
        ),
        encoding="utf-8",
    )
    assert inspect(mesh)["status"] == "FAIL"

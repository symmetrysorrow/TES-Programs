from __future__ import annotations

import json
import math
from pathlib import Path

from scripts.analysis.derive_mortar_control_mesh import derive
from scripts.analysis.evaluate_physical_parity import expected_geometry, mesh_data
from scripts.support.reconcile_project import reconcile_project


ROOT = Path(__file__).resolve().parents[1]


def _write_minimal_mesh(mesh: Path) -> None:
    mesh.mkdir()
    (mesh / "mesh.names").write_text(
        "! ----- names for bodies -----\n"
        "$ abs = 100\n$ TES = 101\n$ Stycast = 102\n$ Membrane_SiNx = 103\n"
        "! ----- names for boundaries -----\n"
        "$ abs__zmin = 1004\n$ TES__zmin = 1104\n$ TES__zmax = 1105\n"
        "$ Stycast__zmin = 1204\n$ Stycast__zmax = 1205\n"
        "$ Membrane_SiNx__zmax = 1305\n",
        encoding="utf-8",
    )
    (mesh / "entities.sif").write_text("! test\n", encoding="utf-8")
    (mesh / "mesh.header").write_text("16 4 6\n2\n303 6\n504 4\n", encoding="utf-8")
    (mesh / "mesh.nodes").write_text(
        "\n".join(
            [
                "1 -1 0 0 0", "2 -1 1 0 0", "3 -1 0 1 0", "4 -1 0 0 1",
                "5 -1 0 0 1", "6 -1 1 0 1", "7 -1 0 1 1", "8 -1 0 0 2",
                "9 -1 0 0 -1", "10 -1 1 0 -1", "11 -1 0 1 -1", "12 -1 0 0 -2",
                "13 -1 0 0 2", "14 -1 1 0 2", "15 -1 0 1 2", "16 -1 0 0 3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (mesh / "mesh.elements").write_text(
        "1 101 504 1 2 3 4\n"
        "2 102 504 5 6 7 8\n"
        "3 100 504 9 10 11 12\n"
        "4 103 504 13 14 15 16\n",
        encoding="utf-8",
    )
    (mesh / "mesh.boundary").write_text(
        "1 1104 1 0 303 1 2 3\n"
        "2 1305 4 0 303 1 2 3\n"
        "3 1105 1 0 303 5 6 7\n"
        "4 1204 2 0 303 5 6 7\n"
        "5 1205 2 0 303 9 10 11\n"
        "6 1004 3 0 303 9 10 11\n",
        encoding="utf-8",
    )


def test_control_derivation_preserves_volume_connectivity_and_duplicates_nodes(tmp_path: Path) -> None:
    source = tmp_path / "source"
    output = tmp_path / "control"
    _write_minimal_mesh(source)

    provenance = derive(source, output)

    source_data = mesh_data(source)
    output_data = mesh_data(output)
    assert provenance["node_count"] == 25
    assert len(output_data[2]) > len(source_data[2])
    assert len(output_data[3]) == len(source_data[3])
    assert len(output_data[4][101]) == len(source_data[4][101])
    assert len(output_data[4][102]) == len(source_data[4][102])
    assert len(output_data[4][100]) == len(source_data[4][100])
    assert len(output_data[4][103]) == len(source_data[4][103])
    assert output_data[6]["TES__zmin"]
    assert output_data[6]["Membrane_SiNx__zmax"]
    assert set(output_data[6]["TES__zmin"]) != set(output_data[6]["Membrane_SiNx__zmax"])


def test_stycast_analytic_geometry_uses_498_um_literal() -> None:
    project = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
    expected, contact = expected_geometry(project)
    stycast_dz = reconcile_project(json.loads(project.read_text(encoding="utf-8")))["parameters"]["Stycast_dz"]
    expected_volume = math.pi * (498.0e-6 / 2.0) ** 2 * stycast_dz
    old_400_um_volume = math.pi * (400.0e-6 / 2.0) ** 2 * stycast_dz
    assert math.isclose(expected["Stycast"], expected_volume, rel_tol=1.0e-15)
    assert math.isclose(contact["TES_Stycast"], expected_volume / stycast_dz, rel_tol=1.0e-15)
    assert not math.isclose(expected["Stycast"], old_400_um_volume, rel_tol=1.0e-3)

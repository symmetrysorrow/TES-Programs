from pathlib import Path

from scripts.support.build_cases import resolve_body_sif_ordinal
from scripts.support.mesh_names import MeshNames, parse_mesh_names


def test_element_body_id_uses_sif_ordinal_not_physical_target() -> None:
    # HeatSolve compares Element%BodyId, which is the one-based SIF Body
    # ordinal, rather than the non-contiguous physical target in mesh.names.
    mesh_names = MeshNames(bodies={"abs": 100, "TES": 101, "Si_1": 102}, boundaries={})
    assert resolve_body_sif_ordinal(mesh_names, "TES") == 2


def test_hybrid_tes_sif_ordinal_is_eight() -> None:
    mesh_names = parse_mesh_names(Path("mesh_hybrid_abs_tet_layers_prism_conformal/mesh.names"))
    assert resolve_body_sif_ordinal(mesh_names, "TES") == 8

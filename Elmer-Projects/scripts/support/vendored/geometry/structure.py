"""Freeform-box StructureSpec assembly.

Vendored 2026-07-14 from Thermal-and-Electoric-Sim: the freeform_boxes subset
of core/geometry/templates.py plus _auto_mesh_sizes from core/project/
builders.py, verbatim except for import paths. The layered_box /
box_with_leads / freeform_paths templates are not used by this project and
were not vendored.
"""
from __future__ import annotations

from scripts.support.vendored.geometry.body_semantics import body_mode_of, body_name_of
from scripts.support.vendored.geometry.primitives import Box
from scripts.support.vendored.geometry.spec import StructureSpec


def _add_all_box_surfaces(
    spec: StructureSpec,
    *,
    box_name: str,
    box_uid: str,
    base_tag: int,
) -> None:
    surfaces = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
    for offset, surface in enumerate(surfaces):
        spec.add_physical_surface(
            name=f"{box_name}__{surface}",
            selector="box_surface",
            tag=base_tag + offset,
            uid=f"{box_uid}__{surface}",
            box_name=box_name,
            surface=surface,
        )


def _auto_reference_box_name(boxes: list[Box]) -> str:
    if not boxes:
        return ""
    additive = [box for box in boxes if body_mode_of(box) == "add"]
    candidates = additive or boxes
    return max(candidates, key=lambda box: float(box.dx) * float(box.dy) * float(box.dz)).name


def build_freeform_boxes_structure(
    *,
    boxes: list[Box],
    bath_surface_name: str = "bath",
    bath_surface_uid: str = "",
    reference_box_name: str = "",
    mesh_min: float = 0.05,
    mesh_max: float = 0.05,
) -> StructureSpec:
    spec = StructureSpec(
        mesh_min=mesh_min,
        mesh_max=mesh_max,
        use_boolean_fragments=True,
    )

    for box in boxes:
        spec.add_box(box)

    grouped: dict[str, list[Box]] = {}
    for box in boxes:
        grouped.setdefault(body_name_of(box), []).append(box)

    volume_tag_base = 100
    face_tag_base = 1000
    for idx, (group_name, members) in enumerate(grouped.items()):
        volume_uid = next((str(box.uid).strip() for box in members if str(box.uid).strip()), "")
        spec.add_physical_volume(
            name=group_name,
            primitive_names=[group_name],
            tag=volume_tag_base + idx,
            uid=volume_uid,
        )
        for box in members:
            _add_all_box_surfaces(
                spec,
                box_name=box.name,
                box_uid=box.uid,
                base_tag=face_tag_base + idx * 100 + members.index(box) * 10,
            )

    reference_box_name = reference_box_name or _auto_reference_box_name(boxes)
    if reference_box_name:
        spec.add_physical_surface(
            name=bath_surface_name,
            selector="all_boundary_surfaces",
            tag=21,
            box_name=reference_box_name,
            uid=bath_surface_uid,
        )

    return spec


def _auto_mesh_sizes(structure: StructureSpec) -> tuple[float, float]:
    if not structure.boxes:
        return (0.05, 0.05)

    spans = []
    feature_sizes = []
    min_x = min_y = min_z = float("inf")
    max_x = max_y = max_z = float("-inf")
    for box in structure.boxes:
        spans.extend([float(box.dx), float(box.dy), float(box.dz)])
        feature_sizes.append(max(float(box.dx), float(box.dy), float(box.dz)))
        min_x = min(min_x, float(box.xmin))
        min_y = min(min_y, float(box.ymin))
        min_z = min(min_z, float(box.zmin))
        max_x = max(max_x, float(box.xmax))
        max_y = max(max_y, float(box.ymax))
        max_z = max(max_z, float(box.zmax))

    if not spans:
        return (0.05, 0.05)

    overall_span = max(max_x - min_x, max_y - min_y, max_z - min_z)
    smallest_feature = min(feature_sizes) if feature_sizes else overall_span
    if overall_span <= 0.0:
        overall_span = max(spans)
    if overall_span <= 0.0:
        overall_span = 0.05
    if smallest_feature <= 0.0:
        smallest_feature = overall_span

    target = min(overall_span / 25.0, smallest_feature / 3.0)
    floor = overall_span / 1000.0
    mesh_max = max(target, floor)
    mesh_min = mesh_max / 2.0
    return (mesh_min, mesh_max)

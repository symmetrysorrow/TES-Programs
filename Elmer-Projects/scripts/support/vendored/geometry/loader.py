"""Project-JSON geometry -> StructureSpec loader.

Replaces the external core/project/io.py + builders.py chain with the subset
this project uses: a `freeform_boxes` geometry given as nested `children`
nodes (groups inherit `group_mode`, e.g. subtraction boxes inside the
substrate mothers). The traversal and box-field mapping reproduce the
upstream behaviour verbatim so that the generated gmsh model is identical.
"""
from __future__ import annotations

from typing import Any

from scripts.support.vendored.geometry.primitives import Box
from scripts.support.vendored.geometry.spec import StructureSpec
from scripts.support.vendored.geometry.structure import (
    _auto_mesh_sizes,
    build_freeform_boxes_structure,
)


def _collect_geometry_leaves(
    nodes: list[dict[str, Any]],
    *,
    parent_group_path: str = "",
    parent_group_mode: str = "add",
) -> list[tuple[dict[str, Any], str, str]]:
    leaves: list[tuple[dict[str, Any], str, str]] = []
    for node in nodes:
        if not isinstance(node, dict):
            continue
        children = node.get("children")
        if isinstance(children, list) and children:
            name = str(node.get("name") or "").strip()
            if name:
                next_group_path = f"{parent_group_path}/{name}" if parent_group_path else name
            else:
                next_group_path = parent_group_path
            next_group_mode = str(node.get("group_mode") or parent_group_mode or "add").strip()
            leaves.extend(
                _collect_geometry_leaves(
                    children,
                    parent_group_path=next_group_path,
                    parent_group_mode=next_group_mode,
                )
            )
        else:
            leaf_mode = str(node.get("group_mode") or parent_group_mode or "add").strip()
            leaves.append((node, parent_group_path, leaf_mode))
    return leaves


def _boxes_from_children(children: list[dict[str, Any]]) -> list[Box]:
    boxes: list[Box] = []
    for node, group_path, group_mode in _collect_geometry_leaves(children):
        shape = str(node.get("shape") or node.get("kind") or "box").strip().lower()
        if shape not in {"box", "wedge"}:
            continue
        boxes.append(
            Box(
                name=str(node.get("name", "")),
                uid=str(node.get("uid", "")),
                x=float(node.get("x", 0.0)),
                y=float(node.get("y", 0.0)),
                z=float(node.get("z", 0.0)),
                dx=float(node.get("dx", 0.0)),
                dy=float(node.get("dy", 0.0)),
                dz=float(node.get("dz", 0.0)),
                kind=shape,
                body_name=str(node.get("body_name") or group_path or node.get("name") or ""),
                body_mode=str(node.get("body_mode") or node.get("group_mode") or group_mode or "add").strip().lower(),
                group_mode=group_mode,
                group_name=group_path,
                ltx=float(node.get("ltx", 0.0)),
                priority=int(float(node.get("priority", 0) or 0)),
            )
        )
    return boxes


def load_structure_spec(model: dict[str, Any]) -> StructureSpec:
    """Build the StructureSpec from a (reconciled) project model dict."""
    geometry = model.get("geometry")
    if not isinstance(geometry, dict):
        raise ValueError("elmer_project.json has no geometry section")
    kind = str(geometry.get("template_kind") or geometry.get("geometry_kind") or "").strip()
    if kind != "freeform_boxes":
        raise ValueError(f"Unsupported geometry template_kind: {kind!r} (only freeform_boxes)")

    children = geometry.get("children")
    if not isinstance(children, list) or not children:
        raise ValueError("geometry.children is empty")

    mesh = model.get("mesh", {})
    structure = build_freeform_boxes_structure(
        boxes=_boxes_from_children(children),
        bath_surface_name=str(geometry.get("bath_surface_name", "bath")),
        bath_surface_uid=str(geometry.get("bath_surface_uid", "")),
        reference_box_name=str(geometry.get("reference_box_name", "")),
        mesh_min=float(mesh.get("mesh_min", 0.05)),
        mesh_max=float(mesh.get("mesh_max", 0.05)),
    )
    if str(mesh.get("mesh_min_mode", "")).strip().lower() == "auto" or str(
        mesh.get("mesh_max_mode", "")
    ).strip().lower() == "auto":
        auto_min, auto_max = _auto_mesh_sizes(structure)
        if str(mesh.get("mesh_min_mode", "")).strip().lower() == "auto":
            structure.mesh_min = auto_min
        if str(mesh.get("mesh_max_mode", "")).strip().lower() == "auto":
            structure.mesh_max = auto_max
    return structure

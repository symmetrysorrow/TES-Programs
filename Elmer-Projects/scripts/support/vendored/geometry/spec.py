"""Vendored 2026-07-14 from Thermal-and-Electoric-Sim core/geometry/spec.py
(import paths adjusted only). See vendored/__init__.py for policy."""
from __future__ import annotations

from dataclasses import dataclass, field

from scripts.support.vendored.geometry.primitives import Box


@dataclass
class PhysicalVolumeSpec:
    name: str
    primitive_names: list[str]
    tag: int
    uid: str = ""


@dataclass
class PhysicalSurfaceSpec:
    name: str
    selector: str
    tag: int
    uid: str = ""
    box_name: str | None = None
    surface: str | None = None


@dataclass
class StructureSpec:
    boxes: list[Box] = field(default_factory=list)
    physical_volumes: list[PhysicalVolumeSpec] = field(default_factory=list)
    physical_surfaces: list[PhysicalSurfaceSpec] = field(default_factory=list)

    mesh_min: float = 0.2
    mesh_max: float = 0.2

    use_boolean_fragments: bool = True

    def add_box(self, box: Box) -> None:
        self.boxes.append(box)

    def add_physical_volume(
        self,
        *,
        name: str,
        primitive_names: list[str],
        tag: int,
        uid: str = "",
    ) -> None:
        self.physical_volumes.append(
            PhysicalVolumeSpec(
                name=name,
                primitive_names=primitive_names,
                tag=tag,
                uid=uid,
            )
        )

    def add_physical_surface(
        self,
        *,
        name: str,
        selector: str,
        tag: int,
        uid: str = "",
        box_name: str | None = None,
        surface: str | None = None,
    ) -> None:
        self.physical_surfaces.append(
            PhysicalSurfaceSpec(
                name=name,
                selector=selector,
                tag=tag,
                uid=uid,
                box_name=box_name,
                surface=surface,
            )
        )

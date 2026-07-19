"""Vendored 2026-07-14 from Thermal-and-Electoric-Sim core/geometry/primitives.py
(import paths adjusted only). See vendored/__init__.py for policy."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, init=False)
class Box:
    name: str
    x: float
    y: float
    z: float
    dx: float
    dy: float
    dz: float
    kind: str = "box"
    uid: str = ""
    body_name: str = ""
    group_name: str = ""
    ltx: float = 0.0
    x_expr: str | None = None
    y_expr: str | None = None
    z_expr: str | None = None
    dx_expr: str | None = None
    dy_expr: str | None = None
    dz_expr: str | None = None
    ltx_expr: str | None = None
    # Priority controls overlap resolution. Higher values win.
    # Semantic names are handled in project configuration, not here.
    priority: int = 0

    def __init__(
        self,
        *,
        name: str,
        x: float,
        y: float,
        z: float,
        dx: float,
        dy: float,
        dz: float,
        kind: str = "box",
        uid: str = "",
        body_mode: str | None = None,
        body_name: str = "",
        group_mode: str = "add",
        group_name: str = "",
        ltx: float = 0.0,
        x_expr: str | None = None,
        y_expr: str | None = None,
        z_expr: str | None = None,
        dx_expr: str | None = None,
        dy_expr: str | None = None,
        dz_expr: str | None = None,
        ltx_expr: str | None = None,
        priority: int = 0,
    ) -> None:
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "y", y)
        object.__setattr__(self, "z", z)
        object.__setattr__(self, "dx", dx)
        object.__setattr__(self, "dy", dy)
        object.__setattr__(self, "dz", dz)
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "uid", uid)
        object.__setattr__(self, "body_name", body_name)
        object.__setattr__(self, "group_name", group_name)
        object.__setattr__(self, "ltx", ltx)
        object.__setattr__(self, "x_expr", x_expr)
        object.__setattr__(self, "y_expr", y_expr)
        object.__setattr__(self, "z_expr", z_expr)
        object.__setattr__(self, "dx_expr", dx_expr)
        object.__setattr__(self, "dy_expr", dy_expr)
        object.__setattr__(self, "dz_expr", dz_expr)
        object.__setattr__(self, "ltx_expr", ltx_expr)
        object.__setattr__(self, "priority", priority)
        resolved_body_mode = str(body_mode).strip() if body_mode is not None else ""
        object.__setattr__(self, "_body_mode", resolved_body_mode or str(group_mode or "add"))
        object.__setattr__(self, "_group_mode", str(group_mode or "add"))

    @property
    def body_mode(self) -> str:
        return str(getattr(self, "_body_mode", self.group_mode or "add"))

    @property
    def group_mode(self) -> str:
        return str(getattr(self, "_group_mode", "add"))

    @property
    def effective_body_name(self) -> str:
        return str(getattr(self, "body_name", "") or "").strip() or str(getattr(self, "group_name", "") or "").strip() or self.name

    @property
    def effective_body_mode(self) -> str:
        return str(getattr(self, "body_mode", "") or "").strip().lower() or str(getattr(self, "group_mode", "") or "").strip().lower() or "add"

    @property
    def xmin(self) -> float:
        return self.x - self.dx / 2.0

    @property
    def x2(self) -> float:
        return self.x + self.dx / 2.0

    @property
    def xmax(self) -> float:
        return self.x + self.dx / 2.0

    @property
    def ymin(self) -> float:
        return self.y - self.dy / 2.0

    @property
    def y2(self) -> float:
        return self.y + self.dy / 2.0

    @property
    def ymax(self) -> float:
        return self.y + self.dy / 2.0

    @property
    def zmin(self) -> float:
        return self.z - self.dz / 2.0

    @property
    def z2(self) -> float:
        return self.z + self.dz / 2.0

    @property
    def zmax(self) -> float:
        return self.z + self.dz / 2.0

    @property
    def volume(self) -> float:
        return self.dx * self.dy * self.dz

    def contains_point(
        self,
        x: float,
        y: float,
        z: float,
        *,
        eps: float = 1.0e-9,
    ) -> bool:
        return (
            self.xmin - eps <= x <= self.xmax + eps
            and self.ymin - eps <= y <= self.ymax + eps
            and self.zmin - eps <= z <= self.zmax + eps
        )

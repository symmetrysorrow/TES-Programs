"""Vendored 2026-07-14 from Thermal-and-Electoric-Sim core/geometry/contact_detection.py
(import paths adjusted only). See vendored/__init__.py for policy."""
from __future__ import annotations

from dataclasses import dataclass

from scripts.support.vendored.geometry.body_semantics import body_mode_of, body_name_of
from scripts.support.vendored.geometry.spec import StructureSpec


@dataclass(frozen=True)
class ContactCandidate:
    name: str
    surface_a: str
    surface_b: str
    box_a: str
    box_b: str
    face_a: str
    face_b: str
    area: float


def _overlap_1d(a0: float, a1: float, b0: float, b1: float) -> float:
    return max(0.0, min(a1, b1) - max(a0, b0))


def _global_bounds(spec: StructureSpec) -> tuple[float, float, float, float, float, float]:
    xmin = min(box.xmin for box in spec.boxes)
    xmax = max(box.xmax for box in spec.boxes)
    ymin = min(box.ymin for box in spec.boxes)
    ymax = max(box.ymax for box in spec.boxes)
    zmin = min(box.zmin for box in spec.boxes)
    zmax = max(box.zmax for box in spec.boxes)
    return xmin, xmax, ymin, ymax, zmin, zmax


def _subtract_rectangle(
    base: tuple[float, float, float, float],
    cut: tuple[float, float, float, float],
    *,
    tol: float,
) -> list[tuple[float, float, float, float]]:
    ax0, ax1, ay0, ay1 = base
    bx0, bx1, by0, by1 = cut
    ix0 = max(ax0, bx0)
    ix1 = min(ax1, bx1)
    iy0 = max(ay0, by0)
    iy1 = min(ay1, by1)
    if ix1 - ix0 <= tol or iy1 - iy0 <= tol:
        return [base]

    pieces: list[tuple[float, float, float, float]] = []
    if ax0 < ix0 - tol:
        pieces.append((ax0, ix0, ay0, ay1))
    if ix1 < ax1 - tol:
        pieces.append((ix1, ax1, ay0, ay1))
    if ay0 < iy0 - tol:
        pieces.append((ix0, ix1, ay0, iy0))
    if iy1 < ay1 - tol:
        pieces.append((ix0, ix1, iy1, ay1))
    return [piece for piece in pieces if (piece[1] - piece[0] > tol and piece[3] - piece[2] > tol)]


def _subtract_cuts_from_rectangles(
    rectangles: list[tuple[float, float, float, float]],
    cuts: list[tuple[float, float, float, float]],
    *,
    tol: float,
) -> list[tuple[float, float, float, float]]:
    remaining = list(rectangles)
    for cut in cuts:
        updated: list[tuple[float, float, float, float]] = []
        for rect in remaining:
            updated.extend(_subtract_rectangle(rect, cut, tol=tol))
        remaining = updated
        if not remaining:
            break
    return remaining


def _rectangles_area(rectangles: list[tuple[float, float, float, float]]) -> float:
    return sum(max(0.0, x1 - x0) * max(0.0, y1 - y0) for x0, x1, y0, y1 in rectangles)


def _contact_area_with_subtracts(
    box_a,
    box_b,
    *,
    face_a: str,
    face_b: str,
    subtract_boxes: list,
    tol: float,
) -> float:
    body_names = {body_name_of(box_a), body_name_of(box_b)}
    relevant_subtract_boxes = [
        sub for sub in subtract_boxes
        if body_name_of(sub) in body_names
    ]
    if face_a in {"zmax", "zmin"} and face_b in {"zmax", "zmin"}:
        plane = float(box_a.zmax if face_a == "zmax" else box_a.zmin)
        base = [(
            max(float(box_a.xmin), float(box_b.xmin)),
            min(float(box_a.xmax), float(box_b.xmax)),
            max(float(box_a.ymin), float(box_b.ymin)),
            min(float(box_a.ymax), float(box_b.ymax)),
        )]
        cuts: list[tuple[float, float, float, float]] = []
        for sub in relevant_subtract_boxes:
            if float(sub.zmin) < plane + tol and float(sub.zmax) > plane - tol:
                cuts.append((
                    float(sub.xmin),
                    float(sub.xmax),
                    float(sub.ymin),
                    float(sub.ymax),
                ))
        return _rectangles_area(_subtract_cuts_from_rectangles(base, cuts, tol=tol))

    if face_a in {"xmax", "xmin"} and face_b in {"xmax", "xmin"}:
        plane = float(box_a.xmax if face_a == "xmax" else box_a.xmin)
        base = [(
            max(float(box_a.ymin), float(box_b.ymin)),
            min(float(box_a.ymax), float(box_b.ymax)),
            max(float(box_a.zmin), float(box_b.zmin)),
            min(float(box_a.zmax), float(box_b.zmax)),
        )]
        cuts = []
        for sub in relevant_subtract_boxes:
            if float(sub.xmin) < plane + tol and float(sub.xmax) > plane - tol:
                cuts.append((
                    float(sub.ymin),
                    float(sub.ymax),
                    float(sub.zmin),
                    float(sub.zmax),
                ))
        return _rectangles_area(_subtract_cuts_from_rectangles(base, cuts, tol=tol))

    plane = float(box_a.ymax if face_a == "ymax" else box_a.ymin)
    base = [(
        max(float(box_a.xmin), float(box_b.xmin)),
        min(float(box_a.xmax), float(box_b.xmax)),
        max(float(box_a.zmin), float(box_b.zmin)),
        min(float(box_a.zmax), float(box_b.zmax)),
    )]
    cuts = []
    for sub in relevant_subtract_boxes:
        if float(sub.ymin) < plane + tol and float(sub.ymax) > plane - tol:
            cuts.append((
                float(sub.xmin),
                float(sub.xmax),
                float(sub.zmin),
                float(sub.zmax),
            ))
    return _rectangles_area(_subtract_cuts_from_rectangles(base, cuts, tol=tol))


def detect_contact_candidates(
    spec: StructureSpec,
    *,
    tol: float = 1.0e-9,
) -> list[ContactCandidate]:
    add_boxes = [box for box in spec.boxes if body_mode_of(box) != "subtract"]
    subtract_boxes = [box for box in spec.boxes if body_mode_of(box) == "subtract"]
    if not add_boxes:
        return []

    xmin, xmax, ymin, ymax, zmin, zmax = _global_bounds(StructureSpec(boxes=add_boxes))
    area_tol = max(tol * tol, 1.0e-30)
    candidates: list[ContactCandidate] = []

    def add_candidate(
        box_a,
        box_b,
        face_a: str,
        face_b: str,
        area: float,
    ) -> None:
        surface_a = f"{box_a.name}__{face_a}"
        surface_b = f"{box_b.name}__{face_b}"
        candidates.append(
            ContactCandidate(
                name=f"contact_{box_a.name}_{face_a}__{box_b.name}_{face_b}",
                surface_a=surface_a,
                surface_b=surface_b,
                box_a=box_a.name,
                box_b=box_b.name,
                face_a=face_a,
                face_b=face_b,
                area=area,
                )
        )

    def maybe_add_contact(
        a,
        b,
        *,
        face_a: str,
        face_b: str,
        area: float,
        plane_value: float,
        lower_bound: float,
        upper_bound: float,
    ) -> None:
        if plane_value <= lower_bound + tol or plane_value >= upper_bound - tol:
            return
        if area <= area_tol:
            return
        add_candidate(a, b, face_a, face_b, area)

    def add_box_pair_contact(a, b) -> None:
        y_overlap = _overlap_1d(a.ymin, a.ymax, b.ymin, b.ymax)
        z_overlap = _overlap_1d(a.zmin, a.zmax, b.zmin, b.zmax)
        x_overlap = _overlap_1d(a.xmin, a.xmax, b.xmin, b.xmax)

        if abs(a.xmax - b.xmin) <= tol and y_overlap > tol and z_overlap > tol:
            area = _contact_area_with_subtracts(a, b, face_a="xmax", face_b="xmin", subtract_boxes=subtract_boxes, tol=tol)
            maybe_add_contact(
                a,
                b,
                face_a="xmax",
                face_b="xmin",
                area=area,
                plane_value=float(a.xmax),
                lower_bound=xmin,
                upper_bound=xmax,
            )
            return
        if abs(b.xmax - a.xmin) <= tol and y_overlap > tol and z_overlap > tol:
            area = _contact_area_with_subtracts(a, b, face_a="xmin", face_b="xmax", subtract_boxes=subtract_boxes, tol=tol)
            maybe_add_contact(
                a,
                b,
                face_a="xmin",
                face_b="xmax",
                area=area,
                plane_value=float(a.xmin),
                lower_bound=xmin,
                upper_bound=xmax,
            )
            return

        if abs(a.ymax - b.ymin) <= tol and x_overlap > tol and z_overlap > tol:
            area = _contact_area_with_subtracts(a, b, face_a="ymax", face_b="ymin", subtract_boxes=subtract_boxes, tol=tol)
            maybe_add_contact(
                a,
                b,
                face_a="ymax",
                face_b="ymin",
                area=area,
                plane_value=float(a.ymax),
                lower_bound=ymin,
                upper_bound=ymax,
            )
            return
        if abs(b.ymax - a.ymin) <= tol and x_overlap > tol and z_overlap > tol:
            area = _contact_area_with_subtracts(a, b, face_a="ymin", face_b="ymax", subtract_boxes=subtract_boxes, tol=tol)
            maybe_add_contact(
                a,
                b,
                face_a="ymin",
                face_b="ymax",
                area=area,
                plane_value=float(a.ymin),
                lower_bound=ymin,
                upper_bound=ymax,
            )
            return

        if abs(a.zmax - b.zmin) <= tol and x_overlap > tol and y_overlap > tol:
            area = _contact_area_with_subtracts(a, b, face_a="zmax", face_b="zmin", subtract_boxes=subtract_boxes, tol=tol)
            maybe_add_contact(
                a,
                b,
                face_a="zmax",
                face_b="zmin",
                area=area,
                plane_value=float(a.zmax),
                lower_bound=zmin,
                upper_bound=zmax,
            )
            return
        if abs(b.zmax - a.zmin) <= tol and x_overlap > tol and y_overlap > tol:
            area = _contact_area_with_subtracts(a, b, face_a="zmin", face_b="zmax", subtract_boxes=subtract_boxes, tol=tol)
            maybe_add_contact(
                a,
                b,
                face_a="zmin",
                face_b="zmax",
                area=area,
                plane_value=float(a.zmin),
                lower_bound=zmin,
                upper_bound=zmax,
            )
            return

    def add_add_to_subtract_contact(a, sub) -> None:
        if body_name_of(a) == body_name_of(sub):
            return
        y_overlap = _overlap_1d(a.ymin, a.ymax, sub.ymin, sub.ymax)
        z_overlap = _overlap_1d(a.zmin, a.zmax, sub.zmin, sub.zmax)
        x_overlap = _overlap_1d(a.xmin, a.xmax, sub.xmin, sub.xmax)

        # A subtract box represents a cavity carved into another body, so its
        # side walls live on the same x/y planes as the inserted add box.
        if abs(a.xmax - sub.xmax) <= tol and y_overlap > tol and z_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="xmax",
                face_b="xmax",
                area=y_overlap * z_overlap,
                plane_value=float(a.xmax),
                lower_bound=xmin,
                upper_bound=xmax,
            )
        if abs(a.xmin - sub.xmin) <= tol and y_overlap > tol and z_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="xmin",
                face_b="xmin",
                area=y_overlap * z_overlap,
                plane_value=float(a.xmin),
                lower_bound=xmin,
                upper_bound=xmax,
            )
        if abs(a.ymax - sub.ymax) <= tol and x_overlap > tol and z_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="ymax",
                face_b="ymax",
                area=x_overlap * z_overlap,
                plane_value=float(a.ymax),
                lower_bound=ymin,
                upper_bound=ymax,
            )
        if abs(a.ymin - sub.ymin) <= tol and x_overlap > tol and z_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="ymin",
                face_b="ymin",
                area=x_overlap * z_overlap,
                plane_value=float(a.ymin),
                lower_bound=ymin,
                upper_bound=ymax,
            )

        if abs(a.xmax - sub.xmin) <= tol and y_overlap > tol and z_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="xmax",
                face_b="xmin",
                area=y_overlap * z_overlap,
                plane_value=float(a.xmax),
                lower_bound=xmin,
                upper_bound=xmax,
            )
            return
        if abs(sub.xmax - a.xmin) <= tol and y_overlap > tol and z_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="xmin",
                face_b="xmax",
                area=y_overlap * z_overlap,
                plane_value=float(a.xmin),
                lower_bound=xmin,
                upper_bound=xmax,
            )
            return

        if abs(a.ymax - sub.ymin) <= tol and x_overlap > tol and z_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="ymax",
                face_b="ymin",
                area=x_overlap * z_overlap,
                plane_value=float(a.ymax),
                lower_bound=ymin,
                upper_bound=ymax,
            )
            return
        if abs(sub.ymax - a.ymin) <= tol and x_overlap > tol and z_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="ymin",
                face_b="ymax",
                area=x_overlap * z_overlap,
                plane_value=float(a.ymin),
                lower_bound=ymin,
                upper_bound=ymax,
            )
            return

        if abs(a.zmax - sub.zmin) <= tol and x_overlap > tol and y_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="zmax",
                face_b="zmin",
                area=x_overlap * y_overlap,
                plane_value=float(a.zmax),
                lower_bound=zmin,
                upper_bound=zmax,
            )
            return
        if abs(sub.zmax - a.zmin) <= tol and x_overlap > tol and y_overlap > tol:
            maybe_add_contact(
                a,
                sub,
                face_a="zmin",
                face_b="zmax",
                area=x_overlap * y_overlap,
                plane_value=float(a.zmin),
                lower_bound=zmin,
                upper_bound=zmax,
            )
            return

    for i, a in enumerate(add_boxes):
        for b in add_boxes[i + 1 :]:
            add_box_pair_contact(a, b)

    for a in add_boxes:
        for sub in subtract_boxes:
            add_add_to_subtract_contact(a, sub)

    return candidates

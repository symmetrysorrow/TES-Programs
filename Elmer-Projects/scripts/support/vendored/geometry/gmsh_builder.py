"""Vendored 2026-07-14 from Thermal-and-Electoric-Sim core/geometry/gmsh_builder.py
(import paths adjusted only). See vendored/__init__.py for policy."""
from __future__ import annotations

from pathlib import Path
from collections import defaultdict
import gmsh
import numpy as np

from scripts.support.vendored.geometry.body_semantics import body_mode_of, body_name_of, uses_body_boolean
from scripts.support.vendored.geometry.contact_detection import ContactCandidate, detect_contact_candidates
from scripts.support.vendored.geometry.primitives import Box
from scripts.support.vendored.geometry.spec import StructureSpec


def _contacting_body_pairs(spec: StructureSpec) -> list[tuple[str, str]]:
    boxes_by_name = {
        str(box.name).strip(): box
        for box in spec.boxes
        if str(box.name).strip()
    }
    pairs: set[tuple[str, str]] = set()
    for candidate in detect_contact_candidates(spec):
        box_a = boxes_by_name.get(str(candidate.box_a).strip())
        box_b = boxes_by_name.get(str(candidate.box_b).strip())
        if box_a is None or box_b is None:
            continue
        body_a = body_name_of(box_a)
        body_b = body_name_of(box_b)
        if not body_a or not body_b or body_a == body_b:
            continue
        pair = tuple(sorted((body_a, body_b)))
        pairs.add(pair)
    return sorted(pairs)


def _contacting_body_components(spec: StructureSpec) -> list[tuple[str, ...]]:
    pairs = _contacting_body_pairs(spec)
    if not pairs:
        return []

    adjacency: dict[str, set[str]] = defaultdict(set)
    for body_a, body_b in pairs:
        adjacency[body_a].add(body_b)
        adjacency[body_b].add(body_a)

    components: list[tuple[str, ...]] = []
    visited: set[str] = set()
    for root in sorted(adjacency):
        if root in visited:
            continue
        stack = [root]
        members: set[str] = set()
        while stack:
            body = stack.pop()
            if body in visited:
                continue
            visited.add(body)
            members.add(body)
            stack.extend(sorted(adjacency.get(body, ()), reverse=True))
        if len(members) > 1:
            components.append(tuple(sorted(members)))
    return sorted(components)


class GmshApiBuilder:
    """
    StructureSpec から Gmsh Python API で直接メッシュを作る。

    BooleanFragments 後の volume を bbox overlap ratio で汎用分類する。
    """

    def __init__(self, *, spec: StructureSpec, verbose: bool = False) -> None:
        self.spec = spec
        self.verbose = bool(verbose)
        self.primitive_tags: dict[str, int] = {}
        self.fragment_sources: dict[int, set[str]] = {}
        self.group_entities: dict[str, list[int]] = {}
        self.build_boxes: list[Box] = []
        self._box_source_names: dict[str, set[str]] = {}

    @staticmethod
    def _box_key(box: Box) -> str:
        uid = str(getattr(box, "uid", "") or "").strip()
        if uid:
            return uid
        return f"id:{id(box)}"

    def build(
        self,
        *,
        model_name: str = "simcore_model",
        conformal_shared_interfaces: bool = False,
        refinement_balls: list[tuple[float, float, float, float, float]] | None = None,
        partition_planes: list[dict] | None = None,
    ) -> None:
        """Build and mesh the geometry.

        refinement_balls: list of (cx, cy, cz, radius, h_fine) tuples.
            Inside each ball of the given radius centred at (cx,cy,cz) the
            target element size is h_fine; outside it falls back to mesh_max.

        partition_planes: list of dicts with keys:
            axis   – "x", "y", or "z"
            value  – coordinate of the plane [metres]
            bodies – (optional) list of body names to embed in; default = all
        Embedded surfaces force conforming mesh faces at the given positions,
        enabling accurate plane_heat_flow() measurements without remeshing
        from a modified CAD model.
        """
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1 if self.verbose else 0)
        gmsh.option.setNumber("General.Verbosity", 5 if self.verbose else 0)
        gmsh.model.add(model_name)
        self._configure_occ_boolean_options()
        self.build_boxes = self._prepare_build_boxes()

        if self._uses_group_boolean_mode():
            self._add_primitives()
            self._build_group_boolean_geometry()
            if conformal_shared_interfaces:
                self._fragment_touching_body_components(_contacting_body_components(self.spec))
        else:
            self._add_boxes()
            self._fragment_if_needed()

        gmsh.model.occ.synchronize()

        self._add_physical_volumes()
        self._add_physical_surfaces()

        gmsh.model.occ.synchronize()

        if partition_planes:
            self._embed_partition_planes(partition_planes)

        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", self.spec.mesh_min)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", self.spec.mesh_max)
        # Frontal (4) is more robust than the default Delaunay (1) for the
        # thin partially overlapping TES/Stycast geometry produced here.
        gmsh.option.setNumber("Mesh.Algorithm3D", 4)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

        field_ids: list[int] = []
        if refinement_balls:
            for fid, (cx, cy, cz, radius, h_fine) in enumerate(refinement_balls, start=1):
                gmsh.model.mesh.field.add("Ball", fid)
                gmsh.model.mesh.field.setNumber(fid, "Radius", radius)
                gmsh.model.mesh.field.setNumber(fid, "VIn", h_fine)
                gmsh.model.mesh.field.setNumber(fid, "VOut", self.spec.mesh_max)
                gmsh.model.mesh.field.setNumber(fid, "XCenter", cx)
                gmsh.model.mesh.field.setNumber(fid, "YCenter", cy)
                gmsh.model.mesh.field.setNumber(fid, "ZCenter", cz)
                field_ids.append(fid)

        next_field_id = len(field_ids) + 1
        for xmin, xmax, ymin, ymax, zmin, zmax, h_fine in self._thin_box_refinements():
            field_id = next_field_id
            next_field_id += 1
            gmsh.model.mesh.field.add("Box", field_id)
            gmsh.model.mesh.field.setNumber(field_id, "VIn", h_fine)
            gmsh.model.mesh.field.setNumber(field_id, "VOut", self.spec.mesh_max)
            gmsh.model.mesh.field.setNumber(field_id, "XMin", xmin)
            gmsh.model.mesh.field.setNumber(field_id, "XMax", xmax)
            gmsh.model.mesh.field.setNumber(field_id, "YMin", ymin)
            gmsh.model.mesh.field.setNumber(field_id, "YMax", ymax)
            gmsh.model.mesh.field.setNumber(field_id, "ZMin", zmin)
            gmsh.model.mesh.field.setNumber(field_id, "ZMax", zmax)
            field_ids.append(field_id)

        shared_interface_refinements = (
            self._shared_interface_refinements() if conformal_shared_interfaces else []
        )
        for xmin, xmax, ymin, ymax, zmin, zmax, h_fine in shared_interface_refinements:
            field_id = next_field_id
            next_field_id += 1
            gmsh.model.mesh.field.add("Box", field_id)
            gmsh.model.mesh.field.setNumber(field_id, "VIn", h_fine)
            gmsh.model.mesh.field.setNumber(field_id, "VOut", self.spec.mesh_max)
            gmsh.model.mesh.field.setNumber(field_id, "XMin", xmin)
            gmsh.model.mesh.field.setNumber(field_id, "XMax", xmax)
            gmsh.model.mesh.field.setNumber(field_id, "YMin", ymin)
            gmsh.model.mesh.field.setNumber(field_id, "YMax", ymax)
            gmsh.model.mesh.field.setNumber(field_id, "ZMin", zmin)
            gmsh.model.mesh.field.setNumber(field_id, "ZMax", zmax)
            field_ids.append(field_id)

        if shared_interface_refinements:
            local_min = min(entry[-1] for entry in shared_interface_refinements)
            gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min(float(self.spec.mesh_min), local_min))

        if field_ids:
            min_id = next_field_id
            gmsh.model.mesh.field.add("Min", min_id)
            gmsh.model.mesh.field.setNumbers(min_id, "FieldsList", field_ids)
            gmsh.model.mesh.field.setAsBackgroundMesh(min_id)

        gmsh.model.mesh.generate(3)

    def write(self, path: str | Path) -> None:
        gmsh.write(str(path))

    def finalize(self) -> None:
        gmsh.finalize()

    def _configure_occ_boolean_options(self) -> None:
        # `Glue = 2` keeps tangent/shared interfaces from producing tiny sliver
        # fragments in partial-face contacts such as TES/Stycast.
        gmsh.option.setNumber("Geometry.OCCBooleanGlue", 2)

    def _uses_group_boolean_mode(self) -> bool:
        if not self.spec.boxes:
            return False
        for box in self.spec.boxes:
            if box.kind.lower().strip() != "box":
                return True
            if uses_body_boolean(box):
                return True
        return False

    def _add_primitives(self) -> None:
        self._configure_occ_boolean_options()
        for box in self.build_boxes:
            tag = self._add_primitive(box)
            self.primitive_tags[self._box_key(box)] = tag

        gmsh.model.occ.synchronize()

    def _add_boxes(self) -> None:
        """
        Compatibility path for the classic non-grouped box workflow.

        This keeps the legacy/simple geometry path alive so existing smoke tests
        and template-driven projects can still generate a mesh without the
        grouped boolean pipeline.
        """
        self._configure_occ_boolean_options()
        for box in self.build_boxes:
            tag = self._add_primitive(box)
            self.primitive_tags[self._box_key(box)] = tag

        gmsh.model.occ.synchronize()

    def _source_names_for_box(self, box: Box) -> set[str]:
        return set(self._box_source_names.get(self._box_key(box), {box.name}))

    def _prepare_build_boxes(self) -> list[Box]:
        boxes = list(self.spec.boxes)
        self._box_source_names = {
            self._box_key(box): {box.name}
            for box in boxes
        }
        if not self._uses_group_boolean_mode():
            boxes = self._split_boxes_for_partial_face_contacts(boxes)
        for box in boxes:
            self._box_source_names.setdefault(self._box_key(box), {box.name})
        return boxes

    def _embed_partition_planes(self, planes: list[dict]) -> None:
        """Embed flat surfaces into volumes to force conforming mesh at those planes.

        Creates one plane surface per (plane, volume) pair, strictly inside each
        volume's bounding box, then embeds it so that gmsh generates conforming
        faces there during mesh.generate().
        """
        EPS = 1e-12   # shrink surface slightly so it is strictly inside the volume

        # Prefer group_entities (populated by group-boolean path); fall back to
        # querying physical groups which are set by _add_physical_volumes() in both
        # the group and non-group code paths.
        all_body_vols: dict[str, list[int]] = dict(self.group_entities)
        if not all_body_vols:
            for pg_dim, pg_tag in gmsh.model.getPhysicalGroups(3):
                pg_name = gmsh.model.getPhysicalName(3, pg_tag)
                vol_tags = list(gmsh.model.getEntitiesForPhysicalGroup(3, pg_tag))
                if pg_name and vol_tags:
                    all_body_vols[pg_name] = vol_tags
        if not all_body_vols:
            # Last resort: embed into every 3-D entity
            all_tags = [t for _d, t in gmsh.model.getEntities(3)]
            all_body_vols = {"__all__": all_tags}

        for spec in planes:
            axis  = str(spec.get("axis", "x")).lower().strip()
            value = float(spec["value"])
            target_bodies = list(spec.get("bodies", sorted(all_body_vols.keys())))

            target_vol_tags: list[int] = []
            for b in target_bodies:
                target_vol_tags.extend(all_body_vols.get(b, []))

            for vt in target_vol_tags:
                try:
                    bnd = gmsh.model.occ.getBoundingBox(3, vt)
                    vxmn, vymn, vzmn, vxmx, vymx, vzmx = bnd

                    # Skip if plane does not intersect this volume
                    if axis == "x" and not (vxmn - EPS < value < vxmx + EPS):
                        continue
                    if axis == "y" and not (vymn - EPS < value < vymx + EPS):
                        continue
                    if axis == "z" and not (vzmn - EPS < value < vzmx + EPS):
                        continue

                    # Build 4 corner points strictly inside the volume
                    if axis == "x":
                        pts = [
                            gmsh.model.occ.addPoint(value, vymn + EPS, vzmn + EPS),
                            gmsh.model.occ.addPoint(value, vymx - EPS, vzmn + EPS),
                            gmsh.model.occ.addPoint(value, vymx - EPS, vzmx - EPS),
                            gmsh.model.occ.addPoint(value, vymn + EPS, vzmx - EPS),
                        ]
                    elif axis == "y":
                        pts = [
                            gmsh.model.occ.addPoint(vxmn + EPS, value, vzmn + EPS),
                            gmsh.model.occ.addPoint(vxmx - EPS, value, vzmn + EPS),
                            gmsh.model.occ.addPoint(vxmx - EPS, value, vzmx - EPS),
                            gmsh.model.occ.addPoint(vxmn + EPS, value, vzmx - EPS),
                        ]
                    else:
                        pts = [
                            gmsh.model.occ.addPoint(vxmn + EPS, vymn + EPS, value),
                            gmsh.model.occ.addPoint(vxmx - EPS, vymn + EPS, value),
                            gmsh.model.occ.addPoint(vxmx - EPS, vymx - EPS, value),
                            gmsh.model.occ.addPoint(vxmn + EPS, vymx - EPS, value),
                        ]

                    lines = [
                        gmsh.model.occ.addLine(pts[0], pts[1]),
                        gmsh.model.occ.addLine(pts[1], pts[2]),
                        gmsh.model.occ.addLine(pts[2], pts[3]),
                        gmsh.model.occ.addLine(pts[3], pts[0]),
                    ]
                    loop = gmsh.model.occ.addCurveLoop(lines)
                    surf = gmsh.model.occ.addPlaneSurface([loop])
                    gmsh.model.occ.synchronize()
                    gmsh.model.mesh.embed(2, [surf], 3, vt)

                except Exception as exc:  # pragma: no cover
                    import warnings
                    warnings.warn(
                        f"partition plane {axis}={value:.4e} in vol {vt}: {exc}",
                        stacklevel=3,
                    )

    def _thin_box_refinements(self) -> list[tuple[float, float, float, float, float, float, float]]:
        if not self.build_boxes:
            return []
        refinements: list[tuple[float, float, float, float, float, float, float]] = []
        global_max = max(float(self.spec.mesh_max), 1.0e-12)
        for box in self.build_boxes:
            if body_mode_of(box) == "subtract":
                continue
            dims = np.asarray([float(box.dx), float(box.dy), float(box.dz)], dtype=float)
            positive_dims = dims[dims > 0.0]
            if positive_dims.size == 0:
                continue
            min_dim = float(np.min(positive_dims))
            max_dim = float(np.max(positive_dims))
            if min_dim >= 0.5 * global_max and max_dim <= 2.0 * global_max:
                continue
            h_fine = min(global_max, max(min_dim * 0.75, global_max / 128.0, float(self.spec.mesh_min)))
            if h_fine >= global_max:
                continue
            refinements.append(
                (
                    float(box.xmin),
                    float(box.xmax),
                    float(box.ymin),
                    float(box.ymax),
                    float(box.zmin),
                    float(box.zmax),
                    float(h_fine),
                )
            )
        return refinements

    def _shared_interface_refinements(self) -> list[tuple[float, float, float, float, float, float, float]]:
        boxes_by_name = {str(box.name): box for box in self.spec.boxes}
        h_fine = max(min(float(self.spec.mesh_min) / 8.0, float(self.spec.mesh_max) / 16.0), 1.0e-9)
        padding = 3.0 * h_fine
        refinements: list[tuple[float, float, float, float, float, float, float]] = []
        seen: set[tuple[float, ...]] = set()

        for candidate in detect_contact_candidates(self.spec):
            box_a = boxes_by_name.get(str(candidate.box_a))
            box_b = boxes_by_name.get(str(candidate.box_b))
            if box_a is None or box_b is None:
                continue
            if "Membrane" not in {body_name_of(box_a), body_name_of(box_b)}:
                continue

            xmin = max(float(box_a.xmin), float(box_b.xmin))
            xmax = min(float(box_a.xmax), float(box_b.xmax))
            ymin = max(float(box_a.ymin), float(box_b.ymin))
            ymax = min(float(box_a.ymax), float(box_b.ymax))
            zmin = max(float(box_a.zmin), float(box_b.zmin))
            zmax = min(float(box_a.zmax), float(box_b.zmax))
            face_a = str(candidate.face_a).strip().lower()
            if face_a.startswith("x"):
                plane = float(box_a.xmin if face_a.endswith("min") else box_a.xmax)
                xmin, xmax = plane - padding, plane + padding
                ymin, ymax = ymin - padding, ymax + padding
                zmin, zmax = zmin - padding, zmax + padding
            elif face_a.startswith("y"):
                plane = float(box_a.ymin if face_a.endswith("min") else box_a.ymax)
                xmin, xmax = xmin - padding, xmax + padding
                ymin, ymax = plane - padding, plane + padding
                zmin, zmax = zmin - padding, zmax + padding
            else:
                # Keep ultra-thin TES contacts on the robust fallback contact
                # path; this refinement targets the lateral membrane exits.
                continue

            if xmax <= xmin or ymax <= ymin or zmax <= zmin:
                continue
            entry = (xmin, xmax, ymin, ymax, zmin, zmax, h_fine)
            key = tuple(round(value, 15) for value in entry)
            if key in seen:
                continue
            seen.add(key)
            refinements.append(entry)
        return refinements

    def _split_boxes_for_partial_face_contacts(self, boxes: list[Box]) -> list[Box]:
        boxes_by_name = {
            str(box.name).strip(): box
            for box in boxes
            if str(box.name).strip()
        }
        split_points: dict[str, dict[str, set[float]]] = defaultdict(
            lambda: {"x": set(), "y": set(), "z": set()}
        )
        tol = 1.0e-12

        for candidate in detect_contact_candidates(StructureSpec(boxes=boxes)):
            box_a = boxes_by_name.get(str(candidate.box_a).strip())
            box_b = boxes_by_name.get(str(candidate.box_b).strip())
            if box_a is None or box_b is None:
                continue
            if box_a.kind.lower().strip() != "box" or box_b.kind.lower().strip() != "box":
                continue
            if body_mode_of(box_a) == "subtract" or body_mode_of(box_b) == "subtract":
                continue
            self._accumulate_contact_split_points(split_points, candidate, box_a, box_b, tol=tol)

        out: list[Box] = []
        new_source_names: dict[str, set[str]] = {}
        for box in boxes:
            if body_mode_of(box) == "subtract" or box.kind.lower().strip() != "box":
                out.append(box)
                new_source_names[self._box_key(box)] = set(self._source_names_for_box(box))
                continue
            cuts = split_points.get(box.name, {})
            parts = self._split_box_by_coordinates(
                box,
                x_cuts=sorted(cuts.get("x", set())),
                y_cuts=sorted(cuts.get("y", set())),
                z_cuts=sorted(cuts.get("z", set())),
                tol=tol,
            )
            out.extend(parts)
            source_names = set(self._source_names_for_box(box))
            for part in parts:
                new_source_names[self._box_key(part)] = set(source_names)

        self._box_source_names = new_source_names
        return out

    def _accumulate_contact_split_points(
        self,
        split_points: dict[str, dict[str, set[float]]],
        candidate: ContactCandidate,
        box_a: Box,
        box_b: Box,
        *,
        tol: float,
    ) -> None:
        face_a = str(candidate.face_a).lower().strip()
        face_b = str(candidate.face_b).lower().strip()
        if face_a.startswith("z") and face_b.startswith("z"):
            overlap_x = (max(box_a.xmin, box_b.xmin), min(box_a.xmax, box_b.xmax))
            overlap_y = (max(box_a.ymin, box_b.ymin), min(box_a.ymax, box_b.ymax))
            self._add_interior_split(split_points[box_a.name]["x"], overlap_x, box_a.xmin, box_a.xmax, tol=tol)
            self._add_interior_split(split_points[box_a.name]["y"], overlap_y, box_a.ymin, box_a.ymax, tol=tol)
            self._add_interior_split(split_points[box_b.name]["x"], overlap_x, box_b.xmin, box_b.xmax, tol=tol)
            self._add_interior_split(split_points[box_b.name]["y"], overlap_y, box_b.ymin, box_b.ymax, tol=tol)
            return
        if face_a.startswith("x") and face_b.startswith("x"):
            overlap_y = (max(box_a.ymin, box_b.ymin), min(box_a.ymax, box_b.ymax))
            overlap_z = (max(box_a.zmin, box_b.zmin), min(box_a.zmax, box_b.zmax))
            self._add_interior_split(split_points[box_a.name]["y"], overlap_y, box_a.ymin, box_a.ymax, tol=tol)
            self._add_interior_split(split_points[box_a.name]["z"], overlap_z, box_a.zmin, box_a.zmax, tol=tol)
            self._add_interior_split(split_points[box_b.name]["y"], overlap_y, box_b.ymin, box_b.ymax, tol=tol)
            self._add_interior_split(split_points[box_b.name]["z"], overlap_z, box_b.zmin, box_b.zmax, tol=tol)
            return
        if face_a.startswith("y") and face_b.startswith("y"):
            overlap_x = (max(box_a.xmin, box_b.xmin), min(box_a.xmax, box_b.xmax))
            overlap_z = (max(box_a.zmin, box_b.zmin), min(box_a.zmax, box_b.zmax))
            self._add_interior_split(split_points[box_a.name]["x"], overlap_x, box_a.xmin, box_a.xmax, tol=tol)
            self._add_interior_split(split_points[box_a.name]["z"], overlap_z, box_a.zmin, box_a.zmax, tol=tol)
            self._add_interior_split(split_points[box_b.name]["x"], overlap_x, box_b.xmin, box_b.xmax, tol=tol)
            self._add_interior_split(split_points[box_b.name]["z"], overlap_z, box_b.zmin, box_b.zmax, tol=tol)

    def _add_interior_split(
        self,
        target: set[float],
        interval: tuple[float, float],
        box_min: float,
        box_max: float,
        *,
        tol: float,
    ) -> None:
        start, end = interval
        if end - start <= tol:
            return
        if box_min + tol < start < box_max - tol:
            target.add(float(start))
        if box_min + tol < end < box_max - tol:
            target.add(float(end))

    def _split_box_by_coordinates(
        self,
        box: Box,
        *,
        x_cuts: list[float],
        y_cuts: list[float],
        z_cuts: list[float],
        tol: float,
    ) -> list[Box]:
        x_edges = [box.xmin] + [value for value in x_cuts if box.xmin + tol < value < box.xmax - tol] + [box.xmax]
        y_edges = [box.ymin] + [value for value in y_cuts if box.ymin + tol < value < box.ymax - tol] + [box.ymax]
        z_edges = [box.zmin] + [value for value in z_cuts if box.zmin + tol < value < box.zmax - tol] + [box.zmax]
        x_edges = sorted(set(float(value) for value in x_edges))
        y_edges = sorted(set(float(value) for value in y_edges))
        z_edges = sorted(set(float(value) for value in z_edges))

        if len(x_edges) == 2 and len(y_edges) == 2 and len(z_edges) == 2:
            return [box]

        parts: list[Box] = []
        for ix, (xmin, xmax) in enumerate(zip(x_edges[:-1], x_edges[1:]), start=1):
            dx = xmax - xmin
            if dx <= tol:
                continue
            for iy, (ymin, ymax) in enumerate(zip(y_edges[:-1], y_edges[1:]), start=1):
                dy = ymax - ymin
                if dy <= tol:
                    continue
                for iz, (zmin, zmax) in enumerate(zip(z_edges[:-1], z_edges[1:]), start=1):
                    dz = zmax - zmin
                    if dz <= tol:
                        continue
                    parts.append(
                        Box(
                            name=f"{box.name}__split_{ix}_{iy}_{iz}",
                            x=(xmin + xmax) / 2.0,
                            y=(ymin + ymax) / 2.0,
                            z=(zmin + zmax) / 2.0,
                            dx=dx,
                            dy=dy,
                            dz=dz,
                            kind=box.kind,
                            body_name=box.body_name,
                            group_name=box.group_name,
                            body_mode=box.body_mode,
                            group_mode=box.group_mode,
                            priority=box.priority,
                        )
                    )
        return parts or [box]

    def _add_primitive(self, box: Box) -> int:
        kind = box.kind.lower().strip()
        if kind == "box":
            return gmsh.model.occ.addBox(box.xmin, box.ymin, box.zmin, box.dx, box.dy, box.dz)
        if kind == "wedge":
            return gmsh.model.occ.addWedge(box.xmin, box.ymin, box.zmin, box.dx, box.dy, box.dz, ltx=box.ltx)
        raise NotImplementedError(f"Unsupported primitive kind: {box.kind}")

    def _build_full_thickness_frame_boxes(self, add_box: Box, sub_box: Box) -> list[Box] | None:
        if add_box.kind.lower().strip() != "box" or sub_box.kind.lower().strip() != "box":
            return None

        eps = 1.0e-12
        if abs(add_box.zmin - sub_box.zmin) > eps or abs(add_box.zmax - sub_box.zmax) > eps:
            return None

        if not (
            add_box.xmin + eps < sub_box.xmin < sub_box.xmax < add_box.xmax - eps
            and add_box.ymin + eps < sub_box.ymin < sub_box.ymax < add_box.ymax - eps
        ):
            return None

        frame_boxes: list[Box] = []

        def make_box(name_suffix: str, xmin: float, xmax: float, ymin: float, ymax: float) -> None:
            dx = xmax - xmin
            dy = ymax - ymin
            if dx <= eps or dy <= eps:
                return
            frame_boxes.append(
                Box(
                    name=f"{add_box.name}__frame_{name_suffix}",
                    x=(xmin + xmax) / 2.0,
                    y=(ymin + ymax) / 2.0,
                    z=add_box.z,
                    dx=dx,
                    dy=dy,
                    dz=add_box.dz,
                    kind="box",
                    body_name=add_box.body_name,
                    group_name=add_box.group_name,
                    body_mode="add",
                    group_mode="add",
                    priority=add_box.priority,
                )
            )

        make_box("left", add_box.xmin, sub_box.xmin, add_box.ymin, add_box.ymax)
        make_box("right", sub_box.xmax, add_box.xmax, add_box.ymin, add_box.ymax)
        make_box("bottom", sub_box.xmin, sub_box.xmax, add_box.ymin, sub_box.ymin)
        make_box("top", sub_box.xmin, sub_box.xmax, sub_box.ymax, add_box.ymax)

        return frame_boxes or None

    def _build_group_boolean_geometry(self) -> None:
        """
        Build grouped solids with per-group union/difference.

        A boolean body is defined by `box.body_name`, with legacy fallback to
        `group_name` and finally `box.name`.
        The body mode comes from `body_mode`, with legacy fallback to
        `group_mode`.
        - `body_mode='add'`: part of the solid
        - `body_mode='subtract'`: cut from the solid

        This is intentionally limited but good enough for TES-style conductor
        blocks and simple voids.
        """
        groups: dict[str, list[Box]] = defaultdict(list)
        for box in self.build_boxes:
            groups[body_name_of(box)].append(box)

        self.group_entities = {}
        for group_name, members in groups.items():
            add_members = [box for box in members if body_mode_of(box) != "subtract"]
            sub_members = [box for box in members if body_mode_of(box) == "subtract"]
            add_tags = [self.primitive_tags[self._box_key(box)] for box in add_members]
            sub_tags = [self.primitive_tags[self._box_key(box)] for box in sub_members]

            if not add_tags:
                continue

            entities: list[tuple[int, int]]
            if len(add_members) == 1 and len(sub_members) == 1:
                frame_boxes = self._build_full_thickness_frame_boxes(add_members[0], sub_members[0])
                if frame_boxes:
                    gmsh.model.occ.remove(
                        [(3, add_tags[0]), (3, sub_tags[0])],
                        recursive=True,
                    )
                    frame_tags = [self._add_primitive(frame_box) for frame_box in frame_boxes]
                    gmsh.model.occ.synchronize()
                    entities = [(3, tag) for tag in frame_tags]
                    sub_tags = []
                else:
                    entities = [(3, tag) for tag in add_tags]
            else:
                entities = [(3, tag) for tag in add_tags]

            if len(entities) > 1:
                # Make overlapping add parts disjoint before cutting holes.
                out_entities, _ = gmsh.model.occ.fragment(entities, [])
                gmsh.model.occ.synchronize()
                entities = [entity for entity in out_entities if entity[0] == 3]

            if sub_tags:
                cut_out, _ = gmsh.model.occ.cut(
                    entities,
                    [(3, tag) for tag in sub_tags],
                    removeObject=True,
                    removeTool=True,
                )
                gmsh.model.occ.synchronize()
                entities = [entity for entity in cut_out if entity[0] == 3]

            final_tags = sorted({tag for dim, tag in entities if dim == 3})
            self.group_entities[group_name] = final_tags
            for tag in final_tags:
                self.fragment_sources.setdefault(tag, set()).add(group_name)
                for box in members:
                    self.fragment_sources[tag].update(self._source_names_for_box(box))

        # The per-body build above already preserves coincident interfaces for
        # contacting solids well enough for meshing. A second fragment pass over
        # touching body chains can collapse ultra-thin layers (for example a TES
        # film sandwiched between SiNx and Stycast) into zero-mass OCC volumes,
        # which later disappear from the mesh entirely.
        #
        # Keep the helper available for future targeted use, but do not run it
        # unconditionally here.
        return

    def _fragment_touching_body_components(self, body_components: list[tuple[str, ...]]) -> None:
        for component in body_components:
            body_tags: dict[str, list[int]] = {
                body_name: [int(tag) for tag in self.group_entities.get(body_name, [])]
                for body_name in component
            }
            if any(not tags for tags in body_tags.values()):
                continue

            entities = [
                (3, tag)
                for body_name in component
                for tag in body_tags[body_name]
            ]
            input_sources = {
                int(tag): set(self.fragment_sources.get(int(tag), set()))
                for _, tag in entities
            }
            input_body_name: dict[int, str] = {}
            for body_name, tags in body_tags.items():
                for tag in tags:
                    input_body_name[int(tag)] = body_name

            _, fragment_map = gmsh.model.occ.fragment(entities, [])
            gmsh.model.occ.synchronize()

            new_tags_by_body: dict[str, list[int]] = defaultdict(list)
            new_sources_by_tag: dict[int, set[str]] = {}
            for (dim, input_tag), result_entities in zip(entities, fragment_map):
                if dim != 3:
                    continue
                source_names = set(input_sources.get(int(input_tag), set()))
                body_name = str(input_body_name.get(int(input_tag), "") or "").strip()
                for out_dim, out_tag in result_entities:
                    if out_dim != 3:
                        continue
                    out_tag_i = int(out_tag)
                    if body_name:
                        new_tags_by_body[body_name].append(out_tag_i)
                    if source_names:
                        new_sources_by_tag.setdefault(out_tag_i, set()).update(source_names)

            for body_name in component:
                if new_tags_by_body.get(body_name):
                    self.group_entities[body_name] = sorted(set(new_tags_by_body[body_name]))

            for tags in body_tags.values():
                for old_tag in tags:
                    self.fragment_sources.pop(int(old_tag), None)
            for new_tag, names in new_sources_by_tag.items():
                self.fragment_sources[int(new_tag)] = set(names)

    def _fragment_if_needed(self) -> None:
        if not self.spec.use_boolean_fragments:
            return

        volumes = [(3, tag) for tag in self.primitive_tags.values()]
        if len(volumes) < 2:
            return

        _, fragment_map = gmsh.model.occ.fragment(volumes, [])
        gmsh.model.occ.synchronize()

        tag_to_name: dict[int, str] = {}
        for box in self.build_boxes:
            primitive_tag = self.primitive_tags.get(self._box_key(box))
            if primitive_tag is None:
                continue
            source_names = self._source_names_for_box(box)
            if source_names:
                tag_to_name[primitive_tag] = sorted(source_names)[0]
        self.fragment_sources = {}

        for (dim, input_tag), result_entities in zip(volumes, fragment_map):
            if dim != 3:
                continue

            source_name = tag_to_name.get(input_tag)
            if source_name is None:
                continue

            for out_dim, out_tag in result_entities:
                if out_dim != 3:
                    continue
                self.fragment_sources.setdefault(out_tag, set()).add(source_name)

    def _add_physical_volumes(self) -> None:
        """
        BooleanFragments 後の volume を primitive に割り当てる。

        汎用ルール:
        1. volume bbox と primitive box の bbox 重なり体積を計算
        2. overlap > 0 の primitive を候補にする
        3. overlap ratio を最優先
        4. ratio が同程度なら priority を優先
        5. さらに同じなら primitive 体積が小さいものを優先
        """

        if self.group_entities:
            for pv in self.spec.physical_volumes:
                selected_tags: list[int] = []
                for primitive_name in pv.primitive_names:
                    selected_tags.extend(self.group_entities.get(primitive_name, []))
                selected_tags = sorted(set(selected_tags))
                if not selected_tags:
                    if self.verbose:
                        print(
                            f"WARNING: physical volume '{pv.name}' has no assigned grouped volumes."
                        )
                    continue
                group = gmsh.model.addPhysicalGroup(3, selected_tags, tag=pv.tag)
                gmsh.model.setPhysicalName(3, group, pv.name)
            return

        volumes = gmsh.model.getEntities(dim=3)

        volume_assignment: dict[int, str] = {}
        unassigned_volumes: list[int] = []

        for dim, vol_tag in volumes:
            candidate_names = self.fragment_sources.get(vol_tag)
            candidates: list[tuple[Box, float, float]] = []

            if candidate_names:
                for box in self.spec.boxes:
                    if box.name in candidate_names:
                        candidates.append((box, 1.0, 1.0))
            else:
                x, y, z = gmsh.model.occ.getCenterOfMass(dim, vol_tag)

                for box in self.spec.boxes:
                    if box.contains_point(x, y, z):
                        candidates.append((box, 1.0, 1.0))

                if not candidates:
                    bbox = gmsh.model.getBoundingBox(dim, vol_tag)
                    vol_bbox_volume = self._bbox_volume(bbox)
                    for box in self.spec.boxes:
                        overlap = self._bbox_overlap_volume(bbox, box)
                        if overlap > 0.0:
                            ratio = overlap / max(vol_bbox_volume, 1.0e-30)
                            candidates.append((box, overlap, ratio))

            if not candidates:
                unassigned_volumes.append(vol_tag)
                continue

            candidates.sort(
                key=lambda item: (
                    -item[2],          # overlap ratio
                    -item[0].priority,
                    item[0].volume,
                    item[0].name,
                )
            )

            volume_assignment[vol_tag] = candidates[0][0].name

        if unassigned_volumes and self.verbose:
            print("WARNING: unassigned volumes:", sorted(unassigned_volumes))

        for pv in self.spec.physical_volumes:
            selected_tags = [
                vol_tag
                for vol_tag, primitive_name in volume_assignment.items()
                if primitive_name in pv.primitive_names
            ]

            selected_tags = sorted(set(selected_tags))

            if not selected_tags:
                if self.verbose:
                    print(
                        f"WARNING: physical volume '{pv.name}' has no assigned volumes."
                    )
                continue

            group = gmsh.model.addPhysicalGroup(
                3,
                selected_tags,
                tag=pv.tag,
            )
            gmsh.model.setPhysicalName(3, group, pv.name)

    def _bbox_volume(
        self,
        bbox: tuple[float, float, float, float, float, float],
    ) -> float:
        xmin, ymin, zmin, xmax, ymax, zmax = bbox
        return (
            max(0.0, xmax - xmin)
            * max(0.0, ymax - ymin)
            * max(0.0, zmax - zmin)
        )

    def _bbox_overlap_volume(
        self,
        bbox: tuple[float, float, float, float, float, float],
        box: Box,
        *,
        eps: float = 1.0e-12,
    ) -> float:
        xmin, ymin, zmin, xmax, ymax, zmax = bbox

        ox = max(0.0, min(xmax, box.xmax) - max(xmin, box.xmin))
        oy = max(0.0, min(ymax, box.ymax) - max(ymin, box.ymin))
        oz = max(0.0, min(zmax, box.zmax) - max(zmin, box.zmin))

        vol = ox * oy * oz
        if vol <= eps:
            return 0.0
        return vol

    def _add_physical_surfaces(self) -> None:
        for ps in self.spec.physical_surfaces:
            if ps.selector not in {"all_boundary_surfaces", "box_surface"}:
                raise NotImplementedError(
                    f"Unsupported surface selector: {ps.selector}"
                )

            if ps.selector == "all_boundary_surfaces":
                outer_volume_tag = self._find_outer_volume_tag(reference_box_name=ps.box_name)
                if outer_volume_tag is None:
                    if self.verbose:
                        print(
                            f"WARNING: physical surface '{ps.name}' has no outer volume."
                        )
                    continue

                surfaces = gmsh.model.getBoundary(
                    [(3, outer_volume_tag)],
                    oriented=False,
                    recursive=False,
                )
                surface_tags = sorted({tag for dim, tag in surfaces if dim == 2})
            else:
                surface_tags = self._find_box_surface_tags(ps.box_name, ps.surface)

            if surface_tags:
                group = gmsh.model.addPhysicalGroup(
                    2,
                    surface_tags,
                    tag=ps.tag,
                )
                gmsh.model.setPhysicalName(2, group, ps.name)

    def _find_outer_volume_tag(self, reference_box_name: str | None = None) -> int | None:
        if not self.spec.boxes:
            return None

        outer_name = reference_box_name or self.spec.boxes[0].name
        volumes = gmsh.model.getEntities(dim=3)
        if not volumes:
            return None

        exact_matches: list[int] = []
        fallback_matches: list[tuple[float, int]] = []

        for _, vol_tag in volumes:
            source_names = self.fragment_sources.get(vol_tag, set())
            if outer_name not in source_names:
                continue

            if source_names == {outer_name}:
                exact_matches.append(vol_tag)
                continue

            bbox = gmsh.model.getBoundingBox(3, vol_tag)
            fallback_matches.append((self._bbox_volume(bbox), vol_tag))

        if exact_matches:
            return max(
                exact_matches,
                key=lambda vol_tag: self._bbox_volume(gmsh.model.getBoundingBox(3, vol_tag)),
            )

        if fallback_matches:
            return max(fallback_matches, key=lambda item: item[0])[1]

        return None

    def _find_box_surface_tags(self, box_name: str | None, surface: str | None) -> list[int]:
        if box_name is None or surface is None:
            return []

        surface = surface.lower().strip()
        target_box = next((box for box in self.spec.boxes if box.name == box_name), None)
        if target_box is None:
            if self.verbose:
                print(f"WARNING: physical surface box '{box_name}' not found.")
            return []

        volumes = gmsh.model.getEntities(dim=3)
        if not volumes:
            return []

        selected_volumes: list[int] = []
        for _, vol_tag in volumes:
            source_names = self.fragment_sources.get(vol_tag, set())
            if box_name in source_names:
                selected_volumes.append(vol_tag)

        if not selected_volumes:
            return []

        # Use getCenterOfMass rather than getBoundingBox for the normal-axis
        # comparison: for a flat face the centroid coordinate is exact, while
        # getBoundingBox pads by ~100 nm which exceeds the TES layer thickness.
        # Use a generous eps for the in-plane axes (those bboxes are fine).
        eps_tight = max(min(target_box.dx, target_box.dy, target_box.dz) * 1.0e-3, 2.0e-9)
        surface_tags: set[int] = set()

        def matches_surface(surf_tag: int) -> bool:
            cx, cy, cz = gmsh.model.occ.getCenterOfMass(2, surf_tag)
            if surface == "xmin":
                return abs(cx - target_box.xmin) <= eps_tight
            if surface == "xmax":
                return abs(cx - target_box.xmax) <= eps_tight
            if surface == "ymin":
                return abs(cy - target_box.ymin) <= eps_tight
            if surface == "ymax":
                return abs(cy - target_box.ymax) <= eps_tight
            if surface == "zmin":
                return abs(cz - target_box.zmin) <= eps_tight
            if surface == "zmax":
                return abs(cz - target_box.zmax) <= eps_tight
            return False

        for vol_tag in selected_volumes:
            surfaces = gmsh.model.getBoundary(
                [(3, vol_tag)],
                oriented=False,
                recursive=False,
            )
            for dim, surf_tag in surfaces:
                if dim != 2:
                    continue
                if matches_surface(surf_tag):
                    surface_tags.add(surf_tag)

        return sorted(surface_tags)

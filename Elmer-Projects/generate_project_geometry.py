import json
import os
import uuid
from pathlib import Path
import sys
from math import pi

import gmsh



ROOT = Path(__file__).resolve().parent
ELMER_PROJECT = ROOT / "elmer_project.json"
sys.path.insert(0, str(ROOT))

from scripts.support.vendored.geometry.body_semantics import body_name_of  # noqa: E402
from scripts.support.vendored.geometry.gmsh_builder import GmshApiBuilder  # noqa: E402
from scripts.support.vendored.geometry.loader import load_structure_spec  # noqa: E402
from scripts.support.vendored.geometry.primitives import Box  # noqa: E402
from scripts.support.vendored.geometry.structure import _add_all_box_surfaces  # noqa: E402
from scripts.support.reconcile_project import reconcile_project  # noqa: E402


_original_set_number = gmsh.option.setNumber

# Seam-rotation angle (radians) for contact-disc imprint tools; any generic
# angle works, it only needs to keep circle mesh nodes off the Stycast rim's.
_CONTACT_DISC_SEAM_ROTATION = 0.082


class TesGmshBuilder(GmshApiBuilder):
    """Project-specific extensions of the vendored builder.

    These used to be monkey patches on the external class; they are ordinary
    overrides now that the builder is vendored.
    """

    def _add_primitive(self, primitive):
        if primitive.kind.lower().strip() == "cylinder":
            return gmsh.model.occ.addCylinder(
                primitive.x,
                primitive.y,
                primitive.zmin,
                0.0,
                0.0,
                primitive.dz,
                primitive.dx / 2.0,
            )
        disc = (getattr(self, "contact_film_discs", None) or {}).get(primitive.name)
        if disc is not None:
            return self._add_contact_film_primitive(primitive, disc)
        return super()._add_primitive(primitive)

    def _add_contact_film_primitive(self, primitive, disc):
        """Build a thin film box pre-split along a contact-disc footprint.

        Any OCC 3D boolean that touches the 0.16 um TES film corrupts the
        solid (the film-height cylindrical wall's edges collapse within the
        OCC tolerance; the codebase already knows the film "can produce
        zero-mass volumes" under fragmentation). So instead of cutting the
        built box, split its FOOTPRINT in 2D (square fragmented by a disk -
        healthy in-plane dimensions only) and extrude both faces in a single
        call: the shared circle edge then extrudes into a single shared
        lateral wall, giving two conformally connected volumes (core disc +
        ring) whose top/bottom faces are exactly the contact disc and its
        complement.

        Returns the core volume tag (the builder tracks one tag per
        primitive); the ring volume is stashed in _extra_primitive_volumes
        and merged into the body's group in _build_group_boolean_geometry.
        """
        cx, cy, radius = disc
        rectangle = gmsh.model.occ.addRectangle(
            primitive.xmin, primitive.ymin, primitive.zmin, primitive.dx, primitive.dy
        )
        disk = gmsh.model.occ.addDisk(cx, cy, primitive.zmin, radius, radius)
        # Rotate the disk's parameterization seam off the +x axis: the Stycast
        # cylinder rims are identical circles in the same planes, and with the
        # default seam both circles mesh onto identical 1D nodes, which
        # ElmerGrid then merges (1e-10) into direct node coupling across the
        # mortar gap. The single-pixel reference mesh has zero shared nodes on
        # these contacts; the rotation (geometry unchanged) restores that.
        if not getattr(self, "merge_mortar_interfaces", False):
            gmsh.model.occ.rotate([(2, disk)], cx, cy, primitive.zmin, 0.0, 0.0, 1.0, _CONTACT_DISC_SEAM_ROTATION)
        # OCCBooleanGlue=2 (the global default here) skips real boolean
        # computation for interfering shapes; disable it for this 2D split.
        _original_set_number("Geometry.OCCBooleanGlue", 0)
        try:
            out, _ = gmsh.model.occ.fragment([(2, rectangle)], [(2, disk)])
        finally:
            _original_set_number("Geometry.OCCBooleanGlue", 2)
        faces = [(dim, tag) for dim, tag in out if dim == 2]
        if len(faces) != 2:
            raise RuntimeError(
                f"contact film {primitive.name!r}: footprint split produced "
                f"{len(faces)} faces (expected disc + ring)"
            )
        extruded = gmsh.model.occ.extrude(faces, 0.0, 0.0, primitive.dz)
        volumes = [tag for dim, tag in extruded if dim == 3]
        if len(volumes) != 2:
            raise RuntimeError(
                f"contact film {primitive.name!r}: extrusion produced "
                f"{len(volumes)} volumes (expected 2)"
            )
        if not hasattr(self, "_extra_primitive_volumes"):
            self._extra_primitive_volumes = {}
        body = body_name_of(primitive)
        self._extra_primitive_volumes.setdefault(body, []).extend(volumes[1:])
        return volumes[0]

    def _build_group_boolean_geometry(self):
        super()._build_group_boolean_geometry()
        # Glue the membrane island bodies to the substrate stack so their shared
        # faces (island<->frame side walls, island bottom<->SiO2_1 sheet) are
        # meshed conformally. Without this the island only touches the frame at a
        # few coincident corner nodes after ElmerGrid node merging, adding a large
        # spurious constriction resistance. abs/TES/Stycast stay out of the
        # fragment: their interfaces are coupled with mortar BCs, and fragmenting
        # the thin overlapping TES/Stycast films can produce zero-mass volumes.
        # "TES"/"Stycast" cover both the single_pixel names and the dual_tes
        # "_L"/"_R" suffixed names (TES_L, TES_R, Stycast_L, Stycast_R).
        stack = tuple(
            name
            for name in self.group_entities
            if name != "abs" and not name.startswith(("TES", "Stycast"))
        )
        if len(stack) > 1:
            # OCCBooleanGlue only handles fully coincident faces; the island faces
            # are subsets of the frame/SiO2_1 faces and need real imprinting, so
            # disable gluing for this fragment only.
            _original_set_number("Geometry.OCCBooleanGlue", 0)
            try:
                self._fragment_touching_body_components([stack])
            finally:
                _original_set_number("Geometry.OCCBooleanGlue", 2)
        self._merge_extra_primitive_volumes()
        self._imprint_contact_discs()

    def _merge_extra_primitive_volumes(self):
        """Fold the ring volumes created by _add_contact_film_primitive into
        their body's group so physical volume/surface assignment sees them."""
        extras = getattr(self, "_extra_primitive_volumes", None) or {}
        for body, extra_tags in extras.items():
            if body not in self.group_entities:
                raise RuntimeError(f"contact film: body {body!r} missing from group entities")
            merged = sorted(set(self.group_entities[body]) | set(extra_tags))
            self.group_entities[body] = merged
            sources: set[str] = set()
            for tag in merged:
                sources |= self.fragment_sources.get(tag, set())
            for tag in merged:
                self.fragment_sources.setdefault(tag, set()).update(sources or {body})

    def _imprint_contact_discs(self):
        """Imprint Stycast contact-disc footprints onto body faces BEFORE
        meshing, so the mesh conforms to the circular mortar contact patches.

        Driven by `self.contact_disc_specs` (set by build(); empty for the
        single_pixel geometry, whose legacy post-mesh code path must stay
        bit-identical): a list of {"body": <group name>, "discs":
        [(x, y, z, radius), ...]} entries. Each body is fragmented with flat
        disk tools lying in its face plane - a pure 2D imprint that splits
        the face into disc + remainder without cutting the volume (the thick
        absorber would survive a 3D cut, but a face imprint adds no interior
        walls). Tool leftovers are removed; group bookkeeping is updated.

        This replaces the intent of the legacy `_fragment_mortar_interfaces`
        (generate_project_geometry.py, post-builder), which never worked: its
        tool cylinders were built at literal y=0 while the geometry sits at
        y=1mm, and it ran after gmsh.model.mesh.generate() where a real
        geometry change would have destroyed the existing mesh.
        """
        specs = getattr(self, "contact_disc_specs", None) or []
        for entry in specs:
            body = entry["body"]
            body_tags = list(self.group_entities.get(body, []))
            if not body_tags:
                raise RuntimeError(f"contact disc imprint: no volumes for body {body!r}")
            objects = [(3, tag) for tag in body_tags]
            tools = []
            for x, y, z, radius in entry["discs"]:
                disk = gmsh.model.occ.addDisk(x, y, z, radius, radius)
                # Same seam rotation as _add_contact_film_primitive: keep the
                # imprinted circle's mesh nodes off the Stycast rim nodes so
                # ElmerGrid does not merge them across the mortar gap.
                if not getattr(self, "merge_mortar_interfaces", False):
                    gmsh.model.occ.rotate([(2, disk)], x, y, z, 0.0, 0.0, 1.0, _CONTACT_DISC_SEAM_ROTATION)
                tools.append((2, disk))
            # OCCBooleanGlue=2 (the builder's global setting) skips the real
            # boolean computation, so the disk tools would not imprint at all;
            # disable it here like the membrane-stack fragment above does.
            _original_set_number("Geometry.OCCBooleanGlue", 0)
            try:
                _, out_map = gmsh.model.occ.fragment(objects, tools)
            finally:
                _original_set_number("Geometry.OCCBooleanGlue", 2)
            gmsh.model.occ.synchronize()
            n_objects = len(objects)
            new_body_tags = sorted(
                {int(tag) for row in out_map[:n_objects] for dim, tag in row if dim == 3}
            )
            # Remove any tool faces that were not absorbed into the body
            # boundary (they would otherwise be meshed as floating surfaces).
            body_faces = {
                int(face_tag)
                for tag in new_body_tags
                for dim, face_tag in gmsh.model.getBoundary(
                    [(3, tag)], oriented=False, recursive=False
                )
                if dim == 2
            }
            leftovers = sorted(
                {
                    int(tag)
                    for row in out_map[n_objects:]
                    for dim, tag in row
                    if dim == 2 and int(tag) not in body_faces
                }
            )
            if leftovers:
                gmsh.model.occ.remove([(2, tag) for tag in leftovers], recursive=True)
                gmsh.model.occ.synchronize()
            sources: set[str] = set()
            for tag in body_tags:
                sources |= self.fragment_sources.pop(tag, set())
            self.group_entities[body] = new_body_tags
            for tag in new_body_tags:
                self.fragment_sources.setdefault(tag, set()).update(sources or {body})

    def _before_mesh_generate(self):
        """Copy meshes across coincident contact faces for Elmer node merging."""
        if not getattr(self, "merge_mortar_interfaces", False):
            return

        def horizontal_faces(body: str, z_target: float) -> list[int]:
            faces: set[int] = set()
            for volume in self.group_entities.get(body, []):
                for dim, face in gmsh.model.getBoundary(
                    [(3, volume)], oriented=False, recursive=False
                ):
                    if dim != 2:
                        continue
                    if abs(gmsh.model.occ.getCenterOfMass(2, face)[2] - z_target) <= 5.0e-8:
                        faces.add(int(face))
            return sorted(faces)

        identity = [
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        ]

        def match_faces(slave: list[int], master: list[int], label: str) -> None:
            unused = set(master)
            if len(slave) != len(master) or not slave:
                raise RuntimeError(
                    f"{label}: contact face count mismatch {len(slave)} vs {len(master)}"
                )
            for slave_face in slave:
                slave_area = abs(gmsh.model.occ.getMass(2, slave_face))
                master_face = min(
                    unused,
                    key=lambda face: abs(abs(gmsh.model.occ.getMass(2, face)) - slave_area),
                )
                master_area = abs(gmsh.model.occ.getMass(2, master_face))
                if abs(master_area - slave_area) > 1.0e-6 * max(master_area, slave_area):
                    raise RuntimeError(
                        f"{label}: unmatched face areas {slave_area} and {master_area}"
                    )
                gmsh.model.mesh.setPeriodic(2, [slave_face], [master_face], identity)
                unused.remove(master_face)

        for slave_body, slave_z, master_body, master_z, target_area, region, label in getattr(
            self, "conformal_contact_pairs", []
        ):
            slave_faces = horizontal_faces(slave_body, slave_z)
            master_faces = horizontal_faces(master_body, master_z)
            if region is not None:
                xmin, xmax, ymin, ymax = region
                slave_faces = [
                    face for face in slave_faces
                    if xmin <= gmsh.model.occ.getCenterOfMass(2, face)[0] <= xmax
                    and ymin <= gmsh.model.occ.getCenterOfMass(2, face)[1] <= ymax
                ]
                master_faces = [
                    face for face in master_faces
                    if xmin <= gmsh.model.occ.getCenterOfMass(2, face)[0] <= xmax
                    and ymin <= gmsh.model.occ.getCenterOfMass(2, face)[1] <= ymax
                ]
            if target_area is not None:
                slave_faces = [
                    face for face in slave_faces
                    if abs(abs(gmsh.model.occ.getMass(2, face)) - target_area)
                    <= 1.0e-6 * target_area
                ]
                master_faces = [
                    face for face in master_faces
                    if abs(abs(gmsh.model.occ.getMass(2, face)) - target_area)
                    <= 1.0e-6 * target_area
                ]
            match_faces(
                slave_faces,
                master_faces,
                label,
            )

    def _thin_box_refinements(self):
        # Optional mesh-convergence hook: set MEMBRANE_REFINE_H (metres) to add a
        # refinement box around the membrane island + SiO2_1 sheet region.
        entries = list(super()._thin_box_refinements())
        h = os.environ.get("MEMBRANE_REFINE_H")
        if h:
            entries.append((-0.6e-3, 0.6e-3, 0.4e-3, 1.6e-3, 1.70e-4, 1.95e-4, float(h)))
            # The builder clamps sizes from below with Mesh.CharacteristicLengthMin;
            # lower it so the requested refinement actually takes effect.
            _original_set_number("Mesh.CharacteristicLengthMin", float(h))

        # Optional mesh-convergence hook: set TES_LOCAL_REFINE_H (metres) to
        # tighten the mesh in a narrow XY column around each TES_<side>
        # contact point, spanning the FULL Z stack (Stycast down through the
        # bath-facing SiO2_2 layer) -- not just the TES/Stycast film pair.
        # A first attempt that only covered TES+Stycast (element counts
        # matched exactly against a full-domain 2x refinement) still left
        # ~1.0% L/R current asymmetry vs 2x's ~0.2%, showing the substrate
        # conduction path's own mesh asymmetry also contributes. Keeping the
        # XY footprint narrow (not the full 3x6mm chip) avoids paying for a
        # full background refinement while still symmetrizing the whole
        # vertical heat path under each TES.
        h_local = os.environ.get("TES_LOCAL_REFINE_H")
        if h_local:
            h_local = float(h_local)
            pad_xy = 0.4e-3
            stack_z: dict[str, list[float]] = {}
            stack_xy: dict[str, tuple[float, float]] = {}
            for box in self.build_boxes:
                name = str(getattr(box, "name", "") or "")
                _base, suffix = _strip_known_suffix(name, ["_L", "_R"])
                if not suffix:
                    continue
                zmin = float(box.z) - float(box.dz) / 2.0
                zmax = float(box.z) + float(box.dz) / 2.0
                lo, hi = stack_z.get(suffix, (zmin, zmax))
                stack_z[suffix] = [min(lo, zmin), max(hi, zmax)]
                if str(_base) == "TES":
                    stack_xy[suffix] = (float(box.x), float(box.y))
            for suffix, (zmin, zmax) in stack_z.items():
                if suffix not in stack_xy:
                    continue
                cx, cy = stack_xy[suffix]
                entries.append((
                    cx - pad_xy, cx + pad_xy,
                    cy - pad_xy, cy + pad_xy,
                    zmin, zmax,
                    h_local,
                ))
            if stack_z:
                _original_set_number("Mesh.CharacteristicLengthMin", min(
                    float(self.spec.mesh_min), h_local
                ))
        return entries


def _set_number_without_optimize(name, value):
    return _original_set_number(name, value)


def _find_box(spec, name: str) -> Box:
    for box in spec.boxes:
        if box.name == name:
            return box
    raise ValueError(f"Box '{name}' not found in structure spec")


# --- multi-side (single_pixel vs. dual_tes) geometry support ---------------
#
# The project JSON's `geometry` tree (injected by build_mesh.py from the
# `geometries` registry) is either a single TES stack (single_pixel: body
# names "TES", "Stycast", ... with no suffix) or two stacks either side of a
# long absorber (dual_tes: "TES_L"/"TES_R", "Stycast_L"/"Stycast_R", ...).
# `build()` below runs the same per-stack processing once per "side" so the
# single-pixel path (side suffix "") is untouched byte-for-byte.

_SIDE_ROLE_ORDER = [
    "TES",
    "Stycast",
    "Membrane_SiNx",
    "SiO2_1",
    "Si_1",
    "SiNx",
    "Si_2",
    "SiO2_2",
    "Membrane_Si1",
]


def _find_geometry_leaf(nodes: list, name: str):
    """Recursively find a leaf (childless) geometry-tree node by name."""
    for node in nodes:
        children = node.get("children")
        if isinstance(children, list) and children:
            found = _find_geometry_leaf(children, name)
            if found is not None:
                return found
        elif node.get("name") == name:
            return node
    return None


def _strip_known_suffix(name: str, suffixes: list[str]) -> tuple[str, str]:
    """Split a body name into (base, suffix) for the first non-empty suffix
    in *suffixes* that it ends with; ("name", "") if none match."""
    for suffix in suffixes:
        if suffix and name.endswith(suffix):
            return name[: -len(suffix)], suffix
    return name, ""


def _rounded_bbox(dim: int, tag: int) -> list[float]:
    return [round(value, 12) for value in gmsh.model.getBoundingBox(dim, tag)]


def _rounded_com(dim: int, tag: int) -> list[float]:
    return [round(value, 12) for value in gmsh.model.occ.getCenterOfMass(dim, tag)]


def _physical_groups_by_name(dim: int) -> dict[str, tuple[int, list[int]]]:
    groups: dict[str, tuple[int, list[int]]] = {}
    for _, physical_tag in gmsh.model.getPhysicalGroups(dim):
        name = gmsh.model.getPhysicalName(dim, physical_tag)
        entities = [int(entity) for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, physical_tag)]
        groups[name] = (physical_tag, entities)
    return groups


def _replace_surface_group(name: str, physical_tag: int, entities: list[int]) -> None:
    existing = {(dim, tag) for dim, tag in gmsh.model.getPhysicalGroups(2)}
    if (2, physical_tag) in existing:
        gmsh.model.removePhysicalGroups([(2, physical_tag)])
    gmsh.model.addPhysicalGroup(2, entities, physical_tag, name=name)


def _replace_volume_group(name: str, physical_tag: int, entities: list[int]) -> None:
    existing = {(dim, tag) for dim, tag in gmsh.model.getPhysicalGroups(3)}
    if (3, physical_tag) in existing:
        gmsh.model.removePhysicalGroups([(3, physical_tag)])
    try:
        gmsh.model.addPhysicalGroup(3, entities, physical_tag, name=name)
    except Exception as exc:
        if "already exists" not in str(exc):
            raise
        gmsh.model.removePhysicalGroups([(3, physical_tag)])
        gmsh.model.addPhysicalGroup(3, entities, physical_tag, name=name)


def _merge_surface_groups(name: str, physical_tag: int, source_names: list[str]) -> None:
    surface_groups = _physical_groups_by_name(2)
    entities: list[int] = []
    removal_tags: list[int] = []
    for source_name in source_names:
        group = surface_groups.get(source_name)
        if group is None:
            continue
        removal_tags.append(group[0])
        entities.extend(group[1])
    if not entities:
        return
    existing = {(dim, tag) for dim, tag in gmsh.model.getPhysicalGroups(2)}
    for tag in removal_tags:
        if (2, tag) in existing and tag != physical_tag:
            gmsh.model.removePhysicalGroups([(2, tag)])
    _replace_surface_group(name, physical_tag, sorted(set(int(entity) for entity in entities)))


def _remove_surface_groups_by_name(names: list[str]) -> None:
    surface_groups = _physical_groups_by_name(2)
    removals = [
        (2, surface_groups[name][0])
        for name in names
        if name in surface_groups
    ]
    if removals:
        gmsh.model.removePhysicalGroups(removals)


def _upward_volumes(surface_tag: int) -> list[int]:
    upward, _ = gmsh.model.getAdjacencies(2, surface_tag)
    return [int(tag) for tag in upward]


def _surface_area(surface_tag: int) -> float:
    return abs(gmsh.model.occ.getMass(2, surface_tag))


def _surface_zmid(surface_tag: int) -> float:
    bbox = gmsh.model.getBoundingBox(2, surface_tag)
    return 0.5 * (bbox[2] + bbox[5])


def _find_surfaces_by_plane(
    *,
    zmid_target: float,
    allowed_upward_volumes: set[int],
    upward_count: int = 1,
    z_tol: float = 1.0e-6,
    flat_tol: float = 5.0e-7,
) -> list[int]:
    matches: list[int] = []
    for _, surface_tag in gmsh.model.getEntities(2):
        upward_volumes = _upward_volumes(surface_tag)
        if len(upward_volumes) != upward_count:
            continue
        if not set(upward_volumes).issubset(allowed_upward_volumes):
            continue
        bbox = gmsh.model.getBoundingBox(2, surface_tag)
        if abs(bbox[5] - bbox[2]) > flat_tol:
            continue
        if abs(0.5 * (bbox[2] + bbox[5]) - zmid_target) > z_tol:
            continue
        matches.append(int(surface_tag))
    return sorted(matches)


def _find_volume_tags_for_boxes(boxes: list[Box]) -> list[int]:
    matches: set[int] = set()
    for _, volume_tag in gmsh.model.getEntities(3):
        x, y, z = gmsh.model.occ.getCenterOfMass(3, volume_tag)
        for box in boxes:
            if box.contains_point(x, y, z, eps=1.0e-9):
                matches.add(int(volume_tag))
                break
    return sorted(matches)


def _find_surface_by_plane_and_area(
    *,
    zmid_target: float,
    area_target: float,
    allowed_upward_volumes: set[int] | None = None,
    upward_count: int = 1,
    z_tol: float = 1.0e-6,
    rel_area_tol: float = 1.5e-1,
) -> int:
    best_tag: int | None = None
    best_score: float | None = None
    fallback_tag: int | None = None
    fallback_score: float | None = None

    for _, surface_tag in gmsh.model.getEntities(2):
        upward_volumes = _upward_volumes(surface_tag)
        if len(upward_volumes) < upward_count:
            continue

        bbox = gmsh.model.getBoundingBox(2, surface_tag)
        zmid = 0.5 * (bbox[2] + bbox[5])
        if abs(zmid - zmid_target) > z_tol:
            continue

        if allowed_upward_volumes is not None:
            overlap = set(upward_volumes).intersection(allowed_upward_volumes)
            if not overlap:
                continue

        area = _surface_area(surface_tag)
        area_error = abs(area - area_target) / area_target
        score = abs(zmid - zmid_target) + area_error
        if fallback_score is None or score < fallback_score:
            fallback_score = score
            fallback_tag = int(surface_tag)

        if area_error > rel_area_tol:
            continue

        if best_score is None or score < best_score:
            best_score = score
            best_tag = int(surface_tag)

    if best_tag is None:
        if fallback_tag is not None:
            return fallback_tag
        raise RuntimeError(
            f"Could not find surface near z={zmid_target:.12g} with area {area_target:.12g}"
        )

    return best_tag


def _split_box_xy_ring(
    *,
    source: Box,
    inner_dx: float,
    inner_dy: float,
    name_prefix: str,
    uid_prefix: str,
) -> list[Box]:
    outer_dx = float(source.dx)
    outer_dy = float(source.dy)
    margin_x = 0.5 * (outer_dx - inner_dx)
    margin_y = 0.5 * (outer_dy - inner_dy)
    if margin_x < 0.0 or margin_y < 0.0:
        raise ValueError(f"Inner box larger than outer box for {name_prefix}")

    common = dict(
        z=source.z,
        dz=source.dz,
        kind=source.kind,
        body_name=source.body_name,
        group_name=source.group_name,
        ltx=source.ltx,
        z_expr=source.z_expr,
        dz_expr=source.dz_expr,
        ltx_expr=source.ltx_expr,
        priority=source.priority,
    )

    parts = [
        Box(
            name=f"{name_prefix}_contact",
            uid=f"{uid_prefix}_contact",
            x=source.x,
            y=source.y,
            dx=inner_dx,
            dy=inner_dy,
            x_expr=source.x_expr,
            y_expr=source.y_expr,
            dx_expr=None,
            dy_expr=None,
            **common,
        ),
        Box(
            name=f"{name_prefix}_free_left",
            uid=f"{uid_prefix}_free_left",
            x=source.x - 0.5 * inner_dx - 0.5 * margin_x,
            y=source.y,
            dx=margin_x,
            dy=outer_dy,
            x_expr=None,
            y_expr=source.y_expr,
            dx_expr=None,
            dy_expr=source.dy_expr,
            **common,
        ),
        Box(
            name=f"{name_prefix}_free_right",
            uid=f"{uid_prefix}_free_right",
            x=source.x + 0.5 * inner_dx + 0.5 * margin_x,
            y=source.y,
            dx=margin_x,
            dy=outer_dy,
            x_expr=None,
            y_expr=source.y_expr,
            dx_expr=None,
            dy_expr=source.dy_expr,
            **common,
        ),
        Box(
            name=f"{name_prefix}_free_bottom",
            uid=f"{uid_prefix}_free_bottom",
            x=source.x,
            y=source.y - 0.5 * inner_dy - 0.5 * margin_y,
            dx=inner_dx,
            dy=margin_y,
            x_expr=source.x_expr,
            y_expr=None,
            dx_expr=None,
            dy_expr=None,
            **common,
        ),
        Box(
            name=f"{name_prefix}_free_top",
            uid=f"{uid_prefix}_free_top",
            x=source.x,
            y=source.y + 0.5 * inner_dy + 0.5 * margin_y,
            dx=inner_dx,
            dy=margin_y,
            x_expr=source.x_expr,
            y_expr=None,
            dx_expr=None,
            dy_expr=None,
            **common,
        ),
    ]
    return parts


def _write_mortar_surface_report(out_path: Path) -> None:
    report: dict[str, object] = {"surface_groups": {}, "volume_groups": {}, "surface_entities": []}

    surface_groups = _physical_groups_by_name(2)
    for name, (physical_tag, entities) in surface_groups.items():
        report["surface_groups"][name] = {
            "physical_tag": physical_tag,
            "entities": entities,
        }

    volume_groups = _physical_groups_by_name(3)
    for name, (physical_tag, entities) in volume_groups.items():
        report["volume_groups"][name] = {
            "physical_tag": physical_tag,
            "entities": entities,
        }

    for _, surface_tag in gmsh.model.getEntities(2):
        upward, downward = gmsh.model.getAdjacencies(2, surface_tag)
        report["surface_entities"].append(
            {
                "tag": surface_tag,
                "bbox": _rounded_bbox(2, surface_tag),
                "center_of_mass": _rounded_com(2, surface_tag),
                "upward_volumes": [int(tag) for tag in upward],
                "downward_edges": [int(tag) for tag in downward],
            }
        )

    out_path.parent.mkdir(exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


def _fragment_mortar_interfaces(enabled: bool) -> None:
    if not enabled:
        return

    volume_groups = _physical_groups_by_name(3)
    required = ("abs", "TES")
    if any(name not in volume_groups for name in required):
        missing = [name for name in required if name not in volume_groups]
        raise RuntimeError(f"Missing volume groups for mortar splitting: {missing}")

    abs_volume = volume_groups["abs"][1][0]
    tes_volume = volume_groups["TES"][1][0]

    tools = [
        (
            "TES",
            101,
            [(3, tes_volume)],
            [
                (
                    3,
                    gmsh.model.occ.addCylinder(
                        0.0,
                        0.0,
                        191.9e-6,
                        0.0,
                        0.0,
                        0.4e-6,
                        249.0e-6,
                    ),
                )
            ],
        ),
        (
            "abs",
            100,
            [(3, abs_volume)],
            [
                (
                    3,
                    gmsh.model.occ.addCylinder(
                        0.0,
                        0.0,
                        212.06e-6,
                        0.0,
                        0.0,
                        700.2e-6,
                        249.0e-6,
                    ),
                )
            ],
        ),
    ]

    for name, physical_tag, objects, tool in tools:
        _, out_dim_tags_map = gmsh.model.occ.fragment(objects, tool)
        gmsh.model.occ.synchronize()
        split_volumes = sorted(
            {
                int(tag)
                for dim, tag in out_dim_tags_map[0]
                if dim == 3
            }
        )
        if not split_volumes:
            raise RuntimeError(f"No fragment volumes found for {name}")
        _replace_volume_group(name, physical_tag, split_volumes)


def _retag_mortar_interface_surfaces(
    enabled: bool,
    membrane_sinx_parts: list[Box],
    *,
    name_prefix: str = "Membrane_SiNx",
    volume_tag: int = 103,
    face_base: int = 1300,
    legacy_contact: bool = True,
    conformal_shared_interfaces: bool = False,
) -> None:
    """Re-tag the Membrane_SiNx (island) volume + zmin/zmax surface groups from
    their transient per-zone pieces into the final semantic groups.

    Parameterized by *name_prefix*/*volume_tag*/*face_base* so this can run
    once per geometry side (single_pixel: defaults, unchanged; dual_tes:
    "Membrane_SiNx_L"/"_R" with that side's tags) without touching the
    single-pixel code path's literal output.
    """
    membrane_sinx_volume_tags = _find_volume_tags_for_boxes(membrane_sinx_parts)
    _replace_volume_group(name_prefix, volume_tag, membrane_sinx_volume_tags)
    volume_groups = _physical_groups_by_name(3)
    membrane_sinx_volumes = set(membrane_sinx_volume_tags)
    if conformal_shared_interfaces:
        def touching_faces_at_z(z_target: float) -> list[int]:
            matches: list[int] = []
            for _, surface_tag in gmsh.model.getEntities(2):
                if not membrane_sinx_volumes.intersection(_upward_volumes(surface_tag)):
                    continue
                bbox = gmsh.model.getBoundingBox(2, surface_tag)
                if abs(bbox[5] - bbox[2]) > 5.0e-7:
                    continue
                if abs(_surface_zmid(surface_tag) - z_target) <= 2.5e-7:
                    matches.append(int(surface_tag))
            return sorted(matches)

        membrane_bottom_surfaces = touching_faces_at_z(191.0e-6)
        membrane_top_surfaces = touching_faces_at_z(192.0e-6)
    else:
        membrane_bottom_surfaces = _find_surfaces_by_plane(
            zmid_target=191.0e-6,
            allowed_upward_volumes=membrane_sinx_volumes,
            z_tol=2.5e-7,
        )
        membrane_top_surfaces = _find_surfaces_by_plane(
            zmid_target=192.0e-6,
            allowed_upward_volumes=membrane_sinx_volumes,
            z_tol=2.5e-7,
        )
    membrane_contact_surface = min(
        membrane_top_surfaces,
        key=lambda tag: abs(_surface_area(tag) - 500.0e-6 * 500.0e-6),
    )
    membrane_free_surfaces = [tag for tag in membrane_top_surfaces if tag != membrane_contact_surface]
    _replace_surface_group(f"{name_prefix}__zmin", face_base + 4, membrane_bottom_surfaces)
    _replace_surface_group(f"{name_prefix}__zmax", face_base + 5, [membrane_contact_surface])
    _replace_surface_group(f"{name_prefix}__zmax_free", face_base + 6, membrane_free_surfaces)
    _remove_surface_groups_by_name(
        [
            f"{name_prefix}_contact__zmin",
            f"{name_prefix}_contact__zmax",
            f"{name_prefix}_free_left__zmin",
            f"{name_prefix}_free_left__zmax",
            f"{name_prefix}_free_right__zmin",
            f"{name_prefix}_free_right__zmax",
            f"{name_prefix}_free_bottom__zmin",
            f"{name_prefix}_free_bottom__zmax",
            f"{name_prefix}_free_top__zmin",
            f"{name_prefix}_free_top__zmax",
        ]
    )

    surface_groups = _physical_groups_by_name(2)
    if not enabled or not legacy_contact:
        return

    tes_volumes = set(volume_groups["TES"][1])
    abs_volumes = set(volume_groups["abs"][1])
    tes_stycast_contact = _find_surface_by_plane_and_area(
        zmid_target=192.16e-6,
        area_target=pi * (249.0e-6**2),
        allowed_upward_volumes=tes_volumes,
    )
    abs_stycast_contact = _find_surface_by_plane_and_area(
        zmid_target=212.16e-6,
        area_target=pi * (249.0e-6**2),
        allowed_upward_volumes=abs_volumes,
    )

    _replace_surface_group(
        "TES__zmax",
        surface_groups.get("TES__zmax", (1105, []))[0],
        [tes_stycast_contact],
    )
    _replace_surface_group(
        "abs__zmin",
        surface_groups.get("abs__zmin", (1004, []))[0],
        [abs_stycast_contact],
    )

    optional_free_groups = (
        ("TES__zmax_free", 1106, 192.16e-6, 500.0e-6 * 500.0e-6 - pi * (249.0e-6**2)),
        ("abs__zmin_free", 1006, 212.16e-6, 1.0e-6 - pi * (249.0e-6**2)),
    )
    for name, tag, zmid_target, area_target in optional_free_groups:
        try:
            free_surface = _find_surface_by_plane_and_area(
                zmid_target=zmid_target,
                area_target=area_target,
                allowed_upward_volumes={
                    "TES__zmax_free": tes_volumes,
                    "abs__zmin_free": abs_volumes,
                }[name],
            )
        except RuntimeError:
            continue
        gmsh.model.addPhysicalGroup(2, [free_surface], tag, name=name)


def _find_body_faces_at_z(volumes: set[int], z_target: float, z_tol: float) -> list[int]:
    """Boundary faces (exactly one parent volume, parent in *volumes*) whose
    center-of-mass z is within *z_tol* of *z_target*. Center-of-mass is exact
    for flat faces, unlike getBoundingBox which pads by ~100 nm (larger than
    the 0.16 um TES film thickness)."""
    matches: list[int] = []
    for _, surface_tag in gmsh.model.getEntities(2):
        upward = _upward_volumes(surface_tag)
        if len(upward) != 1 or upward[0] not in volumes:
            continue
        _, _, cz = gmsh.model.occ.getCenterOfMass(2, surface_tag)
        if abs(cz - z_target) <= z_tol:
            matches.append(int(surface_tag))
    return sorted(matches)


def _retag_contact_surfaces_multi(
    enabled: bool,
    *,
    sides: list[str],
    spec,
    stycast_boxes: dict[str, Box],
    face_bases: dict[str, int],
) -> None:
    """Multi-side (dual_tes) version of the Stycast-contact surface retag.

    The geometry was already split along the contact-disc footprints before
    meshing (TesGmshBuilder._fragment_contact_discs), so here the split faces
    only need relabelling into contact/free physical groups, classified by
    area (the disc and its complement ring share the same center of mass, so
    position cannot distinguish them):
    - TES_<side>__zmax      = the contact disc (area = pi*r^2)
    - TES_<side>__zmax_free = the rest of the TES top face
    - abs__zmin             = both contact discs (one per side; Elmer mortar
                              BCs let both Stycast_<side>__zmax slaves target
                              this single master boundary)
    - abs__zmin_free        = the rest of the absorber bottom face
    """
    if not enabled:
        return
    volume_groups = _physical_groups_by_name(3)
    surface_groups = _physical_groups_by_name(2)

    def split_by_disc_area(face_tags: list[int], disc_area: float, group_label: str) -> tuple[list[int], list[int]]:
        contact = [
            tag for tag in face_tags if abs(_surface_area(tag) - disc_area) / disc_area < 0.05
        ]
        free = [tag for tag in face_tags if tag not in contact]
        if not contact:
            raise RuntimeError(
                f"{group_label}: no face with the contact-disc area {disc_area:.6g} found "
                f"among {[(tag, _surface_area(tag)) for tag in face_tags]}"
            )
        return sorted(contact), sorted(free)

    for suffix in sides:
        tes_name = f"TES{suffix}"
        tes_box = _find_box(spec, tes_name)
        radius = float(stycast_boxes[suffix].dx) / 2.0
        disc_area = pi * radius * radius
        tes_volumes = set(volume_groups[tes_name][1])
        top_faces = _find_body_faces_at_z(tes_volumes, tes_box.zmax, z_tol=tes_box.dz / 4.0)
        contact, free = split_by_disc_area(top_faces, disc_area, f"{tes_name}__zmax")
        _replace_surface_group(
            f"{tes_name}__zmax",
            surface_groups.get(f"{tes_name}__zmax", (face_bases[tes_name] + 5, []))[0],
            contact,
        )
        _replace_surface_group(f"{tes_name}__zmax_free", face_bases[tes_name] + 6, free)

    abs_box = _find_box(spec, "abs")
    abs_volumes = set(volume_groups["abs"][1])
    bottom_faces = _find_body_faces_at_z(abs_volumes, abs_box.zmin, z_tol=1.0e-6)
    # Every side's disc has the same radius here; classify against each side's
    # disc area so differing radii would also work.
    contact_discs: list[int] = []
    for suffix in sides:
        radius = float(stycast_boxes[suffix].dx) / 2.0
        disc_area = pi * radius * radius
        contact_discs += [
            tag
            for tag in bottom_faces
            if tag not in contact_discs
            and abs(_surface_area(tag) - disc_area) / disc_area < 0.05
        ]
    free_faces = [tag for tag in bottom_faces if tag not in contact_discs]
    if len(contact_discs) != len(sides):
        raise RuntimeError(
            f"abs__zmin: expected {len(sides)} contact discs, found {len(contact_discs)} "
            f"among {[(tag, _surface_area(tag)) for tag in bottom_faces]}"
        )
    _replace_surface_group(
        "abs__zmin",
        surface_groups.get("abs__zmin", (face_bases["abs"] + 4, []))[0],
        sorted(contact_discs),
    )
    _replace_surface_group("abs__zmin_free", face_bases["abs"] + 6, sorted(free_faces))


def _identify_contact_surface_groups() -> dict[str, list[int]]:
    surface_groups = _physical_groups_by_name(2)
    contacts: dict[str, list[int]] = {}

    def zmid(bbox: list[float]) -> float:
        return 0.5 * (bbox[2] + bbox[5])

    for name, (_, entities) in surface_groups.items():
        if name == "Membrane_SiNx__zmax":
            contacts["Membrane_SiNx__zmax_surfaces"] = entities
        elif name == "TES__zmin":
            contacts["TES__zmin_surfaces"] = entities
        elif name == "TES__zmax":
            contacts["TES__zmax_surfaces"] = entities
        elif name == "Stycast__zmin":
            contacts["Stycast__zmin_surfaces"] = entities
        elif name == "Stycast__zmax":
            contacts["Stycast__zmax_surfaces"] = entities
        elif name == "abs__zmin":
            contacts["abs__zmin_surfaces"] = entities

    # Surface groups can still contain several geometric faces under one physical tag.
    # Keep a second view grouped by the Z-plane to make mortar-facing splits easier to inspect.
    z_planes: dict[str, list[int]] = {}
    for name, (_, entities) in surface_groups.items():
        if "__" not in name:
            continue
        if not name.endswith(("__zmin", "__zmax")):
            continue
        plane_key = f"{name}@{round(zmid(_rounded_bbox(2, entities[0])), 12)}"
        z_planes[plane_key] = [int(entity) for entity in entities]
    contacts["z_plane_groups"] = z_planes

    return contacts


def set_input_jsons(project_path: Path, geometry_path: Path | None = None) -> None:
    global ELMER_PROJECT
    ELMER_PROJECT = project_path
    if geometry_path is not None and geometry_path != project_path:
        raise ValueError(
            "A separate geometry JSON is no longer supported; "
            "mesh/role_map now live in the single project JSON."
        )


def build(write_mesh: bool = True) -> None:
    raw_project = json.loads(ELMER_PROJECT.read_text(encoding="utf-8"))
    reconciled_project = reconcile_project(raw_project)
    if reconciled_project != raw_project:
        print(
            "NOTE: stale numeric fields in elmer_project.json; using values "
            "derived from expressions (run sync_elmer_parameters.py --fix-numerics to update the file)"
        )
    raw_project = reconciled_project
    spec = load_structure_spec(reconciled_project)
    geometry_children = raw_project.get("geometry", {}).get("children", [])

    box_names = {box.name for box in spec.boxes}
    if "TES" in box_names:
        sides = [""]
    elif "TES_L" in box_names and "TES_R" in box_names:
        sides = ["_L", "_R"]
    else:
        raise RuntimeError(
            "generate_project_geometry.build(): could not find a 'TES' body "
            "(single_pixel) or 'TES_L'+'TES_R' bodies (dual_tes) in the geometry"
        )

    stycast_names = {f"Stycast{suffix}" for suffix in sides}
    spec.boxes = [
        box
        for box in spec.boxes
        if box.name not in stycast_names and box.dx > 0.0 and box.dy > 0.0 and box.dz > 0.0
    ]
    spec.physical_volumes = [
        volume for volume in spec.physical_volumes if volume.name not in stycast_names
    ]
    spec.physical_surfaces = [
        surface
        for surface in spec.physical_surfaces
        if not any(surface.name.startswith(f"{name}__") for name in stycast_names)
    ]

    elmer_overrides = raw_project.get("elmer_overrides", {})

    stycast_boxes: dict[str, Box] = {}
    for suffix in sides:
        if suffix == "":
            stycast_boxes[suffix] = Box(
                name="Stycast",
                uid="d52e6d3377c94ca5b3e5868ba17a9d39",
                x=elmer_overrides.get("stycast_x", 0.0),
                y=elmer_overrides.get("stycast_y", 0.0),
                z=elmer_overrides.get("stycast_z", 0.00020216),
                dx=elmer_overrides.get("stycast_dx", 498.0e-6),
                dy=elmer_overrides.get("stycast_dy", 498.0e-6),
                dz=elmer_overrides.get("stycast_dz", 20.0e-6),
                kind=elmer_overrides.get("stycast_shape", "cylinder"),
                body_name="Stycast",
                group_name="Stycast",
            )
        else:
            # elmer_overrides only carries one (single-pixel) Stycast position;
            # for dual_tes sides, take the resolved leaf straight from the
            # geometry tree (loader.py skips "cylinder"-shaped leaves, so this
            # is the only place their position/uid are available).
            node = _find_geometry_leaf(geometry_children, f"Stycast{suffix}")
            if node is None:
                raise RuntimeError(f"Stycast{suffix} leaf not found in geometry tree")
            stycast_boxes[suffix] = Box(
                name=f"Stycast{suffix}",
                uid=str(node["uid"]),
                x=float(node["x"]),
                y=float(node["y"]),
                z=float(node["z"]),
                dx=float(node["dx"]),
                dy=float(node["dy"]),
                dz=float(node["dz"]),
                kind=str(node.get("shape") or "cylinder"),
                body_name=f"Stycast{suffix}",
                group_name=f"Stycast{suffix}",
            )
    for suffix in sides:
        spec.boxes.insert(2, stycast_boxes[suffix])

    membrane_names = {f"Membrane{suffix}" for suffix in sides}
    spec.boxes = [box for box in spec.boxes if box.name not in membrane_names]
    spec.physical_volumes = [
        volume for volume in spec.physical_volumes if volume.name not in membrane_names
    ]
    spec.physical_surfaces = [
        surface
        for surface in spec.physical_surfaces
        if not any(surface.name.startswith(f"{name}__") for name in membrane_names)
    ]

    expected_tags: dict[str, int] = {"abs": 100}
    membrane_si1_boxes: dict[str, Box] = {}
    membrane_sinx_boxes: dict[str, Box] = {}
    membrane_sinx_parts_by_side: dict[str, list[Box]] = {}

    for side_index, suffix in enumerate(sides):
        tag_base = 101 + 9 * side_index
        for role_index, role in enumerate(_SIDE_ROLE_ORDER):
            expected_tags[f"{role}{suffix}"] = tag_base + role_index

        membrane_si1_source = _find_box(spec, f"Si_1_sub{suffix}")
        membrane_sinx_source = _find_box(spec, f"SiNx_sub{suffix}")
        tes_source = _find_box(spec, f"TES{suffix}")

        if suffix == "":
            membrane_si1_uid = "59fef5ec76f74e6199f79a8e8c90a203"
            membrane_sinx_uid = "a1b85a1f1f3846cf8fb1d0ce6640ef7f"
        else:
            # Synthesized bodies with no counterpart in the JSON tree; a fresh
            # uid only needs to be unique within this run.
            membrane_si1_uid = uuid.uuid4().hex
            membrane_sinx_uid = uuid.uuid4().hex

        membrane_si1 = Box(
            name=f"Membrane_Si1{suffix}",
            uid=membrane_si1_uid,
            x=membrane_si1_source.x,
            y=membrane_si1_source.y,
            z=membrane_si1_source.z,
            dx=membrane_si1_source.dx,
            dy=membrane_si1_source.dy,
            dz=membrane_si1_source.dz,
            kind="box",
            body_name=f"Membrane_Si1{suffix}",
            group_name=f"Membrane_Si1{suffix}",
            ltx=membrane_si1_source.ltx,
            x_expr=membrane_si1_source.x_expr,
            y_expr=membrane_si1_source.y_expr,
            z_expr=membrane_si1_source.z_expr,
            dx_expr=membrane_si1_source.dx_expr,
            dy_expr=membrane_si1_source.dy_expr,
            dz_expr=membrane_si1_source.dz_expr,
            ltx_expr=membrane_si1_source.ltx_expr,
            priority=membrane_si1_source.priority,
        )
        membrane_sinx = Box(
            name=f"Membrane_SiNx{suffix}",
            uid=membrane_sinx_uid,
            x=membrane_sinx_source.x,
            y=membrane_sinx_source.y,
            z=membrane_sinx_source.z,
            dx=membrane_sinx_source.dx,
            dy=membrane_sinx_source.dy,
            dz=membrane_sinx_source.dz,
            kind="box",
            body_name=f"Membrane_SiNx{suffix}",
            group_name=f"Membrane_SiNx{suffix}",
            ltx=membrane_sinx_source.ltx,
            x_expr=membrane_sinx_source.x_expr,
            y_expr=membrane_sinx_source.y_expr,
            z_expr=membrane_sinx_source.z_expr,
            dx_expr=membrane_sinx_source.dx_expr,
            dy_expr=membrane_sinx_source.dy_expr,
            dz_expr=membrane_sinx_source.dz_expr,
            ltx_expr=membrane_sinx_source.ltx_expr,
            priority=membrane_sinx_source.priority,
        )
        membrane_sinx_parts = _split_box_xy_ring(
            source=membrane_sinx,
            inner_dx=tes_source.dx,
            inner_dy=tes_source.dy,
            name_prefix=f"Membrane_SiNx{suffix}",
            uid_prefix=membrane_sinx.uid,
        )
        spec.boxes.extend([membrane_si1, *membrane_sinx_parts])
        membrane_si1_boxes[suffix] = membrane_si1
        membrane_sinx_boxes[suffix] = membrane_sinx
        membrane_sinx_parts_by_side[suffix] = membrane_sinx_parts

    for volume in spec.physical_volumes:
        volume.tag = expected_tags[volume.name]
    for suffix in sides:
        spec.add_physical_volume(
            name=f"Stycast{suffix}",
            primitive_names=[f"Stycast{suffix}"],
            tag=expected_tags[f"Stycast{suffix}"],
            uid=stycast_boxes[suffix].uid,
        )
        spec.add_physical_volume(
            name=f"Membrane_SiNx{suffix}",
            primitive_names=[box.name for box in membrane_sinx_parts_by_side[suffix]],
            tag=expected_tags[f"Membrane_SiNx{suffix}"],
            uid=membrane_sinx_boxes[suffix].uid,
        )
        spec.add_physical_volume(
            name=f"Membrane_Si1{suffix}",
            primitive_names=[f"Membrane_Si1{suffix}"],
            tag=expected_tags[f"Membrane_Si1{suffix}"],
            uid=membrane_si1_boxes[suffix].uid,
        )

    face_bases = {name: 1000 + (tag - 100) * 100 for name, tag in expected_tags.items()}
    face_offsets = {"xmin": 0, "xmax": 1, "ymin": 2, "ymax": 3, "zmin": 4, "zmax": 5}
    for surface in spec.physical_surfaces:
        if surface.name == "bath" or "__" not in surface.name:
            continue
        body, face = surface.name.split("__", 1)
        if not body:
            continue
        base, suf = _strip_known_suffix(body, sides)
        if base.endswith("_sub"):
            parent = base.removesuffix("_sub") + suf
            surface.tag = face_bases[parent] + 10 + face_offsets[face]
        else:
            surface.tag = face_bases[body] + face_offsets[face]

    for side_index, suffix in enumerate(sides):
        stycast = stycast_boxes[suffix]
        membrane_si1 = membrane_si1_boxes[suffix]
        membrane_sinx = membrane_sinx_boxes[suffix]
        name_prefix = f"Membrane_SiNx{suffix}"

        _add_all_box_surfaces(
            spec,
            box_name=f"Stycast{suffix}",
            box_uid=stycast.uid,
            base_tag=face_bases[f"Stycast{suffix}"],
        )
        _add_all_box_surfaces(
            spec,
            box_name=f"Membrane_Si1{suffix}",
            box_uid=membrane_si1.uid,
            base_tag=face_bases[f"Membrane_Si1{suffix}"],
        )

        membrane_sinx_base = face_bases[f"Membrane_SiNx{suffix}"]
        for face, box_name in (
            ("xmin", f"{name_prefix}_free_left"),
            ("xmax", f"{name_prefix}_free_right"),
            ("ymin", f"{name_prefix}_free_bottom"),
            ("ymax", f"{name_prefix}_free_top"),
        ):
            spec.add_physical_surface(
                name=f"{name_prefix}__{face}",
                selector="box_surface",
                tag=membrane_sinx_base + face_offsets[face],
                uid=f"{membrane_sinx.uid}__{face}",
                box_name=box_name,
                surface=face,
            )

        membrane_surface_parts = (
            (f"{name_prefix}_contact", f"{membrane_sinx.uid}_contact"),
            (f"{name_prefix}_free_left", f"{membrane_sinx.uid}_free_left"),
            (f"{name_prefix}_free_right", f"{membrane_sinx.uid}_free_right"),
            (f"{name_prefix}_free_bottom", f"{membrane_sinx.uid}_free_bottom"),
            (f"{name_prefix}_free_top", f"{membrane_sinx.uid}_free_top"),
        )
        # Transient tag range for these zone parts (removed by name before the
        # final write, in _retag_mortar_interface_surfaces): 2300 exactly
        # reproduces the single_pixel literal; for dual_tes it must dodge the
        # face_bases range, which spans up to ~2800 for the 19-body scheme.
        temp_base = 2300 if len(sides) == 1 else 90000 + 100 * side_index
        for offset, (box_name, uid) in enumerate(membrane_surface_parts):
            for face_offset, face in ((4, "zmin"), (5, "zmax")):
                spec.add_physical_surface(
                    name=f"{box_name}__{face}",
                    selector="box_surface",
                    tag=temp_base + offset * 10 + face_offset,
                    uid=f"{uid}__{face}",
                    box_name=box_name,
                    surface=face,
                )

    fragment_mortar_interfaces = bool(
        elmer_overrides.get("fragment_mortar_interfaces", False)
    )

    conformal_shared_interfaces = bool(
        elmer_overrides.get("conformal_shared_interfaces", False)
    )
    conformal_mortar_interfaces = bool(
        elmer_overrides.get("conformal_mortar_interfaces", False)
    )
    builder = TesGmshBuilder(spec=spec, verbose=False)
    if fragment_mortar_interfaces and (len(sides) > 1 or conformal_mortar_interfaces):
        # Multi-side (dual_tes) contact refinement, all BEFORE meshing, with
        # each side's real Stycast position/radius. The single_pixel path
        # keeps the legacy post-builder calls below untouched for
        # bit-identical output.
        # - TES films: built pre-split along the disc footprint (OCC cannot
        #   boolean the 0.16 um film without corrupting it; see
        #   _add_contact_film_primitive).
        # - absorber: its bottom face gets a 2D disc imprint per side
        #   (_imprint_contact_discs; the 700 um body itself is unproblematic).
        builder.contact_film_discs = {
            f"TES{suffix}": (
                stycast_boxes[suffix].x,
                stycast_boxes[suffix].y,
                float(stycast_boxes[suffix].dx) / 2.0,
            )
            for suffix in sides
        }
        abs_box = _find_box(spec, "abs")
        builder.contact_disc_specs = [
            {
                "body": "abs",
                "discs": [
                    (
                        stycast_boxes[suffix].x,
                        stycast_boxes[suffix].y,
                        abs_box.zmin,
                        float(stycast_boxes[suffix].dx) / 2.0,
                    )
                    for suffix in sides
                ],
            }
        ]
        if conformal_mortar_interfaces:
            for suffix in sides:
                tes_box = _find_box(spec, f"TES{suffix}")
                stycast = stycast_boxes[suffix]
                builder.contact_disc_specs.append(
                    {
                        "body": f"Membrane_SiNx{suffix}",
                        "discs": [
                            (
                                stycast.x,
                                stycast.y,
                                tes_box.zmin,
                                float(stycast.dx) / 2.0,
                            )
                        ],
                    }
                )
            builder.merge_mortar_interfaces = True
            builder.conformal_contact_pairs = []
            for suffix in sides:
                tes_box = _find_box(spec, f"TES{suffix}")
                stycast = stycast_boxes[suffix]
                disc_area = pi * (float(stycast.dx) / 2.0) ** 2
                builder.conformal_contact_pairs.extend(
                    [
                        (
                            f"TES{suffix}", tes_box.zmin,
                            f"Membrane_SiNx{suffix}", tes_box.zmin,
                            None,
                            (tes_box.xmin, tes_box.xmax, tes_box.ymin, tes_box.ymax),
                            f"TES{suffix}/membrane",
                        ),
                        (
                            f"Stycast{suffix}", stycast.zmin,
                            f"TES{suffix}", tes_box.zmax,
                            disc_area, None, f"Stycast{suffix}/TES",
                        ),
                        (
                            f"Stycast{suffix}", stycast.zmax,
                            "abs", _find_box(spec, "abs").zmin,
                            disc_area, None, f"Stycast{suffix}/absorber",
                        ),
                    ]
                )

    try:
        gmsh.option.setNumber = _set_number_without_optimize
        builder.build(conformal_shared_interfaces=conformal_shared_interfaces)
        if conformal_mortar_interfaces:
            for body_name, physical_tag in expected_tags.items():
                body_tags = sorted(set(builder.group_entities.get(body_name, [])))
                if not body_tags:
                    body_boxes = [
                        box for box in spec.boxes
                        if body_name_of(box) == body_name
                    ]
                    body_tags = _find_volume_tags_for_boxes(body_boxes)
                if not body_tags:
                    raise RuntimeError(
                        f"Periodic contact meshing lost material volume {body_name!r}"
                    )
                _replace_volume_group(body_name, physical_tag, body_tags)
        if conformal_shared_interfaces:
            # OCC fragmentation changes entity tags. Re-register every material
            # group from its semantic source boxes before writing the mesh.
            for body_name, physical_tag in expected_tags.items():
                body_boxes = [
                    box for box in spec.boxes
                    if body_name_of(box) == body_name
                ]
                body_tags = _find_volume_tags_for_boxes(body_boxes)
                if not body_tags:
                    raise RuntimeError(
                        f"Conformal fragmentation lost material volume {body_name!r}"
                    )
                _replace_volume_group(body_name, physical_tag, body_tags)
        if not conformal_shared_interfaces and not conformal_mortar_interfaces:
            if len(sides) == 1:
                _fragment_mortar_interfaces(fragment_mortar_interfaces)
            for suffix in sides:
                _retag_mortar_interface_surfaces(
                    fragment_mortar_interfaces,
                    membrane_sinx_parts_by_side[suffix],
                    name_prefix=f"Membrane_SiNx{suffix}",
                    volume_tag=expected_tags[f"Membrane_SiNx{suffix}"],
                    face_base=face_bases[f"Membrane_SiNx{suffix}"],
                    legacy_contact=len(sides) == 1,
                )
            if len(sides) > 1:
                _retag_contact_surfaces_multi(
                    fragment_mortar_interfaces,
                    sides=sides,
                    spec=spec,
                    stycast_boxes=stycast_boxes,
                    face_bases=face_bases,
                )
            (ROOT / "generated" / "mortar_surface_groups.json").write_text(
                json.dumps(_identify_contact_surface_groups(), indent=2),
                encoding="utf-8",
            )
            _write_mortar_surface_report(ROOT / "generated" / "mortar_surface_report.json")
        gmsh.write(str(ROOT / "gmsh" / "project_shifted.brep"))
        if write_mesh:
            # ElmerGrid (Elmer 26.1) misreads multi-block msh 4.1 $Nodes sections and
            # scrambles node coordinates; msh 2.2 converts correctly.
            # See docs/mesh_conversion_defect_20260712.md.
            _original_set_number("Mesh.MshFileVersion", 2.2)
            builder.write(ROOT / "gmsh" / "project.msh")
    finally:
        gmsh.option.setNumber = _original_set_number
        builder.finalize()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        project_path = Path(sys.argv[1])
        if not project_path.is_absolute():
            project_path = (ROOT / project_path).resolve()
        geometry_path: Path | None = None
        if len(sys.argv) > 2:
            geometry_path = Path(sys.argv[2])
            if not geometry_path.is_absolute():
                geometry_path = (ROOT / geometry_path).resolve()
        set_input_jsons(project_path, geometry_path)
    build()

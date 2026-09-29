"""Parser for Elmer `mesh.names` files.

ElmerGrid's `-names` output writes two MATC-comment sections, each holding
`$ <name> = <id>` assignment lines: one mapping body names to target body
IDs, one mapping boundary names to target boundary IDs. build_cases.py uses
this to derive bodies and boundary conditions from a case's mesh instead of
hardcoding IDs (see docs/dual_tes_plan.md, Phase B).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import NamedTuple

_ASSIGNMENT_RE = re.compile(r"^\$\s*([A-Za-z_]\w*)\s*=\s*(-?\d+)\s*$")
_BODIES_HEADER = "names for bodies"
_BOUNDARIES_HEADER = "names for boundaries"


class MeshNames(NamedTuple):
    """Name -> target ID maps parsed from a mesh.names file."""

    bodies: dict[str, int]
    boundaries: dict[str, int]


def parse_mesh_names(path: Path) -> MeshNames:
    """Parse a `mesh.names` file into body and boundary name->ID maps.

    Section membership is tracked by the `! ----- names for bodies -----` /
    `! ----- names for boundaries -----` banner comments; everything before
    the first banner is ignored.
    """
    bodies: dict[str, int] = {}
    boundaries: dict[str, int] = {}
    target: dict[str, int] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("!"):
            if _BODIES_HEADER in stripped:
                target = bodies
            elif _BOUNDARIES_HEADER in stripped:
                target = boundaries
            continue
        if target is None:
            continue
        match = _ASSIGNMENT_RE.match(stripped)
        if match:
            target[match.group(1)] = int(match.group(2))
    return MeshNames(bodies, boundaries)

from __future__ import annotations

import math
from pathlib import Path

from scripts.analysis.evaluate_physical_parity import oriented_face_flux, result_values, tetra_gradient


LEFT = [(0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0), (0.0, 0.0, 0.0)]
RIGHT = [(0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0), (0.0, 0.0, 2.0)]
FACE = (0, 1, 2)


def test_two_material_constant_flux_has_opposite_outward_flux() -> None:
    left_temperatures = [1.0, 1.0, 1.0, 2.0]
    right_temperatures = [1.0, 1.0, 1.0, 0.5]
    left = oriented_face_flux(LEFT, left_temperatures, FACE, 2.0)
    right = oriented_face_flux(RIGHT, right_temperatures, FACE, 4.0)
    assert left is not None and right is not None
    assert math.isclose(left["integrated_flux_W"] + right["integrated_flux_W"], 0.0, abs_tol=1.0e-12)
    assert math.isclose(abs(left["integrated_flux_W"]), abs(right["integrated_flux_W"]), rel_tol=1.0e-12)


def test_face_node_order_does_not_change_outward_flux() -> None:
    temperatures = [1.0, 1.0, 1.0, 2.0]
    forward = oriented_face_flux(LEFT, temperatures, (0, 1, 2), 2.0)
    reversed_order = oriented_face_flux(LEFT, temperatures, (2, 1, 0), 2.0)
    assert forward is not None and reversed_order is not None
    assert math.isclose(forward["integrated_flux_W"], reversed_order["integrated_flux_W"], rel_tol=1.0e-12)


def test_zero_temperature_gradient_has_zero_flux() -> None:
    assert tetra_gradient(LEFT, [1.0, 1.0, 1.0, 1.0]) == (0.0, 0.0, 0.0)
    result = oriented_face_flux(LEFT, [1.0, 1.0, 1.0, 1.0], FACE, 123.0)
    assert result is not None
    assert result["integrated_flux_W"] == 0.0


def test_result_parser_uses_mesh_node_column_and_last_field(tmp_path: Path) -> None:
    result = tmp_path / "field.result"
    result.write_text(
        "ASCII 3\n"
        "Time: 1 1 1.0\n"
        "temperature\n"
        "Perm: 3 3\n"
        "1 3\n2 1\n3 2\n"
        "30\n10\n20\n"
        "Time: 2 2 2.0\n"
        "temperature\n"
        "Perm: 3 3\n"
        "1 3\n2 1\n3 2\n"
        "300\n100\n200\n",
        encoding="utf-8",
    )
    assert result_values(result, field_index=0) == {1: 30.0, 2: 10.0, 3: 20.0}
    assert result_values(result) == {1: 300.0, 2: 100.0, 3: 200.0}

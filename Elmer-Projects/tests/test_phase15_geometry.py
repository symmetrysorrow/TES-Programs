import math

from generate_hybrid_prism_geometry import (
    DEFAULT_OUT,
    GLOBAL_STACK_SIZE,
    absorber_local_field_profile,
    configure_local_fields,
    configure_mesh_size_options,
    parse_args,
    stack_local_field_profile,
)


def test_default_cli_preserves_legacy_uniform_mesh_recipe() -> None:
    args = parse_args([])
    assert args.output == DEFAULT_OUT
    assert args.stack_local_size is None
    assert args.absorber_local_size is None
    assert args.absorber_local_radius == 150e-6
    assert args.disable_mesh_size_extend_from_boundary is False
    assert stack_local_field_profile(None, args.stack_local_half_width, 0.0, 1e-3, 0.0, 1e-3) is None


def test_stack25_profile_is_central_and_uses_50um_transition() -> None:
    args = parse_args(["--output", "gmsh/project_hybrid_prism_stack25.msh", "--stack-local-size", "25e-6"])
    profile = stack_local_field_profile(args.stack_local_size, args.stack_local_half_width, 0.0, 1e-3, 2e-6, 20e-6)
    expected = {
        "VIn": 25e-6,
        "VOut": GLOBAL_STACK_SIZE,
        "XMin": -0.4e-3,
        "XMax": 0.4e-3,
        "YMin": 0.6e-3,
        "YMax": 1.4e-3,
        "ZMin": 2e-6 - GLOBAL_STACK_SIZE,
        "ZMax": 20e-6,
        "Thickness": GLOBAL_STACK_SIZE,
    }
    assert profile is not None
    assert profile.keys() == expected.keys()
    for key, value in expected.items():
        assert math.isclose(profile[key], value, rel_tol=0.0, abs_tol=1e-18)


def test_absorber_ball_profile_is_centred_and_uses_50um_transition() -> None:
    args = parse_args(["--absorber-local-size", "16.6666666666667e-6"])
    profile = absorber_local_field_profile(args.absorber_local_size, args.absorber_local_radius, 0.0, 1e-3, 20e-6, 525e-6)
    expected = {
        "VIn": 16.6666666666667e-6,
        "VOut": GLOBAL_STACK_SIZE,
        "XCenter": 0.0,
        "YCenter": 1e-3,
        "ZCenter": 282.5e-6,
        "Radius": 150e-6,
        "Thickness": GLOBAL_STACK_SIZE,
    }
    assert profile is not None
    assert profile.keys() == expected.keys()
    for key, value in expected.items():
        assert math.isclose(profile[key], value, rel_tol=0.0, abs_tol=1e-18)


def test_stack_and_absorber_fields_are_combined_with_min(monkeypatch) -> None:
    calls: list[tuple] = []

    class Field:
        def add(self, field_type: str) -> int:
            calls.append(("add", field_type))
            return len([call for call in calls if call[0] == "add"])

        def setNumber(self, field: int, name: str, value: float) -> None:
            calls.append(("setNumber", field, name, value))

        def setNumbers(self, field: int, name: str, values: list[int]) -> None:
            calls.append(("setNumbers", field, name, values))

        def setAsBackgroundMesh(self, field: int) -> None:
            calls.append(("background", field))

    class Mesh:
        field = Field()

    class Model:
        mesh = Mesh()

    class Gmsh:
        model = Model()

    monkeypatch.setattr("generate_hybrid_prism_geometry.gmsh", Gmsh())
    configure_local_fields({"VIn": 20e-6}, {"VIn": 10e-6})
    assert calls == [
        ("add", "Box"),
        ("setNumber", 1, "VIn", 20e-6),
        ("add", "Ball"),
        ("setNumber", 2, "VIn", 10e-6),
        ("add", "Min"),
        ("setNumbers", 3, "FieldsList", [1, 2]),
        ("background", 3),
    ]


def test_mesh_size_extension_override_is_opt_in(monkeypatch) -> None:
    calls: list[tuple[str, str, int]] = []

    class Option:
        def setNumber(self, name: str, value: int) -> None:
            calls.append(("setNumber", name, value))

    class Gmsh:
        option = Option()

    monkeypatch.setattr("generate_hybrid_prism_geometry.gmsh", Gmsh())
    configure_mesh_size_options(False)
    assert calls == []
    configure_mesh_size_options(True)
    assert calls == [("setNumber", "Mesh.MeshSizeExtendFromBoundary", 0)]

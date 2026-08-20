import math

import pytest

from generate_hybrid_prism_geometry import (
    DEFAULT_OUT,
    GLOBAL_STACK_SIZE,
    assess_mesh_quality,
    absorber_local_field_profile,
    configure_local_fields,
    configure_mesh_size_options,
    inspect_mesh_file,
    parse_args,
    stack_local_field_profile,
)


def test_mesh_quality_contract_accepts_balanced_interfaces_and_good_prisms() -> None:
    assert assess_mesh_quality(2116, 3822, {"Membrane_SiNx": [0.01797], "Membrane_Si1": [0.5]}) == []


def test_mesh_quality_contract_rejects_coarse_mortar_side_and_bad_prisms() -> None:
    reasons = assess_mesh_quality(3710, 174, {"Membrane_SiNx": [0.007744], "Membrane_Si1": [0.058]})
    assert any("element-count ratio" in reason for reason in reasons)
    assert any("Membrane_SiNx" in reason and "minSICN" in reason for reason in reasons)


def test_mesh_quality_contract_ignores_thin_nonmembrane_prism_quality() -> None:
    # SiO2_1 layer splitting can make the global prism minSICN low without
    # degrading either membrane volume or the TES/membrane Mortar interface.
    assert assess_mesh_quality(4274, 7476, {"Membrane_SiNx": [0.015], "Membrane_Si1": [0.058]}) == []


def test_default_cli_preserves_legacy_uniform_mesh_recipe() -> None:
    args = parse_args([])
    assert args.output == DEFAULT_OUT
    assert args.stack_local_size is None
    assert args.absorber_local_size is None
    assert args.absorber_local_radius == 150e-6
    assert args.disable_mesh_size_extend_from_boundary is False
    assert args.mesh_algorithm == 6
    assert args.stycast_layers == 1
    assert args.tes_layers == 1
    assert args.sio2_1_layers == 1
    assert stack_local_field_profile(None, args.stack_local_half_width, 0.0, 1e-3, 0.0, 1e-3) is None


def test_mesh_algorithm_override() -> None:
    assert parse_args(["--mesh-algorithm", "5"]).mesh_algorithm == 5


def test_global_mesh_size_defaults_to_legacy_50um() -> None:
    assert parse_args([]).global_mesh_size == GLOBAL_STACK_SIZE


def test_global_mesh_size_override_raises_the_local_size_ceiling() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--absorber-local-size", "75e-6"])
    args = parse_args(["--global-mesh-size", "100e-6", "--absorber-local-size", "75e-6"])
    assert args.global_mesh_size == 100e-6
    assert args.absorber_local_size == 75e-6


def test_absorber_profile_uses_overridden_global_size() -> None:
    profile = absorber_local_field_profile(
        75e-6, 150e-6, 0.0, 1e-3, 20e-6, 525e-6, global_size=100e-6
    )
    assert profile is not None
    assert profile["VIn"] == 75e-6
    assert profile["VOut"] == 100e-6
    assert profile["Thickness"] == 100e-6
    # At the legacy 50 um background, the same 75 um local size would be a no-op.
    assert absorber_local_field_profile(75e-6, 150e-6, 0.0, 1e-3, 20e-6, 525e-6) is None


def test_inspect_mesh_file_always_finalizes(monkeypatch, tmp_path) -> None:
    calls: list[tuple[str, str] | tuple[str]] = []

    class Gmsh:
        def initialize(self) -> None:
            calls.append(("initialize",))

        def open(self, path: str) -> None:
            calls.append(("open", path))

        def finalize(self) -> None:
            calls.append(("finalize",))

    monkeypatch.setattr("generate_hybrid_prism_geometry.gmsh", Gmsh())
    monkeypatch.setattr("generate_hybrid_prism_geometry.measure_mesh_quality", lambda: (_ for _ in ()).throw(RuntimeError("bad mesh")))
    try:
        inspect_mesh_file(tmp_path / "bad.msh")
    except RuntimeError as error:
        assert str(error) == "bad mesh"
    else:
        raise AssertionError("inspection error was not propagated")
    assert calls[-1] == ("finalize",)


def test_stycast_layer_override() -> None:
    args = parse_args(["--stycast-layers", "8"])
    assert args.stycast_layers == 8


def test_independent_stack_layer_overrides() -> None:
    args = parse_args(["--tes-layers", "2", "--sinx-layers", "4", "--si-2-layers", "3"])
    assert args.tes_layers == 2
    assert args.sinx_layers == 4
    assert args.si_2_layers == 3


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

from scripts.support import run_phase24_gate4_5_nomortar as runner


def test_refinement_spec_is_full_constrained_steady_projection(monkeypatch) -> None:
    monkeypatch.setattr(runner, "MESH", "mesh_test")
    monkeypatch.setattr(runner, "APPLY_MORTAR_BCS", True)
    monkeypatch.setattr(runner, "NONLINEAR_TOLERANCE", 1.0e-8)
    monkeypatch.setattr(runner, "NONLINEAR_MAX_ITERATIONS", 120)

    spec = runner.refinement_spec(
        {
            "template": "pulse",
            "pulse": {"energy": "1332[keV]"},
            "timesteps": [["18[us]", 1]],
            "solver": {"linear_system": "iterative_hypre_flexgmres_boomeramg"},
        },
        name="refined",
        source_case="gate3",
    )

    assert spec["template"] == "steady"
    assert spec["apply_mortar_bcs"] is True
    assert spec["restart_from"] is None
    assert spec["preexisting_restart"] is True
    assert spec["restart_file_base"] == "gate3"
    assert spec["restart_file_path"].endswith("/gate3.result")
    assert spec["output_result"] is True
    assert spec["solver"]["linear_system"] == "mumps"
    assert spec["solver"]["steady_state_convergence_tolerance"] == 1.0e-9
    assert "pulse" not in spec
    assert "timesteps" not in spec

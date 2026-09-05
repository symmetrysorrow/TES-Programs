import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "PoST_Simulations" / "cases" / "tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2"
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))
import target_case_audit  # noqa: E402


def test_target_case_does_not_inherit_generic_physics_values():
    target = json.loads((CASE / "input.json").read_text())
    generic = json.loads((ROOT / "PoST_Simulations" / "input.json").read_text())
    unresolved = [
        "T_c", "R", "R_l", "alpha", "beta", "L", "n", "C_abs", "C_tes",
        "G_abs-abs", "G_abs-tes", "G_tes-bath", "n_abs", "E",
    ]
    assert all(target[name] is None for name in unresolved)
    assert any(target[name] != generic[name] for name in unresolved if target[name] is None)


def test_readiness_is_split_by_capability_and_physical_sources_stay_enabled():
    result = target_case_audit.audit(CASE)
    capabilities = result["capabilities"]
    assert capabilities["operating_point_ready"]["ready"] is False
    assert capabilities["python_stability_ready"]["ready"] is False
    assert capabilities["reduced_noise_ready"]["ready"] is False
    assert capabilities["cpp_parity_ready"]["ready"] is False
    assert capabilities["normalized_comparison_ready"]["ready"] is False
    assert capabilities["absolute_comparison_ready"]["ready"] is False
    assert "T_c" in capabilities["operating_point_ready"]["missing_parameters"]
    assert "E" not in capabilities["reduced_noise_ready"]["required_parameters"]
    assert result["noise_semantics"]["tes_johnson"] == "enabled physical source"
    assert result["noise_semantics"]["load_johnson"] == "enabled physical source"
    assert result["noise_semantics"]["tes_bath_tfn"] == "enabled physical source"
    assert result["noise_semantics"]["tes_absorber_tfn"] == "enabled physical source"
    assert result["noise_semantics"]["post_filter_white_asd_A_rtHz"] == 0.0
    assert result["noise_semantics"]["readout_white_asd_A_rtHz"] == 0.0
    assert result["noise_semantics"]["tes_resistance_fluctuation_model"] == "none"
    assert result["noise_semantics"]["tes_internal_model"] == "none"
    assert "noise_sources" not in json.loads((CASE / "input.json").read_text())


def test_report_names_post_analysis_and_does_not_back_calculate_pre_analysis():
    report = json.loads((CASE / "comparison_summary.json").read_text())
    assert report["comparison_kind"] == "target_case_intrinsic_physical_noise_only"
    assert report["spectrum_semantics"]["experimental_pre_analysis_asd"] is None
    assert "experimental_post_analysis_asd" in report["spectrum_semantics"]
    for row in report["points"]:
        assert "experimental_post_analysis_asd" in row
        assert "experimental_post_analysis_normalized" in row
        assert "simulation_post_analysis_asd" in row
        assert row["simulation_post_analysis_asd"] is None
        assert row["simulation_over_experiment"] is None


def test_identifiability_and_dependency_maps_are_actionable():
    identifiability = json.loads((CASE / "parameter_identifiability.json").read_text())
    dependency = json.loads((CASE / "parameter_dependency.json").read_text())
    for name in ("alpha", "beta", "C_tes", "G_tes-bath", "readout_calibration_A_per_V"):
        entry = identifiability["parameters"][name]
        assert entry["blocks"]
        assert entry["possible_existing_sources"]
        assert entry["observable_needed"]
        assert entry["do_not_use"]
    assert "operating_point" in dependency["dependency_graph"]
    assert "absolute_comparison_ready" in dependency["dependency_graph"]["absolute_asd"]["unlocks"]

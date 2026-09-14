from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))

import Opt_noise as optimizer  # noqa: E402


def test_absolute_level_residual_is_log10_ratio() -> None:
    residual = optimizer.absolute_level_residual(2.0e-12, 1.0e-12)
    assert residual == pytest.approx(np.log10(2.0))


def test_absorber_heat_capacity_prior_is_factor_two_about_material_value() -> None:
    nominal = optimizer.C_ABS_MATERIAL_J_PER_K
    assert optimizer.C_ABS_FIT_MIN_J_PER_K == pytest.approx(0.5 * nominal)
    assert optimizer.C_ABS_FIT_MAX_J_PER_K == pytest.approx(2.0 * nominal)


def test_boundary_diagnostics_uses_log_position_for_log_bound() -> None:
    bound = optimizer.Bound(1.0e-10, 1.0e-8, logarithmic=True)
    candidate = {"L": 1.0e-9}
    result = optimizer.parameter_boundary_diagnostics(
        candidate,
        ["L"],
        {"L": bound},
    )["L"]
    assert result["position_fraction"] == pytest.approx(0.5)
    assert result["nearest_bound"] == "lower"
    assert result["within_1pct_of_bound"] is False


def test_absolute_asd_anchor_is_disabled_by_default() -> None:
    assert optimizer.ABSOLUTE_ASD_WEIGHT_DEFAULT == 0.0



def _stable_stycast_candidate() -> dict:
    return {
        "C_abs": 4.5e-9,
        "C_tes": 8.0e-13,
        "C_stycast": 1.1e-11,
        "G_abs-abs": 5.9e-7,
        "G_tes-stycast": 3.0e-8,
        "G_stycast-abs": 3.0e-8,
        "G_tes-bath": 4.0e-8,
        "R": 0.0175,
        "R_l": 0.004,
        "T_c": 0.25,
        "T_bath": 0.215,
        "alpha": 0.5,
        "beta": 1.0,
        "L": 5.0e-9,
        "n": 3.0,
        "excess_johnson_M": 0.0,
        "post_filter_white_fraction": 1.0e-6,
        "rate": 500_000.0,
        "thermal_link_model": "stycast_node",
    }


def _stable_stycast_rc_candidate() -> dict:
    candidate = _stable_stycast_candidate()
    candidate.update(
        {
            "electrical_link_model": optimizer.ELECTRICAL_LINK_MODEL_RC,
            "R_rc": 1.0e-3,
            "f_rc_Hz": 100_000.0,
        }
    )
    return candidate

def _stable_stycast_hanging_candidate() -> dict:
    candidate = _stable_stycast_candidate()
    candidate.update(
        {
            "thermal_extension": optimizer.THERMAL_EXTENSION_TES_HANGING,
            "C_hanging": candidate["C_tes"],
            "G_tes-hanging": (
                2.0 * np.pi * 10_000.0 * candidate["C_tes"]
            ),
            "electrical_link_model": "rl",
        }
    )
    return candidate


def _stable_stycast_series_candidate() -> dict:
    candidate = _stable_stycast_candidate()
    base_g = candidate["G_tes-stycast"]
    candidate.update(
        {
            "thermal_extension": optimizer.THERMAL_EXTENSION_TES_STYCAST_SERIES,
            "C_series": candidate["C_tes"],
            "G_tes-series": 2.0 * base_g,
            "G_series-stycast": 2.0 * base_g,
            "electrical_link_model": "rl",
        }
    )
    return candidate



def test_source_class_diagnostics_reconstructs_total_model_psd() -> None:
    candidate = _stable_stycast_candidate()
    result = optimizer.post_analysis_source_class_asd(
        candidate,
        np.array([1_000.0, 10_000.0, 100_000.0]),
    )
    assert set(result["class_asd_A_rtHz"]) >= {
        "TES_Johnson",
        "load_Johnson",
        "TES_bath_TFN",
        "TES_Stycast_TFN",
        "Stycast_absorber_TFN",
        "post_filter_white",
    }
    assert result["reconstruction_max_relative_error"] < 1.0e-10
    assert np.all(result["total_asd_A_rtHz"] > 0.0)


def test_source_class_transfer_shape_reports_normalized_shapes_and_delta_p() -> None:
    frequency = np.asarray([1_000.0, 10_000.0, 70_000.0, 150_000.0])
    class_asd = {
        "source_A": np.asarray([2.0, 3.0, 1.0, 1.0]),
        "source_B": np.asarray([1.0, 1.0, 1.0, 1.0]),
    }
    total_asd = np.sqrt(sum(asd**2 for asd in class_asd.values()))
    curves = {
        "frequencies_Hz": frequency,
        "class_asd_A_rtHz": class_asd,
        "total_asd_A_rtHz": total_asd,
    }
    target = np.asarray([1.0, 1.6, 0.5, 0.5])
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
    )

    result = optimizer.source_class_transfer_shape_diagnostics(
        curves,
        target,
        frequency,
        args,
    )

    mid = result["fit_bands"]["5000-15000_Hz"]["sources"]["source_A"]
    high = result["fit_bands"]["40000-100000_Hz"]["sources"]["source_A"]
    assert mid["mean_psd_fraction_p_j"] == pytest.approx(0.9)
    assert mid["mean_asd_scale_sensitivity_delta_p_j"] == pytest.approx(0.1)
    assert high["mean_psd_fraction_p_j"] == pytest.approx(0.5)
    assert high["mean_asd_scale_sensitivity_delta_p_j"] == pytest.approx(-0.3)
    assert result["current_residual_direction"]["5_15kHz_is_deficit"] is True
    assert result["current_residual_direction"]["40_100kHz_is_excess"] is True
    assert result[
        "can_move_5_15kHz_deficit_and_40_100kHz_excess_simultaneously"
    ] is True
    assert {
        "source_class": "source_A",
        "asd_scale_direction": "increase_source_ASD",
        "midband_delta_p_j": pytest.approx(0.1),
        "highband_delta_p_j": pytest.approx(-0.3),
    } in result["source_candidates_for_simultaneous_correct_direction"]


def test_source_class_transfer_shape_is_zero_for_self_target() -> None:
    candidate = _stable_stycast_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 81)
    curves = optimizer.post_analysis_source_class_asd(candidate, fit_freq)
    target = optimizer.normalize_at(
        fit_freq,
        curves["total_asd_A_rtHz"],
        reference_hz=optimizer.ABSOLUTE_ASD_REFERENCE_HZ,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
    )
    result = optimizer.source_class_transfer_shape_diagnostics(
        curves,
        target,
        fit_freq,
        args,
    )
    assert result["production_optimizer_unchanged"] is True
    assert result["best_fit_parameters_held_fixed"] is True
    assert result[
        "can_move_5_15kHz_deficit_and_40_100kHz_excess_simultaneously"
    ] is False
    for band in result["fit_bands"].values():
        assert band["current_mean_log10_model_over_measurement"] == pytest.approx(
            0.0,
            abs=1.0e-12,
        )


def test_eigenmode_diagnostics_reports_stable_seven_state_modes() -> None:
    result = optimizer.eigenmode_diagnostics(_stable_stycast_candidate())
    assert result["state_count"] == 7
    assert result["state_order"] == [
        "I1",
        "TES1",
        "Stycast1",
        "Pb_center",
        "Stycast2",
        "TES2",
        "I2",
    ]
    assert result["stable"] is True
    modes = result["modes_sorted_by_natural_frequency"]
    assert len(modes) == 7
    assert all(mode["real_s_inv"] < 0.0 for mode in modes)
    assert all(mode["natural_frequency_Hz"] > 0.0 for mode in modes)
    for mode in modes:
        assert sum(mode["state_participation"].values()) == pytest.approx(1.0)
        assert mode["dominant_state"] in result["state_order"]
        assert 0.0 <= mode["dominant_state_participation"] <= 1.0


def test_source_ablation_is_diagnostic_and_reports_band_residuals() -> None:
    candidate = _stable_stycast_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 121)
    curves = optimizer.post_analysis_source_class_asd(candidate, fit_freq)
    target = optimizer.normalize_at(
        fit_freq,
        curves["total_asd_A_rtHz"],
        reference_hz=optimizer.ABSOLUTE_ASD_REFERENCE_HZ,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    result = optimizer.source_ablation_diagnostics(
        curves,
        target,
        fit_freq,
        args,
    )
    assert result["diagnostic_only"] is True
    assert result["full_model_shape_score"] == pytest.approx(0.0, abs=1.0e-12)
    assert "TES_Johnson" in result["remove_one_source_class"]
    johnson = result["remove_one_source_class"]["TES_Johnson"]
    assert johnson["shape_score"] > 0.0
    assert "40000-100000_Hz" in johnson["bands"]
    assert "100000-200000_Hz" in johnson["bands"]



def test_johnson_beta_audit_reports_shared_convention() -> None:
    candidate = _stable_stycast_candidate()
    candidate["beta"] = 2.0
    result = optimizer.johnson_beta_audit(candidate)
    assert result["status"] == "shared_production_convention"
    assert result["source_psd_beta_factor"] == pytest.approx(5.0)
    assert result["source_asd_factor_vs_beta0_same_M"] == pytest.approx(
        np.sqrt(5.0)
    )
    assert result["source_voltage_asd_V_rtHz"] > 0.0
    assert result["electrical_corner_Hz"] > 0.0


def test_beta_sweep_keeps_other_parameters_fixed_and_includes_best_beta() -> None:
    candidate = _stable_stycast_candidate()
    candidate["beta"] = 1.7
    fit_freq = np.geomspace(1_000.0, 200_000.0, 61)
    model, _ = optimizer.deterministic_simulated_spectrum(
        candidate.copy(),
        fit_freq,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    result = optimizer.beta_sweep_diagnostics(
        candidate,
        model,
        fit_freq,
        args,
    )
    betas = [row["beta"] for row in result["rows"]]
    assert 1.7 in betas
    matching = next(row for row in result["rows"] if row["beta"] == 1.7)
    assert matching["stable"] is True
    assert matching["shape_score"] == pytest.approx(0.0, abs=1.0e-12)
    assert "40000-100000_Hz" in matching["bands"]
    assert "100000-200000_Hz" in matching["bands"]



def test_johnson_source_scale_diagnostic_recovers_full_model_at_scale_one() -> None:
    candidate = _stable_stycast_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 81)
    curves = optimizer.post_analysis_source_class_asd(candidate, fit_freq)
    target = optimizer.normalize_at(
        fit_freq,
        curves["total_asd_A_rtHz"],
        reference_hz=optimizer.ABSOLUTE_ASD_REFERENCE_HZ,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    result = optimizer.johnson_source_scale_diagnostics(
        curves,
        target,
        fit_freq,
        args,
    )
    scales = [row["TES_Johnson_asd_scale"] for row in result["rows"]]
    assert scales == pytest.approx(np.linspace(0.0, 1.0, 11))
    best = result["best_global_shape_score_row"]
    assert best["TES_Johnson_asd_scale"] == pytest.approx(1.0)
    assert best["shape_score"] == pytest.approx(0.0, abs=1.0e-12)
    assert result["best_scale_by_band_mean_ratio"][
        "100000-200000_Hz"
    ]["TES_Johnson_asd_scale"] == pytest.approx(1.0)


def test_johnson_source_scale_zero_matches_source_ablation_total() -> None:
    candidate = _stable_stycast_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 81)
    curves = optimizer.post_analysis_source_class_asd(candidate, fit_freq)
    target = optimizer.normalize_at(
        fit_freq,
        curves["total_asd_A_rtHz"],
        reference_hz=optimizer.ABSOLUTE_ASD_REFERENCE_HZ,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    scale_result = optimizer.johnson_source_scale_diagnostics(
        curves,
        target,
        fit_freq,
        args,
    )
    ablation_result = optimizer.source_ablation_diagnostics(
        curves,
        target,
        fit_freq,
        args,
    )
    zero = next(
        row for row in scale_result["rows"]
        if row["TES_Johnson_asd_scale"] == 0.0
    )
    assert zero["shape_score"] == pytest.approx(
        ablation_result["remove_one_source_class"]["TES_Johnson"][
            "shape_score"
        ]
    )



def test_required_transfer_diagnostic_is_measurement_over_model() -> None:
    candidate = _stable_stycast_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 81)
    model = np.ones_like(fit_freq)
    target = np.geomspace(1.0, 0.25, len(fit_freq))
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
    )
    summary, curves = optimizer.required_transfer_diagnostics(
        candidate,
        model,
        target,
        fit_freq,
        args,
    )
    np.testing.assert_allclose(
        curves["required_ASD_transfer"],
        target / model,
    )
    np.testing.assert_allclose(
        curves["required_correction_dB"],
        20.0 * np.log10(target / model),
    )
    alias = curves["first_alias_fold_psd_fraction"]
    assert np.all(alias >= 0.0)
    assert np.all(alias <= 1.0)
    assert summary["diagnostic_only"] is True
    assert "40000-100000_Hz" in summary["fit_bands"]
    assert "100000-200000_Hz" in summary["fit_bands"]


def test_required_transfer_unity_target_reports_zero_db_correction() -> None:
    candidate = _stable_stycast_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 41)
    model = 1.0 + 0.2 * np.log10(fit_freq / fit_freq[0])
    target = model.copy()
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
    )
    summary, curves = optimizer.required_transfer_diagnostics(
        candidate,
        model,
        target,
        fit_freq,
        args,
    )
    np.testing.assert_allclose(curves["required_ASD_transfer"], 1.0)
    np.testing.assert_allclose(curves["required_correction_dB"], 0.0)
    for band in summary["fit_bands"].values():
        assert band["geometric_mean_required_ASD_transfer"] == pytest.approx(1.0)
        assert band["mean_required_correction_dB"] == pytest.approx(0.0)



def test_rc_relaxation_adds_two_passive_electrical_states() -> None:
    candidate = _stable_stycast_rc_candidate()
    point = optimizer.tes_operating_point(candidate)
    assert point["valid"] is True
    assert point["stable"] is True
    assert point["electrical_link_model"] == optimizer.ELECTRICAL_LINK_MODEL_RC
    assert point["R_rc_ohm"] == pytest.approx(candidate["R_rc"])
    assert point["f_rc_Hz"] == pytest.approx(candidate["f_rc_Hz"])
    assert point["tau_rc_s"] == pytest.approx(
        1.0 / (2.0 * np.pi * candidate["f_rc_Hz"])
    )
    matrix = optimizer.tes_linearized_matrix(candidate, 0.0)
    assert matrix.shape == (9, 9)


def test_rc_relaxation_keeps_tes_resistance_instantaneous_state_count() -> None:
    result = optimizer.eigenmode_diagnostics(_stable_stycast_rc_candidate())
    assert result["state_count"] == 9
    assert result["state_order"] == [
        "I1",
        "Vrc1",
        "TES1",
        "Stycast1",
        "Pb_center",
        "Stycast2",
        "TES2",
        "Vrc2",
        "I2",
    ]
    assert all("R_TES" not in name for name in result["state_order"])
    assert result["stable"] is True


def test_rc_relaxation_noise_includes_its_johnson_source_and_reconstructs() -> None:
    candidate = _stable_stycast_rc_candidate()
    result = optimizer.post_analysis_source_class_asd(
        candidate,
        np.array([1_000.0, 40_000.0, 100_000.0, 200_000.0]),
    )
    assert "RC_branch_Johnson" in result["class_asd_A_rtHz"]
    assert np.all(result["class_asd_A_rtHz"]["RC_branch_Johnson"] > 0.0)
    assert result["reconstruction_max_relative_error"] < 1.0e-10


def test_rc_relaxation_bounds_are_small_two_parameter_extension() -> None:
    assert optimizer.R_RC_FIT_MIN_OHM > 0.0
    assert optimizer.R_RC_FIT_MAX_OHM == pytest.approx(20.0e-3)
    assert optimizer.F_RC_FIT_MIN_HZ == pytest.approx(5_000.0)
    assert optimizer.F_RC_FIT_MAX_HZ == pytest.approx(500_000.0)



def test_rc_ablation_reports_direct_score_gain() -> None:
    candidate = _stable_stycast_rc_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 81)
    target, _ = optimizer.deterministic_simulated_spectrum(
        candidate.copy(),
        fit_freq,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    result = optimizer.electrical_rc_ablation_diagnostics(
        candidate,
        target,
        fit_freq,
        args,
    )
    assert result["rc_shape_score"] == pytest.approx(0.0, abs=1.0e-12)
    assert result["same_parameters_without_rc_shape_score"] > 0.0
    assert result["score_improvement_from_rc_at_same_parameters"] > 0.0
    assert result["tes_resistance_response"] == "instantaneous alpha/beta in both models"



def test_rc_degeneracy_profile_preserves_exact_best_seed() -> None:
    candidate = _stable_stycast_rc_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 31)
    target, _ = optimizer.deterministic_simulated_spectrum(
        candidate.copy(),
        fit_freq,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    result = optimizer.electrical_rc_degeneracy_diagnostics(
        candidate,
        target,
        fit_freq,
        args,
        fixed_rl_grid_ohm=(candidate["R_l"], 8.0e-3),
        profile_maxfev=20,
    )
    matching = next(
        row for row in result["fixed_R_l_profile"]["rows"]
        if row["R_l_ohm"] == pytest.approx(candidate["R_l"])
    )
    assert matching["shape_score"] == pytest.approx(0.0, abs=1.0e-12)
    assert matching["R_rc_ohm"] == pytest.approx(candidate["R_rc"])
    assert matching["f_rc_Hz"] == pytest.approx(candidate["f_rc_Hz"])
    assert result["tes_resistance_response"] == (
        "instantaneous alpha/beta (unchanged)"
    )


def test_rc_degeneracy_constant_sum_sweep_contains_best_decomposition() -> None:
    candidate = _stable_stycast_rc_candidate()
    fit_freq = np.geomspace(1_000.0, 200_000.0, 31)
    target, _ = optimizer.deterministic_simulated_spectrum(
        candidate.copy(),
        fit_freq,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    result = optimizer.electrical_rc_degeneracy_diagnostics(
        candidate,
        target,
        fit_freq,
        args,
        fixed_rl_grid_ohm=(candidate["R_l"],),
        profile_maxfev=5,
    )
    rows = result["constant_dc_sum_partition_sweep"]["rows"]
    matching = next(
        row for row in rows
        if np.isclose(row["R_l_ohm"], candidate["R_l"])
    )
    assert matching["R_rc_ohm"] == pytest.approx(candidate["R_rc"])
    assert matching["R_l_plus_R_rc_ohm"] == pytest.approx(
        candidate["R_l"] + candidate["R_rc"]
    )
    assert matching["shape_score"] == pytest.approx(0.0, abs=1.0e-12)
    assert result["constant_dc_sum_partition_sweep"][
        "shape_score_span"
    ] >= 0.0



def test_static_effective_series_decode_maps_to_plain_rl_model() -> None:
    original = _stable_stycast_candidate()
    original["R_series_eff"] = original["R_l"]
    bounds = {
        "R_series_eff": optimizer.Bound(
            optimizer.R_SERIES_EFF_FIT_MIN_OHM,
            optimizer.R_SERIES_EFF_FIT_MAX_OHM,
        )
    }
    target_resistance = 20.0e-3
    vector = np.asarray([np.log10(target_resistance)])
    candidate = optimizer.decode(
        vector,
        original,
        ["R_series_eff"],
        bounds,
        sample_rate=500_000.0,
        post_filter_white_asd=0.0,
        electrical_link_model="rl",
    )
    assert candidate["electrical_link_model"] == "rl"
    assert candidate["R_series_eff"] == pytest.approx(target_resistance)
    assert candidate["R_l"] == pytest.approx(target_resistance)
    assert "R_rc" not in candidate
    assert "f_rc_Hz" not in candidate


def test_effective_series_resistance_diagnostic_does_not_claim_decomposition() -> None:
    candidate = _stable_stycast_candidate()
    candidate["R_series_eff"] = 20.5e-3
    candidate["R_l"] = candidate["R_series_eff"]
    result = optimizer.effective_series_resistance_diagnostics(candidate)
    assert result["electrical_fit_model"] == "static_effective_series_resistance"
    assert result["R_series_eff_ohm"] == pytest.approx(20.5e-3)
    assert result["extra_electrical_state_count"] == 0
    assert "not a decomposition" in result["interpretation_guardrail"]


def test_static_effective_series_range_covers_rc_profile_total() -> None:
    assert optimizer.R_SERIES_EFF_FIT_MIN_OHM == pytest.approx(2.5e-3)
    assert optimizer.R_SERIES_EFF_FIT_MAX_OHM == pytest.approx(30.0e-3)
    assert optimizer.R_SERIES_EFF_FIT_MAX_OHM > (
        optimizer.R_L_FIT_MAX_OHM + 8.8e-3
    )



def test_hanging_thermal_body_adds_two_symmetric_states() -> None:
    candidate = _stable_stycast_hanging_candidate()
    point = optimizer.tes_operating_point(candidate)
    assert point["valid"] is True
    assert point["stable"] is True
    assert point["thermal_extension"] == optimizer.THERMAL_EXTENSION_TES_HANGING
    assert point["f_hanging_Hz"] == pytest.approx(10_000.0)
    matrix = optimizer.tes_linearized_matrix(candidate, 0.0)
    assert matrix.shape == (9, 9)

    modes = optimizer.eigenmode_diagnostics(candidate)
    assert modes["state_order"] == [
        "I1",
        "TES1",
        "Hanging1",
        "Stycast1",
        "Pb_center",
        "Stycast2",
        "Hanging2",
        "TES2",
        "I2",
    ]
    assert modes["stable"] is True


def test_hanging_thermal_body_noise_reconstructs_with_tfn_class() -> None:
    candidate = _stable_stycast_hanging_candidate()
    result = optimizer.post_analysis_source_class_asd(
        candidate,
        np.array([1_000.0, 10_000.0, 40_000.0, 100_000.0]),
    )
    assert "TES_hanging_TFN" in result["class_asd_A_rtHz"]
    assert np.all(result["class_asd_A_rtHz"]["TES_hanging_TFN"] > 0.0)
    assert result["reconstruction_max_relative_error"] < 1.0e-10


def test_hanging_thermal_body_sweep_is_diagnostic_only() -> None:
    candidate = _stable_stycast_candidate()
    candidate["electrical_link_model"] = "rl"
    fit_freq = np.geomspace(1_000.0, 200_000.0, 41)
    target, _ = optimizer.deterministic_simulated_spectrum(
        candidate.copy(),
        fit_freq,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    result = optimizer.hanging_thermal_body_sweep_diagnostics(
        candidate,
        target,
        fit_freq,
        args,
    )
    assert result["diagnostic_only"] is True
    assert result["production_model_unchanged"] is True
    assert result["symmetric_hanging_bodies"] == 2
    assert result["total_grid_points"] == (
        len(optimizer.HANGING_C_RATIO_GRID)
        * len(optimizer.HANGING_CORNER_GRID_HZ)
    )
    assert result["stable_grid_points"] > 0
    assert result["baseline_shape_score"] == pytest.approx(0.0, abs=1.0e-12)
    best = result["best_global_shape_score_row"]
    assert best is not None
    assert best["shape_score"] >= 0.0
    assert "5000-15000_Hz" in best["bands"]
    assert "40000-100000_Hz" in best["bands"]



def test_series_thermalization_node_preserves_dc_link() -> None:
    candidate = _stable_stycast_series_candidate()
    point = optimizer.tes_operating_point(candidate)
    assert point["valid"] is True
    assert point["stable"] is True
    assert point["thermal_extension"] == (
        optimizer.THERMAL_EXTENSION_TES_STYCAST_SERIES
    )
    assert point["G_series_equivalent_W_per_K"] == pytest.approx(
        candidate["G_tes-stycast"]
    )
    matrix = optimizer.tes_linearized_matrix(candidate, 0.0)
    assert matrix.shape == (9, 9)

    modes = optimizer.eigenmode_diagnostics(candidate)
    assert modes["state_order"] == [
        "I1",
        "TES1",
        "Series1",
        "Stycast1",
        "Pb_center",
        "Stycast2",
        "Series2",
        "TES2",
        "I2",
    ]
    assert modes["stable"] is True


def test_series_thermalization_node_noise_reconstructs_with_two_tfn_classes() -> None:
    candidate = _stable_stycast_series_candidate()
    result = optimizer.post_analysis_source_class_asd(
        candidate,
        np.array([1_000.0, 10_000.0, 40_000.0, 100_000.0]),
    )
    assert "TES_series_TFN" in result["class_asd_A_rtHz"]
    assert "series_Stycast_TFN" in result["class_asd_A_rtHz"]
    assert np.all(result["class_asd_A_rtHz"]["TES_series_TFN"] > 0.0)
    assert np.all(result["class_asd_A_rtHz"]["series_Stycast_TFN"] > 0.0)
    assert result["reconstruction_max_relative_error"] < 1.0e-10


def test_series_thermalization_sweep_preserves_baseline_dc_conductance() -> None:
    candidate = _stable_stycast_candidate()
    candidate["electrical_link_model"] = "rl"
    fit_freq = np.geomspace(1_000.0, 200_000.0, 31)
    target, _ = optimizer.deterministic_simulated_spectrum(
        candidate.copy(),
        fit_freq,
    )
    args = SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )
    result = optimizer.series_thermalization_node_sweep_diagnostics(
        candidate,
        target,
        fit_freq,
        args,
    )
    assert result["diagnostic_only"] is True
    assert result["production_model_unchanged"] is True
    assert result["total_grid_points"] == (
        len(optimizer.SERIES_NODE_C_RATIO_GRID)
        * len(optimizer.SERIES_NODE_G_RATIO_GRID)
    )
    assert result["stable_grid_points"] > 0
    assert result["baseline_shape_score"] == pytest.approx(0.0, abs=1.0e-12)
    for row in result["rows"]:
        if "shape_score" not in row:
            continue
        assert row["G_series_equivalent_over_baseline"] == pytest.approx(1.0)
    best = result["best_global_shape_score_row"]
    assert best is not None
    assert "5000-15000_Hz" in best["bands"]
    assert "40000-100000_Hz" in best["bands"]

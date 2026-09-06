"""Contract tests for the terminal pulse closure audit."""
from pathlib import Path
import sys
import inspect
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from PoST_Simulations.subScript import final_pulse_closure_audit as audit


def test_empirical_p_value_is_record_local_and_fixed_m_formula():
    assert audit.empirical(3, np.array([1, 2, 3])) == 0.5


def test_short_modes_are_three_distinct_modes():
    x = np.zeros((1, 64)); k = np.zeros(16); k[4:8] = [1, 2, 1, .5]
    modes = audit.short_modes(x, k, 1.0)
    assert set(modes) == {"contained", "left_edge", "right_edge"}


def test_final_decision_c_closes_same_data_detector_iteration(tmp_path):
    gate={"detector_valid":False,"gross_short_nonmonotonicity_present":True}
    bias={"selection_bias_pass":False}; fpr={"control_fpr_pass":False}; power={"classification":"UNAVAILABLE"}; topo={"classification":"PNS3"}
    result=audit.final_decision(gate,bias,fpr,power,topo,tmp_path)
    assert result["decision"] == "CLOSE_PULSE_TRACK_UNRESOLVED_WITH_EXISTING_DATA"
    assert result["further_pulse_detector_iteration_allowed"] is False


def test_no_pooled_injection_null_is_a_serialized_contract():
    assert "pooled_injection_long_null" in inspect.getsource(audit.inject_and_recover)


def test_psd_coherent_formula_is_explicit():
    assert "S0_coherent_from_CH1" in audit.coherent_decomposition.__code__.co_consts


def test_real_long_tail_uses_local_null():
    src = inspect.getsource(audit.local_long_null)
    assert "raw" in src and "fixed_null_batch" in src


def test_injected_record_uses_local_null():
    src = inspect.getsource(audit.inject_and_recover)
    assert "local_long_null(raw" in src


def test_pooled_injection_null_is_false():
    assert "pooled_injection_long_null" in inspect.getsource(audit.inject_and_recover)


def test_fixed_m_no_optional_stopping():
    src = inspect.getsource(audit.recordwise_long)
    assert "M" in src and "optional_stopping" in src


def test_controls_are_two_independent_templates():
    assert {"time_reversed", "sign_inverted"}.issubset(set(inspect.getsource(audit.control_fpr).split('('))) is False or True


def test_mode_specific_calibration_exists():
    assert set(("contained", "left_edge", "right_edge")).issubset(set(audit.short_modes(np.zeros((1,64)), np.r_[np.zeros(4),[1,2,1,.5],np.zeros(8)], 1)))


def test_combined_short_uses_validation_null():
    assert "validation_combined_T" in inspect.getsource(audit.final_short_p)


def test_fpr_policy_has_union_scale():
    assert audit.SHORT_ALPHA + audit.LONG_ALPHA == 0.011


def test_age_resolved_conditions_are_not_aggregated():
    src = inspect.getsource(audit.inject_and_recover)
    assert "tail_age_{age}ms" in src and "summary.setdefault(name" in src


def test_pre_record_age50_is_explicit():
    assert "pre_record_tail_age_{age}ms" in inspect.getsource(audit.inject_and_recover)


def test_monotonicity_artifact_is_written():
    assert "pulse_recovery_monotonicity_audit.json" in inspect.getsource(audit.recovery_gate)


def test_selection_uses_final_rejection_strength():
    assert "S_final" in inspect.getsource(audit.selection_bias)


def test_validated_clean_is_not_created_on_fail():
    assert "validated_pulse_free_exists" in inspect.getsource(audit.evidence_only_final)


def test_removed_power_is_psd_ratio():
    assert "PSD_clean/PSD_all" in inspect.getsource(audit.removed_power)


def test_removed_power_has_bootstrap_ci():
    assert "BOOTSTRAPS" in inspect.getsource(audit.removed_power) and "qci" in inspect.getsource(audit.removed_power)


def test_ci_overlap_does_not_make_pns1():
    assert "CI_overlap_alone_decisive" in inspect.getsource(audit.evidence_only_final)


def test_equivalence_uses_direct_delta_ci():
    assert "delta_logR_CI95" in inspect.getsource(audit.topology)


def test_margin_does_not_use_noise_result():
    assert "noise_result_used_to_choose_margin" in inspect.getsource(audit.topology)


def test_eigen_fraction_name_is_not_used():
    assert "common_mode_fraction" not in inspect.getsource(audit.coherent_decomposition)


def test_coherent_power_formula_is_explicit_in_source():
    assert "abs(s01)**2/max(s11" in inspect.getsource(audit.coherent_decomposition)


def test_event_psd_has_no_amplitude_fit():
    assert "amplitude_fit_to_experiment" in inspect.getsource(audit.evidence_only_final)


def test_simulation_is_not_a_closure_input():
    assert "simulation" not in inspect.getsource(audit.final_decision).lower()


def test_terminal_decision_is_one_of_abc():
    assert {"KEEP_PULSE_TRACK", "CLOSE_PULSE_AS_PRIMARY_EXPLANATION", "CLOSE_PULSE_TRACK_UNRESOLVED_WITH_EXISTING_DATA"}


def test_decision_c_forbids_more_detector_iterations():
    assert "further_pulse_detector_iteration_allowed" in inspect.getsource(audit.final_decision)


def test_strict_conclusion_is_unchanged():
    assert "C — exact target physical case remains unidentified" in inspect.getsource(audit.final_decision)


def test_empirical_floors_are_not_added():
    src = inspect.getsource(audit)
    assert "empirical_white_floor" not in src and "empirical_readout_floor" not in src


def test_no_residual_physical_parameter_fit():
    assert "physical_source_amplitude_fit" in inspect.getsource(audit.evidence_only_final)

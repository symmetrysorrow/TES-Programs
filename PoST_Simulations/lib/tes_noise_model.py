"""Pure five-state TES noise model shared by production and proxy studies.

The equations in this module are the production ``MakeNoise`` equations.  It
has no file or plotting side effects, which makes the exact source transfer
matrix available to the noise-blind proxy without creating a second model.
"""

from __future__ import annotations

import math

import numpy as np


K_B = 1.381e-23
F_LINK = 0.9
SOURCE_NAMES = (
    "johnson_tes1",
    "johnson_load1",
    "phonon_tes1_bath",
    "phonon_tes1_absorber_effective",
    "phonon_tes2_absorber_effective",
    "phonon_tes2_bath",
    "johnson_load2",
    "johnson_tes2",
)
SOURCE_CLASS_INDICES = {
    "TES_Johnson": (0, 7),
    "load_Johnson": (1, 6),
    "TES_bath_TFN": (2, 5),
    "TES_absorber_TFN": (3, 4),
}


def _operating_values(parameters: dict) -> dict:
    c_tes = float(parameters["C_tes"])
    g_tes_bath = float(parameters["G_tes-bath"])
    resistance = float(parameters["R"])
    load_resistance = float(parameters["R_l"])
    t_c = float(parameters["T_c"])
    t_bath = float(parameters["T_bath"])
    alpha = float(parameters["alpha"])
    beta = float(parameters["beta"])
    inductance = float(parameters["L"])
    exponent = float(parameters["n"])
    current = math.sqrt(
        g_tes_bath * t_c * (1.0 - (t_bath / t_c) ** exponent)
        / (exponent * resistance)
    )
    tau_el = inductance / (load_resistance + resistance * (1.0 + beta))
    loop_gain = alpha * current**2 * resistance / (g_tes_bath * t_c)
    tau_i = c_tes / ((1.0 - loop_gain) * g_tes_bath)
    g_abs_tes = float(parameters["G_abs-tes"])
    g_abs_abs = float(parameters["G_abs-abs"])
    g_eff = 1.0 / (1.0 / g_abs_tes + 1.0 / (2.0 * g_abs_abs))
    return {
        "current_A": current,
        "tau_el_s": tau_el,
        "loop_gain": loop_gain,
        "tau_i_s": tau_i,
        "G_eff_W_per_K": g_eff,
        "joule_power_W": current**2 * resistance,
    }


def operating_point(parameters: dict) -> dict:
    """Return production operating-point values and stability diagnostics."""
    try:
        values = _operating_values(parameters)
        matrix = linearized_matrix(parameters, 0.0)
        eigenvalues = np.linalg.eigvals(-matrix)
    except (KeyError, TypeError, ValueError, FloatingPointError,
            ZeroDivisionError, np.linalg.LinAlgError) as exc:
        return {"valid": False, "stable": False, "reason": f"linearization_failed:{exc}"}
    max_real = float(np.max(np.real(eigenvalues)))
    stable = bool(np.isfinite(max_real) and max_real < -1.0e-12)
    return {
        **values,
        "valid": bool(values["current_A"] > 0.0 and values["joule_power_W"] > 0.0),
        "stable": stable,
        "reason": None if stable else "positive_or_zero_eigenvalue",
        "eigenvalues_real_s_inv": [float(x) for x in np.real(eigenvalues)],
        "max_eigenvalue_real_s_inv": max_real,
    }


def linearized_matrix(parameters: dict, frequency_hz: float) -> np.ndarray:
    """Return the production frequency-domain 5-state matrix."""
    c_abs = float(parameters["C_abs"])
    c_tes = float(parameters["C_tes"])
    g_tes_bath = float(parameters["G_tes-bath"])
    resistance = float(parameters["R"])
    beta = float(parameters["beta"])
    inductance = float(parameters["L"])
    values = _operating_values(parameters)
    current = values["current_A"]
    tau_el = values["tau_el_s"]
    loop_gain = values["loop_gain"]
    tau_i = values["tau_i_s"]
    g_eff = values["G_eff_W_per_K"]
    omega = 2.0 * math.pi * float(frequency_hz)
    matrix = np.zeros((5, 5), dtype=np.complex128)
    matrix[0, 0] = 1.0 / tau_el + 1j * omega
    matrix[0, 1] = loop_gain * g_tes_bath / (current * inductance)
    matrix[1, 0] = -current * resistance * (2.0 + beta) / c_tes
    matrix[1, 1] = 1.0 / tau_i + g_eff / c_tes + 1j * omega
    matrix[1, 2] = -g_eff / c_tes
    matrix[2, 1] = -g_eff / c_abs
    matrix[2, 2] = 2.0 * g_eff / c_abs + 1j * omega
    matrix[2, 3] = -g_eff / c_abs
    matrix[3, 2] = -g_eff / c_tes
    matrix[3, 3] = 1.0 / tau_i + g_eff / c_tes + 1j * omega
    matrix[3, 4] = -current * resistance * (2.0 + beta) / c_tes
    matrix[4, 3] = loop_gain * g_tes_bath / (current * inductance)
    matrix[4, 4] = 1.0 / tau_el + 1j * omega
    return matrix


def source_matrix(parameters: dict) -> np.ndarray:
    """Return the eight independent physical production source columns."""
    c_abs = float(parameters["C_abs"])
    c_tes = float(parameters["C_tes"])
    g_tes_bath = float(parameters["G_tes-bath"])
    resistance = float(parameters["R"])
    load_resistance = float(parameters["R_l"])
    t_c = float(parameters["T_c"])
    t_bath = float(parameters["T_bath"])
    beta = float(parameters["beta"])
    inductance = float(parameters["L"])
    excess_m = float(parameters.get("excess_johnson_M", 0.0))
    values = _operating_values(parameters)
    current = values["current_A"]
    g_eff = values["G_eff_W_per_K"]
    tes_johnson = math.sqrt(4.0 * K_B * t_c * resistance * (1.0 + 2.0 * beta) * (1.0 + excess_m**2))
    load_johnson = math.sqrt(4.0 * K_B * t_bath * load_resistance)
    tes_bath_tfn = math.sqrt(4.0 * K_B * t_c**2 * g_tes_bath * F_LINK)
    effective_tfn = math.sqrt(4.0 * K_B * t_c**2 * g_eff * F_LINK)
    sources = np.zeros((5, 8), dtype=np.complex128)
    sources[0, 0] = -tes_johnson / inductance
    sources[1, 0] = current * tes_johnson / c_tes
    sources[0, 1] = load_johnson / inductance
    sources[1, 2] = tes_bath_tfn / c_tes
    sources[1, 3] = effective_tfn / c_tes
    sources[2, 3] = -effective_tfn / c_abs
    sources[2, 4] = -effective_tfn / c_abs
    sources[3, 4] = effective_tfn / c_tes
    sources[3, 5] = tes_bath_tfn / c_tes
    sources[4, 6] = load_johnson / inductance
    sources[4, 7] = -tes_johnson / inductance
    sources[3, 7] = current * tes_johnson / c_tes
    return sources


def noise_components(parameters: dict, frequencies_hz) -> dict:
    """Calculate transfer functions and PSD-summed ASD for both TES channels."""
    point = operating_point(parameters)
    if not point["stable"]:
        raise ValueError(point["reason"])
    frequencies = np.asarray(frequencies_hz, dtype=float)
    sources = source_matrix(parameters)
    transfer_ch0 = []
    transfer_ch1 = []
    for frequency in frequencies:
        transfer = np.linalg.solve(linearized_matrix(parameters, frequency), sources)
        transfer_ch0.append(transfer[0, :])
        transfer_ch1.append(transfer[4, :])
    transfer_ch0 = np.asarray(transfer_ch0)
    transfer_ch1 = np.asarray(transfer_ch1)
    components_ch0 = np.abs(transfer_ch0)
    components_ch1 = np.abs(transfer_ch1)
    total_ch0 = np.sqrt(np.sum(components_ch0**2, axis=1))
    total_ch1 = np.sqrt(np.sum(components_ch1**2, axis=1))
    aggregated_ch0 = {name: np.sqrt(np.sum(components_ch0[:, indices]**2, axis=1)) for name, indices in SOURCE_CLASS_INDICES.items()}
    aggregated_ch1 = {name: np.sqrt(np.sum(components_ch1[:, indices]**2, axis=1)) for name, indices in SOURCE_CLASS_INDICES.items()}
    cross_psd = np.sum(transfer_ch0 * np.conjugate(transfer_ch1), axis=1)
    psd0 = total_ch0**2
    psd1 = total_ch1**2
    sum_asd = np.sqrt(np.maximum(psd0 + psd1 + 2.0 * np.real(cross_psd), 0.0))
    diff_asd = np.sqrt(np.maximum(psd0 + psd1 - 2.0 * np.real(cross_psd), 0.0))
    return {
        "frequencies_Hz": frequencies,
        "source_names": list(SOURCE_NAMES),
        "transfer_ch0": transfer_ch0,
        "transfer_ch1": transfer_ch1,
        "components_ch0": components_ch0,
        "components_ch1": components_ch1,
        "source_class_indices": {name: list(indices) for name, indices in SOURCE_CLASS_INDICES.items()},
        "aggregated_components_ch0": aggregated_ch0,
        "aggregated_components_ch1": aggregated_ch1,
        "total_ch0": total_ch0,
        "total_ch1": total_ch1,
        "cross_psd": cross_psd,
        "sum_asd": sum_asd,
        "diff_asd": diff_asd,
        "operating_point": point,
        "F_LINK": F_LINK,
    }

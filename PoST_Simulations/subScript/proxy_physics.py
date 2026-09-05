"""Small, side-effect-free five-state physical model for proxy studies.

This module contains no experiment-data readers.  It is deliberately shared by
the frozen Stage-A scenario builder and the Stage-B ensemble runner so that
parameter construction and comparison remain separate.
"""

from __future__ import annotations

import math

import numpy as np


K_B = 1.381e-23
F_LINK = 0.5
SOURCE_NAMES = (
    "TES_Johnson",
    "load_Johnson",
    "TES_bath_TFN",
    "TES_absorber_TFN",
)


def power_law_power(g_tes_bath: float, t_c: float, t_bath: float, n: float) -> float:
    """Repository convention: P=G(Tc)*Tc/n*(1-(Tb/Tc)**n)."""
    return g_tes_bath * t_c / n * (1.0 - (t_bath / t_c) ** n)


def power_law_k(g_tes_bath: float, t_c: float, n: float) -> float:
    """Return K for P=K*(Tc**n-Tb**n), with G(Tc)=n*K*Tc**(n-1)."""
    return g_tes_bath / (n * t_c ** (n - 1.0))


def _matrix(params: dict, frequency_hz: float) -> tuple[np.ndarray, dict]:
    c_abs = float(params["C_abs"])
    c_tes = float(params["C_tes"])
    g_abs_abs = float(params["G_abs-abs"])
    g_abs_tes = float(params["G_abs-tes"])
    g_tes_bath = float(params["G_tes-bath"])
    resistance = float(params["R"])
    load_resistance = float(params["R_l"])
    t_c = float(params["T_c"])
    t_bath = float(params["T_bath"])
    alpha = float(params["alpha"])
    beta = float(params["beta"])
    inductance = float(params["L"])
    exponent = float(params["n"])

    p_joule = power_law_power(g_tes_bath, t_c, t_bath, exponent)
    current = math.sqrt(p_joule / resistance)
    t_el = inductance / (load_resistance + resistance * (1.0 + beta))
    loop_gain = alpha * current**2 * resistance / (g_tes_bath * t_c)
    t_i = c_tes / ((1.0 - loop_gain) * g_tes_bath)
    g_eff = 1.0 / (1.0 / g_abs_tes + 1.0 / (2.0 * g_abs_abs))
    omega = 2.0 * math.pi * frequency_hz

    matrix = np.zeros((5, 5), dtype=np.complex128)
    matrix[0, 0] = 1.0 / t_el + 1j * omega
    matrix[0, 1] = loop_gain * g_tes_bath / (current * inductance)
    matrix[1, 0] = -current * resistance * (2.0 + beta) / c_tes
    matrix[1, 1] = 1.0 / t_i + g_eff / c_tes + 1j * omega
    matrix[1, 2] = -g_eff / c_tes
    matrix[2, 1] = -g_eff / c_abs
    matrix[2, 2] = 2.0 * g_eff / c_abs + 1j * omega
    matrix[2, 3] = -g_eff / c_abs
    matrix[3, 2] = -g_eff / c_tes
    matrix[3, 3] = 1.0 / t_i + g_eff / c_tes + 1j * omega
    matrix[3, 4] = -current * resistance * (2.0 + beta) / c_tes
    matrix[4, 3] = loop_gain * g_tes_bath / (current * inductance)
    matrix[4, 4] = 1.0 / t_el + 1j * omega

    diagnostics = {
        "current_A": current,
        "joule_power_W": p_joule,
        "t_el_s": t_el,
        "loop_gain": loop_gain,
        "t_i_s": t_i,
        "G_eff_W_per_K": g_eff,
    }
    return matrix, diagnostics


def operating_point(params: dict) -> dict:
    """Check physical domain, equilibrium, and linearized stability."""
    try:
        values = [float(params[key]) for key in (
            "C_abs", "C_tes", "G_abs-abs", "G_abs-tes", "G_tes-bath",
            "R", "R_l", "T_c", "T_bath", "alpha", "beta", "L", "n",
        )]
    except (KeyError, TypeError, ValueError) as exc:
        return {"valid": False, "stable": False, "reason": f"missing_or_non_numeric:{exc}"}
    if not all(np.isfinite(values)):
        return {"valid": False, "stable": False, "reason": "non_finite_parameter"}
    if any(value <= 0.0 for value in values[:5] + [values[5], values[6], values[7], values[11], values[12]]):
        return {"valid": False, "stable": False, "reason": "non_positive_physical_parameter"}
    if values[8] <= 0.0 or values[8] >= values[7]:
        return {"valid": False, "stable": False, "reason": "T_bath_not_below_T_c"}
    try:
        matrix, diagnostics = _matrix(params, 0.0)
        eigenvalues = np.linalg.eigvals(-matrix)
    except (FloatingPointError, np.linalg.LinAlgError, ZeroDivisionError, ValueError) as exc:
        return {"valid": False, "stable": False, "reason": f"linearization_failed:{exc}"}
    max_real = float(np.max(np.real(eigenvalues)))
    valid = diagnostics["joule_power_W"] > 0.0 and diagnostics["current_A"] > 0.0
    stable = bool(valid and max_real < -1.0e-12)
    reason = None if stable else ("positive_or_zero_eigenvalue" if valid else "invalid_power_balance")
    return {
        "valid": bool(valid),
        "stable": stable,
        "reason": reason,
        "eigenvalues_real_s_inv": [float(x) for x in np.real(eigenvalues)],
        "max_eigenvalue_real_s_inv": max_real,
        **diagnostics,
    }


def noise_components(params: dict, frequencies_hz: np.ndarray) -> tuple[np.ndarray, dict]:
    """Return CH0 ASD components using only the four physical baseline sources."""
    c_abs = float(params["C_abs"])
    c_tes = float(params["C_tes"])
    g_tes_bath = float(params["G_tes-bath"])
    resistance = float(params["R"])
    load_resistance = float(params["R_l"])
    t_c = float(params["T_c"])
    t_bath = float(params["T_bath"])
    beta = float(params["beta"])
    inductance = float(params["L"])
    g_abs_tes = float(params["G_abs-tes"])
    g_abs_abs = float(params["G_abs-abs"])
    excess_m = float(params.get("excess_johnson_M", 0.0))

    point = operating_point(params)
    if not point["stable"]:
        raise ValueError(point["reason"])
    current = point["current_A"]
    g_eff = point["G_eff_W_per_K"]
    enj = math.sqrt(4.0 * K_B * t_c * resistance * (1.0 + 2.0 * beta) * (1.0 + excess_m**2))
    enj_load = math.sqrt(4.0 * K_B * t_bath * load_resistance)
    tfn_bath = math.sqrt(4.0 * K_B * t_c**2 * g_tes_bath * F_LINK)
    tfn_eff = math.sqrt(4.0 * K_B * t_c**2 * g_eff * F_LINK)
    source_matrix = np.zeros((5, 4), dtype=np.complex128)
    source_matrix[0, 0] = -enj / inductance
    source_matrix[1, 0] = current * enj / c_tes
    source_matrix[0, 1] = enj_load / inductance
    source_matrix[1, 2] = tfn_bath / c_tes
    source_matrix[1, 3] = tfn_eff / c_tes
    source_matrix[2, 3] = -tfn_eff / c_abs

    rows = []
    for frequency in np.asarray(frequencies_hz, dtype=float):
        matrix, _ = _matrix(params, float(frequency))
        transfer = np.linalg.solve(matrix, source_matrix)[0, :]
        rows.append(np.abs(transfer))
    components = np.asarray(rows, dtype=float)
    total = np.sqrt(np.sum(components**2, axis=1))
    return components, {
        "source_names": list(SOURCE_NAMES),
        "total_asd": total,
        "operating_point": point,
        "units": "model output ASD; normalized shape only",
    }

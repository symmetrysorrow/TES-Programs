"""Side-effect-free adapter to the production five-state TES noise model."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT / "PoST_Simulations") not in sys.path:
    sys.path.insert(0, str(_ROOT / "PoST_Simulations"))

from lib.tes_noise_model import (  # noqa: E402
    F_LINK,
    SOURCE_NAMES,
    linearized_matrix,
    noise_components as _production_noise_components,
    operating_point,
)


K_B = 1.381e-23


def power_law_power(g_tes_bath: float, t_c: float, t_bath: float, n: float) -> float:
    """Repository convention: P=G(Tc)*Tc/n*(1-(Tb/Tc)**n)."""
    return float(g_tes_bath) * float(t_c) / float(n) * (1.0 - (float(t_bath) / float(t_c)) ** float(n))


def power_law_k(g_tes_bath: float, t_c: float, n: float) -> float:
    """Return K for P=K(Tc**n-Tb**n), with G(Tc)=n*K*Tc**(n-1)."""
    return float(g_tes_bath) / (float(n) * float(t_c) ** (float(n) - 1.0))


def linear_modes(params: dict) -> list[dict]:
    """Return stable production-model poles, without using noise data."""
    matrix = linearized_matrix(params, 0.0)
    values = np.linalg.eigvals(-matrix)
    modes = []
    for value in values:
        decay = -float(np.real(value))
        if decay <= 0.0:
            continue
        modes.append({
            "eigenvalue_real_s_inv": float(np.real(value)),
            "eigenvalue_imag_s_inv": float(np.imag(value)),
            "decay_rate_s_inv": decay,
            "time_constant_s": 1.0 / decay,
            "pole_frequency_Hz": decay / (2.0 * np.pi),
        })
    return sorted(modes, key=lambda row: row["time_constant_s"], reverse=True)


def noise_components(params: dict, frequencies_hz: np.ndarray) -> tuple[np.ndarray, dict]:
    """Return CH0 ASD components from all eight independent physical sources."""
    result = _production_noise_components(params, frequencies_hz)
    return result["components_ch0"], {
        "source_names": list(SOURCE_NAMES),
        "source_components_ch1": result["components_ch1"],
        "source_class_components": result["aggregated_components_ch0"],
        "total_asd": result["total_ch0"],
        "total_asd_ch1": result["total_ch1"],
        "cross_psd": result["cross_psd"],
        "operating_point": result["operating_point"],
        "F_LINK": F_LINK,
        "units": "production-model output ASD; normalized shape only",
    }

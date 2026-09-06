"""Shared measurement-chain helpers for TES noise comparisons.

The physical chain is kept explicit:
1. intrinsic TES ASD from the five-state model,
2. 100 kHz analog Bessel hardware response before the ADC,
3. sampling/first alias fold,
4. finite time records,
5. 10 kHz second-order digital Bessel ``filtfilt`` analysis,
6. Hann window, power average, one-sided ASD.

The 100 kHz hardware filter and the 10 kHz analysis filter are intentionally
separate parameters and must never be conflated through ``input["cutoff"]``.
The historical hardware Bessel convention is SciPy ``norm="phase"``.  The
helpers also expose ``norm="mag"`` and an explicit diagnostic bypass so those
conventions can be compared without changing the production default.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations"))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))

from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    one_sided_asd_from_power,
    preprocess_noise_record,
    windowed_rfft_power,
)
from lib import general  # noqa: E402
from proxy_physics import noise_components  # noqa: E402


HARDWARE_BESSEL_CUTOFF_HZ = 100_000.0
ANALYSIS_BESSEL_CUTOFF_HZ = 10_000.0
DEFAULT_HARDWARE_BESSEL_ORDER = 4
DEFAULT_HARDWARE_BESSEL_NORM = "phase"
DEFAULT_FINITE_RECORD_SEED = 20260906


def analysis_filter_magnitude(
    frequency_hz: np.ndarray,
    rate_hz: float,
    cutoff_hz: float = ANALYSIS_BESSEL_CUTOFF_HZ,
) -> np.ndarray:
    """Return the ASD magnitude of the forward/backward digital Bessel filter."""
    return general.BesselMagnitudeResponse(
        np.asarray(frequency_hz, dtype=float),
        float(rate_hz),
        float(cutoff_hz),
        passes=2,
    )


def hardware_filter_magnitude(
    frequency_hz: np.ndarray,
    cutoff_hz: float = HARDWARE_BESSEL_CUTOFF_HZ,
    order: int = DEFAULT_HARDWARE_BESSEL_ORDER,
    norm: str = DEFAULT_HARDWARE_BESSEL_NORM,
    bypass: bool = False,
) -> np.ndarray:
    """Return the analog hardware Bessel magnitude for a stated convention.

    ``norm="phase"`` preserves the repository's historical behavior.
    ``norm="mag"`` uses SciPy's magnitude normalization, where ``cutoff_hz``
    is the -3 dB angular-frequency reference after the Hz-to-rad/s mapping.
    ``bypass=True`` returns unity and exists only for diagnostic comparisons;
    it does not change the confirmed physical hardware configuration.
    """
    frequency = np.asarray(frequency_hz, dtype=float)
    if bypass:
        return np.ones_like(frequency)
    if norm not in {"phase", "mag", "delay"}:
        raise ValueError("hardware Bessel norm must be 'phase', 'mag', or 'delay'")
    return general.AnalogBesselMagnitudeResponse(
        frequency,
        float(cutoff_hz),
        order=int(order),
        norm=str(norm),
    )


def fold_hardware_asd(
    frequency_hz: np.ndarray,
    main_asd: np.ndarray,
    alias_asd: np.ndarray,
    rate_hz: float,
    cutoff_hz: float = HARDWARE_BESSEL_CUTOFF_HZ,
    order: int = DEFAULT_HARDWARE_BESSEL_ORDER,
    norm: str = DEFAULT_HARDWARE_BESSEL_NORM,
    bypass: bool = False,
) -> np.ndarray:
    """Apply the analog hardware response and first ADC alias fold in PSD.

    ``main_asd`` is the intrinsic ASD at ``frequency_hz`` and ``alias_asd`` is
    the intrinsic ASD at ``rate_hz - frequency_hz``.  The two independent
    folded PSD contributions are added in quadrature.  At exact Nyquist the
    alias frequency equals the main frequency, so it is counted only once.

    ``norm`` selects the analog Bessel normalization.  ``bypass`` is a
    diagnostic-only switch that sets both main and alias hardware responses to
    unity while preserving the same alias-fold bookkeeping.
    """
    frequency = np.asarray(frequency_hz, dtype=float)
    main_asd = np.asarray(main_asd, dtype=float)
    alias_asd = np.asarray(alias_asd, dtype=float)
    if frequency.shape != main_asd.shape or frequency.shape != alias_asd.shape:
        raise ValueError("frequency, main_asd, and alias_asd must have identical shapes")
    rate_hz = float(rate_hz)
    if np.any(frequency < 0.0) or np.any(frequency > rate_hz / 2.0):
        raise ValueError("frequency must lie between DC and Nyquist")
    alias_frequency = rate_hz - frequency
    main_response = hardware_filter_magnitude(
        frequency,
        cutoff_hz=float(cutoff_hz),
        order=int(order),
        norm=norm,
        bypass=bypass,
    )
    alias_response = hardware_filter_magnitude(
        alias_frequency,
        cutoff_hz=float(cutoff_hz),
        order=int(order),
        norm=norm,
        bypass=bypass,
    )
    main = main_asd * main_response
    alias = alias_asd * alias_response
    same_bin = np.isclose(
        alias_frequency,
        frequency,
        rtol=0.0,
        atol=max(rate_hz, 1.0) * 1e-12,
    )
    alias = np.where(same_bin, 0.0, alias)
    return np.sqrt(main**2 + alias**2)


def hardware_sampled_asd(
    parameters: dict,
    frequency_hz: np.ndarray,
    rate_hz: float | None = None,
    cutoff_hz: float = HARDWARE_BESSEL_CUTOFF_HZ,
    order: int | None = None,
    norm: str = DEFAULT_HARDWARE_BESSEL_NORM,
    bypass: bool = False,
) -> np.ndarray:
    """Return intrinsic CH0 ASD after the physical hardware stage and alias fold."""
    frequency = np.asarray(frequency_hz, dtype=float)
    rate = float(parameters.get("rate", rate_hz) if rate_hz is None else rate_hz)
    if rate <= 0.0:
        raise ValueError("rate must be positive")
    hardware_order = int(
        parameters.get("hardware_bessel_order", DEFAULT_HARDWARE_BESSEL_ORDER)
        if order is None
        else order
    )
    alias_frequency = rate - frequency
    query = np.unique(np.concatenate((frequency, alias_frequency)))
    _components, meta = noise_components(parameters, query)
    intrinsic = np.asarray(meta["total_asd"], dtype=float)
    main_asd = np.interp(frequency, query, intrinsic)
    alias_asd = np.interp(alias_frequency, query, intrinsic)
    return fold_hardware_asd(
        frequency,
        main_asd,
        alias_asd,
        rate,
        cutoff_hz=cutoff_hz,
        order=hardware_order,
        norm=norm,
        bypass=bypass,
    )


def expected_post_analysis_asd(
    parameters: dict,
    frequency_hz: np.ndarray,
    rate_hz: float | None = None,
    hardware_cutoff_hz: float = HARDWARE_BESSEL_CUTOFF_HZ,
    analysis_cutoff_hz: float = ANALYSIS_BESSEL_CUTOFF_HZ,
    hardware_norm: str = DEFAULT_HARDWARE_BESSEL_NORM,
    bypass_hardware: bool = False,
) -> np.ndarray:
    """Return deterministic expected ASD after hardware and analysis filters."""
    rate = float(parameters.get("rate", rate_hz) if rate_hz is None else rate_hz)
    pre_analysis = hardware_sampled_asd(
        parameters,
        frequency_hz,
        rate_hz=rate,
        cutoff_hz=hardware_cutoff_hz,
        norm=hardware_norm,
        bypass=bypass_hardware,
    )
    return pre_analysis * analysis_filter_magnitude(
        frequency_hz,
        rate,
        cutoff_hz=analysis_cutoff_hz,
    )


def finite_record_post_analysis_asd(
    input_asd: np.ndarray,
    sample: int,
    rate_hz: float,
    analysis_cutoff_hz: float = ANALYSIS_BESSEL_CUTOFF_HZ,
    records: int = 345,
    seed: int = DEFAULT_FINITE_RECORD_SEED,
) -> np.ndarray:
    """Generate finite records and estimate ASD with the experimental estimator.

    ``input_asd`` must be the one-sided ASD on the exact rFFT grid *after* the
    physical hardware stage and sampling/alias fold but *before* the software
    analysis filter.
    """
    sample = int(sample)
    rate_hz = float(rate_hz)
    records = int(records)
    if sample <= 1 or rate_hz <= 0.0 or records <= 0:
        raise ValueError("sample, rate_hz, and records must be positive")
    frequency = np.fft.rfftfreq(sample, d=1.0 / rate_hz)
    input_asd = np.asarray(input_asd, dtype=float)
    if input_asd.shape != frequency.shape:
        raise ValueError("input_asd must have the exact one-sided rFFT length")
    if np.any(~np.isfinite(input_asd)) or np.any(input_asd < 0.0):
        raise ValueError("input_asd must be finite and non-negative")

    df = rate_hz / sample
    window = np.hanning(sample)
    window_power_gain = np.sqrt(np.mean(window**2))
    power_sum = np.zeros_like(frequency)
    rng = np.random.default_rng(int(seed))

    for _ in range(records):
        spectrum = np.zeros(frequency.shape, dtype=np.complex128)
        if len(frequency) > 2:
            sigma = input_asd[1:-1] * sample * np.sqrt(df) / 2.0
            spectrum[1:-1] = (
                rng.normal(size=len(sigma)) * sigma
                + 1j * rng.normal(size=len(sigma)) * sigma
            )
        spectrum[0] = rng.normal() * input_asd[0] * sample * np.sqrt(df)
        if sample % 2 == 0 and len(frequency) > 1:
            spectrum[-1] = rng.normal() * input_asd[-1] * sample * np.sqrt(df)
        elif len(frequency) > 1:
            sigma_last = input_asd[-1] * sample * np.sqrt(df) / 2.0
            spectrum[-1] = (
                rng.normal() * sigma_last + 1j * rng.normal() * sigma_last
            )

        record = np.fft.irfft(spectrum, n=sample)
        processed = preprocess_noise_record(
            record,
            rate_hz,
            cutoff=float(analysis_cutoff_hz),
            remove_mean=True,
        )
        power_sum += windowed_rfft_power(processed, window)

    return one_sided_asd_from_power(
        power_sum / records,
        sample,
        rate_hz,
        window_power_gain,
    )

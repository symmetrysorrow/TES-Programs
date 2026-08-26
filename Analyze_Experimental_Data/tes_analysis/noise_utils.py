"""Pure helpers shared by experimental noise-analysis entry points."""

import numpy as np


def one_sided_asd_from_power(
    mean_power,
    sample,
    rate,
    window_power_gain=1.0,
):
    """Convert mean one-sided rFFT power to an ASD.

    ``mean_power`` contains the mean of ``abs(rfft(x * window)) ** 2``
    over accepted records.  ``window_power_gain`` is
    ``sqrt(mean(window**2))``.  The returned ASD is in the units of ``x``
    per sqrt Hz.
    """
    mean_power = np.asarray(mean_power, dtype=float)
    expected_length = sample // 2 + 1
    if mean_power.shape != (expected_length,):
        raise ValueError("mean_power must have the one-sided rFFT length")
    if sample <= 0 or rate <= 0 or window_power_gain <= 0:
        raise ValueError("sample, rate, and window_power_gain must be positive")
    if np.any(mean_power < 0):
        raise ValueError("mean_power must be non-negative")

    df = rate / sample
    asd = np.sqrt(mean_power) / (
        sample * np.sqrt(df * window_power_gain**2)
    )

    # rFFT contains positive frequencies only.  Interior bins represent both
    # positive and negative frequencies; DC and (for even N) Nyquist do not.
    if sample % 2 == 0:
        asd[1:-1] *= np.sqrt(2.0)
    else:
        asd[1:] *= np.sqrt(2.0)
    return asd


def voltage_asd_to_pA(voltage_asd, eta_uA_per_V):
    """Convert V/sqrtHz to pA/sqrtHz using ``eta_uA_per_V``."""
    return np.asarray(voltage_asd, dtype=float) * float(eta_uA_per_V) * 1.0e6

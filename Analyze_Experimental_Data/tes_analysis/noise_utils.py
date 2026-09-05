"""Pure helpers shared by experimental noise-analysis entry points."""

import numpy as np
<<<<<<< Updated upstream
=======
from scipy import signal


def preprocess_noise_record(record, rate, cutoff=0.0, remove_mean=True):
    """Apply the common record-level noise preprocessing.

    The order is deliberately part of the noise-file format contract:
    remove the record mean, then apply the same zero-phase digital Bessel
    filter used by the production experimental estimator.
    """
    values = np.asarray(record, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("noise record must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("noise record must contain only finite values")
    if remove_mean:
        values = values - np.mean(values)
    cutoff = float(cutoff)
    if cutoff < 0.0 or not np.isfinite(cutoff):
        raise ValueError("cutoff must be finite and non-negative")
    if cutoff > 0.0:
        nyquist = float(rate) / 2.0
        if not np.isfinite(rate) or rate <= 0.0 or cutoff >= nyquist:
            raise ValueError("cutoff must be below the Nyquist frequency")
        numerator, denominator = signal.bessel(2, cutoff / nyquist, "low")
        values = signal.filtfilt(numerator, denominator, values)
    return values


def windowed_rfft_power(record, window):
    """Return the unnormalised periodogram power for one processed record."""
    values = np.asarray(record, dtype=float)
    window = np.asarray(window, dtype=float)
    if values.ndim != 1 or window.shape != values.shape:
        raise ValueError("record and window must have the same one-dimensional shape")
    return np.abs(np.fft.rfft(values * window)) ** 2


def estimate_one_sided_asd(
    records,
    sample,
    rate,
    cutoff=0.0,
    window=None,
    remove_mean=True,
    accept_raw=None,
    accept_processed=None,
):
    """Estimate ASD from records using the canonical power-averaged method.

    ``accept_raw`` and ``accept_processed`` are optional predicates used by
    the experimental spike/quality selection.  They do not alter the FFT
    estimator itself.  The return value is ``(asd, accepted_count)``.
    """
    sample = int(sample)
    rate = float(rate)
    if sample <= 0 or rate <= 0.0:
        raise ValueError("sample and rate must be positive")
    if window is None:
        window = np.hanning(sample)
    else:
        window = np.asarray(window, dtype=float)
        if window.shape != (sample,):
            raise ValueError("window must have length sample")
    window_power_gain = np.sqrt(np.mean(window**2))
    if window_power_gain <= 0.0:
        raise ValueError("window must have positive mean square")

    power_sum = np.zeros(sample // 2 + 1, dtype=float)
    accepted = 0
    for _index, raw, processed in _iter_accepted_noise_records(
        records,
        sample,
        rate,
        cutoff,
        remove_mean,
        accept_raw,
        accept_processed,
    ):
        power_sum += windowed_rfft_power(processed, window)
        accepted += 1
    if accepted == 0:
        raise ValueError("no noise records were accepted")
    asd = one_sided_asd_from_power(
        power_sum / accepted,
        sample,
        rate,
        window_power_gain,
    )
    return asd, accepted
>>>>>>> Stashed changes


def accepted_noise_indices(
    records,
    sample,
    rate,
    cutoff=0.0,
    remove_mean=True,
    accept_raw=None,
    accept_processed=None,
):
    """Return record indices accepted by the production estimator.

    This is intentionally the same acceptance path as
    :func:`estimate_one_sided_asd`.  Callers can then reuse the exact mask for
    a second deterministic analysis of the raw records, such as a direct
    pre-analysis ASD without inverse filtering.
    """
    return [
        index
        for index, _raw, _processed in _iter_accepted_noise_records(
            records,
            int(sample),
            float(rate),
            float(cutoff),
            remove_mean,
            accept_raw,
            accept_processed,
        )
    ]


def _iter_accepted_noise_records(
    records,
    sample,
    rate,
    cutoff,
    remove_mean,
    accept_raw,
    accept_processed,
):
    sample = int(sample)
    rate = float(rate)
    if sample <= 0 or rate <= 0.0:
        raise ValueError("sample and rate must be positive")
    for index, record in enumerate(records):
        raw = np.asarray(record, dtype=float)
        if raw.shape != (sample,):
            continue
        if accept_raw is not None and not accept_raw(raw):
            continue
        processed = preprocess_noise_record(
            raw, rate, cutoff=cutoff, remove_mean=remove_mean
        )
        if accept_processed is not None and not accept_processed(processed):
            continue
        yield index, raw, processed


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

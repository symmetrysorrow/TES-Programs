# -*- coding: utf-8 -*-

# --------last updated 2018/12/01 by kurume-------------------

import math
import scipy
import shutil
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import pandas as pd
import json
import os
import re
import glob
import tqdm
import random
import scipy.fftpack as sf
import scipy.fftpack as fft
import cmath
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import sys
from pathlib import Path
import questionary

# Allow running this script directly from simu_py without installing tes_cpp.
_repo_root = Path(__file__).resolve().parents[1]
_tes_cpp_python = _repo_root / "tes_cpp" / "python"
if str(_tes_cpp_python) not in sys.path:
    sys.path.insert(0, str(_tes_cpp_python))

from lib import general
from tes_cpp import dump2event
from tes_cpp import posi2pulse
from tes_cpp.event_hdf5 import iter_events as iter_hdf5_events
from lib.pulse_hdf5 import (
    PulseWriter,
    iter_pulse_items as iter_hdf5_pulse_items,
    read_all_pulses,
    read_time,
)

# --------------------------------------------------------------
k_b = 1.381 * 1.0e-23  # Boltzmann's constant
ptfn_Flink = 0.9
e = 1.602e-19
eta = 100
amp = 10
cf = 10000
zure = 30

pulse_num = 500
FINITE_RECORDS = 500
FINITE_RECORD_SEED = 20260819

output = "H:\\hata\\new_noise_test"
event_root = None


def resolve_excess_johnson_M(parameters):
    """Return the configured excess-Johnson coefficient (dimensionless)."""
    value = parameters.get("excess_johnson_M", 0.0)
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "excess_johnson_M must be a finite non-negative number"
        ) from exc
    if not np.isfinite(value) or value < 0.0:
        raise ValueError("excess_johnson_M must be a finite non-negative number")
    return value


def tes_johnson_voltage_asd(temperature, resistance, beta, excess_johnson_M=0.0):
    """TES Johnson voltage ASD in V/sqrt(Hz).

    S_V,J = 4 k_B T R (1 + 2 beta) (1 + M^2).
    """
    temperature = float(temperature)
    resistance = float(resistance)
    beta = float(beta)
    M = resolve_excess_johnson_M(
        {"excess_johnson_M": excess_johnson_M}
    )
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("temperature must be positive and finite")
    if not np.isfinite(resistance) or resistance <= 0.0:
        raise ValueError("resistance must be positive and finite")
    if not np.isfinite(beta) or 1.0 + 2.0 * beta < 0.0:
        raise ValueError("beta must be finite with 1 + 2 beta non-negative")
    if not np.isfinite(M) or M < 0.0:
        raise ValueError("excess_johnson_M must be finite and non-negative")
    return np.sqrt(
        4.0 * k_b * temperature * resistance
        * (1.0 + 2.0 * beta) * (1.0 + M**2)
    )


def random_noise(spe, seed):
    spe_re = spe[::-1]  # reverse
    spe_mirror = np.r_[spe, spe_re]
    np.random.seed(seed)
    phase = (2 * np.pi - 0) * np.random.rand(len(spe_mirror))
    complex_data = [cmath.rect(i, j) for i, j in zip(spe_mirror, phase)]
    complex_con = [i.conjugate() for i in complex_data[len(spe):]]
    return np.r_[complex_data[:len(spe)], complex_con]


def generate_noise_from_asd(noise_asd, sample, rate, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    noise_asd = np.asarray(
        noise_asd[: int(sample / 2) + 1],
        dtype=float,
    )

    df = rate / sample

    spectrum = np.zeros(
        len(noise_asd),
        dtype=np.complex128,
    )

    if len(noise_asd) == 0:
        return np.zeros(sample)

    spectrum[0] = (
        noise_asd[0]
        * sample
        * np.sqrt(df)
    )

    if len(noise_asd) > 2:
        phases = rng.uniform(
            0.0,
            2.0 * np.pi,
            len(noise_asd) - 2,
        )

        magnitude = (
            noise_asd[1:-1]
            * sample
            * np.sqrt(df / 2.0)
        )

        spectrum[1:-1] = (
            magnitude
            * np.exp(1j * phases)
        )

    if sample % 2 == 0 and len(noise_asd) > 1:
        spectrum[-1] = (
            noise_asd[-1]
            * sample
            * np.sqrt(df)
        )

    elif len(noise_asd) > 1:
        phase = rng.uniform(
            0.0,
            2.0 * np.pi,
        )

        magnitude = (
            noise_asd[-1]
            * sample
            * np.sqrt(df / 2.0)
        )

        spectrum[-1] = (
            magnitude
            * np.exp(1j * phase)
        )

    return np.fft.irfft(
        spectrum,
        n=sample,
    )


def asd_from_rfft(noise_fft, sample, rate):
    noise_fft = np.asarray(noise_fft)

    df = rate / sample

    amp_dens = np.zeros(
        len(noise_fft),
        dtype=float,
    )

    if len(noise_fft) == 0:
        return amp_dens

    amp_dens[0] = (
        np.abs(noise_fft[0])
        / (
            sample
            * np.sqrt(df)
        )
    )

    if len(noise_fft) > 2:
        amp_dens[1:-1] = (
            np.sqrt(2.0)
            * np.abs(noise_fft[1:-1])
            / (
                sample
                * np.sqrt(df)
            )
        )

    if len(noise_fft) > 1:
        if sample % 2 == 0:
            amp_dens[-1] = (
                np.abs(noise_fft[-1])
                / (
                    sample
                    * np.sqrt(df)
                )
            )
        else:
            amp_dens[-1] = (
                np.sqrt(2.0)
                * np.abs(noise_fft[-1])
                / (
                    sample
                    * np.sqrt(df)
                )
            )

    return amp_dens


def make_noise_time_from_asd(
    noise_spe_dens,
    sample,
    rate,
    rng=None,
):
    noise_spe_dens = np.asarray(
        noise_spe_dens,
        dtype=float,
    )

    return generate_noise_from_asd(
        noise_spe_dens,
        sample,
        rate,
        rng=rng,
    )


def MakePulse():
    with open(
        f"{output}/input.json",
        "r",
    ) as f:
        para = json.load(f)

    # Generate a pulse at every one-based absorber block.
    positions = list(
        range(
            1,
            int(para["n_abs"]) + 1,
        )
    )

    input_path = (
        f"{output}/input.json"
    )

    posi2pulse(
        input_path,
        positions,
        output_path=f"{output}/pulses.h5",
    )

    settling_time = int(para.get("SettlingTime", 0))
    if settling_time == 0:
        # A zero value means that the settling window has not been selected
        # yet.  Use the weakest CH1 reference pulse across all absorber
        # positions, then take the sample where that waveform peaks.
        weakest = None
        for position, pulse in iter_hdf5_pulse_items(f"{output}/pulses.h5"):
            waveform = np.asarray(pulse["ch1"], dtype=float)
            if waveform.size == 0 or not np.all(np.isfinite(waveform)):
                continue
            peak = float(np.max(waveform))
            if weakest is None or peak < weakest[0]:
                weakest = (peak, int(position), int(np.argmax(waveform)))

        if weakest is None:
            raise ValueError("Could not determine SettlingTime from CH1 reference pulses")

        _, settling_position, settling_time = weakest
        print(
            f"SettlingTime was 0; selected the weakest CH1 waveform at "
            f"position {settling_position}: sample {settling_time}"
        )

    para["SettlingTime"] = settling_time

    with open(
        f"{output}/input.json",
        "w",
    ) as f:
        json.dump(
            para,
            f,
            indent=4,
        )


def LoadPulses():
    """Load the shared-time, position-keyed HDF5 pulse schema."""

    return read_all_pulses(
        f"{output}/pulses.h5"
    )


def _event_data_root():
    """Return the root containing external dumpall.dat/event.h5 data."""

    if event_root is None:
        return Path(output)

    return Path(event_root).expanduser()


def IterJsonObjectItems(
    path,
    chunk_size=1024 * 1024,
):
    """Yield top-level JSON object items without loading the entire file."""

    decoder = json.JSONDecoder()

    with open(
        path,
        "r",
        encoding="utf-8",
    ) as file:

        buffer = ""

        def read_more():
            nonlocal buffer

            chunk = file.read(
                chunk_size
            )

            if not chunk:
                raise ValueError(
                    f"unexpected end of JSON input: {path}"
                )

            buffer += chunk

        def parse_value():
            nonlocal buffer

            while True:
                buffer = buffer.lstrip()

                try:
                    value, end = (
                        decoder.raw_decode(
                            buffer
                        )
                    )

                except json.JSONDecodeError:
                    read_more()
                    continue

                buffer = buffer[end:]

                return value

        while not buffer:
            read_more()

        buffer = buffer.lstrip()

        if not buffer.startswith("{"):
            raise ValueError(
                f"expected a JSON object at the root of: {path}"
            )

        buffer = buffer[1:]

        while True:

            while True:
                buffer = buffer.lstrip()

                if buffer:
                    break

                read_more()

            if buffer.startswith("}"):
                return

            if buffer.startswith(","):
                buffer = buffer[1:]
                continue

            key = parse_value()

            if not isinstance(
                key,
                str,
            ):
                raise ValueError(
                    f"expected a string key in: {path}"
                )

            while True:
                buffer = buffer.lstrip()

                if buffer:
                    break

                read_more()

            if not buffer.startswith(":"):
                raise ValueError(
                    f"expected ':' after key {key!r} in: {path}"
                )

            buffer = buffer[1:]

            yield key, parse_value()


def IterPulseItems(
    path,
    chunk_size=1024 * 1024,
):
    """Yield pulse records from HDF5 (or legacy JSON during migration)."""

    if Path(path).suffix.lower() in {
        ".h5",
        ".hdf5",
    }:
        yield from iter_hdf5_pulse_items(
            path
        )

        return

    decoder = json.JSONDecoder()

    with open(
        path,
        "r",
        encoding="utf-8",
    ) as file:

        buffer = ""

        def read_more():
            nonlocal buffer

            chunk = file.read(
                chunk_size
            )

            if not chunk:
                raise ValueError(
                    f"unexpected end of JSON input: {path}"
                )

            buffer += chunk

        def parse_value():
            nonlocal buffer

            while True:
                buffer = buffer.lstrip()

                try:
                    value, end = (
                        decoder.raw_decode(
                            buffer
                        )
                    )

                except json.JSONDecodeError:
                    read_more()
                    continue

                buffer = buffer[end:]

                return value

        def next_non_whitespace():
            nonlocal buffer

            while True:
                buffer = buffer.lstrip()

                if buffer:
                    return

                read_more()

        next_non_whitespace()

        if not buffer.startswith("{"):
            raise ValueError(
                f"expected a JSON object at the root of: {path}"
            )

        buffer = buffer[1:]

        while True:
            next_non_whitespace()

            if buffer.startswith("}"):
                return

            if buffer.startswith(","):
                buffer = buffer[1:]
                continue

            key = parse_value()

            if not isinstance(
                key,
                str,
            ):
                raise ValueError(
                    f"expected a string key in: {path}"
                )

            next_non_whitespace()

            if not buffer.startswith(":"):
                raise ValueError(
                    f"expected ':' after key {key!r} in: {path}"
                )

            buffer = buffer[1:]

            if key != "pulses":
                parse_value()
                continue

            next_non_whitespace()

            if not buffer.startswith("{"):
                raise ValueError(
                    f"expected pulses to be an object in: {path}"
                )

            buffer = buffer[1:]

            while True:
                next_non_whitespace()

                if buffer.startswith("}"):
                    buffer = buffer[1:]
                    break

                if buffer.startswith(","):
                    buffer = buffer[1:]
                    continue

                pulse_id = parse_value()

                if not isinstance(
                    pulse_id,
                    str,
                ):
                    raise ValueError(
                        f"expected a pulse ID key in: {path}"
                    )

                next_non_whitespace()

                if not buffer.startswith(":"):
                    raise ValueError(
                        f"expected ':' after pulse ID {pulse_id!r} in: {path}"
                    )

                buffer = buffer[1:]

                yield (
                    pulse_id,
                    parse_value(),
                )

            return


def FitRatios():
    with open(
        f"{output}/input.json",
        "r",
        encoding="utf-8",
    ) as f:
        para = json.load(f)

    pulses_by_position = (
        LoadPulses()
    )

    center_position = (
        int(para["n_abs"]) + 1
    ) / 2

    pitch = (
        float(para["length"])
        / int(para["n_abs"])
    )

    rows = []

    for (
        position,
        pulse,
    ) in pulses_by_position.items():

        ch0_max = max(
            pulse["ch0"]
        )

        ch1_max = max(
            pulse["ch1"]
        )

        if ch0_max == 0:
            raise ValueError(
                f"CH0 maximum is zero at position {pulse['position']}"
            )

        rows.append(
            {
                "position": position,
                "x_mm": (
                    position
                    - center_position
                )
                * pitch,
                "ch1_ch0_max_ratio": (
                    ch1_max
                    / ch0_max
                ),
            }
        )

    pd.DataFrame(
        rows
    ).to_csv(
        f"{output}/ratios.csv",
        index=False,
    )


def MakeNoise():
    """
    Calculate detector noise with a reduced five-state thermal model.

    State vector:
        0 : I_TES1
        1 : T_TES1
        2 : T_abs
        3 : T_TES2
        4 : I_TES2

    The absorber node represents the absorber center.

    The thermal path from either TES to the absorber center contains:
        1) TES-Pb interface thermal resistance
        2) one half of the Pb end-to-end thermal resistance

    Therefore:

        1 / G_eff
            =
        1 / G_abs_tes
            +
        1 / (2 * G_abs_abs)

    and the reduced thermal network is:

        TES1 -- G_eff -- Absorber -- G_eff -- TES2

    There is no independent absorber-absorber TFN source in this
    one-absorber-node representation.

    Complex source-to-output transfer functions for both TESs are retained
    so that correlated two-channel detector noise can be generated.
    """

    with open(
        f"{output}/input.json",
        "r",
    ) as f:
        para = json.load(f)

    C_abs = para[
        "C_abs"
    ]

    C_tes = para[
        "C_tes"
    ]

    G_abs_abs = float(
        para["G_abs-abs"]
    )

    G_abs_tes = para[
        "G_abs-tes"
    ]

    G_tes_bath = para[
        "G_tes-bath"
    ]

    R = para[
        "R"
    ]

    R_l = para[
        "R_l"
    ]

    T_c = para[
        "T_c"
    ]

    T_bath = para[
        "T_bath"
    ]

    a = para[
        "alpha"
    ]

    b = para[
        "beta"
    ]

    excess_johnson_M = resolve_excess_johnson_M(para)

    L = para[
        "L"
    ]

    n = para[
        "n"
    ]

    rate = int(
        para["rate"]
    )

    samples = int(
        para["samples"]
    )

    frequency = np.arange(
        0,
        rate,
        rate / samples,
    )

    # --------------------------------------------------------
    # TES operating point
    # --------------------------------------------------------

    I = np.sqrt(
        (
            G_tes_bath
            * T_c
            * (
                1
                - (
                    T_bath
                    / T_c
                ) ** n
            )
        )
        / (
            n
            * R
        )
    )

    t_el = (
        L
        / (
            R_l
            + R
            * (
                1
                + b
            )
        )
    )

    L_I = (
        a
        * I**2
        * R
        / (
            G_tes_bath
            * T_c
        )
    )

    t_I = (
        C_tes
        / (
            (
                1
                - L_I
            )
            * G_tes_bath
        )
    )

    # --------------------------------------------------------
    # Effective TES-to-absorber-center conductance
    # --------------------------------------------------------

    G_eff = 1.0 / (
        1.0 / G_abs_tes
        + 1.0 / (
            2.0
            * G_abs_abs
        )
    )

    print(
        "G_abs_tes =",
        G_abs_tes,
    )

    print(
        "G_abs_abs =",
        G_abs_abs,
    )

    print(
        "G_eff     =",
        G_eff,
    )

    # --------------------------------------------------------
    # Thermal fluctuation noise amplitudes
    #
    # ptfn_Flink is intentionally left unchanged here.
    # --------------------------------------------------------

    ptfn_tes_bath = np.sqrt(
        4
        * k_b
        * T_c**2
        * G_tes_bath
        * ptfn_Flink
    )

    ptfn_eff = np.sqrt(
        4
        * k_b
        * T_c**2
        * G_eff
        * ptfn_Flink
    )

    # --------------------------------------------------------
    # Johnson noise amplitudes
    # --------------------------------------------------------

    enj = tes_johnson_voltage_asd(
        T_c,
        R,
        b,
        excess_johnson_M,
    )

    print(f"TES excess Johnson M = {excess_johnson_M:g}")

    enj_R = np.sqrt(
        4
        * k_b
        * T_bath
        * R_l
    )

    source_names = [
        "johnson_tes1",
        "johnson_load1",
        "phonon_tes1_bath",
        "phonon_tes1_absorber_effective",
        "phonon_tes2_absorber_effective",
        "phonon_tes2_bath",
        "johnson_load2",
        "johnson_tes2",
    ]

    # --------------------------------------------------------
    # Noise source matrix
    #
    # Every column is one independent physical noise source.
    # --------------------------------------------------------

    def matrix_N():

        X = np.zeros(
            (
                5,
                len(source_names),
            ),
            dtype=np.complex128,
        )

        # TES1 Johnson noise
        X[0, 0] = (
            -enj
            / L
        )

        X[1, 0] = (
            I
            * enj
            / C_tes
        )

        # Load1 Johnson noise
        X[0, 1] = (
            enj_R
            / L
        )

        # TES1-bath TFN
        X[1, 2] = (
            ptfn_tes_bath
            / C_tes
        )

        # TES1-absorber effective TFN
        #
        # Same physical heat-current fluctuation is injected
        # with opposite signs into the two connected nodes.
        X[1, 3] = (
            +ptfn_eff
            / C_tes
        )

        X[2, 3] = (
            -ptfn_eff
            / C_abs
        )

        # Absorber-TES2 effective TFN
        X[2, 4] = (
            -ptfn_eff
            / C_abs
        )

        X[3, 4] = (
            +ptfn_eff
            / C_tes
        )

        # TES2-bath TFN
        X[3, 5] = (
            ptfn_tes_bath
            / C_tes
        )

        # Load2 Johnson noise
        X[4, 6] = (
            enj_R
            / L
        )

        # TES2 Johnson noise
        X[4, 7] = (
            -enj
            / L
        )

        X[3, 7] = (
            I
            * enj
            / C_tes
        )

        return X

    # --------------------------------------------------------
    # Linearized electrothermal matrix
    # --------------------------------------------------------

    def matrix_M(
        omega,
    ):

        X = np.zeros(
            (
                5,
                5,
            ),
            dtype=np.complex128,
        )

        # TES1 electrical
        X[0, 0] = (
            1
            / t_el
            + 1j
            * omega
        )

        X[0, 1] = (
            L_I
            * G_tes_bath
            / (
                I
                * L
            )
        )

        # TES1 thermal
        X[1, 0] = (
            -I
            * R
            * (
                2
                + b
            )
            / C_tes
        )

        X[1, 1] = (
            1
            / t_I
            + G_eff
            / C_tes
            + 1j
            * omega
        )

        X[1, 2] = (
            -G_eff
            / C_tes
        )

        # Absorber thermal
        X[2, 1] = (
            -G_eff
            / C_abs
        )

        X[2, 2] = (
            2
            * G_eff
            / C_abs
            + 1j
            * omega
        )

        X[2, 3] = (
            -G_eff
            / C_abs
        )

        # TES2 thermal
        X[3, 2] = (
            -G_eff
            / C_tes
        )

        X[3, 3] = (
            1
            / t_I
            + G_eff
            / C_tes
            + 1j
            * omega
        )

        X[3, 4] = (
            -I
            * R
            * (
                2
                + b
            )
            / C_tes
        )

        # TES2 electrical
        X[4, 3] = (
            L_I
            * G_tes_bath
            / (
                I
                * L
            )
        )

        X[4, 4] = (
            1
            / t_el
            + 1j
            * omega
        )

        return X

    # --------------------------------------------------------
    # Transfer functions
    # --------------------------------------------------------

    N = matrix_N()

    omega_array = (
        frequency
        * 2
        * math.pi
    )

    transfer_ch0 = []
    transfer_ch1 = []

    for omg in tqdm.tqdm(
        omega_array
    ):

        H = np.linalg.solve(
            matrix_M(omg),
            N,
        )

        transfer_ch0.append(
            H[0, :]
        )

        transfer_ch1.append(
            H[4, :]
        )

    # Shape:
    #   frequency x source
    transfer_ch0 = np.asarray(
        transfer_ch0
    )

    transfer_ch1 = np.asarray(
        transfer_ch1
    )

    # --------------------------------------------------------
    # Individual source ASD
    # --------------------------------------------------------

    asd_ch0_components = np.abs(
        transfer_ch0
    )

    asd_ch1_components = np.abs(
        transfer_ch1
    )

    # --------------------------------------------------------
    # Total ASD
    #
    # Independent physical noise sources add in PSD:
    #
    #   ASD_total = sqrt(sum(ASD_i^2))
    #
    # --------------------------------------------------------

    noise_total_ch0 = np.sqrt(
        np.sum(
            asd_ch0_components**2,
            axis=1,
        )
    )

    noise_total_ch1 = np.sqrt(
        np.sum(
            asd_ch1_components**2,
            axis=1,
        )
    )

    # --------------------------------------------------------
    # CH0-CH1 cross PSD
    # --------------------------------------------------------

    cross_psd = np.sum(
        transfer_ch0
        * np.conjugate(
            transfer_ch1
        ),
        axis=1,
    )

    psd_ch0 = (
        noise_total_ch0**2
    )

    psd_ch1 = (
        noise_total_ch1**2
    )

    # --------------------------------------------------------
    # Sum / difference channel ASD
    # --------------------------------------------------------

    noise_sum = np.sqrt(
        np.maximum(
            psd_ch0
            + psd_ch1
            + 2.0
            * np.real(
                cross_psd
            ),
            0.0,
        )
    )

    noise_diff = np.sqrt(
        np.maximum(
            psd_ch0
            + psd_ch1
            - 2.0
            * np.real(
                cross_psd
            ),
            0.0,
        )
    )

    # --------------------------------------------------------
    # Save noise model
    # --------------------------------------------------------

    import h5py

    with h5py.File(
        f"{output}/noise.h5",
        "w",
    ) as f:

        f.attrs[
            "input_json"
        ] = json.dumps(
            para,
            separators=(
                ",",
                ":",
            ),
        )

        f.attrs[
            "noise_model"
        ] = (
            "five_state_effective_conductance"
        )

        f.attrs["excess_johnson_M"] = excess_johnson_M
        f.attrs["johnson_model"] = "standard_tes_johnson_with_excess_factor"
        f.attrs["johnson_formula"] = (
            "S_V_J = 4*k_B*T*R*(1+2*beta)*(1+M^2)"
        )

        f.attrs[
            "G_eff_W_per_K"
        ] = G_eff

        f.create_dataset(
            "frequency",
            data=frequency,
            compression="gzip",
        )

        # Compatibility with existing LoadNoise()
        f.create_dataset(
            "total",
            data=noise_total_ch0,
            compression="gzip",
        )

        f.create_dataset(
            "total_ch0",
            data=noise_total_ch0,
            compression="gzip",
        )

        f.create_dataset(
            "total_ch1",
            data=noise_total_ch1,
            compression="gzip",
        )

        f.create_dataset(
            "sum_asd",
            data=noise_sum,
            compression="gzip",
        )

        f.create_dataset(
            "diff_asd",
            data=noise_diff,
            compression="gzip",
        )

        f.create_dataset(
            "cross_psd",
            data=cross_psd,
            compression="gzip",
        )

        # Complex transfer functions for correlated
        # two-channel noise generation.
        f.create_dataset(
            "transfer_ch0",
            data=transfer_ch0,
            compression="gzip",
        )

        f.create_dataset(
            "transfer_ch1",
            data=transfer_ch1,
            compression="gzip",
        )

        group0 = f.create_group(
            "components_ch0"
        )

        group1 = f.create_group(
            "components_ch1"
        )

        for (
            idx,
            name,
        ) in enumerate(
            source_names
        ):

            group0.create_dataset(
                name,
                data=(
                    asd_ch0_components[
                        :,
                        idx,
                    ]
                ),
                compression="gzip",
            )

            group1.create_dataset(
                name,
                data=(
                    asd_ch1_components[
                        :,
                        idx,
                    ]
                ),
                compression="gzip",
            )

    # --------------------------------------------------------
    # Plot CH0 noise
    # --------------------------------------------------------

    plt.figure(
        figsize=(8, 8)
    )

    for (
        idx,
        name,
    ) in enumerate(
        source_names
    ):

        plt.plot(
            frequency,
            asd_ch0_components[
                :,
                idx,
            ],
            linewidth=1.5,
            label=name,
        )

    plt.plot(
        frequency,
        noise_total_ch0,
        linewidth=3,
        label="Total Noise (CH0)",
    )

    plt.xlabel(
        "Frequency [Hz]",
        fontsize=20,
    )

    plt.ylabel(
        "Noise Spectral Density [A/rtHz]",
        fontsize=20,
    )

    plt.ylim(
        1e-13,
        1e-9,
    )

    plt.xlim(
        1,
        1e6,
    )

    plt.loglog()
    plt.grid()

    plt.legend(
        loc="best",
        fancybox=True,
        fontsize=9,
        ncol=2,
    )

    plt.tight_layout()

    #plt.savefig(
    #    f"{output}/noise_all.png",
    #    dpi=700,
    #)

    plt.cla()


def LoadNoise():
    """Read the CH0 total noise spectral density from noise.h5."""

    import h5py

    with h5py.File(
        f"{output}/noise.h5",
        "r",
    ) as f:
        return f[
            "total"
        ][:]


def _simulation_pre_analysis_asd(frequency, rate, para):
    """Build the ASD entering the finite-record analysis."""
    import h5py

    with h5py.File(f"{output}/noise.h5", "r") as file:
        source_frequency = file["frequency"][:]
        detector_asd = file["total"][:]

    hardware_order = int(para.get("hardware_bessel_order", 4))
    hardware_main = general.AnalogBesselMagnitudeResponse(
        frequency,
        100000.0,
        order=hardware_order,
    )
    detector_main = (
        np.interp(frequency, source_frequency, detector_asd)
        * hardware_main
    )

    # The analog 100 kHz Bessel is before the ADC, so the part above Nyquist
    # aliases into the measured band.
    alias_frequency = rate - frequency
    hardware_alias = general.AnalogBesselMagnitudeResponse(
        alias_frequency,
        100000.0,
        order=hardware_order,
    )
    detector_alias = (
        np.interp(alias_frequency, source_frequency, detector_asd)
        * hardware_alias
    )

    white_asd = float(
        para.get(
            "post_filter_white_asd_A_rtHz",
            para.get("readout_white_asd_A_rtHz", 0.0),
        )
    )
    if white_asd < 0:
        raise ValueError("post/readout white ASD must be non-negative")

    return np.sqrt(
        detector_main**2
        + detector_alias**2
        + white_asd**2
    )


def finite_record_simulation_spectrum(
    input_asd,
    sample,
    rate,
    cutoff,
    records=FINITE_RECORDS,
    seed=FINITE_RECORD_SEED,
):
    """Estimate the post-analysis ASD using finite records."""
    frequency = np.fft.rfftfreq(sample, d=1 / rate)
    input_asd = np.asarray(input_asd, dtype=float)
    if input_asd.shape != frequency.shape:
        raise ValueError("input ASD must have the one-sided rFFT length")

    df = rate / sample
    window = np.hanning(sample)
    window_power_gain = np.sqrt(np.mean(window**2))
    numerator, denominator = scipy.signal.bessel(
        2,
        cutoff / (rate / 2),
        "low",
    )
    rng = np.random.default_rng(seed)
    # Average power spectra before taking the square root.  Averaging the
    # magnitudes directly estimates the mean of a Rayleigh-distributed
    # random variable, rather than the ASD, and also misses the rFFT
    # normalization/one-sided factors handled by ``asd_from_rfft``.
    power_sum = np.zeros(len(frequency))

    for _ in range(records):
        spectrum = np.zeros(len(frequency), dtype=np.complex128)
        sigma = input_asd[1:-1] * sample * np.sqrt(df) / 2.0
        spectrum[1:-1] = (
            rng.normal(size=len(sigma)) * sigma
            + 1j * rng.normal(size=len(sigma)) * sigma
        )
        spectrum[0] = rng.normal() * input_asd[0] * sample * np.sqrt(df)
        spectrum[-1] = rng.normal() * input_asd[-1] * sample * np.sqrt(df)

        noise = np.fft.irfft(spectrum, n=sample)
        noise = noise - np.mean(noise)
        noise = scipy.signal.filtfilt(numerator, denominator, noise)
        noise_fft = np.fft.rfft(noise * window) / window_power_gain
        power_sum += np.abs(noise_fft) ** 2

    mean_amplitude = np.sqrt(power_sum / records)
    return asd_from_rfft(mean_amplitude, sample, rate)


def LoadNoiseTransfers():
    """Read complex source-to-CH0/CH1 transfer functions."""

    import h5py

    with h5py.File(
        f"{output}/noise.h5",
        "r",
    ) as f:

        transfer_ch0 = f[
            "transfer_ch0"
        ][:]

        transfer_ch1 = f[
            "transfer_ch1"
        ][:]

    return (
        transfer_ch0,
        transfer_ch1,
    )


def generate_correlated_noise_pair_from_transfer(
    transfer_ch0,
    transfer_ch1,
    sample,
    rate,
    rng=None,
):
    """
    Generate one correlated CH0/CH1 detector-noise realization.

    Each source column receives one independent random Fourier realization.
    The same source realization is propagated through both TES transfer
    functions, preserving the calculated CH0-CH1 cross spectrum.
    """

    if rng is None:
        rng = np.random.default_rng()

    n_rfft = (
        sample
        // 2
        + 1
    )

    H0 = np.asarray(
        transfer_ch0[
            :n_rfft
        ],
        dtype=np.complex128,
    )

    H1 = np.asarray(
        transfer_ch1[
            :n_rfft
        ],
        dtype=np.complex128,
    )

    if H0.shape != H1.shape:
        raise ValueError(
            "CH0 and CH1 transfer-function shapes do not match"
        )

    df = (
        rate
        / sample
    )

    n_source = (
        H0.shape[1]
    )

    z = np.zeros(
        (
            n_rfft,
            n_source,
        ),
        dtype=np.complex128,
    )

    stop = (
        n_rfft - 1
        if sample % 2 == 0
        else n_rfft
    )

    if stop > 1:

        phase = rng.uniform(
            0.0,
            2.0 * np.pi,
            size=(
                stop - 1,
                n_source,
            ),
        )

        z[
            1:stop
        ] = (
            sample
            * np.sqrt(
                df
                / 2.0
            )
            * np.exp(
                1j
                * phase
            )
        )

    spectrum0 = np.sum(
        H0 * z,
        axis=1,
    )

    spectrum1 = np.sum(
        H1 * z,
        axis=1,
    )

    # Remove random DC component.
    spectrum0[0] = 0.0
    spectrum1[0] = 0.0

    # Nyquist component must be purely real.
    # Here it is simply set to zero.
    if sample % 2 == 0:
        spectrum0[-1] = 0.0
        spectrum1[-1] = 0.0

    noise0 = np.fft.irfft(
        spectrum0,
        n=sample,
    )

    noise1 = np.fft.irfft(
        spectrum1,
        n=sample,
    )

    return (
        noise0,
        noise1,
    )


def _ShowNoiseSpectrum():

    with open(
        f"{output}/input.json",
        "r",
    ) as f:
        para = json.load(f)

    sample = int(
        para["samples"]
    )

    rate = para[
        "rate"
    ]

    freq = np.fft.rfftfreq(
        sample,
        d=1 / rate,
    )
    input_asd = _simulation_pre_analysis_asd(freq, rate, para)
    amp_dens = (
        finite_record_simulation_spectrum(
            input_asd,
            sample,
            rate,
            para["cutoff"],
        )
        * 1e6
    )

    plt.plot(
        freq,
        amp_dens,
    )

    plt.xlabel(
        "Frequency [Hz]",
        fontsize=20,
    )

    plt.ylabel(
        "Amplitude [uA/rtHz]",
        fontsize=20,
    )

    plt.loglog()

    plt.savefig(
        f"{output}/noise_total-bessel100k.png",
        dpi=350,
    )

    plt.clf()

    np.savetxt(
        f"{output}/noise_total-bessel100k.dat",
        amp_dens,
    )

def ShowSamples():

    with open(
        f"{output}/input.json",
        "r",
    ) as f:
        para = json.load(f)

    pulses_by_position = (
        LoadPulses()
    )

    def pulse_channel(
        position,
        channel,
    ):

        try:
            values = (
                pulses_by_position[
                    int(position)
                ][
                    f"ch{channel}"
                ]
            )

        except KeyError as error:
            raise ValueError(
                f"pulses.h5 does not contain position {position}"
            ) from error

        return np.asarray(
            values,
            dtype=float,
        )

    sample = int(
        para["samples"]
    )

    rate = para[
        "rate"
    ]

    # Simulated noise frequency domain
    noise_spe_dens = (
        LoadNoise()
    )

    noise_time = (
        make_noise_time_from_asd(
            noise_spe_dens,
            sample,
            rate,
        )
    )

    print(
        f"spe len:{len(noise_spe_dens)} "
        f"time len:{len(noise_time)}"
    )

    time = np.linspace(
        0,
        para["samples"]
        / rate,
        int(
            para["samples"]
        ),
    )

    frequency_pulse = np.arange(
        0,
        rate,
        rate
        / int(
            para["samples"]
        ),
    )

    pulse_cmap_name = (
        "viridis"
    )

    # --------------------------------------------------------
    # Pulse time domain
    # --------------------------------------------------------

    fig, ax = plt.subplots()

    cmap = plt.get_cmap(
        pulse_cmap_name
    )

    norm = Normalize(
        vmin=0.0,
        vmax=float(
            para["length"]
        ),
    )

    sm = cm.ScalarMappable(
        norm=norm,
        cmap=cmap,
    )

    sm.set_array([])

    for i in para[
        "position"
    ]:

        data = pulse_channel(
            i,
            0,
        )

        distance = (
            float(i)
            / float(
                para["n_abs"]
            )
            * float(
                para["length"]
            )
        )

        ax.plot(
            time * 1e3,
            data * 1e6,
            color=cmap(
                norm(
                    distance
                )
            ),
            linewidth=1.5,
        )

    ax.set_xlabel(
        "Time [ms]",
        fontsize=20,
    )

    ax.set_ylabel(
        "Current [uA]",
        fontsize=20,
    )

    plt.xlim(
        -0.1,
        5,
    )

    ax.grid()

    fig.tight_layout()

    cbar = fig.colorbar(
        sm,
        ax=ax,
        ticks=np.arange(
            0,
            float(
                para["length"]
            )
            + 0.1,
            2,
        ),
    )

    cbar.set_label(
        "Distance [mm]",
        fontsize=20,
    )

    cbar.ax.tick_params(
        labelsize=12
    )

    cbar.ax.invert_yaxis()

    #fig.savefig(
    #    f"{output}/checkpulse_pulses.png",
    #    dpi=350,
    #    transparent=True,
    #)

    plt.close(
        fig
    )

    # --------------------------------------------------------
    # Noise time domain
    # --------------------------------------------------------

    plt.plot(
        time * 1e3,
        noise_time * 1e6,
        linewidth=1.5,
    )

    plt.xlabel(
        "Time [ms]",
        fontsize=20,
    )

    plt.ylabel(
        "Current [uA]",
        fontsize=20,
    )

    plt.grid()
    plt.tight_layout()

    #plt.savefig(
    #    f"{output}/checkpulse_noise.png",
    #    dpi=350,
    #    transparent=True,
    #)

    plt.clf()

    # --------------------------------------------------------
    # Pulse with noise
    # --------------------------------------------------------

    fig, ax = plt.subplots()

    cmap = plt.get_cmap(
        pulse_cmap_name
    )

    norm = Normalize(
        vmin=0.0,
        vmax=float(
            para["length"]
        ),
    )

    sm = cm.ScalarMappable(
        norm=norm,
        cmap=cmap,
    )

    sm.set_array([])

    for i in para[
        "position"
    ]:

        data = pulse_channel(
            i,
            0,
        )

        distance = (
            float(i)
            / float(
                para["n_abs"]
            )
            * float(
                para["length"]
            )
        )

        data += (
            make_noise_time_from_asd(
                noise_spe_dens,
                sample,
                rate,
            )
        )

        ax.plot(
            time * 1e3,
            data * 1e6,
            color=cmap(
                norm(
                    distance
                )
            ),
            linewidth=1.5,
        )

    ax.set_xlabel(
        "Time [ms]",
        fontsize=20,
    )

    ax.set_ylabel(
        "Current [uA]",
        fontsize=20,
    )

    ax.set_xlim(
        -0.1,
        5,
    )

    ax.grid()
    fig.tight_layout()

    cbar = fig.colorbar(
        sm,
        ax=ax,
        ticks=np.arange(
            0,
            float(
                para["length"]
            )
            + 0.1,
            2,
        ),
    )

    cbar.set_label(
        "Distance [mm]",
        fontsize=20,
    )

    cbar.ax.tick_params(
        labelsize=12
    )

    cbar.ax.invert_yaxis()

    #fig.savefig(
    #    f"{output}/checkpulse_pulses_with_noise.png",
    #    dpi=350,
    #    transparent=True,
    #)

    plt.close(
        fig
    )

    os.makedirs(
        f"{output}/checkpulse_sn_ratio",
        exist_ok=True,
    )

    # --------------------------------------------------------
    # Position ratio
    # --------------------------------------------------------

    data = np.loadtxt(
        f"{output}/ratios.csv",
        delimiter=",",
        skiprows=1,
    )

    x = data[:, 1]
    y = data[:, 2]

    plt.scatter(
        x,
        y,
        c=x,
        cmap="coolwarm",
        s=50,
    )

    plt.xlabel(
        "Position[mm]"
    )

    plt.ylabel(
        "CH1/CH0[-]"
    )

    plt.tight_layout()
    plt.grid(True)

    #plt.savefig(
    #    f"{output}/checkpulse_max_ratio.png",
    #    dpi=350,
    #)

    plt.cla()

    # --------------------------------------------------------
    # S/N ratio
    # --------------------------------------------------------

    for i in para[
        "position"
    ]:

        pulse = pulse_channel(
            i,
            0,
        )

        pulse_noise = (
            pulse
            + make_noise_time_from_asd(
                noise_spe_dens,
                sample,
                rate,
            )
        )

        pulse_rfft = np.fft.rfft(
            pulse
        )

        pulse_noise_rfft = (
            np.fft.rfft(
                pulse_noise
            )
        )

        plt.plot(
            frequency_pulse[
                : len(
                    pulse_noise_rfft
                )
            ],
            np.abs(
                pulse_noise_rfft
            ),
            label="posi+noise",
            linestyle="--",
        )

        plt.plot(
            frequency_pulse[
                : len(
                    pulse_rfft
                )
            ],
            np.abs(
                pulse_rfft
            ),
            label="pulse",
        )

        plt.plot(
            frequency_pulse,
            noise_spe_dens,
            color="black",
            label="noise",
        )

        plt.xlim(
            1,
            5e4,
        )

        plt.loglog()
        plt.legend()

        plt.xlabel(
            "Frequency [Hz]"
        )

        plt.ylabel(
            "Amplitude [A/rtHz]"
        )

        #plt.savefig(
        #    f"{output}/checkpulse_sn_ratio/"
        #    f"sn_ratio_position_{i}.png",
        #    dpi=350,
        #)

        plt.cla()

    _ShowNoiseSpectrum()


def _legacy_MultiPulse():

    pulse_num = 300

    with open(
        f"{output}/input.json",
        "r",
    ) as f:
        para = json.load(f)

    sample = int(
        para["samples"]
    )

    rate = para[
        "rate"
    ]

    def AddPulse(
        noise_spe_dens,
        data,
    ):

        noise_time = (
            make_noise_time_from_asd(
                noise_spe_dens,
                sample,
                rate,
            )
        )

        return (
            data
            + noise_time[
                : len(data)
            ]
        )

    def Process(
        pulse,
        output,
        noise_spe_dens,
        ch,
        posi,
        k,
        para,
    ):

        noised_pulse = AddPulse(
            noise_spe_dens,
            pulse,
        )

        np.savetxt(
            f'{output}/{para["E"]}keV_{posi}/'
            f'pulse_noise/CH{ch}/CH{ch}_{k}.dat',
            noised_pulse,
        )

    noise_spe_dens = (
        LoadNoise()
    )

    for posi in tqdm.tqdm(
        para["position"]
    ):

        noise_path = (
            f'{output}/{para["E"]}keV_{posi}/pulse_noise'
        )

        if os.path.exists(
            noise_path
        ):
            shutil.rmtree(
                noise_path
            )

        for ch in [
            0,
            1,
        ]:

            os.makedirs(
                f'{output}/{para["E"]}keV_{posi}/'
                f'pulse_noise/CH{ch}',
                exist_ok=True,
            )

            pulse = np.loadtxt(
                f'{output}/{para["E"]}keV_{posi}/'
                f'pulse/CH{ch}/CH{ch}_1.dat'
            )

            with concurrent.futures.ThreadPoolExecutor() as executor:

                futures = [
                    executor.submit(
                        Process,
                        pulse,
                        output,
                        noise_spe_dens,
                        ch,
                        posi,
                        k,
                        para,
                    )
                    for k in range(
                        pulse_num
                    )
                ]

                for future in futures:
                    try:
                        future.result()

                    except Exception as e:
                        print(
                            f"error:{e}"
                        )


def Pulse_Noise():
    """
    Create noisy single-site pulses with correlated CH0/CH1 detector noise.
    """

    pulse_num = 500

    with open(
        f"{output}/input.json",
        "r",
        encoding="utf-8",
    ) as f:
        para = json.load(f)

    pulses_by_position = (
        LoadPulses()
    )

    sample = int(
        para["samples"]
    )

    rate = para[
        "rate"
    ]

    (
        transfer_ch0,
        transfer_ch1,
    ) = LoadNoiseTransfers()

    rng = (
        np.random.default_rng()
    )

    for posi in tqdm.tqdm(
        para["position"]
    ):

        try:
            pulse = (
                pulses_by_position[
                    int(posi)
                ]
            )

        except KeyError as error:
            raise ValueError(
                f"pulses.h5 does not contain position {posi}"
            ) from error

        ch0 = np.asarray(
            pulse["ch0"],
            dtype=float,
        )

        ch1 = np.asarray(
            pulse["ch1"],
            dtype=float,
        )

        position_directory = (
            f"{output}/{posi}"
        )

        os.makedirs(
            position_directory,
            exist_ok=True,
        )

        path = (
            f"{position_directory}/pulse_noise.h5"
        )

        with PulseWriter(
            path,
            pulse["time"],
            para,
        ) as f:

            for k in range(
                pulse_num
            ):

                (
                    noise0,
                    noise1,
                ) = generate_correlated_noise_pair_from_transfer(
                    transfer_ch0,
                    transfer_ch1,
                    sample,
                    rate,
                    rng=rng,
                )

                f.append(
                    k,
                    ch0 + noise0,
                    ch1 + noise1,
                )


def Dump2Event():
    """Convert external dumpall.dat files into event.h5 files.

    ``input.json`` supplies the input energy and position list, while the
    dump/event files themselves can live under ``event_root``. If no external
    root is configured, the existing ``output`` layout is used.
    """

    with open(
        f"{output}/input.json",
        "r",
        encoding="utf-8",
    ) as f:
        para = json.load(f)

    input_energy_mev = (
        float(
            para["E"]
        )
        / 1000.0
    )

    data_root = _event_data_root()

    for posi in tqdm.tqdm(
        para["position"],
        desc="dumpall to event",
    ):

        position_directory = data_root / str(posi)

        dump_path = (
            position_directory / "dumpall.dat"
        )

        event_path = (
            position_directory / "event.h5"
        )

        if not os.path.isfile(
            dump_path
        ):
            raise FileNotFoundError(
                f"dumpall.dat was not found: {dump_path}"
            )

        dump2event(
            dump_path,
            event_path,
            input_energy=input_energy_mev,
            save_all=True,
        )


def Pulse_Ms():
    """
    Synthesize one CH0/CH1 pulse for each event in every position folder.

    For each input.json position, read event.h5 from ``event_root`` and
    write pulse_MS.h5 under ``output``. event.h5 stores event IDs as its outer
    keys; all deposits in one event are summed into a single pulse.
    """

    with open(
        f"{output}/input.json",
        "r",
        encoding="utf-8",
    ) as f:
        para = json.load(f)

    reference_pulses = (
        LoadPulses()
    )

    n_abs = int(
        para["n_abs"]
    )

    length_mm = float(
        para["length"]
    )

    reference_energy_kev = float(
        para["E"]
    )

    event_data_root = _event_data_root()

    if (
        length_mm <= 0
        or reference_energy_kev <= 0
    ):
        raise ValueError(
            "input.json length and E must be positive"
        )

    if (
        set(
            range(
                1,
                n_abs + 1,
            )
        )
        - reference_pulses.keys()
    ):
        raise ValueError(
            "pulses.h5 does not contain all absorber positions"
        )

    def block_from_x_deposit(
        x_deposit_cm,
    ):

        x_mm = (
            float(
                x_deposit_cm
            )
            * 10.0
        )

        fraction = (
            x_mm
            + length_mm
            / 2.0
        ) / length_mm

        if (
            fraction < 0.0
            or fraction > 1.0
        ):
            return None

        return min(
            n_abs,
            int(
                fraction
                * n_abs
            )
            + 1,
        )

    time = (
        reference_pulses[
            1
        ][
            "time"
        ]
    )

    sample_count = len(
        time
    )

    for posi in tqdm.tqdm(
        para["position"],
        desc="MS pulse positions",
    ):

        position_directory = Path(output) / str(posi)
        event_directory = event_data_root / str(posi)

        position_directory.mkdir(parents=True, exist_ok=True)

        event_path = (
            event_directory / "event.h5"
        )

        output_path = (
            position_directory / "pulse_MS.h5"
        )

        if not os.path.isfile(
            event_path
        ):
            raise FileNotFoundError(
                f"event.h5 was not found: {event_path}"
            )

        with PulseWriter(
            output_path,
            time,
            para,
        ) as f:

            for (
                index,
                (
                    event_id,
                    event,
                ),
            ) in enumerate(
                tqdm.tqdm(
                    iter_hdf5_events(
                        event_path
                    ),
                    desc=f"MS pulses {posi}",
                    leave=False,
                )
            ):

                energy_by_block = np.zeros(
                    n_abs,
                    dtype=float,
                )

                if not isinstance(
                    event,
                    dict,
                ):
                    raise ValueError(
                        f"event.h5 event {event_id} is not an object"
                    )

                for particle in event.values():

                    positions = particle.get(
                        "x_deposit",
                        [],
                    )

                    energies = particle.get(
                        "E_deposit",
                        [],
                    )

                    if (
                        len(positions)
                        != len(energies)
                    ):
                        raise ValueError(
                            f"event.h5 event {event_id} has mismatched "
                            "x_deposit and E_deposit lengths"
                        )

                    for (
                        x_deposit,
                        e_deposit_mev,
                    ) in zip(
                        positions,
                        energies,
                    ):

                        block = block_from_x_deposit(
                            x_deposit
                        )

                        if block is not None:
                            energy_by_block[
                                block - 1
                            ] += float(
                                e_deposit_mev
                            )

                ch0 = np.zeros(
                    sample_count,
                    dtype=float,
                )

                ch1 = np.zeros(
                    sample_count,
                    dtype=float,
                )

                for (
                    block,
                    energy_mev,
                ) in enumerate(
                    energy_by_block,
                    start=1,
                ):

                    if energy_mev == 0.0:
                        continue

                    scale = (
                        energy_mev
                        * 1000.0
                        / reference_energy_kev
                    )

                    reference = (
                        reference_pulses[
                            block
                        ]
                    )

                    ch0 += (
                        scale
                        * np.asarray(
                            reference[
                                "ch0"
                            ],
                            dtype=float,
                        )
                    )

                    ch1 += (
                        scale
                        * np.asarray(
                            reference[
                                "ch1"
                            ],
                            dtype=float,
                        )
                    )

                f.append(
                    event_id,
                    ch0,
                    ch1,
                )


def _legacy_MS_Noise():

    with open(
        f"{output}/input.json",
        "r",
    ) as f:
        para = json.load(f)

    sample = int(
        para["samples"]
    )

    rate = para[
        "rate"
    ]

    def AddPulse(
        noise_spe_dens,
        data,
    ):

        noise_time = (
            make_noise_time_from_asd(
                noise_spe_dens,
                sample,
                rate,
            )
        )

        return (
            data
            + noise_time[
                : len(data)
            ]
        )

    noise_spe_dens = (
        LoadNoise()
    )

    def Process(
        file,
        output,
        noise_spe_dens,
        ch,
        num,
        posi,
        para,
    ):

        pulse = np.loadtxt(
            file
        )

        noised_pulse = AddPulse(
            noise_spe_dens,
            pulse,
        )

        np.savetxt(
            f'{output}/{para["E"]}keV_{posi}/'
            f'pulse_noise_ms_test/CH{ch}/CH{ch}_{num}.dat',
            noised_pulse,
        )

    for posi in tqdm.tqdm(
        para["position"]
    ):

        for ch in [
            0,
            1,
        ]:

            os.makedirs(
                f'{output}/{para["E"]}keV_{posi}/'
                f'pulse_noise_ms_test/CH{ch}',
                exist_ok=True,
            )

            file_pattern = (
                f'{output}/{para["E"]}keV_{posi}/'
                f'pulse_ms/CH{ch}/CH{ch}_*.dat'
            )

            file_list = glob.glob(
                file_pattern
            )

            numbers = []

            for file in file_list:

                match = re.search(
                    r"CH\d+_(\d+).dat",
                    file,
                )

                if match:
                    numbers.append(
                        int(
                            match.group(
                                1
                            )
                        )
                    )

            with concurrent.futures.ThreadPoolExecutor() as executor:

                futures = [
                    executor.submit(
                        Process,
                        file,
                        output,
                        noise_spe_dens,
                        ch,
                        number,
                        posi,
                        para,
                    )
                    for (
                        file,
                        number,
                    ) in zip(
                        file_list,
                        numbers,
                    )
                ]

                for future in futures:
                    future.result()


def Pulse_MS_Noise():
    """
    Add correlated detector noise to multi-site CH0/CH1 pulses.
    """

    with open(
        f"{output}/input.json",
        "r",
        encoding="utf-8",
    ) as f:
        para = json.load(f)

    sample = int(
        para["samples"]
    )

    rate = para[
        "rate"
    ]

    (
        transfer_ch0,
        transfer_ch1,
    ) = LoadNoiseTransfers()

    rng = (
        np.random.default_rng()
    )

    time = (
        LoadPulses()[
            1
        ][
            "time"
        ]
    )

    for posi in tqdm.tqdm(
        para["position"],
        desc="MS noise positions",
    ):

        position_directory = (
            f"{output}/{posi}"
        )

        source_path = (
            f"{position_directory}/pulse_MS.h5"
        )

        output_path = (
            f"{position_directory}/pulse_MS_noise.h5"
        )

        if not os.path.isfile(
            source_path
        ):
            raise FileNotFoundError(
                f"pulse_MS.h5 was not found: {source_path}"
            )

        with PulseWriter(
            output_path,
            time,
            para,
        ) as f:

            for (
                event_id,
                pulse,
            ) in tqdm.tqdm(
                IterPulseItems(
                    source_path
                ),
                desc=f"MS noise {posi}",
                leave=False,
            ):

                ch0 = np.asarray(
                    pulse["ch0"],
                    dtype=float,
                )

                ch1 = np.asarray(
                    pulse["ch1"],
                    dtype=float,
                )

                if (
                    len(ch0) != sample
                    or len(ch1) != sample
                ):
                    raise ValueError(
                        f"pulse {event_id} in {source_path} "
                        f"does not match samples={sample}"
                    )

                (
                    noise0,
                    noise1,
                ) = generate_correlated_noise_pair_from_transfer(
                    transfer_ch0,
                    transfer_ch1,
                    sample,
                    rate,
                    rng=rng,
                )

                f.append(
                    event_id,
                    ch0 + noise0,
                    ch1 + noise1,
                )


def _run_selected_actions(selected_actions):
    """Run selected processing steps in dependency-safe order."""

    actions = {
        "make_pulse": ("MakePulse", MakePulse),
        "fit_ratios": ("FitRatios", FitRatios),
        "make_noise": ("MakeNoise", MakeNoise),
        "show_noise_spectrum": ("_ShowNoiseSpectrum", _ShowNoiseSpectrum),
        "show_samples": ("ShowSamples", ShowSamples),
        "pulse_noise": ("Pulse_Noise", Pulse_Noise),
        "dump_to_event": ("Dump2Event", Dump2Event),
        "pulse_ms": ("Pulse_Ms", Pulse_Ms),
        "pulse_ms_noise": ("Pulse_MS_Noise", Pulse_MS_Noise),
    }

    # The menu is not executed in the order in which the user moves the
    # cursor.  This keeps the common pulse -> ratio/noise -> derived-data
    # pipeline reproducible when several items are selected.
    execution_order = (
        "make_pulse",
        "fit_ratios",
        "make_noise",
        "show_noise_spectrum",
        "show_samples",
        "pulse_noise",
        "dump_to_event",
        "pulse_ms",
        "pulse_ms_noise",
    )

    selected_actions = set(selected_actions)
    for action_key in execution_order:
        if action_key not in selected_actions:
            continue

        action_name, action = actions[action_key]
        print(f"\n=== {action_name} ===")
        action()


def _run_interactive():
    """Ask which simulation steps should be executed."""

    global event_root

    choices = [
        questionary.Choice(
            "パルスを生成 (MakePulse) -> pulses.h5",
            value="make_pulse",
        ),
        questionary.Choice(
            "パルス比を計算 (FitRatios) -> ratios.csv",
            value="fit_ratios",
        ),
        questionary.Choice(
            "ノイズを生成 (MakeNoise) -> noise.h5",
            value="make_noise",
        ),
        questionary.Choice(
            "ノイズスペクトルを出力 (_ShowNoiseSpectrum)",
            value="show_noise_spectrum",
        ),
        questionary.Choice(
            "サンプル波形を表示 (ShowSamples)",
            value="show_samples",
        ),
        questionary.Choice(
            "単一サイト・ノイズパルスを生成 (Pulse_Noise)",
            value="pulse_noise",
        ),
        questionary.Choice(
            "外部dumpall.dat を event.h5 に変換 (Dump2Event)",
            value="dump_to_event",
        ),
        questionary.Choice(
            "外部event.h5 から複数相互作用パルスを生成 (Pulse_Ms)",
            value="pulse_ms",
        ),
        questionary.Choice(
            "複数相互作用 + ノイズを生成 (Pulse_MS_Noise)",
            value="pulse_ms_noise",
        ),
    ]

    print(f"出力先: {output}")
    selected_actions = questionary.checkbox(
        "実行する処理を選択してください（Space: ON/OFF、Enter: 実行）",
        choices=choices,
    ).ask()

    if selected_actions is None:
        print("キャンセルしました。")
        return

    if not selected_actions:
        print("処理が選択されていません。")
        return

    external_data_actions = {
        "dump_to_event",
        "pulse_ms",
    }
    if external_data_actions.intersection(selected_actions) and event_root is None:
        data_root = questionary.path(
            "dumpall.dat / event.h5 のルートを指定してください",
            default=output,
        ).ask()
        if data_root is None:
            print("キャンセルしました。")
            return
        event_root = str(Path(data_root).expanduser())

    _run_selected_actions(selected_actions)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--noise-only",
        action="store_true",
        help="Generate noise.h5 and its filtered ASD without pulse diagnostics.",
    )
    parser.add_argument(
        "--output",
        default=output,
        help="Directory containing input.json and receiving simulation outputs.",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Show the interactive action menu even when other arguments are given.",
    )
    parser.add_argument(
        "--event-root",
        default=None,
        help=(
            "Root directory containing position/<id>/dumpall.dat and "
            "event.h5. pulse_MS outputs remain under --output."
        ),
    )
    script_args = parser.parse_args()
    output = script_args.output
    event_root = script_args.event_root

    if script_args.interactive or len(sys.argv) == 1:
        if script_args.noise_only:
            parser.error("--noise-only cannot be combined with interactive mode.")
        _run_interactive()
    else:
        # Preserve the previous command-line behavior for existing scripts.
        MakeNoise()

        if script_args.noise_only:
            _ShowNoiseSpectrum()
        else:
            ShowSamples()

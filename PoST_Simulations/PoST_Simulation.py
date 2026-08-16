# -*- coding: utf-8 -*-

# --------last updated 2018/12/01 by kurume-------------------

import math
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

# Allow running this script directly from simu_py without installing tes_cpp.
_repo_root = Path(__file__).resolve().parents[1]
_tes_cpp_python = _repo_root / "tes_cpp" / "python"
if str(_tes_cpp_python) not in sys.path:
    sys.path.insert(0, str(_tes_cpp_python))

from lib import general
from tes_cpp import dump2event
from tes_cpp import posi2pulse
from tes_cpp.event_hdf5 import iter_events as iter_hdf5_events
from lib.pulse_hdf5 import PulseWriter, iter_pulse_items as iter_hdf5_pulse_items, read_all_pulses, read_time
# --------------------------------------------------------------
k_b = 1.381 * 1.0e-23  # Boltzmann's constant
ptfn_Flink = 0.5
e = 1.602e-19 
eta = 100
amp = 10
cf = 10000
zure = 30 

pulse_num=500

output="H:\\hata\\662_142_136_300split"

def random_noise(spe, seed):
    spe_re = spe[::-1]  # reverce
    spe_mirror = np.r_[spe, spe_re]
    np.random.seed(seed)
    phase = (2 * np.pi - 0) * np.random.rand(len(spe_mirror))  # random phase
    complex = [cmath.rect(i, j) for i, j in zip(spe_mirror, phase)]
    complex_con = [i.conjugate() for i in complex[len(spe) :]]  # conjugate
    return np.r_[complex[: len(spe)], complex_con]

def generate_noise_from_asd(noise_asd, sample, rate, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    noise_asd = np.asarray(noise_asd[: int(sample / 2) + 1], dtype=float)
    df = rate / sample

    spectrum = np.zeros(len(noise_asd), dtype=np.complex128)
    if len(noise_asd) == 0:
        return np.zeros(sample)

    spectrum[0] = noise_asd[0] * sample * np.sqrt(df)

    if len(noise_asd) > 2:
        phases = rng.uniform(0.0, 2.0 * np.pi, len(noise_asd) - 2)
        magnitude = noise_asd[1:-1] * sample * np.sqrt(df / 2.0)
        spectrum[1:-1] = magnitude * np.exp(1j * phases)

    if sample % 2 == 0 and len(noise_asd) > 1:
        spectrum[-1] = noise_asd[-1] * sample * np.sqrt(df)
    elif len(noise_asd) > 1:
        phase = rng.uniform(0.0, 2.0 * np.pi)
        magnitude = noise_asd[-1] * sample * np.sqrt(df / 2.0)
        spectrum[-1] = magnitude * np.exp(1j * phase)

    return np.fft.irfft(spectrum, n=sample)

def asd_from_rfft(noise_fft, sample, rate):
    noise_fft = np.asarray(noise_fft)
    df = rate / sample
    amp_dens = np.zeros(len(noise_fft), dtype=float)

    if len(noise_fft) == 0:
        return amp_dens

    amp_dens[0] = np.abs(noise_fft[0]) / (sample * np.sqrt(df))

    if len(noise_fft) > 2:
        amp_dens[1:-1] = (
            np.sqrt(2.0) * np.abs(noise_fft[1:-1]) / (sample * np.sqrt(df))
        )

    if len(noise_fft) > 1:
        if sample % 2 == 0:
            amp_dens[-1] = np.abs(noise_fft[-1]) / (sample * np.sqrt(df))
        else:
            amp_dens[-1] = (
                np.sqrt(2.0) * np.abs(noise_fft[-1]) / (sample * np.sqrt(df))
            )

    return amp_dens


def make_noise_time_from_asd(noise_spe_dens, sample, rate, rng=None):
    noise_spe_dens = np.asarray(noise_spe_dens, dtype=float)
    return generate_noise_from_asd(noise_spe_dens, sample, rate, rng=rng)


def MakePulse():
    with open(f"{output}/input.json", "r") as f:
        para = json.load(f)

    # Generate a pulse at every one-based absorber block.  Do not restrict this
    # to input.json's ``position`` array.
    positions = list(range(1, int(para["n_abs"]) + 1))
    input_path = f"{output}/input.json"

    # The native generator is wrapped in the tes_cpp package and writes the
    # shared-time pulse schema as compressed HDF5.
    posi2pulse(input_path, positions, output_path=f"{output}/pulses.h5")

    # Keep the settling-time metadata without reading the large JSON file back.
    reference_pulse = posi2pulse(input_path, [1])[0]
    SettlingTime = int(np.argmax(reference_pulse.ch1))

    para["SettlingTime"]=SettlingTime
    with open(f"{output}/input.json","w") as f:
        json.dump(para,f,indent=4)

def LoadPulses():
    """Load the shared-time, position-keyed HDF5 pulse schema."""
    return read_all_pulses(f"{output}/pulses.h5")

def IterJsonObjectItems(path, chunk_size=1024 * 1024):
    """Yield top-level JSON object items without loading the entire file."""
    decoder = json.JSONDecoder()
    with open(path, "r", encoding="utf-8") as file:
        buffer = ""

        def read_more():
            nonlocal buffer
            chunk = file.read(chunk_size)
            if not chunk:
                raise ValueError(f"unexpected end of JSON input: {path}")
            buffer += chunk

        def parse_value():
            nonlocal buffer
            while True:
                buffer = buffer.lstrip()
                try:
                    value, end = decoder.raw_decode(buffer)
                except json.JSONDecodeError:
                    read_more()
                    continue
                buffer = buffer[end:]
                return value

        while not buffer:
            read_more()
        buffer = buffer.lstrip()
        if not buffer.startswith("{"):
            raise ValueError(f"expected a JSON object at the root of: {path}")
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
            if not isinstance(key, str):
                raise ValueError(f"expected a string key in: {path}")

            while True:
                buffer = buffer.lstrip()
                if buffer:
                    break
                read_more()
            if not buffer.startswith(":"):
                raise ValueError(f"expected ':' after key {key!r} in: {path}")
            buffer = buffer[1:]
            yield key, parse_value()

def IterPulseItems(path, chunk_size=1024 * 1024):
    """Yield pulse records from HDF5 (or legacy JSON during migration)."""
    if Path(path).suffix.lower() in {".h5", ".hdf5"}:
        yield from iter_hdf5_pulse_items(path)
        return
    """Yield a shared-time JSON document's position/event-keyed pulses."""
    decoder = json.JSONDecoder()
    with open(path, "r", encoding="utf-8") as file:
        buffer = ""

        def read_more():
            nonlocal buffer
            chunk = file.read(chunk_size)
            if not chunk:
                raise ValueError(f"unexpected end of JSON input: {path}")
            buffer += chunk

        def parse_value():
            nonlocal buffer
            while True:
                buffer = buffer.lstrip()
                try:
                    value, end = decoder.raw_decode(buffer)
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
            raise ValueError(f"expected a JSON object at the root of: {path}")
        buffer = buffer[1:]

        while True:
            next_non_whitespace()
            if buffer.startswith("}"):
                return
            if buffer.startswith(","):
                buffer = buffer[1:]
                continue

            key = parse_value()
            if not isinstance(key, str):
                raise ValueError(f"expected a string key in: {path}")
            next_non_whitespace()
            if not buffer.startswith(":"):
                raise ValueError(f"expected ':' after key {key!r} in: {path}")
            buffer = buffer[1:]

            if key != "pulses":
                parse_value()
                continue

            next_non_whitespace()
            if not buffer.startswith("{"):
                raise ValueError(f"expected pulses to be an object in: {path}")
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
                if not isinstance(pulse_id, str):
                    raise ValueError(f"expected a pulse ID key in: {path}")
                next_non_whitespace()
                if not buffer.startswith(":"):
                    raise ValueError(
                        f"expected ':' after pulse ID {pulse_id!r} in: {path}"
                    )
                buffer = buffer[1:]
                yield pulse_id, parse_value()
            return

def FitRatios():
    with open(f"{output}/input.json", "r", encoding="utf-8") as f:
        para = json.load(f)
    pulses_by_position = LoadPulses()

    # Absorber blocks are one-based.  With an odd n_abs, the central block has
    # x = 0; x is expressed in mm using the absorber-block pitch.
    center_position = (int(para["n_abs"]) + 1) / 2
    pitch = float(para["length"]) / int(para["n_abs"])
    rows = []
    for position, pulse in pulses_by_position.items():
        ch0_max = max(pulse["ch0"])
        ch1_max = max(pulse["ch1"])
        if ch0_max == 0:
            raise ValueError(f"CH0 maximum is zero at position {pulse['position']}")
        rows.append({
            "position": position,
            "x_mm": (position - center_position) * pitch,
            "ch1_ch0_max_ratio": ch1_max / ch0_max,
        })

    pd.DataFrame(rows).to_csv(f"{output}/ratios.csv", index=False)

def MakeNoise():
        # add omega at diagonal
    def add_omega(M, n_abs, omega):
        omega_list = np.full(n_abs + 4, omega * 1.0j, dtype=np.complex128)
        omega_diag = np.diag(omega_list)
        return M + omega_diag
    
    with open (f"{output}/input.json", "r") as f:
        para = json.load(f)
    n_abs = para["n_abs"]  # absorber pixel
    C_abs = para["C_abs"]# heat capacity per 1-pixel
    C_tes = para["C_tes"]  # heat capacity (TES)
    G_abs_abs = float(para["G_abs-abs"])
    G_abs_tes = para["G_abs-tes"]  # thermal conductivity (absorber-TES)
    G_tes_bath = para["G_tes-bath"]  # thermal conductivity (TES-bath)
    R = para["R"]  # R_TES
    R_l = para["R_l"]  # R_load (shunt)
    T_c = para["T_c"]  # T_c
    T_bath = para["T_bath"]  # T_bath
    a = para["alpha"]  # alpha
    b = para["beta"]  # beta
    L = para["L"]  # Indactance
    n = para["n"]  # dimensionless constant (dominant thermal transport mechanism )
    E = para["E"]  # energy
    length = para["length"]  # length
    rate = int(para["rate"])  # sample rate
    samples = int(para["samples"])  # samples

    time = np.linspace(0, samples / rate, samples)
    frequency = np.arange(0, rate, rate / samples)

    I = np.sqrt((G_tes_bath * T_c * (1 - ((T_bath / T_c) ** n))) / (n * R))  # I_tes

    t_el = L / (R_l + R * (1 + b))  # tau_electron
    L_I = (a * (I**2) * R) / (G_tes_bath * T_c)  # Loop gain
    t_I = C_tes / ((1 - L_I) * G_tes_bath)  # tau_?

    # ----- Noises -----------------------------------------
    # Thremal Fluction Noise
    ptfn_tes_bath = np.sqrt(
        4 * k_b * T_c**2 * G_tes_bath * ptfn_Flink
    )  # Phonon Noise (tes-bath) [W/√Hz]
    ptfn_abs_tes = np.sqrt(
        4 * k_b * T_c**2 * G_abs_tes * ptfn_Flink
    )  # Phonon Noise (abs-tes) [W/√Hz]
    ptfn_abs_abs = np.sqrt(
        4 * k_b * T_c**2 * G_abs_abs * ptfn_Flink
    )  # Phonon Noise (abs-abs) [W/√Hz]

    # Johnson Noise
    enj = np.sqrt(4 * k_b * T_c * R * (1 + 2 * b + b**2))  # at TES
    enj_R = np.sqrt(4 * k_b * T_bath * R_l)  # at R_l

    # noise sources matrix N
    def matrix_N(n_abs):
        X = np.zeros((5, 9), dtype=np.complex128)  # initialize matrix
        X[0,0] = -enj / L#
        X[0,1] = enj_R / L#

        X[1,0] = I * enj / C_tes#
        X[1, 2] = ptfn_tes_bath / C_tes#
        X[1, 3] = ptfn_abs_tes / C_tes#
        X[1,4]=ptfn_abs_abs/C_tes

        X[2,4]=-ptfn_abs_tes/C_abs#
        X[2,5]=-2*ptfn_abs_abs/C_abs
        X[2,6]=-ptfn_abs_tes/C_abs#

        X[3,4]=ptfn_abs_abs/C_tes
        X[3, 5] =ptfn_abs_tes / C_tes#
        X[3, 6] = ptfn_tes_bath / C_tes#
        X[3,8]=I * enj / C_tes#

        X[4, 7] = enj_R / L#
        X[4, 8] = -enj / L#
        return X

    # --------------------------------------------------

    # matrix M without omega
    def matrix_M(n_abs, omega):
        X = np.zeros((5, 5), dtype=np.complex128)  # initialize matrix
        
        X[0, 0] = 1 / t_el + omega * 1.0j
        X[0, 1] = L_I * G_tes_bath / (I * L)

        X[1, 0] = -I * R * (2 + b) / C_tes
        X[1, 1] = 1 / t_I + (G_abs_tes / C_tes) + omega * 1.0j
        X[1, 2] = -G_abs_tes / C_tes

        X[2,1]=-G_abs_tes/C_abs
        X[2,2]=2*G_abs_abs/C_abs+omega*1.0j
        X[2,3]=-G_abs_tes/C_abs

        X[3, 2] = -G_abs_tes / C_tes
        X[3, 3] = 1 / t_I + (G_abs_tes / C_tes) + omega * 1.0j
        X[3, 4] = -I * R * (2 + b) / C_tes

        X[4, 3] = L_I * G_tes_bath / (I * L)
        X[4, 4] = 1 / t_el + omega * 1.0j

        return X


    N = matrix_N(1)

    omega = frequency * 2 * math.pi
    noise = []

    cnt = 0
    for omg in tqdm.tqdm(omega):
        M_inv = np.linalg.inv(matrix_M(n_abs, omg))
        noise_out = np.abs(M_inv[0,:]@N)
        noise.append(noise_out)
    noise = np.array(noise).T
    noise_total = np.sqrt(np.sum(noise ** 2, axis=0))
    components = {
        "johnson_tes1": noise[0, :].tolist(),
        "johnson_load1": noise[1, :].tolist(),
        "phonon_tes1_bath": noise[2, :].tolist(),
        "phonon_tes1_absorber": noise[3, :].tolist(),
        "phonon_absorber_absorber": noise[4, :].tolist(),
        "phonon_tes2_absorber": noise[5, :].tolist(),
        "phonon_tes2_bath": noise[6, :].tolist(),
        "johnson_load2": noise[7, :].tolist(),
        "johnson_tes2": noise[8, :].tolist(),
    }
    import h5py
    with h5py.File(f"{output}/noise.h5", "w") as f:
        f.attrs["input_json"] = json.dumps(para, separators=(",", ":"))
        f.create_dataset("frequency", data=frequency, compression="gzip")
        component_group = f.create_group("components")
        for name, values in components.items():
            component_group.create_dataset(name, data=values, compression="gzip")
        f.create_dataset("total", data=noise_total, compression="gzip")

    # --- grough Noise Spectral Density--------------------------------------------
    plt.figure(figsize=(8, 8))
    plt.plot(
        frequency,
        noise[0,:],
        color="red",
        linewidth=2,
        linestyle=(0, (5, 1)),
        label="Johnson Noise (TES1)",
    )
    plt.plot(
        frequency,
        noise[-1,:],
        color="orange",
        linewidth=2,
        linestyle=(0, (5, 1)),
        label="Johnson Noise (TES2)",
    )
    plt.plot(
        frequency,
        noise[1,:],
        color="lawngreen",
        linewidth=2,
        linestyle=(0, (5, 5)),
        label="Johnson Noise (Load1)",
    )
    plt.plot(
        frequency,
        noise[7,:],
        color="greenyellow",
        linewidth=2,
        linestyle=(0, (5, 5)),
        label="Johnson Noise (Load2)",
    )
    plt.plot(
        frequency,
        noise[2,:],
        color="blue",
        linewidth=2,
        linestyle=(0, (3, 5, 1, 5)),
        label="Phonon Noise (TES1-Bath)",
    )
    plt.plot(
        frequency,
        noise[6,:],
        color="royalblue",
        linewidth=2,
        linestyle=(0, (3, 5, 1, 5)),
        label="Phonon Noise (TES2-Bath)",
    )
    plt.plot(
        frequency,
        noise[3,:],
        color="magenta",
        linewidth=2,
        linestyle=(0, (3, 1, 1, 1)),
        label="Phonon Noise (TES1-Absorber)",
    )
    plt.plot(
        frequency,
        noise[5,:],
        color="pink",
        linewidth=2,
        linestyle=(0, (3, 1, 1, 1)),
        label="Phonon Noise (TES2-Absorber)",
    )

    plt.plot(
        frequency,
        noise[4,:],
        color="dodgerblue",
        linewidth=2,
        linestyle=(0, (3, 1, 1, 1, 1, 1)),
        label="Phonon Noise sum (Absorber-Absorber)",
    )

    plt.plot(
        frequency, noise_total, color="black", linewidth=3, label="Total Noise"
    )

    plt.xlabel("Frequency [Hz]", fontsize=20)
    plt.ylabel("Noise Spectral Density [A/rtHz]", fontsize=20)
    plt.ylim(10e-14, 10e-10)
    plt.xlim(10e-1, 10e5)
    plt.loglog()
    plt.grid()
    plt.legend(loc="best", fancybox=True, fontsize=10, ncol=2)
    plt.tight_layout()
    plt.savefig(f"{output}/noise_all.png", dpi=700)
    #plt.show()
    plt.cla()

def LoadNoise():
    """Read the total noise spectral density from noise.h5."""
    import h5py
    with h5py.File(f"{output}/noise.h5", "r") as f:
        return f["total"][:]

def _ShowNoiseSpectrum():
    with open(f'{output}/input.json', "r") as f:
        para = json.load(f)
    noise_spe_dens = LoadNoise()
    sample=int(para['samples'])
    rate = para["rate"]
    noise_spe_dens = noise_spe_dens[: int(sample / 2) + 1]
    power_model = np.zeros(len(noise_spe_dens))
    for i in tqdm.tqdm(range(100)):
        noise_time = generate_noise_from_asd(noise_spe_dens, sample, rate)
        noise_time = general.Bessel(noise_time, rate, 100000)
        noise_time = general.Bessel(noise_time, rate, para["cutoff"])
        noise_fft = np.fft.rfft(noise_time)
        power_model += np.abs(noise_fft) ** 2
    power_model /= 100
    amp_dens = np.sqrt(power_model)
    amp_dens = asd_from_rfft(amp_dens, sample, rate) * eta * 1e+6
    freq = np.fft.rfftfreq(sample, d=1 / rate)
    if len(freq) > 1:
        freq = freq[:-1]
        amp_dens = amp_dens[:-1]
    plt.plot(freq, amp_dens)
    plt.xlabel("Frequency [Hz]", fontsize=20)
    plt.ylabel("Amplitude [uA/rtHz]", fontsize=20)
    plt.loglog()
    plt.savefig(f"{output}/noise_total-bessel100k.png", dpi=350)
    plt.clf()
    np.savetxt(f"{output}/noise_total-bessel100k.dat", amp_dens)

    if False:
        pre=np.loadtxt(f"{output}/noise_total-bessel100k.dat")
        plt.plot(pre,label="pre")
        plt.plot(amp_dens,label="post")
        plt.xlabel("Frequency [Hz]", fontsize=20)
        plt.ylabel("Amplitude [uA/rtHz]", fontsize=20)
        plt.loglog()
        plt.legend()
        plt.show()

def ShowSamples():
    with open(f'{output}/input.json', "r") as f:
        para = json.load(f)
    pulses_by_position = LoadPulses()

    def pulse_channel(position, channel):
        try:
            values = pulses_by_position[int(position)][f"ch{channel}"]
        except KeyError as error:
            raise ValueError(
                f"pulses.h5 does not contain position {position}"
            ) from error
        # A new ndarray is required because the pulse-with-noise plot modifies it.
        return np.asarray(values, dtype=float)

    sample=int(para['samples'])
    rate = para["rate"]

    # simulated noise frequency domain
    noise_spe_dens = LoadNoise()
    noise_time = make_noise_time_from_asd(noise_spe_dens, sample, rate)
    print(f"spe len:{len(noise_spe_dens)} time len:{len(noise_time)}")

    time = np.linspace(0, para['samples'] / rate, int(para['samples']))
    frequency_pulse = np.arange(0, rate, rate / int(para['samples']))
    # Use a non red/blue palette so the position gradient is easier to read.
    pulse_cmap_name = "viridis"

    # --- Pulse time domain
    fig, ax = plt.subplots()
    cmap = plt.get_cmap(pulse_cmap_name)
    norm = Normalize(vmin=0.0, vmax=float(para["length"]))
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    for i in para["position"]:
        data = pulse_channel(i, 0)
        distance = float(i) / float(para["n_abs"]) * float(para["length"])
        ax.plot(
            time * 1e3,
            data * 1e6,
            color=cmap(norm(distance)),
            linewidth=1.5,
        )
    ax.set_xlabel("Time [ms]", fontsize=20)
    ax.set_ylabel("Current [uA]", fontsize=20)
    plt.xlim(-0.1, 5)
    #plt.ylim(0,2)
    ax.grid()
    fig.tight_layout()
    cbar = fig.colorbar(
        sm,
        ax=ax,
        ticks=np.arange(0, float(para["length"]) + 0.1, 2),
    )
    cbar.set_label("Distance [mm]", fontsize=20)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.invert_yaxis()
    
    fig.savefig(f"{output}/checkpulse_pulses.png", dpi=350, transparent=True)
    plt.close(fig)

    # ---- noise time domain --------
    plt.plot(time * 1e3,noise_time* 1e6,linewidth=1.5)
    plt.xlabel("Time [ms]", fontsize=20)
    plt.ylabel("Current [uA]", fontsize=20)
    plt.grid()
    plt.tight_layout()
    #plt.legend(fontsize=12, loc='upper right')
    plt.savefig(f"{output}/checkpulse_noise.png", dpi=350, transparent=True)
    plt.clf()

    # ---- pulse with noise time domain -------

    fig, ax = plt.subplots()
    cmap = plt.get_cmap(pulse_cmap_name)
    norm = Normalize(vmin=0.0, vmax=float(para["length"]))
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    for i in para["position"]:
        data = pulse_channel(i, 0)
        distance = float(i) / float(para["n_abs"]) * float(para["length"])

        data += make_noise_time_from_asd(noise_spe_dens, sample, rate)
        #data = general.Bessel(data,para['rate'],para['cutoff'])

        ax.plot(
            time * 1e3,
            data * 1e6,
            color=cmap(norm(distance)),
            linewidth=1.5,
        )
    ax.set_xlabel("Time [ms]", fontsize=20)
    ax.set_ylabel("Current [uA]", fontsize=20)
    ax.set_xlim(-0.1, 5)
    ax.grid()
    fig.tight_layout()
    cbar = fig.colorbar(
        sm,
        ax=ax,
        ticks=np.arange(0, float(para["length"]) + 0.1, 2),
    )
    cbar.set_label("Distance [mm]", fontsize=20)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.invert_yaxis()
    fig.savefig(f"{output}/checkpulse_pulses_with_noise.png", dpi=350, transparent=True)
    plt.close(fig)

    os.makedirs(f"{output}/checkpulse_sn_ratio", exist_ok=True)

    # ---- Show Ratio---
    data=np.loadtxt(f"{output}/ratios.csv", delimiter=',', skiprows=1)
    # 1列目をX軸、2列目をY軸
    x = data[:, 1]
    y = data[:, 2]

    plt.scatter(x, y, c=x, cmap='coolwarm', s=50)
    plt.xlabel('Position[mm]')
    plt.ylabel('CH1/CH0[-]')
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(f"{output}/checkpulse_max_ratio.png", dpi=350)
    plt.cla()

    #SN ratio
    for i in para["position"]:
        pulse = pulse_channel(i, 0)
        pulse_noise = pulse + make_noise_time_from_asd(noise_spe_dens, sample, rate)
        pulse_rfft= np.fft.rfft(pulse)
        pulse_noise_rfft= np.fft.rfft(pulse_noise)
        plt.plot(frequency_pulse[:len(pulse_noise_rfft)],np.abs(pulse_noise_rfft),label=f"posi+noise",linestyle='--')
        plt.plot(frequency_pulse[:len(pulse_rfft)],np.abs(pulse_rfft),label=f"pulse")
        plt.plot(frequency_pulse,noise_spe_dens,color='black',label="noise")
        plt.xlim(1,5e4)
        plt.loglog()
        plt.legend()
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Amplitude [A/rtHz]")
        plt.savefig(f"{output}/checkpulse_sn_ratio/sn_ratio_position_{i}.png", dpi=350)
        plt.cla()

    _ShowNoiseSpectrum()

def _legacy_MultiPulse():

    pulse_num=300
    with open(f"{output}/input.json", "r") as f:
        para = json.load(f)
    sample = int(para["samples"])
    rate = para["rate"]

    def AddPulse(noise_spe_dens,data):
        noise_time = make_noise_time_from_asd(noise_spe_dens, sample, rate)
        return data + noise_time[:len(data)]

    def Process(pulse,output,noise_spe_dens,ch,posi,k,para):
        noised_pulse=AddPulse(noise_spe_dens,pulse)
        np.savetxt(f"{output}/{para["E"]}keV_{posi}/pulse_noise/CH{ch}/CH{ch}_{k}.dat",noised_pulse)

    noise_spe_dens = LoadNoise()

    for posi in tqdm.tqdm(para["position"]):
        noise_path=f"{output}/{para["E"]}keV_{posi}/pulse_noise"
        if os.path.exists(noise_path):  # ディレクトリが存在するか確認
            shutil.rmtree(noise_path)

        for ch in [0,1]:
            os.makedirs(f"{output}/{para["E"]}keV_{posi}/pulse_noise/CH{ch}", exist_ok=True)
            pulse=np.loadtxt(f"{output}/{para["E"]}keV_{posi}/pulse/CH{ch}/CH{ch}_1.dat")

            with concurrent.futures.ThreadPoolExecutor() as executor:
                # zipでfile_listとnumbersを組み合わせ、process関数を並行して実行
                futures = [executor.submit(Process, pulse,output,noise_spe_dens,ch,posi,k,para) for k in range(pulse_num)]

                # 結果を待機して処理が終了したら次に進む
                for future in futures:
                    try:
                        future.result()  # 処理結果が必要な場合、ここで結果を取得
                    except Exception as e:
                        print(f"error:{e}")

def Pulse_Noise():
    """Create noisy pulse ensembles from pulses.h5.

    One JSON file is written per requested position. Its schema stores the
    input, one shared time array, and noise realizations keyed by index.
    """
    pulse_num = 500
    with open(f"{output}/input.json", "r", encoding="utf-8") as f:
        para = json.load(f)
    pulses_by_position = LoadPulses()

    sample = int(para["samples"])
    rate = para["rate"]
    noise_spe_dens = LoadNoise()

    def add_noise(data):
        noise_time = make_noise_time_from_asd(noise_spe_dens, sample, rate)
        return data + noise_time[:len(data)]

    for posi in tqdm.tqdm(para["position"]):
        try:
            pulse = pulses_by_position[int(posi)]
        except KeyError as error:
            raise ValueError(
                f"pulses.h5 does not contain position {posi}"
            ) from error

        ch0 = np.asarray(pulse["ch0"], dtype=float)
        ch1 = np.asarray(pulse["ch1"], dtype=float)
        position_directory = f"{output}/{posi}"
        os.makedirs(position_directory, exist_ok=True)
        path = f"{position_directory}/pulse_noise.h5"
        with PulseWriter(path, pulse["time"], para) as f:
            for k in range(pulse_num):
                f.append(k, add_noise(ch0), add_noise(ch1))

def Dump2Event():
    """Convert each position folder's dumpall.dat into event.h5."""
    with open(f"{output}/input.json", "r", encoding="utf-8") as f:
        para = json.load(f)

    # dump2event expects MeV, while input.json E is in keV.
    input_energy_mev = float(para["E"]) / 1000.0
    for posi in tqdm.tqdm(para["position"], desc="dumpall to event"):
        position_directory = f"{output}/{posi}"
        dump_path = f"{position_directory}/dumpall.dat"
        event_path = f"{position_directory}/event.h5"
        if not os.path.isfile(dump_path):
            raise FileNotFoundError(f"dumpall.dat was not found: {dump_path}")
        dump2event(
            dump_path,
            event_path,
            input_energy=input_energy_mev,
            save_all=True,
        )

def Pulse_Ms():
    """Synthesize one CH0/CH1 pulse for each event in every position folder.

    For each input.json position, read its event.h5 and write pulse_MS.h5.
    event.h5 stores event IDs as its outer keys; all deposits in one event
    are summed into a single pulse.
    """
    with open(f"{output}/input.json", "r", encoding="utf-8") as f:
        para = json.load(f)
    reference_pulses = LoadPulses()
    n_abs = int(para["n_abs"])
    length_mm = float(para["length"])
    reference_energy_kev = float(para["E"])
    if length_mm <= 0 or reference_energy_kev <= 0:
        raise ValueError("input.json length and E must be positive")
    if set(range(1, n_abs + 1)) - reference_pulses.keys():
        raise ValueError("pulses.h5 does not contain all absorber positions")

    # PHITS coordinates are cm, while input.json length is mm.  Return None
    # for deposits outside the absorber span instead of assigning them to an
    # edge block.
    def block_from_x_deposit(x_deposit_cm):
        x_mm = float(x_deposit_cm) * 10.0
        fraction = (x_mm + length_mm / 2.0) / length_mm
        if fraction < 0.0 or fraction > 1.0:
            return None
        return min(n_abs, int(fraction * n_abs) + 1)

    time = reference_pulses[1]["time"]
    sample_count = len(time)
    for posi in tqdm.tqdm(para["position"], desc="MS pulse positions"):
        position_directory = f"{output}/{posi}"
        event_path = f"{position_directory}/event.h5"
        output_path = f"{position_directory}/pulse_MS.h5"
        if not os.path.isfile(event_path):
            raise FileNotFoundError(f"event.h5 was not found: {event_path}")

        with PulseWriter(output_path, time, para) as f:
            for index, (event_id, event) in enumerate(
                tqdm.tqdm(iter_hdf5_events(event_path), desc=f"MS pulses {posi}", leave=False)
            ):
                energy_by_block = np.zeros(n_abs, dtype=float)
                if not isinstance(event, dict):
                    raise ValueError(f"event.h5 event {event_id} is not an object")
                for particle in event.values():
                    positions = particle.get("x_deposit", [])
                    energies = particle.get("E_deposit", [])
                    if len(positions) != len(energies):
                        raise ValueError(
                            f"event.h5 event {event_id} has mismatched "
                            "x_deposit and E_deposit lengths"
                        )
                    for x_deposit, e_deposit_mev in zip(positions, energies):
                        block = block_from_x_deposit(x_deposit)
                        if block is not None:
                            energy_by_block[block - 1] += float(e_deposit_mev)

                ch0 = np.zeros(sample_count, dtype=float)
                ch1 = np.zeros(sample_count, dtype=float)
                for block, energy_mev in enumerate(energy_by_block, start=1):
                    if energy_mev == 0.0:
                        continue
                    scale = energy_mev * 1000.0 / reference_energy_kev
                    reference = reference_pulses[block]
                    ch0 += scale * np.asarray(reference["ch0"], dtype=float)
                    ch1 += scale * np.asarray(reference["ch1"], dtype=float)

                f.append(event_id, ch0, ch1)

def _legacy_MS_Noise():
    with open(f"{output}/input.json", "r") as f:
        para = json.load(f)
    sample = int(para["samples"])
    rate = para["rate"]

    def AddPulse(noise_spe_dens,data):
        noise_time = make_noise_time_from_asd(noise_spe_dens, sample, rate)
        return data + noise_time[:len(data)]
    
    noise_spe_dens = LoadNoise()

    def Process(file,output,noise_spe_dens,ch,num,posi,para):
        pulse=np.loadtxt(file)
        noised_pulse=AddPulse(noise_spe_dens,pulse)
        np.savetxt(f"{output}/{para["E"]}keV_{posi}/pulse_noise_ms_test/CH{ch}/CH{ch}_{num}.dat",noised_pulse)

    for posi in tqdm.tqdm(para["position"]):
        for ch in [0,1]:
            os.makedirs(f"{output}/{para["E"]}keV_{posi}/pulse_noise_ms_test/CH{ch}", exist_ok=True)

            file_pattern = f"{output}/{para["E"]}keV_{posi}/pulse_ms/CH{ch}/CH{ch}_*.dat"
            file_list = glob.glob(file_pattern)
            numbers = []

            for file in file_list:
            # ファイル名から数字部分を抽出 (例: CH0_1.dat -> 1)
                match = re.search(r'CH\d+_(\d+).dat', file)
                if match:
                # 数字をリストに追加
                    numbers.append(int(match.group(1)))

            with concurrent.futures.ThreadPoolExecutor() as executor:
            # zipでfile_listとnumbersを組み合わせ、process関数を並行して実行
                futures = [executor.submit(Process, file, output,noise_spe_dens,ch,number,posi,para) for file, number in zip(file_list, numbers)]

                # 結果を待機して処理が終了したら次に進む
                for future in futures:
                    future.result()  # 処理結果が必要な場合、ここで結果を取得

def Pulse_MS_Noise():
    """Add independent random noise to every pulse in each pulse_MS.h5."""
    with open(f"{output}/input.json", "r", encoding="utf-8") as f:
        para = json.load(f)
    sample = int(para["samples"])
    rate = para["rate"]
    noise_spe_dens = LoadNoise()
    time = LoadPulses()[1]["time"]

    def add_noise(data):
        noise_time = make_noise_time_from_asd(noise_spe_dens, sample, rate)
        return data + noise_time[:len(data)]

    for posi in tqdm.tqdm(para["position"], desc="MS noise positions"):
        position_directory = f"{output}/{posi}"
        source_path = f"{position_directory}/pulse_MS.h5"
        output_path = f"{position_directory}/pulse_MS_noise.h5"
        if not os.path.isfile(source_path):
            raise FileNotFoundError(f"pulse_MS.h5 was not found: {source_path}")

        with PulseWriter(output_path, time, para) as f:
            for index, (event_id, pulse) in enumerate(
                tqdm.tqdm(IterPulseItems(source_path), desc=f"MS noise {posi}", leave=False)
            ):
                ch0 = np.asarray(pulse["ch0"], dtype=float)
                ch1 = np.asarray(pulse["ch1"], dtype=float)
                if len(ch0) != sample or len(ch1) != sample:
                    raise ValueError(
                        f"pulse {event_id} in {source_path} does not match samples={sample}"
                    )
                f.append(event_id, add_noise(ch0), add_noise(ch1))

#MakePulse()
#FitRatios()
#MakeNoise()
ShowSamples()
#Pulse_Noise()
#Dump2Event()
#Pulse_Ms()
#Pulse_MS_Noise()

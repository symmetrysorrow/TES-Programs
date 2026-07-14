# -*- coding: utf-8 -*-

# --------last updated 2018/12/01 by kurume-------------------

import math
import ctypes
import shutil
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import matplotlib.cm as cm
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
import general
import tqdm
import random

# --------------------------------------------------------------
k_b = 1.381 * 1.0e-23  # Boltzmann's constant
ptfn_Flink = 0.5
e = 1.602e-19 
eta = 100
amp = 10
cf = 10000
zure = 30 

pulse_num=500

output="H:/hata2025/1332_215_195-trial"
ENERGY_CANDIDATES = [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300]

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


def MakeNoiseAndPulse():
    with open(f'{output}/input.json', "r") as f:
        para = json.load(f)
    noise_spe_dens = np.loadtxt(f"{output}/noise_total.dat")
    sample=int(para['samples'])
    rate = para["rate"]
    noise_spe_dens = noise_spe_dens[: int(sample / 2) + 1]
    power_model = np.zeros(len(noise_spe_dens))
    pulse=np.loadtxt(f"{output}/1332keV_17/Pulse/CH0/CH0_1.dat")
    pulse_start=20000
    for i in tqdm.tqdm(range(100)):
        random_coeff=np.random.rand()
        noise_time = generate_noise_from_asd(noise_spe_dens, sample, rate)
        noise_time+=random_coeff*pulse[pulse_start:pulse_start+sample]
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
    plt.plot(amp_dens,label="with pulse tail")
    plt.xlabel("Frequency [Hz]", fontsize=20)
    plt.ylabel("Amplitude [uA/rtHz]", fontsize=20)
    plt.loglog()
    plt.savefig(f"{output}/noise_total-bessel100k-with-pulseTail.png", dpi=350)
    amp_pre=np.loadtxt(f"{output}/noise_total-bessel100k.dat")
    plt.plot(amp_pre,label="pre")
    plt.legend()
    plt.show()
    plt.clf()
    np.savetxt(f"{output}/noise_total-bessel100k-with-pulseTail.dat",amp_dens)

def ShowPulse():
    with open(f"{output}/input.json", "r") as f:
        para = json.load(f)
    time=general.GetTime(para["rate"],para["samples"])
    for ene in ENERGY_CANDIDATES:
        Pulse=np.loadtxt(f"{output}/pulse/CH1_ene={ene}.dat")
        plt.plot(time,Pulse,label=f"{ene}keV")
    plt.legend()
    plt.show()


MakeNoiseAndPulse()

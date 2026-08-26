# 2022/09/12
# それぞれのデータの周波数スペクトルの平均（モデルノイズ）を出力


import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
import os
from natsort import natsorted
import glob
import shutil
import libs.getpara as gp
import pandas as pd
import json
import tqdm

try:
    from tes_analysis.noise_utils import (
        one_sided_asd_from_power,
        voltage_asd_to_pA,
    )
except ImportError:  # Package import (e.g. from the repository root).
    from .tes_analysis.noise_utils import (
        one_sided_asd_from_power,
        voltage_asd_to_pA,
    )


# 実行
def main():
    set = gp.loadJson()
    if "eta_uA_per_V" not in set["Config"]:
        if "eta" in set["Config"]:
            # Accept old settings files once, while using the explicit unit
            # name for all subsequent calculations and saved settings.
            set["Config"]["eta_uA_per_V"] = float(set["Config"]["eta"])
        else:
            eta_uA_per_V = input("eta [uA/V]:")
            set["Config"]["eta_uA_per_V"] = float(eta_uA_per_V)
        jsn = json.dumps(set, indent=4)
        with open("setting.json", "w") as file:
            file.write(jsn)
    os.chdir(set["Config"]["path"])
    ch = set["Config"]["channel"]
    rate, samples = set["Config"]["rate"], set["Config"]["samples"]
    eta_uA_per_V = float(set["Config"]["eta_uA_per_V"])
    time = gp.data_time(rate, samples)
    fq = np.arange(0, rate, rate / samples)
    output = f'CH{set["Config"]["channel"]}_noise/output/{set["Config"]["output"]}'

    window = np.hanning(samples)
    window_power_gain = np.sqrt(np.mean(window**2))
    power_model = np.zeros(samples // 2 + 1)
    accepted = 0
    noise = natsorted(glob.glob(f"CH{ch}_noise/rawdata/CH{ch}_*.dat"))

    for i in tqdm.tqdm(noise):
        try:
            data = gp.loadbi(i, "binary")
            base, data_ba = gp.baseline(data, set["Config"]["presamples"], 1000, 500)
            # Remove each record's DC component before filtering and FFT so all
            # noise-generation paths use the same AC-only analysis.
            data = data - np.mean(data)
            if set["main"]["cutoff"] > 0:
                data = gp.BesselFilter(data, rate, set["main"]["cutoff"])
            peak = np.max(data_ba)
            if (base <= -3 and base >= 3) or peak >= float(set["Config"]["threshold"]):
                print("error")
                continue
            else:
                data_fft = np.fft.rfft(data * window)
            power_model += np.abs(data_fft) ** 2
            accepted += 1

        except FileNotFoundError:
            continue
    
    if accepted == 0:
        raise RuntimeError("No accepted noise records found")
    amp_dens = one_sided_asd_from_power(
        power_model / accepted,
        samples,
        rate,
        window_power_gain,
    )
    amp_dens = voltage_asd_to_pA(amp_dens, eta_uA_per_V)
    print(amp_dens)
    np.savetxt(f"{output}/modelnoise.txt", amp_dens)

    # スペクトルをグラフ化
    plt.plot(fq[: int(samples / 2) + 1], amp_dens, linestyle="-", linewidth=0.7)
    plt.loglog()
    plt.xlabel("Frequency[Hz]")
    plt.ylabel("Intensity[pA/Hz$^{1/2}$]")
    plt.grid()
    plt.savefig(f"{output}/modelnoise.png")
    plt.show()


if __name__ == "__main__":
    main()

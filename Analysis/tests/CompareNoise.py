import numpy as np
from tes_analysis import analysis_utils as gn
import matplotlib.pyplot as plt
import questionary
import glob
import os

def CompareNoise():
    target_freq = 1000  # 正規化したい周波数[Hz]

    sim_path = "h:/hata2025/1332_200_180"
    sim_noise_data_be = gn.LoadTxt(f"{sim_path}/noise_total-bessel100k.dat")
    para = gn.LoadJson(f"{sim_path}/input.json")
    rate_sim = para["rate"]
    fq_sim_be = gn.GetFreq(rate_sim/2, len(sim_noise_data_be))

    # シミュレーション側も同様に100Hz正規化
    idx_100_sim = (np.abs(fq_sim_be - target_freq)).argmin()
    sim_noise_data_be /= sim_noise_data_be[idx_100_sim]

    plt.plot(fq_sim_be, sim_noise_data_be, label="New Simulation")

    print("old simulation")
    sim_path = "h:/hata2025/1332_215_195-old"
    sim_noise_data_be = gn.LoadTxt(f"{sim_path}/noise_total-bessel100k.dat")
    para = gn.LoadJson(f"{sim_path}/input.json")
    rate_sim = para["rate"]
    fq_sim_be = gn.GetFreq(rate_sim/2, len(sim_noise_data_be))

    # シミュレーション側も同様に100Hz正規化
    idx_100_sim = (np.abs(fq_sim_be - target_freq)).argmin()
    sim_noise_data_be /= sim_noise_data_be[idx_100_sim]

    #plt.plot(fq_sim_be, sim_noise_data_be, label="Old Simulation")

    print("experimental data")
    path = "G:/TSURUTA/20230616_post/room1-ch2-3_180mK_570uA_100kHz_g10"
    ch = "CH0_noise"
    noise_data = gn.LoadTxt(f"{path}/{ch}/modelnoise.txt")
    para = gn.LoadJson(f"{path}/PulseConfig.json")
    rate_exp = para["Readout"]["Rate"]
    fq_exp = gn.GetFreq(rate_exp/2, len(noise_data))

    # --- ここで100Hzで正規化 ---
    idx_100 = (np.abs(fq_exp - target_freq)).argmin()  # 最も近い周波数のインデックス
    noise_data /= noise_data[idx_100]

    plt.plot(fq_exp, noise_data, label="Experimental Data",linestyle="--")

    
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized Amplitude Density [1/√Hz]")
    plt.legend()
    plt.loglog()
    plt.show()

CompareNoise()

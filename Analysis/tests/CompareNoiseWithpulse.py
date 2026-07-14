import numpy as np
from tes_analysis import analysis_utils as gn
import matplotlib.pyplot as plt
import questionary
import glob
import os

def CompareNoise():
    target_freq = 1000  # 正規化したい周波数[Hz]

    sim_path="h:/hata2025/1332_215_195-trial"
    #sim_noise_data_be = gn.LoadTxt(f"{sim_path}/noise_total-bessel100k-with-pulsetail.dat")
    sim_noise_data_be = gn.LoadTxt(f"{sim_path}/noise_total-bessel100k-with-pulseTail.dat")
    para = gn.LoadJson(f"{sim_path}/input.json")
    rate_sim = para["rate"]
    fq_sim_be = gn.GetFreq(rate_sim/2, len(sim_noise_data_be))

    # シミュレーション側も同様に100Hz正規化
    idx_100_sim = (np.abs(fq_sim_be - target_freq)).argmin()
    sim_noise_data_be /= sim_noise_data_be[idx_100_sim]

    plt.plot(fq_sim_be, sim_noise_data_be, label="Sim : Noise+Pulse")

    print("old simulation")
    sim_path = "h:/hata2025/1332_215_195-trial"
    sim_noise_data_be = gn.LoadTxt(f"{sim_path}/noise_total-bessel100k.dat")
    para = gn.LoadJson(f"{sim_path}/input.json")
    rate_sim = para["rate"]
    fq_sim_be = gn.GetFreq(rate_sim/2, len(sim_noise_data_be))

    # シミュレーション側も同様に100Hz正規化
    idx_100_sim = (np.abs(fq_sim_be - target_freq)).argmin()
    sim_noise_data_be /= sim_noise_data_be[idx_100_sim]

    plt.plot(fq_sim_be, sim_noise_data_be, label="Sim : Noise Only")

    print("experimental data")
    path="G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2"
    ch = "CH0_noise"
    noise_data = gn.LoadTxt(f"{path}/{ch}/modelnoise.txt")
    para = gn.LoadJson(f"{path}/PulseConfig.json")
    rate_exp = para["Readout"]["Rate"]
    fq_exp = gn.GetFreq(rate_exp/2, len(noise_data))

    # --- ここで100Hzで正規化 ---
    idx_100 = (np.abs(fq_exp - target_freq)).argmin()  # 最も近い周波数のインデックス
    noise_data /= noise_data[idx_100]

    #plt.plot(fq_exp, noise_data, label="FilteredExperimental Data",linestyle=":")

    noise_data_pre=gn.LoadTxt(f"{path}/{ch}/modelnoise-pre.txt")
    noise_data_pre /= noise_data_pre[idx_100]
    plt.plot(fq_exp, noise_data_pre, label="Experimental Data",linestyle="--",color="green")

    
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized Amplitude Density [1/√Hz]")
    plt.legend()
    plt.loglog()
    plt.xlim(10,1e6)
    plt.savefig("D:/desktop/CompareNoiseWithPulse.png",transparent=True)
    plt.show()

CompareNoise()

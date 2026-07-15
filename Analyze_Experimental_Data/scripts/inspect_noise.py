from tes_analysis import analysis_utils as gn
import matplotlib.pyplot as plt
import os

sim_path = "h:/hata2025/200_180"
exp_path="G:/TSURUTA/20230616_post/room1-ch2-3_180mK_570uA_100kHz_g10"
exp_path="G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2"

noise_key=input("Enter noise key: ")

noise=gn.LoadBin(f"{exp_path}/CH0_noise/rawdata/CH0_{noise_key}.dat")

para=gn.LoadJson(f"{exp_path}/PulseConfig.json")
time=gn.GetTime(para["Readout"]["Rate"],para["Readout"]["Sample"])
plt.plot(time,noise,label="Experimental Data",linestyle="--",color="green")
noise=gn.Bessel(noise,para["Readout"]["Rate"],para["Analysis"]["CutoffFrequency"])
plt.plot(time,noise,label="Filtered Experimental Data",linestyle=":")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude [V]")
plt.legend()
plt.show()

import numpy as np
from tes_analysis import analysis_utils as gn
import matplotlib.pyplot as plt
import glob
import questionary
import re
import os
import pandas as pd

path=input("Enter experimental data path: ")

folders = glob.glob(os.path.join(path, "CH*_pulse"))
chs = [re.search(r'CH(.*)_pulse', os.path.basename(f)).group(1) for f in folders]

para=gn.LoadJson(f"{path}/PulseConfig.json")

time=gn.GetTime(para["Readout"]["Rate"], para["Readout"]["Sample"])

ch=questionary.select("Select Channel to plot pulses:", choices=chs).ask()

SelectedKeys=gn.LoadTxt(f"{path}/SelectedKeys.txt")

Presample=int(para["Readout"]["Rate"]/10)

keys=[]

for key in SelectedKeys:
    pulse=gn.LoadBin(f"{path}/CH{ch}_pulse/rawdata/CH{ch}_{int(key)}.dat")
    pulse=pulse.copy()
    pulse=gn.Bessel(pulse, para["Readout"]["Rate"], para["Analysis"]["CutoffFrequency"])
    
    pulse-=np.mean(pulse[0:Presample])
    if pulse[Presample]<0.13 or pulse[Presample*2]<0.1 or np.max(pulse[Presample*2:])>0.4:
        continue
    plt.plot(time,pulse,color="gray",alpha=0.5)
    keys.append(int(key))

plt.xlabel("Sample")
plt.ylabel("Amplitude")
plt.show()

np.savetxt(f"{path}/SelectedKeys_fromScatter.txt", keys)

pandas_df0=pd.read_csv(f"{path}/CH0_pulse/output.csv")
pandas_df1=pd.read_csv(f"{path}/CH1_pulse/output.csv")

peaks_0 = pandas_df0.loc[pandas_df0["key"].isin(keys), "Peak"].to_numpy()
peaks_1 = pandas_df1.loc[pandas_df1["key"].isin(keys), "Peak"].to_numpy()

plt.scatter(peaks_0, peaks_1)
plt.xlabel("CH0 Peak")
plt.ylabel("CH1 Peak")
plt.show()

peaks_0=np.array(peaks_0)
peaks_0_min=np.min(peaks_0)
peaks_0/=peaks_0_min

plt.hist(peaks_0, bins=50)
plt.show()

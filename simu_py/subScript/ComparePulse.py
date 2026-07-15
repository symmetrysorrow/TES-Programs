import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import general as gn
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path

amp=10
eta=100

def detect_rise_index(pulse, frac=0.1):
    peak = np.max(pulse)
    return np.argmax(pulse > frac * peak)

def align_pulse_to_presample(pulse, rise_idx, presample):
    N = len(pulse)
    aligned = np.zeros_like(pulse)

    shift = presample - rise_idx

    if shift > 0:
        ncopy = N - shift
        aligned[shift:shift+ncopy] = pulse[:ncopy]
    else:
        ncopy = N + shift
        aligned[:ncopy] = pulse[-shift:-shift+ncopy]

    return aligned

exp_path="G:/TSURUTA/20230616_post/room1-ch2-3_180mK_570uA_100kHz_g10"
selected_keys=gn.LoadTxt(f"{exp_path}/SelectedKeys_fromScatter.txt")
para=gn.LoadJson(f"{exp_path}/PulseConfig.json")

Presample=int(para["Readout"]["PreSample"])

t=gn.GetTime(para["Readout"]["Rate"], para["Readout"]["Sample"])

for key in selected_keys:
    pulse = gn.LoadBin(f"{exp_path}/CH0_pulse/rawdata/CH0_{int(key)}.dat")
    pulse = pulse.copy()
    pulse-=np.mean(pulse[0:Presample])
    pulse = gn.Bessel(pulse, para["Readout"]["Rate"], para["Analysis"]["CutoffFrequency"])
    rise_idx = detect_rise_index(pulse)
    pulse = align_pulse_to_presample(pulse, rise_idx, Presample)
    plt.plot(t,pulse*eta/amp/1e6/48, color="gray", alpha=0.5)

#sim_path="h:/hata2025/200_180"
#para_sim=gn.LoadJson(f"{sim_path}/input.json")
#if False:
#    t=gn.GetTime(para_sim["rate"], para_sim["samples"])
#    t*=1
#cnt=0
#for i in para_sim["position"]:
#    pulse_sim = np.loadtxt(f"{sim_path}/{para_sim["E"]}keV_{i}/pulse/CH0/CH0_1.dat")
#    pulse_sim = pulse_sim.copy()
#    #pulse_sim = gn.Bessel(pulse_sim, para_sim["rate"], 100000)
#    zeros=np.zeros(int(Presample/1))
#    pulse_sim=np.concatenate([zeros, pulse_sim])[:int(para_sim["samples"])]
#    plt.plot(t,pulse_sim, color=cm.hsv((float(cnt)) / float(len(para_sim["position"]))))
#    cnt+=1

def load_post_center(path: Path) -> np.ndarray:
    """Load COMSOL text table while ignoring comment/header lines."""

    data = np.loadtxt(path, comments="%", encoding="utf-8")
    if data.ndim != 2 or data.shape[1] < 5:
        raise ValueError(
            "Expected at least 5 numeric columns: time, temperatures, and currents."
        )
    return data

data=load_post_center(r"D:\Github\simu_py\PoST_186.txt")

time_ms = data[:, 0]
Current_TES1=data[:, 3]
Current_TS2=data[:, 4]

# 0~100msのゼロデータを最初に追加

Current_TES1-=np.mean(Current_TES1[0])    
Current_TS2-=np.mean(Current_TS2[0]) 

prepend_time = np.arange(0, 100, 1)
time_ms+=100
time_ms = np.concatenate([prepend_time, time_ms])
Current_TES1 = np.concatenate([np.zeros(len(prepend_time)), Current_TES1])
Current_TS2 = np.concatenate([np.zeros(len(prepend_time)), Current_TS2])

   

plt.plot(time_ms/1e3,Current_TES1, label="TES1")
plt.plot(time_ms/1e3,Current_TS2, label="TES2")

plt.xlabel("Time[s]")
plt.ylabel("Current[A]")
#plt.xlim(0.098,0.1025)
plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
import os
from natsort import natsorted
import glob
import shutil

import pandas as pd
from tkinter import filedialog

eta = 98



def main():
    set = gp.loadJson()
    os.chdir('H:/Matsumi/data/20230512/room2-2_140mK_870uA_gain10_trig0.4_500kHz/CH0_pulse/filter')
    ch = set['Config']['channel']
    rate,samples= 1e6,1e5
    fq = np.arange(0,rate,rate/samples)



    five_hund = np.loadtxt("modelnoise_500kHz.txt")
    hund = np.loadtxt("modelnoise_100kHz.txt")
    ten = np.loadtxt("modelnoise_10kHz.txt")

    five_hund = np.sqrt(five_hund[:int(samples/2)+1])*int(eta)*1e+6*np.sqrt(1/rate/samples)
    hund = np.sqrt(hund[:int(samples/2)+1])*int(eta)*1e+6*np.sqrt(1/rate/samples)
    ten = np.sqrt(ten[:int(samples/2)+1])*int(eta)*1e+6*np.sqrt(1/rate/samples)

    plt.plot(fq[:int(samples/2)+1],five_hund,label="500 kHz")
    plt.plot(fq[:int(samples/2)+1],hund,label="100 kHz")
    plt.plot(fq[:int(samples/2)+1],ten,label="10 kHz")
    plt.xlabel('Frequency[Hz]')
    plt.ylabel('Intensity[pA/kHz$^{1/2}$]')
    plt.ylim(10,1e5)
    plt.legend()
    plt.grid()
    plt.loglog()
    plt.savefig(f'H:/Matsumi/data/20230512/room2-2_140mK_870uA_gain10_trig0.4_500kHz/CH0_pulse/filter/noise_matome.png')
    
    plt.show()


if __name__ == '__main__':
    fld = filedialog.askdirectory()
    print(fld)
    
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import os
import pandas as pd
from  scipy import signal


rate = 1e6
samples = 1e5
cmap = cm.get_cmap("hsv")


def main():
    
    hund = np.loadtxt("H:/Matsumi/data/20230510/room2-2_140mK_870uA_gain10_trig0.1_10kHz/compare/modelnoise_100kHz.txt")
    hund_pyt = np.loadtxt("H:/Matsumi/data/20230510/room2-2_140mK_870uA_gain10_trig0.1_10kHz/compare/modelnoise_100kHz_python.txt")
    time = np.arange(0,1/rate*samples,1/rate)
    fq = np.arange(0,rate,rate/samples)

    plt.plot(fq[:int(samples/2)+1],hund[:int(samples/2)+1],label = "SRS")
    plt.plot(fq[:int(samples/2)+1],hund_pyt[:int(samples/2)+1],label = "python")
    plt.xlabel("Freaquency(Hz)")
    plt.ylabel('Intensity[pA/kHz$^{1/2}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('model noise spectrum (100kHz)')
    plt.grid()
    plt.legend()
    plt.show()
    plt.cla()

if __name__ == "__main__":
    main()
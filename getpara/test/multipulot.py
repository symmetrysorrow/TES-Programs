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

def BesselFilter(x,rate,fs):
    #fn = rate/2
    #ws = fs/fn
    ws = fs/rate*2
    b,a = signal.bessel(2,ws,"low")
    y = signal.filtfilt(b,a,x)
    return y


def graugh_spe(fq,data):
    plt.plot(fq[:int(samples/2)+1],data[:int(samples/2)+1],label = "10kHz")

    plt.xlabel("Freaquency(Hz)")
    plt.ylabel('Intensity[pA/kHz$^{1/2}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Average pulse spectrum')
    plt.grid()
    plt.legend()
    plt.show()
    plt.cla()

def main():
    
    ten = np.loadtxt("H:/Matsumi/data/20230510/room2-2_140mK_870uA_gain10_trig0.1_10kHz/compare/average_pulse_10kHz.txt")
    hund = np.loadtxt("H:/Matsumi/data/20230512/room2-2_140mK_870uA_gain10_trig0.4_500kHz/CH0_pulse/output_500kHz/select/average_pulse.txt")
    time = np.arange(0,1/rate*samples,1/rate)
    fq = np.arange(0,rate,rate/samples)


    cutoff = [10,20,40,60,80,90,100,200,250]
    cnt = 0
    plt.plot(time,ten,label = f"10 kHz (raw)")
    for i in cutoff:
        filt = BesselFilter(hund,1e6,i*1000)
        plt.plot(time,filt,'--',label = f"{i} kHz",color = cmap((float(cnt))/float(len(cutoff))))
        cnt +=1
    plt.plot(time,hund,label =f"500 kHz (raw)")
    plt.title('Average pulse')
    plt.xlabel("time [s]")
    plt.ylabel('Volt [V]')
    plt.grid()
    plt.legend()
    plt.show()

    f_100k = np.fft.fft(hund)
    F_100k =np.abs(f_100k)
    f_10k = np.fft.fft(ten)
    F_10k =np.abs(f_10k)

    cnt = 0
    plt.plot(fq[:int(samples/2)+1],F_10k[:int(samples/2)+1],label = f"10 kHz (raw)")
    for i in cutoff:
        filt = BesselFilter(hund,1e6,i*1000)
        f = np.fft.fft(filt)
        F = np.abs(f)
        plt.plot(fq[:int(samples/2)+1],F[:int(samples/2)+1],'--',label = f"{i} kHz",color = cmap((float(cnt))/float(len(cutoff))))
        cnt +=1
    
    plt.plot(fq[:int(samples/2)+1],F_100k[:int(samples/2)+1],label =f"500 kHz (raw)")
    plt.xlabel("Freaquency(Hz)")
    plt.ylabel('Intensity[pA/kHz$^{1/2}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Average pulse spectrum')
    plt.grid()
    plt.legend()
    plt.show()

    """
    hund_10k = BesselFilter(hund,1e6,1e4)

    f_100k = np.fft.fft(hund)
    F_100k =np.abs(f_100k)

    f_10k = np.fft.fft(hund_10k)
    F_10k =np.abs(f_10k)
    plt.xlabel("Freaquency(Hz)")
    plt.ylabel('Intensity[pA/kHz$^{1/2}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Average pulse spectrum')
    plt.grid()
    plt.legend()
    plt.show()
    plt.cla()
    """



    

    """
    plt.plot(fq[:int(samples/2)+1],ten[:int(samples/2)+1],label = "10kHz")
    plt.plot(fq[:int(samples/2)+1],hund[:int(samples/2)+1],label = "100kHz")
    plt.xlabel("Freaquency(Hz)")
    plt.ylabel('Intensity[pA/kHz$^{1/2}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Average pulse spectrum')
    plt.grid()
    plt.legend()
    plt.show()
    """


if __name__ == "__main__":
    main()
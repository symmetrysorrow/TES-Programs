
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
import cmath

# eta[uA/V]

#--------create random noise-----------
def random_noise(spe):
    spe_re = spe[::-1] #reverce
    spe_mirror = np.r_[spe,spe_re]
    phase = (2*np.pi) * np.random.rand(len(spe_mirror)) #random phase
    complex = [cmath.rect(i,j) for i,j  in zip(spe_mirror,phase)] 
    complex_con = [i.conjugate() for i in complex[len(spe):]] #conjugate
    return np.r_[complex[:len(spe)],complex_con]
#--------------------------------------


def reverse_complex(abs):
    random_phase = np.random.uniform(0, 2 * np.pi,len(abs))
    return abs * np.exp(1j * random_phase)


#実行
def main():
    set = gp.loadJson()
    if not 'eta' in set['Config']:
        eta = input("eta: ")
        set['Config']['eta'] = float(eta)
        jsn = json.dumps(set,indent=4)
        with open("setting.json", 'w') as file:
            file.write(jsn) 
    os.chdir(set["Config"]["path"])
    ch = set['Config']['channel']
    rate,samples= set["Config"]["rate"],set["Config"]["samples"]
    eta = int(set['Config']['eta'])*1e-6 # [A/V]
    time = gp.data_time(rate,samples)
    fq = np.arange(0,rate,rate/samples)
    output = f'CH{set["Config"]["channel"]}_noise/output/{set["Config"]["output"]}'

    # rawdata
    model = np.loadtxt(f'{output}/modelnoise.txt')
    rawdata = gp.loadbi(f'CH{ch}_noise/rawdata/CH0_1.dat')
    
    base,rawdata = gp.baseline(rawdata,set['Config']['presamples'],1000,500)
    rawdata_fft = np.fft.fft(rawdata)
    freq = np.fft.fftfreq(samples,d=1/rate)


    amp_scale = np.abs(rawdata_fft/(samples/2)) #/(samples/2)
    df = rate/samples
    power = amp_scale**2/df
    noise_spe_dens = np.sqrt(power)*eta
    noise_spe_dens_V = np.sqrt(power)


    amp = np.abs(rawdata_fft)**2


    noise_rd =  random_noise(model)#*(samples/2)
    noise_time = np.fft.ifft(noise_rd,samples).real
    
    base,noise_time_ifft = gp.baseline(noise_time,set['Config']['presamples'],1000,500)


    plt.plot(time,rawdata)
    plt.title('rawdata')
    plt.xlabel('time [s]')
    plt.ylabel('Volt [V]')
    plt.show()
    plt.cla()

    plt.plot(freq[1:int(samples/2)],noise_spe_dens[1:int(samples/2)]*1e12)
    plt.grid()
    plt.title('fft')
    plt.xlabel('Frequency[Hz]')
    plt.ylabel('Intensity[pA/Hz$^{1/2}$]')
    plt.loglog()
    plt.show()

    plt.plot(freq[1:int(samples/2)],noise_spe_dens_V[1:int(samples/2)])
    plt.grid()
    plt.title('fft')
    plt.xlabel('Frequency[Hz]')
    plt.ylabel('Intensity[V/Hz$^{1/2}$]')
    plt.loglog()
    plt.show()

    plt.plot(time,noise_time_ifft)
    plt.title('ifft')
    plt.xlabel('time [s]')
    plt.ylabel('Volt [V]')
    plt.show()
    plt.cla()


if __name__ == '__main__':
    main()
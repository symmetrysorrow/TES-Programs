import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import shutil
from natsort import natsorted
import libs.getpara as gp
import libs.fft_spectrum as sp

os.chdir('/Users/matsumi/teion/data/20220902/CH1_150mK_710uA_trig0.1_gain10_10kHz') 

rate = int(250000)
samples = int(50000)
presamples = int(5100)
threshold = 0.02
time = np.arange(0,1/rate*samples,1/rate)

path = natsorted(glob.glob('CH1/rawdata/CH1_*.dat'))


data_array = []
for num in path:
    print(os.path.basename(num))
    data = gp.loadtxt(num)
    error = gp.check_samples(data,samples)
    base,data = gp.baseline(data,presamples)
    peak_index,peak_time,peak_av,peak = gp.peak(data,presamples,time)
    rise = gp.risetime(data,peak_av,peak_index,rate)
    decay = gp.decaytime(data,peak_av,peak_index,rate)
    event = gp.silicon_event(decay)
    noise = gp.noise_rm(data,peak)
    #gp.graugh_save(num,data,time)

    data_column = [error,base,peak_av,peak_time,rise,decay,event,noise]
    data_array.append(data_column)
                #DataFrameの作成
df = pd.DataFrame(data_array,\
columns=["error","base","peak","peaktime","rise","decay","event",'noise'],\
index=path)
output_path = ('CH1/output')
gp.output(output_path,df)
print('end')
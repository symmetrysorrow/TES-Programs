import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import shutil
from natsort import natsorted
import libs.getpara as gp
import libs.fft_spectrum as sp
import json
import pprint
import sys
import re
import matplotlib.cm as cm


num = 8083
ch1 = gp.loadtxt("F:/実験/Nakano/修論保管用/両チャンネル_Bi吸収体_コリメータ実験/20210830/145mK_ch1-897uA_ch2-800uA_gain10_trig0.05_10kHz/CH1/rawdata/CH1_{}.dat".format(num))
ch2 = gp.loadtxt("F:/実験/Nakano/修論保管用/両チャンネル_Bi吸収体_コリメータ実験/20210830/145mK_ch1-897uA_ch2-800uA_gain10_trig0.05_10kHz/CH2/rawdata/CH2_{}.dat".format(num))

plt.plot(ch1)
plt.plot(ch2)
plt.show()
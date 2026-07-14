import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.fftpack as fft
from scipy.optimize import curve_fit
import getpara as gp
import glob
from natsort import natsorted
import os
import sys
import re
import json
import matplotlib.cm as cm


R_SH = 3.9e-3  # shant resistance


# one degree func
def func(x, a, b):
    return a * x + b


def main():
    ax = sys.argv
    path = input("path: ")
    os.chdir(path)

        
    files = natsorted(glob.glob("calibration/IV_*.txt"))
    
    print(files)
    I_bias = []
    V_out = []

    cnt = 0
    for i in files:
        temp = re.sub(r"\D", "", i)
        data = np.loadtxt(i)
        I_bias = data[0]
        V_out = data[1]
        plt.plot(
            I_bias,
            V_out,
            marker="o",
            linewidth=1,
            markersize=4,
            label=f"{temp} mK",
            color=cm.hsv((float(cnt)) / float(len(files))),
        )
        cnt += 1
    plt.title(f"I-V ")
    plt.xlabel("I_bias[uA]")
    plt.ylabel("V_out[V]")
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"calibration/IV.pdf")
    plt.show()

    cnt = 0
    for i in files:
        temp = re.sub(r"\D", "", i)
        data = np.loadtxt(i)
        I_bias = data[0]
        V_out = data[1]

        popt, cov = curve_fit(func, I_bias[:10], V_out[:10])
        if cnt == 0:
            eta = 1 / popt[0]
        print(eta)
        I_tes = eta * V_out
        I_sh = I_bias - I_tes
        V_tes = I_sh * R_SH
        R_tes = V_tes[1:] / I_tes[1:]
        R_tes = np.append(0.0, R_tes)

        plt.plot(
            I_bias,
            R_tes*1000,
            marker="o",
            linewidth=1,
            markersize=4,
            label=f"{temp} mK",
            color=cm.hsv((float(cnt)) / float(len(files))),
        )
        cnt += 1

    plt.title(f"I-R ")
    plt.xlabel("I_bias[uA]")
    plt.ylabel("R_tes[m$\Omega$]")
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"calibration/IR_calibration.pdf")
    plt.show()

    
if __name__ == "__main__":
    main()

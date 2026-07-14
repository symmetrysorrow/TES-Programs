
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
import json

#-----------パラメタの設定----------------------
R_SH = 3.9e-3

# set offset to zero
def offset(data):

    data = data - data[0]

    if np.mean(data[:10]) < 0:
        data = data * -1
    return data

# one degree func
def func(x,a,b):
    return a * x + b

def main():

    if not os.path.exists('calibration'):
        os.mkdir('calibration')

    path = input('path: ')
    os.chdir(path)
    temp = input('temperature [mK]: ')


    files = natsorted(glob.glob(f"{temp}mK/*.dat"))
    I_bias = [] 
    V_out = []
    for i in files:
        data = gp.loadtxt(i)
        V_out.append(np.mean(data))
        I = os.path.splitext(os.path.basename(i))[0][10:-2]
        I_bias.append(int(I))
    
    V_out = np.array(V_out)
    I_bias = np.array(I_bias)
    V_out = offset(V_out)


    
    print(I_bias[-10:-1])

    calib = '0'
    while calib == '0':
        # Fitting and get eta
        popt, cov = curve_fit(func,I_bias[:10],V_out[:10])
        eta = 1 / popt[0]
        print(f"eta (uA/V): {eta}")

        popt2, cov2 = curve_fit(func,I_bias[-10:-1],V_out[-10:-1])

        x_fit = np.arange(0,I_bias[-1],1)

        y_fit = func(x_fit,*popt2)
        y_fit -= y_fit[0]

        plt.plot(I_bias,V_out,marker = "o",c = "red",linewidth = 1,markersize = 6)
        plt.plot(x_fit,func(x_fit,*popt))
        plt.plot(x_fit,y_fit)
        plt.title(f'I-V at {temp}')
        plt.xlabel("I_bias[uA]")
        plt.ylabel("V_out[V]")
        plt.grid(True)
        plt.show()
        
        start = int(input('start:'))
        stop = int(input('stop:'))

        mode = int(input("super:[0], normal[1], delete[2], else[3]:"))
        print(len(V_out))
        if mode == 0:
            diff =  func(start,*popt) - V_out[np.where(I_bias==start)][0]
            V_out[np.where(I_bias==start)[0][0]:np.where(I_bias==stop)[0][0] + 1 ] += diff
        elif mode == 1:
            
            diff =  y_fit[stop-1] - V_out[np.where(I_bias==stop)][0]
            V_out[np.where(I_bias==start)[0][0]:np.where(I_bias==stop)[0][0] + 1 ] += diff
        elif mode == 2:
            print(f"delete")
            V_out = np.delete(V_out,slice(np.where(I_bias==start)[0][0],np.where(I_bias==stop)[0][0] + 1 ))
            I_bias = np.delete(I_bias,slice(np.where(I_bias==start)[0][0],np.where(I_bias==stop)[0][0] + 1 )) 
        elif mode == 3:
            
            diff =  V_out[np.where(I_bias==start)[0]-1] - V_out[np.where(I_bias==start)][0]
            print(f"expect: {diff}")
            change = float(input('change:'))
            V_out[np.where(I_bias==start)[0][0]:np.where(I_bias==stop)[0][0] + 1 ] += change
        else:
            exit()
        

        plt.plot(I_bias,V_out,marker = "o",c = "red",linewidth = 1,markersize = 6)
        plt.plot(x_fit,func(x_fit,*popt))
        plt.plot(x_fit,y_fit)
        plt.title(f'I-V at {temp}')
        plt.xlabel("I_bias[uA]")
        plt.ylabel("V_out[V]")
        plt.grid(True)
        plt.show()
        
        calib =  input('continue? yes[0],no[1]:')

    I_tes = eta * V_out
    I_sh = I_bias - I_tes
    V_tes = I_sh * R_SH
    R_tes = V_tes[1:] / I_tes[1:]
    R_tes =  np.append(0.0,R_tes)
    

    IVR = [I_bias,V_out,R_tes]
    np.savetxt(f"calibration/IV_{temp}.txt",IVR)

    plt.plot(I_bias,V_out,marker = "o",c = "red",linewidth = 1,markersize = 6)
    plt.title(f'I-V at {temp}')
    plt.xlabel("I_bias[uA]")
    plt.ylabel("V_out[V]")
    plt.grid(True)
    plt.savefig(f"calibration/IV_{temp}.pdf")
    plt.show()

    # I-R graugh
    plt.plot(I_bias,R_tes,marker = "o",c = "red",linewidth = 1,markersize = 6)
    plt.title(f'I-R at {temp}')
    plt.xlabel("I_bias[uA]",fontsize = 16)
    plt.ylabel("R_tes[$\Omega$]",fontsize = 16)
    plt.grid(True)
    plt.savefig(f"calibration/IR_{temp}.pdf")
    plt.cla()


if __name__ == '__main__':
    main()

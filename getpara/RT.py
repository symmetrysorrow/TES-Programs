
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import scipy.fftpack as fft
from scipy.optimize import curve_fit
import libs.getpara as gp
import glob
from natsort import natsorted
import os
import sys
import re
import shutil

R_SH = 3.9	# shant resistance[mohm]
cmap = cm.get_cmap("hsv")

# one degree func
def func(x,a,b):
    return a * x + b



def main():
    ax = sys.argv
    ax.pop(0)
    
    os.chdir(ax[0])


    if not os.path.exists('output'):
        os.mkdir('output')
    if not os.path.exists('rawdata'):
        files = natsorted(glob.glob("*.dat"))
        os.mkdir('rawdata')
        for i in files:
            shutil.move(i,"rawdata")
    else:
        files = natsorted(glob.glob("rawdata/*.dat"))

    
    I_bias = []
    V_out = []
    T = []
    for i in files:

        data = gp.loadtxt(i)
        name = os.path.splitext(os.path.basename(i))[0]
        V_out.append(np.mean(data))
        t = re.sub(r"\D", "", name.split('_')[1])
        i_bias = re.sub(r"\D", "", name.split('_')[2])
        T.append(int(t))
        I_bias.append(float(i_bias))
    
    low_temp = T.count(np.min(T))

    # Fitting and get eta
    popt, cov = curve_fit(func,I_bias[:low_temp],V_out[:low_temp])
    eta = 1 / popt[0]
    print(f"eta (uA/V): {eta}")

    T = natsorted(set(T[low_temp:]))
    I_bias_2 = natsorted(set(I_bias[low_temp:]))


    V_out = []
    
    for i in I_bias_2:
        V = []
        for t in T:
            data = gp.loadtxt(f'rawdata/CH1_{t}mK_{int(i)}uA.dat')
            V.append(np.mean(data))
        V_out.append(V)

    V_out = np.array(V_out)
    cnt = 0
    for i in V_out:
        if cnt > 0:
            V_out_base  = i - V_out[0]
            R = R_SH*(I_bias_2[cnt]/(eta*V_out_base)-1)
            plt.title('R-T')
            plt.plot(T,R,marker = 'o',linewidth = 1,label =f'{I_bias_2[cnt]}uA',markersize = 4)
            plt.xlabel('Temperature[mK]',fontsize = 16)
            plt.ylabel('Resistance[m$\Omega$]',fontsize = 16)
            plt.grid(True)
            plt.legend(loc ='best',fancybox = True,shadow = True)
            np.savetxt(f'output/{int(I_bias_2[cnt])}uA.txt',[T,R])
        cnt += 1

    plt.savefig(f"output/RT.png")
    plt.show()



    """

        V_out = offset(V_out)

        # Fitting and get eta
        popt, cov = curve_fit(func,I_bias[:10],V_out[:10])
        eta = 1 / popt[0]
        print(f"eta (uA/V): {eta}")

        
        

        I_tes = eta * V_out
        I_sh = I_bias - I_tes
        V_tes = I_sh * R_SH
        R_tes = V_tes[1:] / I_tes[1:]
        R_tes =  np.append(0.0,R_tes)

        
        np.savetxt(f"rawdata/IV_{t}.txt",[I_bias,V_out])

        # I-V graugh
        plt.plot(I_bias,V_out,marker = "o",c = "red",linewidth = 1,markersize = 6)
        plt.title('I-V at {t}')
        plt.xlabel("I_bias[uA]")
        plt.ylabel("V_out[V]")
        plt.grid(True)
        plt.savefig(f"output/IV_{t}.png")
        #plt.show()
        plt.cla()
        
        # I-R graugh
        plt.plot(I_bias,R_tes,marker = "o",c = "red",linewidth = 1,markersize = 6)
        plt.title('I-R at {t}')
        plt.xlabel("I_bias[uA]",fontsize = 16)
        plt.ylabel("R_tes[$\Omega$]",fontsize = 16)
        plt.grid(True)
        plt.savefig(f"output/IR_{t}.png")
        plt.cla()

    for t in temps:
        files = natsorted(glob.glob(os.path.join(t,"*.dat")))
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
        plt.plot(I_bias,V_out,marker = "o",linewidth = 1,markersize = 6,label=f'{t}mK')

    plt.title('I-V')
    plt.xlabel("I_bias[uA]")
    plt.ylabel("V_out[V]")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"output/IV_matome.png")
    plt.show()
    """
    
        

if __name__ == '__main__':
    main()
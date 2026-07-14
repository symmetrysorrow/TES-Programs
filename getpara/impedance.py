

import glob
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import libs.getpara as gp
import sys
import re
import natsort

RB = 1.0e4
RSH = 3.9e-3

eta = 98 #[uA/V]
rate = 1e6
samples = 1e6
time = gp.data_time(rate,samples)
p0 = [0, 0.47, 1, 0]

def sin(x,a,b,c,d):
    return a + b * np.sin(2 * np.pi / c * x + d)

def main():
    ax = sys.argv
    ax.pop(0)
    
    os.chdir(ax[0])
    dir = glob.glob("*")
    for i in dir:
        txt = os.path.splitext(i)[0].split("_")
        T_c = re.sub(r"\D", "", txt[0])
        I_bias = re.sub(r"\D", "", txt[1])
        freqs =  natsort.natsorted(glob.glob(f"{i}/*Hz.txt"))[1:]
        Re = []
        Im = []

        for freq in freqs:
            print(freq)
            data = np.loadtxt(freq).T    
            V_out1 = data[0]
            V_out2 = data[1]
            fq = int(re.sub(r"\D", "", os.path.basename(freq)))
            length = int(samples/fq)     
            x_fit = time[:length]
            y_fit_1 = V_out1[:length]
            y_fit_2 = V_out2[:length]
			
            try:   
            	popt1,pcov = curve_fit(sin,x_fit,y_fit_1,p0,maxfev = 100000)
            	popt2,pcov = curve_fit(sin,x_fit,y_fit_2,p0,maxfev = 100000)
            except:
                popt1 = [0,0,0,0]
                popt2 = [0,0,0,0]
            
            
            fit1 = sin(time,*popt1)
            fit2 = sin(time,*popt2)
            amp = popt1[1]/popt2[1]
            phase = popt1[3]-popt2[3]
            if amp < 0:
                amp = np.abs(amp)
                phase += np.pi
            
            Z_a =  RSH/RB*amp/eta/1e-6
            Z_phi = phase
            Re_Z = Z_a*np.cos(Z_phi)
            Im_Z = Z_a*np.sin(Z_phi)
            Re.append(Re_Z)
            Im.append(Im_Z)
            
        plt.figure(figsize=(7, 7))
        plt.plot(Re,Im,'o',color='blue',markersize=5)
        plt.xlim(-0.06,0.06)
        plt.ylim(-0.06,0.06)
        plt.xlabel('ReZ[$\Omega$]')
        plt.ylabel('ImZ[$\Omega$]')
        plt.grid()
        plt.savefig("ReZ.pdf",format='pdf')
        plt.show()
        
		



        
        
            
            
            
			
    


if __name__ == '__main__':
    main()
    
	
    
            
            
            
            
            

            
            
        
        

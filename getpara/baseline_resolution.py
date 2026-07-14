

import numpy as np
from scipy.optimize import curve_fit
import scipy.optimize as opt
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import libs.getpara as gp
from natsort import natsorted
import glob



p0 = [800,0,0.05]
pulse_fmin = -1
pulse_fmax = 1
bins = 500


def fit(x,a,b,c):
	return a*x**2+b*x+c

def gausse2(x,A,mu,sigma):
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))

def FWHW(sigma):
	return 2*sigma*(2*np.log(2))**(1/2)
def gauss(p,x,y):
	residual=y-(p[0]*np.exp(-(x-p[2])**2/(2*p[1]**2)))
	return residual

def main():
	set = gp.loadJson()
	os.chdir(set['Config']["path"])
	df = pd.read_csv((f'CH{set["Config"]["channel"]}_pulse/output/output.csv'),index_col=0)
	rate,samples,ch = set["Config"]["rate"],set["Config"]["samples"],set["Config"]["channel"]
	time = gp.data_time(rate,samples)
	template = np.loadtxt(f'CH{set["Config"]["channel"]}_pulse/output/opt_template.txt')

	
	cali = np.loadtxt(f'CH{set["Config"]["channel"]}_pulse/output/calibration_curve_last.txt')

	noise = []
	"""
	with open(f"CH{ch}_noize/random_noise.txt",'r',encoding='latin-1') as f:
		for row in f.read().splitlines():
			noise.append(row)
	"""
        	
	
	path = natsorted(glob.glob(f'CH{ch}_noize/rawdata/CH{ch}_*.dat'))
	
	A = []

	for num in path:
		print(os.path.basename(num))
		try:
			data = gp.loadbi(num)
			
			if np.average(data) < 1 and np.average(data) > -1:
				A.append(np.sum(data*template))
				
		except FileNotFoundError:
			print("error")
			continue
	
	print(len(A))
	pulseheight_f = np.array(A)

	pulseheight_fc = fit(pulseheight_f,*cali)


	plt.hist(pulseheight_fc, bins =bins)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.xlabel('Energy[keV]',fontsize=16, fontname='serif')
	plt.ylabel('Count',fontsize=16, fontname='serif')
	plt.savefig("base_spectrum.pdf",format = 'pdf')
	plt.show()

	
	hist, bin_edges = np.histogram(pulseheight_fc, bins=bins)
	bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
	p0 = [0]*3
	p0[0]= np.max(hist)										#hist内の最大値を返す
	p0[1] = (pulse_fmax - pulse_fmin)/6
	p0[2] = bin_centers[np.argmax(hist)]					#argmaxは最大値の添字を返す
	print ('initial v	alue',p0)

	
	

	output=opt.leastsq(gauss,p0,args=(bin_centers,hist))
	a = output[0][0]
	sigma = output[0][1]
	mu = output[0][2]
	print  (a, sigma, mu)
	x = np.arange(pulse_fmin,pulse_fmax,0.001)
	#y = a/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))
	y = a*np.exp(-(x-mu)**2/(2*sigma**2))
	E = 2*sigma*np.sqrt(2*np.log(2))

	print ('Energy resolution :',E,'keV')

	plt.hist(pulseheight_fc, bins =bins)
	#plt.ylim(0,300)
	plt.plot(x, y, color='red',linewidth=1,linestyle='-')
	plt.xlabel('Energy[keV]')
	plt.ylabel('Count[-]')
	plt.title('dE= '+'{:.3f}'.format(E)+'keV')
	plt.savefig("base_spectrum_fit.pdf",format = 'pdf')
	plt.show()

	    




	

if __name__ == "__main__":
	main()


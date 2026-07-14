

#----------- モジュールのインポート --------------
import math
import glob
import numpy as np
import scipy.optimize as opt
import random
from array import array
import os
import os.path
import matplotlib.pyplot as plt
import getpara as gp
import pandas as pd
from scipy.optimize import curve_fit
import sys
import json
import plt_config
12.16,12.28
#-------------- ガウスフィッティング ----------------
BINS = 128

# K_x  = -77keV

def gausse(x,A,mu,sigma):
	return A*np.exp(-(x-mu)**2/(2.0*sigma**2))

def main():
	ax = sys.argv
	ax.pop(0)
	if len(ax) < 1:
		print("add parameter you want to see as histgrum!")
		sys.exit()
	para = ax[0]
	print(ax)

	set = gp.loadJson()
	ch,path = set["Config"]["channel"],set["Config"]["path"]
	output = f'CH{set["Config"]["channel"]}_pulse/output/{set["Config"]["output"]}'
	
	os.chdir(path)
	df = pd.read_csv((f'{output}/output.csv'),index_col=0)
	df_sel = gp.select_condition(df,set['select'])

	if os.path.exists(f'{output}/{set["select"]["output"]}/gfit.json'):
		with open(f'{output}/{set["select"]["output"]}/gfit.json') as f:
			jsn = json.load(f)
	else:
		jsn = {}

		
	plt.hist(df_sel[para], bins = BINS)
	plt.xlabel('Channel',fontsize = 14)
	plt.ylabel('Count',fontsize = 14)
	plt.tick_params(axis='both', which='both', direction='in',\
					bottom=True, top=True, left=True, right=True)
	plt.grid(True, which='major', color='black', linestyle='-', linewidth=0.2)
	plt.show()

	if '-a' not in ax:
		energy = float(input('input Energy(keV): '))
		min = float(input('min: '))
		max = float(input('max: '))
		
		bins = int(input("bins: "))
		pulseheight = df_sel[(df_sel[para] > min)&(df_sel[para] < max)][para]
		hist, bin_edges = np.histogram(pulseheight, bins=bins)
		p0 = [np.max(hist),(max-min)/6,bin_edges[np.argmax(hist)]]
		
		params,cov = curve_fit(gausse,bin_edges[:-1],hist,p0=p0,maxfev=10000)
		a,mu,sigma=params[0],params[1],np.abs(params[2])
		
		x_fit = np.linspace(min*0.8,max*1.2,1000)
		y_fit = gausse(x_fit,*params)
		fwhw = gp.FWHW(sigma)
		dE = energy*fwhw/mu
		
		fit_para = {energy:
				{
				"A":a,
				"mu":mu,
				"sigma":sigma,
				"fwhw":fwhw,
				"dE":dE,
				"min":min,
				"max":max,
				"bins":bins
				}
			}
		
		jsn.update(fit_para)
		
		
		print(f'\nA: {a}\nmu: {mu}\nsigma: {sigma}')
		print(f'\nFWHW: {fwhw}\nEnergy resoltion: {dE} keV')

		plt.hist(pulseheight, bins =bins,alpha=0.8)
		plt.plot(x_fit, y_fit, color='red',linewidth=1,linestyle='-')

		#plt.title(f'dE = {dE:.3f}keV @{int(energy)} keV')
		plt.xlabel('Channel',fontsize = 14)
		plt.ylabel('Count',fontsize = 14)
		plt.tick_params(axis='both', which='both', direction='in',\
						bottom=True, top=True, left=True, right=True)
		plt.grid(True, which='major', color='black', linestyle='-', linewidth=0.2)
		plt.tight_layout()
		plt.savefig(f'{output}/{set["select"]["output"]}/E_{int(energy)}keV.png')
		plt.show()

	else:
		pulseheight = df_sel[para]
		hist, bin_edges = np.histogram(pulseheight, bins=BINS)
		x_fit = np.linspace(bin_edges[0],bin_edges[-1],1000)
		plt.hist(df_sel[para], bins = BINS,alpha=0.8)
		for i in jsn:
			print(i)
			popt = jsn[i]
			p0 = [popt["A"],popt["mu"],popt["sigma"]]
	
			min = np.argmin(bin_edges<popt["min"])
			max = np.argmin(bin_edges<popt['max'])
			print(min,max)
			params,cov = curve_fit(gausse,bin_edges[:-1][min:max],hist[min:max],p0=p0,maxfev=10000)
			y_fit = gausse(x_fit,*params)
			plt.plot(x_fit, y_fit, color='tab:red',linewidth=1,linestyle='-')

		#plt.title(f'dE = {dE:.3f}keV @{int(energy)} keV')
		plt.xlabel('Channel')
		plt.ylabel('Count')
		plt.tick_params(axis='both', which='both', direction='in',\
						bottom=True, top=True, left=True, right=True)
		plt.grid(True, which='major', color='black', linestyle='-', linewidth=0.2)
		plt.tight_layout()
		plt.savefig(f'{output}/{set["select"]["output"]}/histgrum_fit.png')
		plt.show()



	output_jsn = json.dumps(jsn,indent=4)
	with open(f'{output}/{set["select"]["output"]}/gfit.json', 'w') as file:
		file.write(output_jsn)



if __name__ == '__main__':
	main()

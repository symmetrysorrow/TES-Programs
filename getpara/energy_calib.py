# -*- coding: utf-8 -*-
#Calibration
# 11-6D1

import glob
import numpy as np
import scipy.optimize as opt
import scipy.fftpack as fft
from array import array
import os
import os.path
import time
import matplotlib.pyplot as plt
import libs.getpara as gp
import pandas as pd
from scipy.optimize import curve_fit
import json
import sys
import plt_config

p0 = [0,0]


#----------- パルスハイトの読み込み ----------------
def main():

	ax = sys.argv
	ax.pop(0)
	if len(ax) != 1:
		print("add parameter you want to see as histgrum!")
		sys.exit()
	para = ax[0]

	set = gp.loadJson()
	path,ch = set["Config"]["path"],set["Config"]["channel"]
	output = f'CH{ch}_pulse/output/{set["Config"]["output"]}'
	output_f = f'{output}/{set["select"]["output"]}'
	os.chdir(path)
	df = pd.read_csv((f'{output}/output.csv'),index_col=0)
	with open(f'{output}/{set["select"]["output"]}/gfit.json') as f:
			fit_para = json.load(f)
	
	height = []#[0.0]
	energy = [] #[0.0]
	for i in fit_para:
		energy.append(float(i))
		height.append(fit_para[i]["mu"])
	
	
	popt,pcov = curve_fit(gp.multi_func,energy,height,p0)
	x_fit = np.arange(0,np.max(energy)*1.2,0.1)
	y_fit = gp.multi_func(x_fit,*tuple(popt))

	if len(p0) == 3:
		fit_para = {"a":popt[0],"b":popt[1],"c":popt[2]}
	elif len(p0) == 2:
		fit_para = {"a":popt[1],"b":popt[0]}
	
	output_jsn = json.dumps(fit_para,indent=4)
	with open(f'{output}/{set["select"]["output"]}/energy_calib.json', 'w') as file:
		file.write(output_jsn)


	plt.plot(energy,height,'o',markersize=6)
	plt.plot(x_fit,y_fit,color='red',linewidth=1,linestyle='-')
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.ylabel('Pulseheight[V]')
	plt.xlabel('Energy[keV]')
	plt.grid()
	plt.tight_layout()
	#plt.legend()
	plt.savefig(f"{output}/{set['select']['output']}/energy_calibration.png",dpi = 350)
	plt.show()
	
	y = float(input("height: "))
	x = (y-fit_para['b'])/fit_para['a']

	print(x)

	




#----------------- 2次曲線でfitting --------------------
'''
p0 = [0.2,0.0,0.0]
def fit(p,x,y):
	residual=y-(p[0]*x**2+p[1]*x+p[2])
	return residual

output=opt.leastsq(fit,p0,args=(ch,energy))
a = output[0][0]
b = output[0][1]
c = output[0][2]
print ('a=',a,'b=',b,'c=',c)
x = np.arange(0,100,0.01)
y = a*x**2+b*x+c
'''


'''
#直線で校正
a = (energy[1]-energy[0])/(ch[1]-ch[0])
x = np.arange(0,18,0.01)
y = a*x
'''

'''
#校正曲線の表示
plt.xlim(0,700)#ここの単位はkeV
plt.ylim(0,10.0)#ここの単位はV（パルス高）
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Pulseheight[V]',fontsize=14, fontname='serif')
plt.xlabel('Energy[keV]',fontsize=14, fontname='serif')
plt.plot(energy,ch,'ro',color='blue',markersize=6)
plt.plot(y,x,color='red',linewidth=1,linestyle='-')
#plt.grid()
plt.savefig("cal_curve.pdf",format = 'pdf')
plt.show()
'''

if __name__ == '__main__':
	main()


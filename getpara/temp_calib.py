
# Matsumi
# Temperature Calibration

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import libs.getpara as gp
import json
import tqdm

# Give textfile name as argment when you run this python file.
# ex)
# python .\temp_com_curve.py baseline.txt pulseheight_opt.txt


# pO is initialize parameter to fitting polynomial.
# Make sure the (length of p0 - 1) is the degree of the polynomial.
p0=[0.01,0.01,0.01,0.01,0.01,0.01]


def func(X, *params):
    Y = np.zeros_like(X)
    for i, param in enumerate(params):
        Y = Y + np.array(param * X ** i)
    return Y


def Calibration(x,params):
	array = np.zeros(len(params))
	for i,param in enumerate(params):
		term = param * x ** i
		array[i] = term
		sum = np.sum(array)
	return sum
		

def PlotFitting(x,y,x_fit,y_fit,path):
	plt.plot(x,y,'o',color='blue',markersize=3,label='a')
	plt.plot(x_fit,y_fit,color='red',linewidth=1.0,linestyle='-')
	plt.xlabel('baseline [V]',fontsize = 16)
	plt.ylabel('pulseheight [V]',fontsize = 16)
	plt.grid()
	plt.savefig(f'{path}/fitting.png',dpi=350)
	plt.show()
	plt.cla()

def PlotCalibration(x1,y1,x2,y2,path):
	plt.plot(x1,y1,'o',color='tab:blue',markersize=0.7,label='a')
	plt.plot(x2,y2,'o',color='tab:red',markersize=0.7,label='a')
	plt.xlabel('baseline [V]',fontsize = 16)
	plt.ylabel('pulseheight [V]',fontsize = 16)
	plt.grid()
	plt.savefig(f'{path}/temp_calibration.png',dpi=350)
	plt.show()
	plt.cla()



def main():
	
	set = gp.loadJson()
	ch,path = set["Config"]["channel"],set["Config"]["path"]
	output = f'CH{set["Config"]["channel"]}_pulse/output/{set["Config"]["output"]}'
	output_sel = f'CH{set["Config"]["channel"]}_pulse/output/{set["Config"]["output"]}/{set["select"]["output"]}'

	if (not 'base->' in set['select'])or(not 'base-<' in set['select']):
		base_min = float(input("base > : "))
		base_max = float(input("base < : "))
		set['select']['base->'] = base_min
		set['select']['base-<'] = base_max
		jsn = json.dumps(set,indent=4)
		with open("setting.json", 'w') as file:
			file.write(jsn)
	
	os.chdir(path) 
	df = pd.read_csv((f'{output}/output.csv'),index_col=0)
	x,y = gp.extruct(df,"base","height_opt")


	df_sel = gp.select_condition(df,set)
	
	picked = gp.pickSamples(df_sel,"base","height_opt") 
	baseline = df_sel.loc[picked,"base"]
	pulseheight = df_sel.loc[picked,"height_opt"]
	
		
	print(f"Selected {len(picked)} Samples.")
	
	#Fitting
	popt,pcov = curve_fit(func,baseline,pulseheight,p0)
	x_fit = np.linspace(set['select']['base->'],set['select']['base-<'],100000)
	fitted = func(x_fit,*tuple(popt))
	PlotFitting(baseline,pulseheight,x_fit,fitted,output_sel)        

	#Calibration
	st = np.mean(pulseheight)
	pulseheight_cal = np.zeros(len(df_sel))

	for index,row in tqdm.tqdm(df_sel.iterrows()):
		df.at[index,"height_opt_temp"] = row['height_opt']/Calibration(row['base'],popt)*st

	PlotCalibration(df["base"],df["height_opt_temp"],df.loc[picked,"base"],df.loc[picked,"height_opt_temp"],output_sel)
	df.to_csv(f'{output}/output.csv')

if __name__ == '__main__':
	main()

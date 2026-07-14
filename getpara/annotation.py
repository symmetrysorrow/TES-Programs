import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import getpara as gp
import sys
import shutil
import pprint
import json
import glob
import re
from tkinter import filedialog
import tkinter as tk

BUF = 100

def main():
	
	setting = gp.loadJson()
	config = setting["Config"]

	os.chdir(config["path"])
	
	root  = tk.Tk()

	root.withdraw()

	output = f'CH{config["channel"]}_pulse/output/{config["output"]}'
	df = pd.read_csv((f'{output}/output.csv'),index_col=0)
	rate,samples,presamples,threshold,ch = config["rate"],config["samples"],config["presamples"],config["threshold"],config["channel"]
	time = gp.data_time(rate,samples)
	df_sel = gp.select_condition(df,setting["select"])

	start = 0
	i = 1
	while i != 0:
		os.makedirs(f"{output}/img",exist_ok=True)
		buf = df_sel.iloc[start:start+BUF].index.values
		print("loading...")
		for num in buf:
			path = f'CH{ch}_pulse/rawdata/CH{ch}_{num}.dat'
			analysis = setting["main"]
			data = gp.loadbi(path,config["type"])
			if analysis['cutoff'] > 0:
				data = gp.BesselFilter(data,config["rate"],fs = analysis['cutoff'])
			base,data = gp.baseline(data,presamples,1000,500)
			
			plt.plot(time,data)
			plt.title(num)
			plt.xlabel("time (s)")
			plt.ylabel("volt (V)")
			plt.savefig(f'{output}/img/{num}.png')
			plt.cla()
		
		fle = filedialog.askopenfilenames(initialdir=f"{output}/img")
		for f in fle:
			num =  int(re.findall(r'\d+', os.path.basename(f))[0])    
			df.at[num,"error"] = 0

		shutil.rmtree(f"{output}/img")
		try:
			i = int(input('finish? [0]'))
		except:
			start  += BUF
			continue
		

	df.to_csv(f'{output}/output.csv')
   
	
		
		
	

	"""
	num = 1
	while num != 0:
				
				try:
					picked.remove(f'CH{ch}_pulse/rawdata\\CH{ch}_{num}.dat')
					os.remove(f'{output_f}/img/CH{ch}_{num}.png')
					np.savetxt(f'{output_f}/selected_index.txt',picked,fmt="%s")
					index = f"CH{ch}_pulse/rawdata\CH{ch}_{num}.dat"
					df.at[index,"quality"] = 0
				except:
					print("Not exist file")
	"""
	

	"""
	# delete noise data
	fle = filedialog.askopenfilenames(filetypes=[('画像ファイル','*.png')],initialdir=f"{output_f}/img")
	for f in fle:
		num =  re.findall(r'\d+', os.path.basename(f))[1]
		picked.remove(f'CH{ch}_pulse/rawdata\\CH{ch}_{num}.dat')
		os.remove(f'{output_f}/img/CH{ch}_{num}.png')
		np.savetxt(f'{output_f}/selected_index.txt',picked,fmt="%s")
		index = f"CH{ch}_pulse/rawdata\CH{ch}_{num}.dat"
		df.at[index,"quality"] = 0
	"""

if __name__=='__main__':
	main()

	
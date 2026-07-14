import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import libs.plt_config
import os
import libs.getpara as gp

import shutil
import pprint
import json
import re
from tkinter import filedialog 
import platform
import tkinter as tk
import tqdm


def overlap(df_0,df_1):
	#index_0,index_1 = df_0.index.values,df_1.index.values
	##ch1 = re.findall(r'\d+', index_1[0])[0]
	#num_0 = [re.findall(r'\d+', i)[2] for i in index_0]
	#num_1 = [re.findall(r'\d+', i)[2] for i in index_1]
	df_comp_0 = df_0[df_0.index.isin(df_1.index)]
	df_comp_1 = df_1[df_1.index.isin(df_0.index)]
	
	
	return df_comp_0,df_comp_1


def df_number(df):
	index= df.index.values
	num = [re.findall(r'\d+', i)[2] for i in index]
	df['number'] = num
	return df


def main():
	ch0 = 0
	ch1 = 1
	para = "height"
	setting = gp.loadJson()
	config = setting["Config"]
	rate,samples,presamples,threshold,ch = config["rate"],config["samples"],config["presamples"],config["threshold"],config["channel"]

	os.chdir(config["path"])

	output_0 = f'CH{ch0}_pulse/output/{config["output"]}'
	output_1 = f'CH{ch1}_pulse/output/{config["output"]}'
	
	df_0 =  pd.read_csv((f'{output_0}/output.csv'),index_col=0)
	df_1 =  pd.read_csv((f'{output_1}/output.csv'),index_col=0)


	df_0_clear = gp.select_condition(df_0,setting["select"])
	df_1_clear = gp.select_condition(df_1,setting["select"])


	df_0_over,df_1_over = overlap(df_0_clear,df_1_clear)
	print(df_0_clear)
	x,y = df_0_over[para],df_1_over[para]

	fle = filedialog.askopenfilename()

	selected = gp.loadIndex(fle)
	x_sel,y_sel = df_0_over[df_0_over.index.isin(selected)][para],df_1_over[df_1_over.index.isin(selected)][para]

	plt.scatter(x,y,s=2,alpha=0.7)
	plt.scatter(x_sel,y_sel,s=4)
	plt.xlabel(f'channel {ch0} [V]')
	plt.ylabel(f'channel {ch1} [V]')
	plt.grid()
	#plt.savefig(f'{output_select}/pulse_height_{ch0}_{ch1}selected_index.png')
	plt.show()
	plt.cla()
	
if __name__ == "__main__":
	main()
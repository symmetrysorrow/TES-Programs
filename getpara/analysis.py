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
import tkinter as tk
from tkinter import filedialog
import tqdm
import platform
import plt_config
import natsort
import glob
import matplotlib.cm as cm
from scipy.optimize import curve_fit

arrival_threshold = 0.01

def arrival_time(data,point,x,w):
	
	fit_data = data[point-x:point-x+w]
	fit_range = np.arange(point-x,point-x+w,1)
	
	popt,ccpov = curve_fit(gp.multi_func,fit_range,fit_data,p0=[0,0])

	#arrival = -popt[1]/popt[0]
	return popt#arrival

def arrival_time(data,presamples,x,w):
	diff = np.argmax(gp.diff(data[presamples:])) + presamples
	fit_data = data[diff-x:diff-x+w]
	fit_range = np.arange(diff-x,diff-x+w,1)
	popt,ccpov = curve_fit(gp.multi_func,fit_range,fit_data,p0=[0,0])
	arrival = -popt[0]/popt[1]
	return arrival,popt


def main():
	setting = gp.loadJson()
	config = setting['Config']
	path = config['path']
	output = config['output']
	ch = config['channel']
	para = setting['main']
	time = gp.data_time(rate=config['rate'],samples=config['samples'])
	presamples = int(config['presamples'])
	rate = config['rate']
	os.chdir(path)

	'''
	df_0 = pd.read_csv(f'CH0_pulse/output/{output}/output.csv')
	df_1 = pd.read_csv(f'CH1_pulse/output/{output}/output.csv')

	df_0_clear = gp.select_condition(df_0, setting["select"])
	df_1_clear = gp.select_condition(df_1, setting["select"])

	df_0_over, df_1_over = gp.overlap(df_0_clear, df_1_clear)
	x, y = df_0_over['height'].values, df_1_over['height'].values
	'''
	
	blocks_path_0 = natsort.natsorted(glob.glob(f'CH0_pulse/output/{output}/blocks/average_pulse_*.txt'))
	blocks_path_1 = natsort.natsorted(glob.glob(f'CH1_pulse/output/{output}/blocks/average_pulse_*.txt'))

	block = np.arange(1,len(blocks_path_0)+1,1)
	
	k = 1
	arrival_diffs = []

	

	for i,j in zip(blocks_path_0,blocks_path_1):
		ch_0 = np.loadtxt(i)
		ch_1 = np.loadtxt(j)
		height_0 = gp.peak(ch_0,config["presamples"],para['peak_max'],para['peak_x'],para['peak_w'])
		height_1 = gp.peak(ch_1,config["presamples"],para['peak_max'],para['peak_x'],para['peak_w'])
		ch_1 = ch_1 *1.6967260839791022

		ratio = height_0[1]/height_1[1]
		
		
		
		print(f'block: {k}')
		print(f"height ratio (ch0/ch1): {ratio}")
		#print(arrival_diff/config['rate'])
	
		
		diff_0 = np.argmax(gp.diff(ch_0[presamples:])) + presamples
		diff_1 = np.argmax(gp.diff(ch_1[presamples:])) + presamples

		
		x_fit = np.arange(presamples-10,presamples+30,1)
		arrival_0,popt_0 = arrival_time(ch_0,diff_0,2,5)
		arrival_1,popt_1 = arrival_time(ch_1,diff_1,2,5)
		y_fit_0 = gp.multi_func(x_fit,*popt_0)
		y_fit_1 = gp.multi_func(x_fit,*popt_1)
		arrival_0 = arrival_0/rate
		arrival_1 = arrival_1/rate
		'''
		arrival_0 = gp.arrival_time_threshold(ch_0,0.01)
		arrival_1 = gp.arrival_time_threshold(ch_1,0.01)
		'''

		arrival_diff = (arrival_0 - arrival_1)
		arrival_diffs.append(arrival_diff/config['rate'])

		
		plt.plot(time*1e3,ch_0,markersize = 2.0,label = 'CH0')
		plt.plot(time*1e3,ch_1,markersize = 1.0,label = 'CH1')
		
		'''
		plt.plot(time*1e3,np.linspace(0.01,0.01,config['samples']),"--",c='black')
		plt.scatter(time[arrival_0]*1e3,ch_0[arrival_0],c='red',label='arrival time')
		plt.scatter(time[arrival_1]*1e3,ch_1[arrival_1],c='red',label='arrival time')
		'''
		plt.scatter(time[diff_0-2:diff_0-2+5]*1e3,ch_0[diff_0-2:diff_0-2+5],c = 'tab:green',label='fitting sample')
		plt.scatter(time[diff_0]*1e3,ch_0[diff_0],c='black',label = 'diff max')
		plt.scatter(arrival_0*1e3,0,c='red',label='arrival time')
		plt.plot(x_fit/rate*1e3,y_fit_0,'--',label = 'fitting')

		plt.scatter(time[diff_1-2:diff_1-2+5]*1e3,ch_1[diff_1-2:diff_1-2+5],c = 'tab:green',label='fitting sample')
		plt.scatter(time[diff_1]*1e3,ch_1[diff_1],c='black',label = 'diff max')
		plt.scatter(arrival_1*1e3,0,c='red',label='arrival time')
		plt.plot(x_fit/rate*1e3,y_fit_1,'--',label = 'fitting')
		
		
		plt.title(f"block{k}")
		plt.xlabel("time [ms]")
		plt.ylabel("volt [V]")
		gp.graugh_condition(setting["graugh"])
		plt.grid()
		#plt.legend(fontsize=12)
		plt.tight_layout()
		plt.savefig(f'CH0_pulse/output/{output}/blocks/block{k}.png')
		#plt.show()
		plt.cla()
		k+=1
	
	arrival_diffs = np.array(arrival_diffs)*1e3
	plt.scatter(block,arrival_diffs)
	plt.ylabel("arrival time [ms]")
	plt.xlabel("block")
	#plt.ylim(0-1,setting['length']+1)
	plt.xticks(block)
	plt.grid()
	plt.tight_layout()
	plt.savefig(f'CH0_pulse/output/{output}/blocks/arrival_time.png')
	plt.show()

	
	a = arrival_diffs[:9]-arrival_diffs[5]
	a = np.r_[a,arrival_diffs[9:11]]

	plt.scatter(block,a)
	plt.ylabel("arrival time [ms]")
	plt.xlabel("block")
	#plt.ylim(0-1,setting['length']+1)
	plt.xticks(block)
	plt.grid()
	plt.tight_layout()
	plt.savefig(f'CH0_pulse/output/{output}/blocks/arrival_time.png')
	plt.show()
	
	

	b = 0
	for i,j in zip(blocks_path_0,blocks_path_1):
		block_0 = gp.loadbi(i,"text")
		block_1 = gp.loadbi(j,"text")
		height_0 = gp.peak(block_0,config["presamples"],para['peak_max'],para['peak_x'],para['peak_w'])
		height_1 = gp.peak(block_1,config["presamples"],para['peak_max'],para['peak_x'],para['peak_w'])
		ratio = height_0[1]/height_1[1]
		plt.scatter(b+1,ratio,c=cm.hsv(float(b)/float(len(blocks_path_0))))
		b += 1
		
	
	#plt.plot(x_fit, y_fit, "--")
	# plt.title("pixel vs pulseheight1/ pulseheight2")
	plt.ylabel("pulseheight0/ pulseheight1")
	plt.xlabel("block")
	#plt.ylim(0-1,setting['length']+1)
	plt.xticks(block)
	plt.grid()
	plt.tight_layout()
	plt.savefig(f'CH0_pulse/output/{output}/blocks/PH_PH_ratio.png')
	plt.show()
	


if __name__ == "__main__":
	main()
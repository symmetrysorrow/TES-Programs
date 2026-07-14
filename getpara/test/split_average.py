import matplotlib.pyplot as plt
import getpara as gp
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import os
import tqdm
import glob

# parameter 
para = "height" 
n_block = 11
a_ini = 1.5 # y = ax + 1
resolution = 0.00001 # a decrement



def main():

	setting = gp.loadJson()
	
	config = setting["Config"]
	analysis = setting["main"]
	select = setting["select"]
	ch = config["channel"]
	presamples = config["presamples"]
	time = gp.data_time(config['rate'],config['samples'])
	output = f'CH{ch}_pulse/output/{config["output"]}/blocks'
	os.chdir(config["path"])

	blocks_path = glob.glob(f"CH{ch}_pulse/output/{config['output']}/blocks/block_*.txt")

	
	
	for block in blocks_path:
		print(block)
		idx = gp.loadIndex(block)

		array = []
		n = 0
		for num in idx:
			path = f'CH{ch}_pulse/rawdata/CH{ch}_{num}.dat'
			data = gp.loadbi(path,config["type"])
			base,data = gp.baseline(data,presamples,analysis["base_x"],analysis["base_w"])
			array.append(data)
		av = np.mean(array,axis=0)	
		n+=1		
		np.savetxt(f'{output}/average_pulse_ch{ch}_{os.path.splitext(os.path.basename(block))[0]}.txt',av)
		plt.plot(time,av)
		

	gp.graugh_condition(setting["graugh"])
	plt.xlabel("time(s)")
	plt.ylabel("volt(V)")
	plt.title("average pulse")
	plt.savefig(f'{output}/block_average_pulses_ch{ch}.png')
	plt.show()
	plt.cla()

main()
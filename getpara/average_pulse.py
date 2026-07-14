import getpara as gp 
from tkinter import filedialog
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import shutil

def main():
	

	setting = gp.loadJson()
	config = setting["Config"]
	_main = setting["main"]
	presamples = config["presamples"]
	time = gp.data_time(config['rate'],config['samples'])
	os.chdir(config["path"])

	path = filedialog.askopenfilename()
	idx = gp.loadIndex(path)
	setting["select"] = idx
	output = os.path.dirname(path)


	ch = input("cahnnel: ")
	df_path = f'CH{ch}_pulse/output/{config["output"]}/output.csv'
	df = pd.read_csv(df_path,index_col=0)
	df = gp.select_condition(df,setting["select"])

	# create output dir
	output_select = f'{output}'
	if not os.path.exists(f"{output_select}/img"):
		os.makedirs(f"{output_select}/img",exist_ok=True)
	else:
		shutil.rmtree(f"{output_select}/img")
		os.mkdir(f"{output_select}/img")

	
	ax = ["rise","height"]

	x,y = gp.extruct(df,*ax)
	df_sel = df.loc[idx]
	
	x_sel,y_sel =gp.extruct(df_sel,*ax)

	plt.scatter(x,y,s=2,alpha=0.7)
	plt.scatter(x_sel,y_sel,s=4)
	plt.xlabel(gp.ax_unit[ax[0]])
	plt.ylabel(gp.ax_unit[ax[1]])
	plt.title(f"{ax[0]} vs {ax[1]}")
	plt.grid()
	plt.savefig(f'{output}/{ax[0]} vs {ax[1]}_sel_ch{ch}.png')
	plt.show()
	plt.cla()

	array = []
	for num in idx:
		path = f'CH{ch}_pulse/rawdata/CH{ch}_{num}.dat'
		data = gp.loadbi(path,config["type"])
		base,data = gp.baseline(data,presamples,_main["base_x"],_main["base_w"])
		# triggerd channel?
		analysis = setting["main"]
		
		if analysis['cutoff'] > 0:
			data = gp.BesselFilter(data,config["rate"],fs = analysis['cutoff'])

		array.append(data)

	av = np.mean(array,axis=0)

	
	np.savetxt(f'{output}/selected_average_pulse.txt',av)

	plt.plot(time,av)
	gp.graugh_condition(setting["graugh"])
	plt.xlabel("time(s)")
	plt.ylabel("volt(V)")
	plt.title("average pulse")
	plt.savefig(f'{output}/selected_average_pulse.png')
	plt.show()
	plt.cla()



if __name__ == "__main__":
	main()
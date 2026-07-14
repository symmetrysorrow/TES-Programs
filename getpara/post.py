import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plt_config
import os
import getpara as gp
import sys
import shutil
import pprint
import json
import re
from tkinter import filedialog
import platform
import tkinter as tk
import tqdm




def df_number(df):
	index = df.index.values
	num = [re.findall(r"\d+", i)[2] for i in index]
	df["number"] = num
	return df


def main():

	root  = tk.Tk()
	root.withdraw()

	ch0 = 0
	ch1 = 1
	para = sys.argv[1]
	setting = gp.loadJson()
	config = setting["Config"]
	rate, samples, presamples, threshold, ch = (
		config["rate"],
		config["samples"],
		config["presamples"],
		config["threshold"],
		config["channel"],
	)
	time = gp.data_time(rate, samples)
	os.chdir(config["path"])

	output_0 = f'CH{ch0}_pulse/output/{config["output"]}'
	output_1 = f'CH{ch1}_pulse/output/{config["output"]}'

	df_0 = pd.read_csv((f"{output_0}/output.csv"), index_col=0)
	df_1 = pd.read_csv((f"{output_1}/output.csv"), index_col=0)

	try:
		trig_ch = np.loadtxt("channel.txt")
	except:
		trig_ch = np.zeros(len(df_0))

	df_0_clear = gp.select_condition(df_0, setting["select"])
	df_1_clear = gp.select_condition(df_1, setting["select"])

	df_0_over, df_1_over = gp.overlap(df_0_clear, df_1_clear)
	x, y = df_0_over[para], df_1_over[para]

	plt.scatter(x, y, s=0.4)
	plt.xlabel(f"channel {ch0} [V]")
	plt.ylabel(f"channel {ch1} [V]")
	plt.title(para)
	plt.grid()
	plt.savefig(f"{output_0}/pulse_height_CH{ch0}_CH{ch1}.png")
	plt.show()
	plt.cla()

	# ----- choose creating selected output or no ---------------------------
	output_select = input("output: ")
	if output_select == "":
		print("Exit")
		exit()
	else:
		setting["select"]["output"] = output_select

	# create output dir
	output_select = f"{output_0}/{output_select}"
	os.makedirs(f"{output_select}", exist_ok=True)
	

	# pick samples from graugh

	picked = gp.pickSamples_PoST(df_0_over,df_1_over,para)
	print(f"Selected {len((picked))} samples.")

	# graugh picked samples
	x_sel, y_sel = df_0_over.loc[picked][para], df_1_over.loc[picked][para]
	np.savetxt(f"{output_select}/selected_index.txt", picked, fmt="%s")

	plt.scatter(x, y, s=2, alpha=0.7)
	plt.scatter(x_sel, y_sel, s=4)
	plt.xlabel(f"channel {ch0} [V]")
	plt.ylabel(f"channel {ch1} [V]")
	plt.grid()
	plt.savefig(f"{output_select}/pulse_height_{ch0}_{ch1}selected.png")
	#plt.show()
	plt.cla()

	a = '0'#input("create averagepulse?[0]")
	
	if a != '0':
		print("exit")
		exit()

	for ch in [ch0, ch1]:
		if not os.path.exists(f"{output_select}/CH{ch}/img"):
			os.makedirs(f"{output_select}/CH{ch}/img", exist_ok=True)
		else:
			shutil.rmtree(f"{output_select}/CH{ch}/img")
			os.mkdir(f"{output_select}/CH{ch}/img")
		
	for ch in [ch0, ch1]:
		# ----- create onetime pulse ---------------------------------------
		for num in tqdm.tqdm(picked):
			trig = trig_ch[int(num) - 1]
			if trig == int(ch):
				analysis = setting["main"]
			else:
				analysis = setting["main2"]
			path = f"CH{ch}_pulse/rawdata/CH{ch}_{num}.dat"
			data = gp.loadbi(path, config["type"])
			base, data = gp.baseline(
				data, presamples, analysis["base_x"], analysis["base_w"]
			)
			if analysis["cutoff"] > 0:
				data = gp.BesselFilter(data, rate, analysis["cutoff"])
			# save onetime rawdata figures
			plt.plot(time, data)
			gp.graugh_condition(setting["graugh"])
			plt.title(f"CH{ch} {num}")
			plt.xlabel("time(s)")
			plt.ylabel("volt(V)")
			plt.savefig(f"{output_select}/CH{ch}/img/{num}.png")
			plt.cla()

		# ----- delete noise data ---------------------------------------------

		
		try:
			print(output_select)
			root = tk.Tk()
			fle = filedialog.askopenfilenames(initialdir=f"{output_select}/CH{ch}/img")
			root.withdraw()
		except:
			fle = []
		

		for f in fle:
			num = int(re.findall(r"\d+", os.path.basename(f))[0])
			picked.remove(num)
			os.remove(f"{output_select}/CH{ch}/img/{num}.png")
			np.savetxt(f"{output_select}/selected_index.txt", picked, fmt="%s")
			if ch == ch0:
				df_0.at[num, "error"] = 0
			else:
				df_1.at[num, "error"] = 0

	for ch in [ch0, ch1]:
		# ----- create average pulse ---------------------------------------
		print("Creating Average Pulse...")
		array = []

		for num in tqdm.tqdm(picked):
			trig = trig_ch[int(num) - 1]
			if trig == int(ch):
				analysis = setting["main"]
			else:
				analysis = setting["main2"]
			path = f"CH{ch}_pulse/rawdata/CH{ch}_{num}.dat"
			data = gp.loadbi(path, config["type"])
			base, data = gp.baseline(
				data, presamples, analysis["base_x"], analysis["base_w"]
			)
			if analysis["cutoff"] > 0:
				data = gp.BesselFilter(data, rate, analysis["cutoff"])
			array.append(data)
		av = np.mean(array, axis=0)

		# plot averate pulse
		plt.plot(time, av)
		gp.graugh_condition(setting["graugh"])
		plt.xlabel("time(s)")
		plt.ylabel("volt(V)")
		plt.title("average pulse")
		plt.savefig(f"{output_select}/CH{ch}/average_pulse_{ch}.png")
		plt.show()

		# log scale
		plt.cla()
		plt.plot(time, av)
		plt.xlabel("time(s)")
		plt.ylabel("volt(V)")
		plt.title("average pulse")
		plt.yscale("log")
		plt.savefig(f"{output_select}/CH{ch}/average_pulse_{ch}_log.png")
		plt.show()

		np.savetxt(f"{output_select}/CH{ch}/average_pulse_{ch}.txt", av)

	# a = df.loc['CH0_pulse/rawdata\CH0_47388.dat']


if __name__ == "__main__":
	main()

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

# main.py　の後　output.csv　を用いてパラメータの相関をみる
# 
# 実行例--------------
# ex) python selectdata.py rise height
# exit[1]: 1 (終了)
# exit[1]: Enter (continue)
# 囲んだ所のパルスを抽出、平均パルスを作成
#
# setting.json　を用いてプロットするデータの条件を絞ることができる
# ex) "base->":0
# ex) "rise-<":0.01
# ex) "samples-=":100000


ax_unit = {
	"base":'base[V]',
	"height":'pulse height[V]',
	"peak_index":'peak index',
	"height_opt":'pulse height opt',
	"height_opt_temp":'pulse height opt temp',
	'rise':'rise [ms]',
	'decay':'decay [ms]',
	'area':'area',
	'rise_fit':'rise_fit[s]',
	'tau_rise':'tau_rise[s]',
	'tau_decay':'tau_decay[s]',
	'height_fit':'height_fit[V]',
	'rSquared':'rSquared'
}


def main():

	ax = sys.argv
	ax.pop(0)
	if len(ax) == 0:
		print("axis needed!!")
		exit()

	root  = tk.Tk()
	root.withdraw()
	
	#----- Initialize ------------------------------------------
	setting = gp.loadJson()
	config = setting["Config"]
	os.chdir(config["path"])
	

	output = f'CH{config["channel"]}_pulse/output/{config["output"]}'
	df = pd.read_csv((f'{output}/output.csv'),index_col=0)
	rate,samples,presamples,threshold,ch = config["rate"],config["samples"],config["presamples"],config["threshold"],config["channel"]
	time = gp.data_time(rate,samples)
	
	try:
		trig_ch = np.loadtxt('channel.txt')
	except:
		trig_ch = np.zeros(len(df))
	

	#----- select condition -------------------------------------------
	df_clear = gp.select_condition(df,setting["select"])
	print(f'Pulse : {len(df_clear)} samples')


	#----- plot data -------------------------------------------------
	# time vs ax
	if len(ax) == 1:
		y = gp.extruct(df_clear,ax)
		x = np.arange(len(y[0]))
		plt.scatter(x,y[0],s=0.4)
		plt.xlabel("data number")
		plt.ylabel(ax_unit[ax[0]])
		plt.savefig(f'{output}/{ax[0]}.png')
		plt.show()

	# ax vs ax
	elif len(ax) == 2:
		x,y = gp.extruct(df_clear,*ax)
		plt.scatter(x*1e3,y,s=2,alpha=0.7)
		plt.xlabel(ax_unit[ax[0]])
		plt.ylabel(ax_unit[ax[1]])
		plt.title(f"{ax[0]} vs {ax[1]}")
		plt.grid()
		plt.tight_layout()
		plt.savefig(f'{output}/{ax[0]} vs {ax[1]}.png')
		plt.show()
		plt.cla()
		
		#----- choose creating selected output or no ---------------------------
		output_select = input('output: ')
		if output_select == "":
			print("Exit")
			exit()
		else:
			setting['select']['output'] = output_select

		# create output dir
		output_select = f'{output}/{output_select}'
		if not os.path.exists(f"{output_select}/img"):
			os.makedirs(f"{output_select}/img",exist_ok=True)
		else:
			shutil.rmtree(f"{output_select}/img")
			os.mkdir(f"{output_select}/img")

		#---------------------------------------------------------------------

		# pick samples from graugh
		x,y = gp.extruct(df_clear,*ax)
		picked = gp.pickSamples(df_clear,*ax)
		print(f"Selected {len((picked))} samples.")
		np.savetxt(f"{output_select}/selected_index.txt", picked, fmt="%s")

		# graugh picked samples
		x_sel,y_sel =gp.extruct(df_clear.loc[picked],*ax)
		plt.scatter(x,y,s=2,alpha=0.5)
		plt.scatter(x_sel,y_sel,s=4)
		plt.xlabel(ax_unit[ax[0]])
		plt.ylabel(ax_unit[ax[1]])
		plt.title(f"{ax[0]} vs {ax[1]}")
		plt.grid()
		plt.savefig(f'{output_select}/{ax[0]}_{ax[1]}_selected.png')
		plt.show()
		plt.cla()

		a = input("image write?[0]")
	
		if a != '0':
			print("exit")
			exit()

		# one sumple
		if len(picked) == 1:
			path = f'CH{ch}_pulse/rawdata/CH{ch}_{picked[0]}.dat'
			picked_data = gp.loadbi(path,config["type"])
			gp.graugh(path,picked_data,time)
			gp.graugh_condition(setting["graugh"])
			plt.show()

		# multi sumples
		else:
			for num in tqdm.tqdm(picked):
				path = f'CH{ch}_pulse/rawdata/CH{ch}_{num}.dat'
				data = gp.loadbi(path,config["type"])

				# triggerd channel?
				trig = trig_ch[int(num)-1]
				if trig == int(ch):
					analysis = setting["main"]
				else:
					analysis = setting["main2"]

				
				base,data = gp.baseline(data,presamples,analysis['base_x'],analysis['base_w'])
				if analysis['cutoff'] > 0:
					data = gp.BesselFilter(data,rate,fs = analysis['cutoff'])

				# save onetime rawdata figures
				plt.plot(time*1e3,data)
				gp.graugh_condition(setting["graugh"])
				plt.title(f'CH{ch} {num}')
				plt.xlabel("time(ms)")
				plt.ylabel("volt(V)")
				plt.savefig(f'{output_select}/img/{num}.png')
				plt.cla()

			np.savetxt(f'{output_select}/selected_index.txt',picked,fmt="%s")


			#----- delete noise data ---------------------------------------------
			try:
				fle = filedialog.askopenfilenames(initialdir=f"{output_select}/img")
				root.withdraw()
			except:
				fle = []

			for f in fle:
				num =  int(re.findall(r'\d+', os.path.basename(f))[0])              
				picked.remove(num)
				os.remove(f'{output_select}/img/{num}.png')
				np.savetxt(f'{output_select}/selected_index.txt',picked,fmt="%s")
				df.at[num,"error"] = 0
			#-------------------------------------------------------------------
			

			#----- create average pulse ---------------------------------------
			print('Creating Average Pulse...')
			array = []
			for num in picked:
				path = f'CH{ch}_pulse/rawdata/CH{ch}_{num}.dat'
				data = gp.loadbi(path,config["type"])
				base,data = gp.baseline(data,presamples,analysis['base_x'],analysis['base_w'])
				if analysis['cutoff'] > 0:
					data = gp.BesselFilter(data,rate,analysis['cutoff'])
				array.append(data)
			av = np.mean(array,axis=0)

			# plot averate pulse
			plt.plot(time*1e3,av)
			gp.graugh_condition(setting["graugh"])
			plt.xlabel("time(ms)")
			plt.ylabel("volt(V)")
			plt.title("average pulse")
			plt.savefig(f'{output_select}/average_pulse.png')
			plt.show()

			# log scale
			plt.cla()
			plt.plot(time*1e3,av)
			plt.xlabel("time(ms)")
			plt.ylabel("volt(V)")
			plt.title("average pulse")
			plt.yscale('log')
			plt.savefig(f'{output_select}/average_pulse_log.png')
			plt.show()

			np.savetxt(f'{output_select}/average_pulse.txt',av)

	
	gp.saveJson(setting,path=output_select)
	df.to_csv(f'{output}/output.csv')

#実行
if __name__=='__main__':
	main()
	print('end')

	

	

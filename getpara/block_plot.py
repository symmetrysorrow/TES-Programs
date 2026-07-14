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


def main():
	#----- Initialize ------------------------------------------
	ax = sys.argv
	ax.pop(0)
	if len(ax) == 0:
		print("axis needed!!")
		exit()

	setting = gp.loadJson()
	config = setting["Config"]
	os.chdir(config["path"])
	

	output = f'CH{config["channel"]}_pulse/output/{config["output"]}'
	df = pd.read_csv((f'{output}/output.csv'),index_col=0)
	rate,samples,presamples,threshold,ch = config["rate"],config["samples"],config["presamples"],config["threshold"],config["channel"]
	time = gp.data_time(rate,samples)
	
	#----- select condition -------------------------------------------
	df_clear = gp.select_condition(df,setting["select"])
	print(f'Pulse : {len(df_clear)} samples')

	
	blocks = natsort.natsorted(glob.glob(f"{output}/blocks/block_*.txt"))
	print(blocks)

	for i in blocks:
		idx = gp.loadIndex(i)
		df_sel = df_clear[df_clear.index.isin(idx)]
		x,y = gp.extruct(df_sel,*ax)
		num = int(gp.num(os.path.basename(i))[0])
		plt.scatter(x*1e3,y,s=2,c=cm.hsv((float(num-1))/float(len(blocks))))


	#plt.title(f"{ax[0]} vs {ax[1]}")
	plt.grid()
	plt.xlabel("rise [ms]")
	plt.ylabel("height [V]")
	plt.tight_layout()
	plt.savefig(f'{output}/blocks/{ax[0]}_{ax[1]}_block.pdf')
	plt.show()
	plt.cla()

		
	
main()
	
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import libs.plt_config
import os
import libs.getpara as gp
import sys
import shutil
import pprint
import json
import re
from tkinter import filedialog
import platform
import tkinter as tk
import tqdm



def main():
	ch0 = 0
	ch1 = 1
	para = "area"
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
	area = df_0[para].values + df_1[para].values
	print(area)
	plt.hist(area,bins=1024)
	plt.show()
	
    
if __name__ == '__main__':
	main()
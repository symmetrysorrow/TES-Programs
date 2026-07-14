import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import os
import pandas as pd
import libs.getpara as gp
import json

BIN = 80

def main():

	ax = sys.argv
	ax.pop(0)
	if len(ax) != 1:
		print("add parameter you want to see as histgrum!")
		sys.exit()
	para = ax[0]

	# Load Setting 
	set = gp.loadJson()
	path,ch = set["Config"]["path"],set["Config"]["channel"]
	output = f'CH{set["Config"]["channel"]}_pulse/output/{set["Config"]["output"]}'
	os.chdir(path)

	# Load data and transform histgrum
	df = pd.read_csv((f'{output}/output.csv'),index_col=0)
	df_sel = gp.select_condition(df,set['select'])
	print(df_sel)
	data = df_sel[para]
	hist,bins = np.histogram(data,bins=BIN)
	

	plt.bar(bins[:-1],hist,width = bins[1]-bins[0])
	plt.xlabel('Pulse Height [ch]')
	plt.ylabel('Counts')
	plt.tick_params(axis='both', which='both', direction='in',
					bottom=True, top=True, left=True, right=True)
	plt.grid(True, which='major', color='black', linestyle='-', linewidth=0.2)
	plt.grid(True, which='minor', color='black', linestyle=':', linewidth=0.1)
	plt.savefig(f'{output}/{set["select"]["output"]}/histgrum.png')
	plt.show()

	


if __name__ == "__main__":
	main()
	print('end')
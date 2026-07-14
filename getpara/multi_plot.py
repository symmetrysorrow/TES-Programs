

import matplotlib.pyplot as plt
from tkinter import filedialog
import getpara as gp
import numpy as np
import os
import natsort
import matplotlib.cm as cm
import tkinter as tk
import itertools
import plt_config


def main():
	setting = gp.loadJson()
	config = setting["Config"]
	_main = setting["main"]
	presamples = config["presamples"]
	
	cnt = "1"
	root  = tk.Tk()
	paths = []
	while cnt == "1":
		paths.append(natsort.natsorted(filedialog.askopenfilenames()))
		#root.withdraw()
		cnt = input("continue [1]: ")
	
	paths = list(itertools.chain.from_iterable(paths))
	print(paths)

	cnt=0
	n = len(paths)

	for i in paths:
		data1 = np.loadtxt(i)

		time = gp.data_time(config["rate"],config["samples"])
		if n == 2:
			if cnt == 0:
				plt.plot(time*1000,data1,color = "tab:blue" )
			else:
				plt.plot(time*1000,data1,color = "tab:orange" )

		else:
			plt.plot(time*1000,data1,color = cm.hsv((float(cnt))/float(n)))

		cnt += 1
		
	plt.xlabel("time [ms]")
	plt.ylabel("volt [V]")
	plt.tight_layout()
	plt.grid()
	#plt.savefig(f"{os.path.dirname(paths)}/double_plot.png")
	plt.show()

main()
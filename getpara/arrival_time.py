import matplotlib.pyplot as plt
from tkinter import filedialog
import getpara as gp
import numpy as np

def main():
	path = "G:/TSURUTA/20230703/180mK_550uA_Difftrig5E-5_100kHz_gain10/CH0_pulse/output/run03/anl01"
	rate = 500000
	samples = 100000
	time = np.linspace(0,1/rate*samples,samples)
	print(time)
	data1 = np.loadtxt(f"{path}/average_pulse.txt")
	data2 = np.loadtxt(f"{path}/selected_average_pulse.txt")

	data2 = gp.BesselFilter(data2,rate,1000)

	data1_diff = np.diff(data1)
	data1_diff_diff = np.diff(data1_diff)
	data2_diff = np.diff(data2)
	data2_diff_diff = np.diff(data2_diff)


	plt.plot(time,data1,"o",markersize = 1.0)
	plt.plot(time[:-1],data1_diff *10,"o",markersize = 1.0)
	plt.plot(time[:-2],data1_diff_diff *100,"o",markersize = 1.0)
	plt.plot(time[:-1],data2_diff *100,"o",markersize = 1.0)
	plt.plot(time[:-2],data2_diff_diff *1000,"o",markersize = 1.0)
	plt.plot(time,data2,"o",markersize = 1.0)
	plt.show()
	

main()
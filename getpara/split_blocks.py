import matplotlib.pyplot as plt
import getpara as gp
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import os
import tqdm
import glob
import plt_config
import sys

# parameter 



# liner function
def func(x,a):
	return a*x


def main():
	ax = sys.argv
	ax.pop(0)
	if len(ax) != 1:
		print("add parameter you want to see as histgrum!")
		sys.exit()
	para = ax[0]


	setting = gp.loadJson()
	a_ini = 1.4 	# y = ax + 1
	resolution = 0.00001	# a decrement
	
	config = setting["Config"]
	select = setting["select"]
	os.chdir(config["path"])

	ch = config["channel"]

	output_0 = f'CH{0}_pulse/output/{config["output"]}'
	output_1 = f'CH{1}_pulse/output/{config["output"]}'

	os.makedirs(f"CH{ch}_pulse/output/{config['output']}/blocks",exist_ok=True)

	df_0 = pd.read_csv((f'{output_0}/output.csv'),index_col=0)
	df_1 = pd.read_csv((f'{output_1}/output.csv'),index_col=0)
	
	df_0 = gp.select_condition(df_0,select=select)
	df_1 = gp.select_condition(df_1,select=select)
	
	df_0_lap,df_1_lap = gp.overlap(df_0,df_1)


	n_block = int(input("block: ")) 
	length = len(df_0_lap)
	block = int(length/n_block)

	# each blocks have same samples
	if length%n_block != 0:
		df_0_lap = df_0_lap[:block*n_block]
		df_1_lap = df_1_lap[:block*n_block]

	print(f"length: {len(df_0_lap)}")
	print(f"1-block: {block}")

	a_array = np.arange(a_ini,0.0,resolution*-1)
	x_line = np.arange(0,df_0_lap[para].max()*1.1,0.001)

	y_line_down = func(x_line,a_ini)
	plt.scatter(df_0_lap[para],df_1_lap[para],s=2)
	plt.plot(x_line,y_line_down,"--",linewidth=0.5,markersize=0.7,color="black",alpha=0.7)
	plt.show()
	plt.cla()
	#plt.plot(x_line,y_line_down,"--",linewidth=0.5,markersize=0.7,color="black",alpha=0.7)

	# ---check first data -------
	print("check first data .......")
	for i in a_array:
		d_up = df_1_lap[para].values - func(df_0_lap[para],i).values 
		d_down = df_1_lap[para].values - func(df_0_lap[para],a_ini).values 
		
		idx_up = (d_up > 0).nonzero()[0]
		idx_down = (d_down < 0).nonzero()[0]

		span = list(set(idx_up) & set(idx_down))
		indexs = df_0_lap.iloc[span].index
		
		if len(indexs)==1:
			y_line_down = func(x_line,i)
			#plt.plot(x_line,y_line_down,"--",linewidth=0.5,markersize=0.7,color="black",alpha=0.7)
			a_ini = i+resolution
			break

	n = n_block + 1

	a_array = np.arange(a_ini,0.0,resolution*-1)

	# decrease a per resolution and sepalate blocks
	for i in tqdm.tqdm(a_array):
		d_up = df_1_lap[para].values - func(df_0_lap[para],i).values 
		d_down = df_1_lap[para].values - func(df_0_lap[para],a_ini).values 
		
		idx_up = (d_up >= 0).nonzero()[0]
		idx_down = (d_down <= 0).nonzero()[0]

		span = list(set(idx_up) & set(idx_down))
		indexs = df_0_lap.iloc[span].index
	
		if len(indexs)>=block:
			if len(indexs)!=block:
				print('resolution recquired!!!')
				exit()
			a_ini = i - resolution
			n -= 1

			y_line_down = func(x_line,a_ini)
			np.savetxt(f"CH{ch}_pulse/output/{config['output']}/blocks/block_{n}.txt",indexs)

			df_sel_0,df_sel_1 = df_0_lap.loc[indexs],df_1_lap.loc[indexs]
			x,y = df_sel_0[para],df_sel_1[para]
			#plt.scatter(x,y,s=2,color = cm.hsv((float(n))/float(n_block)))

			if n==n_block:
				plt.scatter(x,y,s=2,color = "tab:red")
			elif n==2:
				plt.scatter(x,y,s=2,color = "tab:blue")
			else:
				plt.scatter(x,y,s=2,color = "0.5")
				

			#plt.plot(x_line,y_line_down,"--",linewidth=0.5,markersize=0.7,color="black",alpha=0.7)

			df_0_lap.drop(indexs, inplace = True)
			df_1_lap.drop(indexs, inplace = True)
			#plt.show()

			x,y = df_0_lap[para],df_1_lap[para]
			#plt.scatter(x,y,s=2,color = cm.hsv((float(n))/float(n_block)))
			#plt.show()

		if n == 1:
			break

	plt.grid()
	plt.xlabel(f"Pulse Height (CH0) [V]")
	plt.ylabel(f"Pulse Height (CH1) [V]")
	plt.tight_layout()
	plt.savefig(f"CH{ch}_pulse/output/{config['output']}/blocks/blocks.png")
	plt.show()
	plt.cla()

	

if __name__ == "__main__":
	main()









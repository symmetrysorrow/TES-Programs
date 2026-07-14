import numpy as np
import pandas as pd
import scipy.fftpack as fft
import libs.getpara as gp
import libs.fft_spectrum as sp
import matplotlib.pyplot as plt
import shutil
import os
import pprint
import tqdm
import scipy


# output.csvは必ず閉じておくこと（ファイルを保存できなくなる）
# ---Parameter--------------------------------------------------




def lowpass(F,fq,cf):
	F[(fq>cf)] = 0
	return F

def main():
	# read setting
	set = gp.loadJson()
	path = set["Config"]["path"]
	pprint.pprint(set)
	os.chdir(path)
	output = f'CH{set["Config"]["channel"]}_pulse/output/{set["Config"]["output"]}'
	output_noise = f'CH{set["Config"]["channel"]}_noise/output/{set["Config"]["output"]}'
	output_sel = os.path.join(output,set['select']['output'])
	ch,rate,samples = set["Config"]["channel"],set["Config"]["rate"],set["Config"]["samples"]
	time = gp.data_time(rate,samples)
	fq = np.arange(0,rate,rate/samples)

	
	
	pulse_av = np.loadtxt(f"{output_sel}/average_pulse.txt")
	noise_spe = np.loadtxt(f"{output_noise}/modelnoise.txt")


	#outputファイルを読み込み
	df  = pd.read_csv(f"{output}/output.csv",index_col=0)
	df_sel = gp.select_condition(df,set)

	eta = set['Config']['eta']
	cf = set['main']['cutoff']

	#Fourier transform
	#pulse_av = sp.hanning(pulse_av,samples)[0]
	F = fft.fft(pulse_av)
	scipy.fftpack
	amp = np.abs(F)
	amp_spe = np.sqrt(amp)*int(eta)*1e+6*np.sqrt(1/rate/samples)
	np.savetxt(f'{output_sel}/averagepulse_spectrum.txt',amp_spe)
	sp.graugh_spe(fq[:int(samples/2)+1],amp_spe[:int(samples/2)+1])
	plt.savefig(f'{output_sel}/averagepulse_spectrum.png')
	plt.show()

	#LowPassFilter
	F2 = lowpass(F,fq,cf)
	amp_filt = np.abs(F2)
	amp_spe_filt = np.sqrt(amp_filt)*int(eta)*1e+6*np.sqrt(1/rate/samples)
	sp.graugh_spe(fq[:int(samples/2)+1],amp_spe_filt[:int(samples/2)+1])
	plt.show()

	filt = fft.ifft(F2/noise_spe)
	filt = filt.real
	np.savetxt(f'{output_sel}/opt_template.txt',filt)
	plt.plot(time,filt)
	plt.savefig(f'{output_sel}/opt_template.png')
	plt.show()

	# optimal filter
	for index,row in tqdm.tqdm(df_sel.iterrows()):
		path = f'CH{ch}_pulse/rawdata/CH{ch}_{index}.dat'
		data = gp.loadbi(path)
		base = df.at[index,"base"]
		data = data-base
		if set['main']['cutoff'] > 0:
			data = gp.BesselFilter(data,rate,fs = set['main']['cutoff'])
		df.at[index,"height_opt"] = np.sum(data*filt)
	df.to_csv(f"{output}/output.csv")
	

if __name__ == "__main__":
	main()
	print("end")

	
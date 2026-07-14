 #  coding: utf-8

#------------------------------------------------------
#  モデルノイズの作成　
#------------------------------------------------------
#last updated Tsuruta 2019.06.15----------------
#最適フィルタをかけるのに必要なモデルノイズを作成します。
#ファイルの数をわざわざ打ち込まなくても良いようにしました。
#python3で使用可

import math
import glob
import numpy as np
import scipy.optimize as opt
import scipy.fftpack as fft
import matplotlib.pyplot as plt
from array import array
import time
import os
import getpara as gp


n = 100000			#データ点数


#----------- ファイルの読み込み ----------------
#globでファイルを読み込む
'''
data = []
for file in glob.glob('*.awd'):
	data_sub = []
	data.append(data_sub)
	del data_sub[:]
#	for line in np.loadtxt(file,comments='#'):
#		data_sub.append(line)
	f = open(file,'r')
	file = f.readlines()
	for i in range(8,n+8,1):
		data_sub.append(float(file[i].strip()))
'''

#filelistから読み込む(filelistが必要)
'''
flug = os.path.exists('./filelist.txt')
if flug:
	data = []
	f = open('filelist.txt','r')
	for file in f:
		data_sub = []
		data.append(data_sub)
		del data_sub[:]
		f = open(file.strip(),'r')
		file = f.readlines()
		for i in range(8,n+8,1):
			data_sub.append(float(file[i].strip()))
	f.close
else:
	print ('filelist.txtがありません')
	exit()
'''

#--------ファイル名からデータ読み込み---------



filename_1 = 'CH0_'
filename_2 = '.dat'

number_file = []
data = []
'''
for i in range(0,m,1):
    filename = filename_1 + str(int(i + 1)) + filename_2
    data_sub = np.loadtxt(filename,comments='#',skiprows = 6)
    data.append(data_sub)
    print (i)
'''
#-----------正常なファイルのみを読み込む-------------------
number_file = []
erorr = []												#erorrが起こった場合、格納する。

for file in glob.glob('*.dat'):
	file_i = file.replace(filename_1,'')
	file_i = file_i.replace(filename_2,'')
	number_file.append(int(file_i))
	number_file.sort()

print ('The Number of file =',len(number_file))			#ファイルの数。つまり信号の数。

data = []

quasi_pulse = []
quasi_pulse_sub = []
for i in range(len(number_file)):
	try:
		filename = filename_1 + str(int(number_file[i])) + filename_2
		data_sub = gp.loadbi(filename,"binary")
	except:
		erorr.append(number_file[i])
		print ('data',str(number_file[i]),'is erorr')
	else:
		if len(data_sub) == n:
			data.append(data_sub)
			quasi_pulse_sub.append(number_file[i])
			print ('Finish storaging of data',str(number_file[i]))
		else:
			erorr.append(number_file[i])
			print ('data',str(number_file[i]),'is erorr')

print ('The Number of readable =',len(data))
#data_subとして読み込んだ後，判定を行い正常なものだけがdataに格納される。
#データ抜けがあるものなどはここで弾かれるはず。

#----------- パラメータの設定 --------------------
time_inc = float(input('time increament[us]:'))*10**-6

eta = 109					# 変換係数 [uA/V]  : 超伝導領域での傾き

#----------- データの選別 -----------------
#dataが正常に読み込めた全データで，ノイズとみなされたものがdata2に入り
#data3にアレイ化されたものが入る。
min = []
max = []
for i in range(len(data)):
	min.append(np.amin(data[i]))
	max.append(np.amax(data[i]))

pulseheight = np.array(max)-np.array(min)
print (pulseheight)
pulsemax = float(input('ノイズと判定する波高の最大値[V]:'))#この波高値以下をノイズと判定
data2 = []
data3=[]
for i in range(len(data)):
	if pulseheight[i] < pulsemax:
		data2.append(data[i])

data3 = np.array(data2)
mp = len(data3)   #指定した最大値より小さいノイズの数
print (mp )#要するに，ノイズと判定されたデータの数が表示される。

#---------- fft ----------
model_n = np.array([0]*n)
for i in range(mp):
	x = fft.fft(data3[i])
	spec = [c.real ** 2 + c.imag ** 2 for c in x]
	model_n = model_n + np.array(spec)
#データに抜け等があると（そのようなデータは弾いているはずだが）
#operands could not be broadcast together with shapes(49999,)(50000,)
#のようなエラーが出る。
	print ('Fourier transform',i)
model_n = model_n/mp
model_f = model_n[:n//2+1]			#領域を半分にする([:n/2+1]は0からn/2+1の要素という意味)
amp_spe = np.sqrt(model_f)*eta*1e+6*np.sqrt(time_inc/n)

#---------- 周波数リスト作成 ----------本来はこっち
'''
frq_inc = 1/(time_inc*n)
print (frq_inc)
freqlist = []
for i in range(n/2+1):
	freqlist.append(frq_inc/1000*i)
freqlist = np.array(freqlist)
mf = len(freqlist)
'''
#----------周波数リスト作成（黒岩方式）-------
frq_inc = 1/(time_inc*n)				                            #周波数間隔[Hz]
freqlist = np.arange(0,(n/2+1)*frq_inc,frq_inc)/1000		        #周波数リスト[kHz]
mf = len(freqlist)
#
#---------- ファイルに出力 ----------

f = open("modelnoise.txt","w") #モデルノイズ
for x in range(len(model_n)):
	f.write(str(model_n[x])+'\n')
f.flush()
f.close()

f = open("noise_intensity.txt","w") #ノイズ強度
for x in range(mf):
	f.write(str(amp_spe[x])+'\n')
f.flush()
f.close()

f = open("noise_spe.txt","w")
f.write('#Timeincreament[us]'+'\n'+str(time_inc)+'\t'+'0'+'\n')
f.write('#Frequency[kHz]'+'\t'+'Intensity[pA/kHz^{1/2}]'+'\n')
for x in range(mf):
	f.write(str(freqlist[x])+'\t'+ str(amp_spe[x])+'\n')
f.flush()
f.close()

f = open("The_number_of_noise.txt","w")
f.write(str(len(data3)))
f.flush()
f.close()

f = open("The_number_of_erorr.txt","w")
f.write(str(len(erorr)))
f.flush()
f.close()

f = open("input_parameter.txt","w")
f.write(str(len(data))+'\n'+str(n)+'\n'+str(time_inc)+'\n'+str(pulsemax*1e+3)+'\n')
f.write(filename_1+'\n'+filename_2+'\n')
f.flush()
f.close()
#---------- ノイズスペクトルグラフ化 ----------

plt.plot(freqlist,amp_spe,linestyle='-',linewidth=0.7)
#plt.xlim(1,1/(2*time_inc)/1000)
#plt.ylim(1e+1,1e+4)
#plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency[kHz]')
plt.ylabel('Intensity[pA/kHz$^{1/2}$]')
plt.grid()
plt.show()

exit()

# -*- coding: utf-8 -*-
#ベースの分解能の算出
#optimal_filter.txtとcalibration_curve.txtをrawdataフォルダの中に入れておく
import math
import os
import glob
import numpy as np
import scipy.optimize as opt
import scipy.fftpack as fft
import matplotlib.pyplot as plt
from array import array


amp = 10

#------------- パラメータの読み込み --------------------
os.chdir('./rawdata')


f = open('input_parameter.txt','r')
fdata = f.readlines()
m = int(fdata[0].strip())
n = int(fdata[1].strip())
time_inc = float(fdata[2].strip())
pulsemax = float(fdata[3].strip())*1e-3
d_id = fdata[4].strip()
ext = fdata[5].strip()

filt = np.loadtxt('optimal_filter.txt',comments='#')
#cali = np.loadtxt('calibration_curve.txt',comments='#')
f = open('calibration_curve_last.txt','r')
fdata = f.readlines()
cali = float(fdata[0].strip())


#os.chdir('..')


#------------- ファイルの読み込み ----------------
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
	print 'filelist.txtがありません'
	exit()
'''

#ファイル名から読み込む
A = []
model_n = 0
mp = 0
for i in range(30000):   #最初は試しに少ない数で実行してみる
	data = []
	f = open(d_id+str(i+1)+ext,'r')			#ファイル名からファイルを開く
	file = f.readlines()
	skip = 6								#ファイルの先頭から読み飛ばす行数
	bool = ''
	print (f)
	for j in range(skip,n+skip,1):
		try:
			data.append(float(file[j].strip())/amp)
			bool = 'y'
		except:
			bool = 'n'
			break                                                  #エラーが出た場合脱出
	f.close()

	if bool =='n':
		continue

	data = -np.array(data)


#----------- 波高値の計算 -----------------
	min = np.amin(data)
	max = np.amax(data)
	pulseheight = max - min

	if pulseheight < pulsemax and pulseheight != 0:
		mp += 1


#----------- 最適フィルタ適用 --------------
		A.append(np.sum(data*filt))
		#for j in range(m):
			#A[i] += plot[i][j]*filt.real[j]


#print len(data)
#print m


pulseheight_f = np.array(A)
#print pulseheight_f

#print pulseheight_f

#エネルギー校正
#pulseheight_fc = cali[0]*pulseheight_f**2+cali[1]*pulseheight_f+cali[2]
pulseheight_fc = cali*pulseheight_f

os.chdir('./output')


#--------------- グラフ化 ------------------
plt.hist(pulseheight_fc, bins =300)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Energy[eV]',fontsize=16, fontname='serif')
plt.ylabel('Count',fontsize=16, fontname='serif')
plt.savefig("base_spectrum.pdf",format = 'pdf')
plt.show()

#ファイルに出力
f = open("base_spectrum.txt","w")
for x in range(len(pulseheight_fc)):
	f.write(str(pulseheight_fc[x])+'\n')
f.flush()
f.close()


#----------- フィッティング ------------
pulse_fmin = float(input('minimum pulseheight:'))
pulse_fmax = float(input('maximum pulseheight:'))
ap = len(pulseheight_fc)
pulseheight_f = []
for i in range(ap):
	if pulse_fmin <= pulseheight_fc[i] <= pulse_fmax:
		pulseheight_f.append(pulseheight_fc[i])

hist, bin_edges = np.histogram(pulseheight_f, bins=40)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

#フィッティングの初期値
p0 = [0]*3
p0[0]= np.max(hist)										#hist内の最大値を返す
p0[1] = (pulse_fmax - pulse_fmin)/6
p0[2] = bin_centers[np.argmax(hist)]					#argmaxは最大値の添字を返す
print ('initial value',p0)

def gauss(p,x,y):
	residual=y-(p[0]*np.exp(-(x-p[2])**2/(2*p[1]**2)))
	return residual

output=opt.leastsq(gauss,p0,args=(bin_centers,hist))
a = output[0][0]
sigma = output[0][1]
mu = output[0][2]
print  (a, sigma, mu)
x = np.arange(pulse_fmin,pulse_fmax,0.001)
#y = a/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))
y = a*np.exp(-(x-mu)**2/(2*sigma**2))
E = 2*sigma*np.sqrt(2*np.log(2))

print ('Energy resolution :',E,'keV')

plt.hist(pulseheight_f, bins =40)
#plt.ylim(0,300)
plt.plot(x, y, color='red',linewidth=1,linestyle='-')
plt.xlabel('Energy[keV]')
plt.ylabel('Count[-]')
plt.title('dE= '+'{:.3f}'.format(E)+'keV')
plt.savefig("base_spectrum_fit.pdf",format = 'pdf')
plt.show()

#---------------- ファイルに出力 -----------------

f = open('baseline_energy_resolution.txt','w')
f.write(str(E)+'[keV]'+'\n')
f.flush()
f.close()
#---------------------------------------------

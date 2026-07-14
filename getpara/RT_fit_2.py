# -*- coding: utf-8 -*-

#last updated 2018/08/29 by kurume-------------------

import math
import os
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit


#-------------変更するparameters---------------------------------

RN_s = 0		#fittng時のパラメータの取り得る範囲。常伝導抵抗値RN[mohm],転移温度Tc[mK],T1[mK],T2[mK]
RN_f = 'inf'					#例えば、パラメタRNを0[mohm]以上で動かしたい場合は、RN_s = 0.0,RN_f = 'inf'。
Tc_s = 0
Tc_f = 'inf'			#T1、T2についても同様。
T1_s = 0.0					#注意として、全てのパラメタで下限は0.0である。
T1_f = 'inf'				#また、計算量を減らすために、RN、T_cについては、目星をつけた値にすることを推奨する。
T2_s = 0.0					#さらに、Tcについては下のList of Hot Bathの最大値と最小値の範囲内におさめること。
T2_f = 'inf'

para0 = [0.1,0.1,0.1,0.1]		#fitting時のパラメタ列の初期値。[RN,Tc,T1,T2]である。0.0から始めることは推奨しない。

sup_R = 90				#転移幅を求めるためのに変更する。転移幅をΔT[mK]とすると、
inf_R = 10					#T[mK]のとき、RNのinf_R[%]、T+ΔT[mK]のとき、RNのsup_R[%]

alpha_percent = 20			#常伝導抵抗値のalpha_percent[%]におけるalpha

#-------------- RT_01.pyで作ったデータを読み込む-----------------


path = input('path: ')
os.chdir(path)

T = np.loadtxt('output/10uA.txt')[0]
R_1 = np.loadtxt('output/10uA.txt')[1]
R_2 = np.loadtxt('output/20uA.txt')[1]


#--------------- R-T曲線をフィッティングする-----------------------

def fit_func(para,x,y):                  #fit関数の定義。paraはパラメタ列。
	RN = para[0]                         #パラメタ列の1番目は常伝導抵抗値(mohm)
	Tc = para[1]                         #パラメタ列の2番目は相転移温度(mK)
	T1 = para[2]                         #パラメタ列の3番目はT1(mK)
	T2 = para[3]                         #パラメタ列の4番目はT2(mK)
	if RN_f =='inf':
		RN_func = np.exp(RN)+RN_s
	else:
		RN_func = ((RN_f-RN_s)/math.pi)*math.atan(RN)+(RN_f+RN_s)/2

	if Tc_f == 'inf':
		Tc_func = np.exp(Tc)+ Tc_s
	else:
		Tc_func = ((Tc_f-Tc_s)/math.pi)*math.atan(Tc)+(Tc_f+Tc_s)/2

	if T1_f =='inf':
		T1_func = np.exp(T1)+T1_s
	else:
		T1_func = ((T1_f-T1_s)/math.pi)*math.atan(T1)+(T1_f+T1_s)/2

	if T2_f =='inf':
		T2_func = np.exp(T2)+T2_s
	else:
		T2_func = ((T2_f-T2_s)/math.pi)*math.atan(T2)+(T2_f+T2_s)/2

	residual = y - (RN_func/((1+np.exp(-(x-Tc_func)/T1_func))*(1+np.exp(-(x-Tc_func)/T2_func))))
	return(residual)

result_1 = opt.leastsq(fit_func,para0,args = (T,R_1),)
result_2 = opt.leastsq(fit_func,para0,args = (T,R_2))

RN_fit_1 = result_1[0][0]
Tc_fit_1 = result_1[0][1]
T1_fit_1 = result_1[0][2]
T2_fit_1 = result_1[0][3]

RN_fit_2 = result_2[0][0]
Tc_fit_2 = result_2[0][1]
T1_fit_2 = result_2[0][2]
T2_fit_2 = result_2[0][3]

if RN_f == 'inf':
	a_1 =  np.exp(RN_fit_1)+RN_s
	print ('RN_fit_10uA = ',a_1,'[mohm]')
	a_2 =  np.exp(RN_fit_2)+RN_s
	print ('RN_fit_20uA = ',a_2,'[mohm]')
else:
	a_1 = ((RN_f-RN_s)/math.pi)*math.atan(RN_fit_1)+(RN_f+RN_s)/2
	print ('RN_fit_10uA = ',a_1,'[mohm]')
	a_2 = ((RN_f-RN_s)/math.pi)*math.atan(RN_fit_2)+(RN_f+RN_s)/2
	print ('RN_fit_20uA = ',a_2,'[mohm]')

if Tc_f == 'inf':
	b_1 = np.exp(Tc_fit_1)+Tc_s
	print ('Tc_fit_10uA = ',b_1,'[mK]')
	b_2 = np.exp(Tc_fit_2)+Tc_s
	print ('Tc_fit_20uA = ',b_2,'[mK]')
else:
	b_1 = ((Tc_f-Tc_s)/math.pi)*math.atan(Tc_fit_1)+(Tc_f+Tc_s)/2
	print ('Tc_fit_10uA = ',b_1,'[mK]')
	b_2 = ((Tc_f-Tc_s)/math.pi)*math.atan(Tc_fit_2)+(Tc_f+Tc_s)/2
	print ('Tc_fit_20uA = ',b_2,'[mK]')

if T1_f == 'inf':
	c_1 = np.exp(T1_fit_1)+T1_s
	print ('T1_fit_10uA = ',c_1,'[mK]')
	c_2 = np.exp(T1_fit_2)+T1_s
	print ('T1_fit_20uA = ',c_2,'[mK]')
else:
	c_1 = ((T1_f-T1_s)/math.pi)*math.atan(T1_fit_1)+(T1_f+T1_s)/2
	print ('T1_fit_10uA = ',c_1,'[mK]')
	c_2 = ((T1_f-T1_s)/math.pi)*math.atan(T1_fit_2)+(T1_f+T1_s)/2
	print ('T1_fit_20uA = ',c_2,'[mK]')

if T2_f == 'inf':
	d_1 = np.exp(T2_fit_1)+T2_s
	print ('T2_fit_10uA = ',d_1,'[mK]')
	d_2 = np.exp(T2_fit_2)+T2_s
	print ('T2_fit_20uA = ',d_2,'[mK]')
else:
	d_1 = ((T2_f-T2_s)/math.pi)*math.atan(T2_fit_1)+(T2_f+T2_s)/2
	print ('T2_fit_10uA = ',d_1,'[mK]')
	d_2 = ((T2_f-T2_s)/math.pi)*math.atan(T2_fit_2)+(T2_f+T2_s)/2
	print ('T2_fit_20uA = ',d_2,'[mK]')

x = np.arange(T[0]-2.0,T[-1]+2.0,0.1)
f_1 = a_1/((1+np.exp(-(x-b_1)/c_1))*(1+np.exp(-(x-b_1)/d_1)))
f_2 = a_2/((1+np.exp(-(x-b_2)/c_2))*(1+np.exp(-(x-b_2)/d_2)))

#--------------------転移幅を求める-----------------------------

RN_sup_1 = a_1*sup_R*0.01
RN_inf_1 = a_1*inf_R*0.01

T_inf_1 = []
T_sup_1 = []

for i in range(len(x)):
	if f_1[i] > RN_inf_1:
		T_inf_1.append(x[i])
		break

for i in range(len(x)):
	if f_1[i] > RN_sup_1:
		T_sup_1.append(x[i])
		break

print(T_sup_1)
print(T_inf_1)

#deltaT_1 = T_sup_1[0] - T_inf_1[0]
#print 'deltaT_1 = ',deltaT_1,'[mK]'

RN_sup_2 = a_2*sup_R*0.01
RN_inf_2 = a_2*inf_R*0.01

T_inf_2 = []
T_sup_2 = []

for i in range(len(x)):
	if f_2[i] > RN_inf_2:
		T_inf_2.append(x[i])
		break

for i in range(len(x)):
	if f_2[i] > RN_sup_2:
		T_sup_2.append(x[i])
		break

#deltaT_2 = T_sup_2[0] - T_inf_2[0]
#print ('deltaT_2 = ',deltaT_2,'[mK]')

#-----------------alphaを求める---------------------------------

y_1 = np.linspace(T_inf_1[0],T_sup_1[0],1000)
y_2 = np.linspace(T_inf_2[0],T_sup_2[0],1000)

alpha_1 = []
R_alpha_1 = []
T_alpha_1 = []
diff_R_1 =[]
z_1 = []

for i in range(len(y_1)-1):
	T_alpha_1.append(y_1[i])
	R_alpha_1.append(a_1/((1+np.exp(-(y_1[i]-b_1)/c_1))*(1+np.exp(-(y_1[i]-b_1)/d_1))))
	diff_R_1.append(((a_1/((1+np.exp(-(y_1[i+1]-b_1)/c_1))*(1+np.exp(-(y_1[i+1]-b_1)/d_1))))-R_alpha_1[i])/(y_1[i+1]-y_1[i]))
	alpha_1.append((T_alpha_1[i]*diff_R_1[i])/R_alpha_1[i])
	z_1.append((100.0*R_alpha_1[i])/a_1)

alpha_2 = []
R_alpha_2 = []
T_alpha_2 = []
diff_R_2 =[]
z_2 = []

for i in range(len(y_2)-1):
	T_alpha_2.append(y_2[i])
	R_alpha_2.append(a_2/((1+np.exp(-(y_2[i]-b_2)/c_2))*(1+np.exp(-(y_2[i]-b_2)/d_2))))
	diff_R_2.append(((a_2/((1+np.exp(-(y_2[i+1]-b_2)/c_2))*(1+np.exp(-(y_2[i+1]-b_2)/d_2))))-R_alpha_2[i])/(y_2[i+1]-y_2[i]))
	alpha_2.append((T_alpha_2[i]*diff_R_2[i])/R_alpha_2[i])
	z_2.append((100.0*R_alpha_2[i])/a_2)

#-----------------ファイルに出力----------------------------------

f = open("output/RT_data.txt","w")
f.write('#10[uA]'+'\t'+'20[uA]'+'\n')
f.write('#RN[mohm]'+'\t'+str(a_1)+'\t'+str(a_2)+'\n')
f.write('#Tc[mK]'+'\t'+str(b_1)+'\t'+str(b_2)+'\n')
f.write('#T1[mK]'+'\t'+str(c_1)+'\t'+str(c_2)+'\n')
f.write('#T2[mK]'+'\t'+str(d_1)+'\t'+str(d_2)+'\n')
f.flush()
f.close()

#-----------------グラフとして出力----------------------------------

plt.title('RT_10[uA]')
plt.xlim(T[0]-2.0,T[-1]+2.0)
plt.xlabel('Temperature[mK]',fontsize = 16)
plt.ylabel('Resistance[m$\Omega$]',fontsize = 16)
plt.plot(T,R_1,'o',color ='lightsalmon',label ='10uA',markersize = 6)
plt.plot(x,f_1,color ='black',label = 'fitting',linestyle ='-',linewidth = 1)
plt.grid(True)
plt.legend(loc ='best',fancybox = True,shadow = True)
plt.savefig("output/RTfit_10uA.png",format = 'png')
plt.show()

plt.title('RT_20[uA]')
plt.xlim(T[0]-2.0,T[-1]+2.0)
plt.xlabel('Temperature[mK]',fontsize = 16)
plt.ylabel('Resistance[m$\Omega$]',fontsize = 16)
plt.plot(T,R_2,'o',color ='lightsalmon',label ='20uA',markersize = 6)
plt.plot(x,f_2,color ='black',label = 'fitting',linestyle ='-',linewidth = 1)
plt.grid(True)
plt.legend(loc ='best',fancybox = True,shadow = True)
plt.savefig("output/RTfit_20uA.png",format = 'png')
plt.show()

#----------重ねて表示-------------
plt.title('ch1:(RT_10[uA], RT_20[uA])')
plt.xlim(T[0]-2.0,T[-1]+2.0)
plt.xlabel('Temperature[mK]',fontsize = 16)
plt.ylabel('Resistance[m$\Omega$]',fontsize = 16)
plt.plot(T,R_1,'o',color ='red',label ='10uA',markersize = 5)
plt.plot(x,f_1,color ='red',label = 'fitting',linestyle ='-',linewidth = 2)
plt.plot(T,R_2,'o',color ='blue',label ='20uA',markersize = 5)
plt.plot(x,f_2,color ='blue',label = 'fitting',linestyle ='-',linewidth = 2)
plt.grid(True)
plt.legend(loc ='best',fancybox = True,shadow = True)
plt.savefig("output/RTfit_All.png",format = 'png')
plt.show()


plt.title('alpha_10[uA]')
plt.xlim(inf_R,sup_R)
plt.xlabel('Bias Point[%]',fontsize = 16)
plt.ylabel('alpha[-]',fontsize = 16)
plt.plot(z_1,alpha_1,color ='lightsalmon',label ='10uA',markersize = 6)
plt.grid(True)
plt.legend(loc ='best',fancybox = True,shadow = True)
plt.savefig("output/alpha_10uA.png",format = 'png')
plt.show()

plt.title('alpha_20[uA]')
plt.xlim(inf_R,sup_R)
plt.xlabel('Bias Point[%]',fontsize = 16)
plt.ylabel('alpha[-]',fontsize = 16)
plt.plot(z_2,alpha_2,color ='lightsalmon',label ='20uA',markersize = 6)
plt.grid(True)
plt.legend(loc ='best',fancybox = True,shadow = True)
plt.savefig("output/alpha_20uA.png",format = 'png')
plt.show()

plt.title('ch1:(alpha_10[uA],alpha_20[uA])')
plt.xlim(inf_R,sup_R)
plt.xlabel('Bias Point[%]',fontsize = 16)
plt.ylabel('alpha[-]',fontsize = 16)
plt.plot(z_1,alpha_1,color ='red',label ='10uA',markersize = 6)
plt.plot(z_2,alpha_2,color ='blue',label ='20uA',markersize = 6)
plt.grid(True)
plt.legend(loc ='best',fancybox = True,shadow = True)
plt.savefig("output/alpha_All.png",format = 'png')
plt.show()

# -*- coding: utf-8 -*-

import math
import glob
from sys import argv
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib as mpl
import os.path

#---------------------------変更するparameters-------------------------

R_n = 0.11					#常伝導抵抗値[ohm]
R_sh = 3.9e-3				#シャント抵抗値[ohm]
PoR_n = 50.0				#Percentage of R_n[-]
							#常伝導抵抗値の何％のジュール発熱を得たいのか。0 < PofR_n < 100を満たす任意の有理数。
T_bath_list = [110,115,120,125,130,135,140]					#使用する熱浴温度のリスト[mK]

G_s = 0.0					#fittng時のパラメータの取り得る範囲。G[W/K],T_c[mK],n[-]
G_f = 'inf'					#例えば、パラメタGを0[W/K]以上で動かしたい場合は、G_s = 0.0,G_f = 'inf'。
T_c_s = 120				#同様に、パラメタT_cを150~160[mK]で動かしたい場合は、T_c_s = 150.0,T_c_f = 160.0。
T_c_f = 160			#nについても同様。
#n_s = 4.0	#nは3,4,5でグラフを重ねるプログラム
#n_f = 4.0

para0 = [1.0,1.0,1.0]		#fitting時のパラメタ列の初期値。[G,T_c,n]であるが、T_cとnにおいては、0で割ることは許されていない。

#------------------------------------------------------------------

T_n = len(T_bath_list)					#熱浴リスト内の温度数
R_set = R_n*PoR_n*1.0e-2				#常伝導抵抗値のPoR_n％は、R_set[ohm]
print( 'List of Hot Bath =',T_bath_list,'[mK]')
print( 'The Number of Hot Bath =',T_n)
print( 'Operating Point =',PoR_n,'%')
print( 'Resistance Value at Operating Point =',R_set,'[ohm]')

filename_1 = 'IV_'            		#ファイル名の始め
filename_2 = 'mK.txt'               	#ファイル名の後ろその１

I_b = []
V_out = []
R_tes = []
eta = []

for i in range(T_n):
	filename_a = filename_1 + str(int(T_bath_list[i])) + filename_2
	data_a = np.loadtxt(filename_a,comments='#',skiprows = 1)
	data_sub1 = data_a[:,0]			#各熱浴温度におけるバイアス電流を抜き出す。
	I_b.append(data_sub1)
	data_sub2 = data_a[:,1]
	V_out.append(data_sub2)			#各熱浴温度における出力電圧を抜き出す。
	data_sub3 = data_a[:,2]
	R_tes.append(data_sub3)			#各熱浴温度におけるTESの抵抗値を抜き出す。
	data_sub4 = data_a[:,3]
	eta.append(data_sub4[0])				#各熱浴温度におけるetaを抜き出す。

#print (R_tes)

def getNearestValue(list, num):						#関数の定義。listの中から、ある値（num）に一番近い値を返す。
	idx = np.abs(np.asarray(list) - num).argmin()	#ただし、最も近い値が複数存在すれば、始めの値を返す。
	return list[idx]

R_tes_Percent = []
num_R_tes_Percent = []

for i in range(T_n):
	R_tes_list = R_tes[i].tolist()								#ndarrayをlistに変換
	R_tes_Percent.append(getNearestValue(R_tes_list, R_set))	#R_setに最も近い抵抗値をndarray化
	a = R_tes_list.index(getNearestValue(R_tes_list, R_set))	#R_setに最も近い抵抗値が何番目にあるのか
	num_R_tes_Percent.append(a)

#print (R_tes_Percent)
#print (num_R_tes_Percent)

P_j = []

for i in range(T_n):
	j = num_R_tes_Percent[i]
	P_j.append(eta[i]*V_out[i][j]*R_sh*(I_b[i][j]-eta[i]*V_out[i][j]))		#P_j[pW]

#print(P_j)

T_bath_list_K = np.array(T_bath_list)*1.0e-3		#熱浴温度を[mK]から[K]に変換。
P_j_W = np.array(P_j)*1.0e-12						#ジュール発熱を[pW]から[W]に変換。

#print(T_bath_list_K)
#print(P_j_W)

def fit_func(para,x,y):               #fit関数の定義。paraはパラメタ列。
	G = para[0]
	T_c = para[1]
	n = para[2]
	if G_f == 'inf':
		G_func = abs(G)+ G_s
	else:
		G_func = ((G_f-G_s)/math.pi)*math.atan(G)+(G_f+G_s)/2
	if T_c_f == 'inf':
		T_c_func = abs(T_c)+ T_c_s*1.0e-3
	else:
		T_c_func = ((T_c_f*1.0e-3-T_c_s*1.0e-3)/math.pi)*math.atan(T_c)+(T_c_f*1.0e-3+T_c_s*1.0e-3)/2
	if 3 == 'inf':
		n_func = abs(n)+ 3
	else:
		n_func = ((3-3)/math.pi)*math.atan(n)+(3+3)/2
	residual = y-(G_func*T_c_func*(1.0/n_func))*(1.0-(np.array(x)/T_c_func)**n_func)
	return(residual)

result = opt.leastsq(fit_func,para0,args = (T_bath_list_K,P_j_W))	#パラメタ列の初期値から始めて、(Ib_1,Vout_1)に対して,
																	#residualが最小になるようなa,bを返す
G_fit = result[0][0]
T_c_fit = result[0][1]
n_fit = result[0][2]

#print(G_fit)
#print(n_fit)
#print(T_c_fit)

if G_f == 'inf':
	a = abs(G_fit)+ G_s
	print( 'G_fit = ',a,'[W/K]')
else:
	a = ((G_f-G_s)/math.pi)*math.atan(G_fit)+(G_f+G_s)/2
	print( 'G_fit = ',a,'[W/K]')
if T_c_f == 'inf':
	b = abs(T_c_fit)+ T_c_s
	print( 'T_c_fit = ',b,'[mK]')
else:
	b = ((T_c_f-T_c_s)/math.pi)*math.atan(T_c_fit)+(T_c_f+T_c_s)/2
	print( 'T_c_fit = ',b,'[mK]')
if 3 == 'inf':
	c = abs(n_fit)+ 3
	print( 'n_fit = ',c,'[-]')
else:
	c = ((3-3)/math.pi)*math.atan(n_fit)+(3+3)/2
	print( 'n_fit = ',c,'[-]')

x0 = np.arange(np.array(T_bath_list)[0]-2.0,np.array(T_bath_list)[T_n-1]+2.0,0.00001)
f0 = (a*(b*1.0e-3)*(1.0/c))*(1.0-(np.array(x0)*1.0e-3/(b*1.0e-3))**c)

#--------------------------------------------------------------------------------------------------------------------------------
def fit_func(para,x,y):               #fit関数の定義。paraはパラメタ列。
	G = para[0]
	T_c = para[1]
	n = para[2]
	if G_f == 'inf':
		G_func = abs(G)+ G_s
	else:
		G_func = ((G_f-G_s)/math.pi)*math.atan(G)+(G_f+G_s)/2
	if T_c_f == 'inf':
		T_c_func = abs(T_c)+ T_c_s*1.0e-3
	else:
		T_c_func = ((T_c_f*1.0e-3-T_c_s*1.0e-3)/math.pi)*math.atan(T_c)+(T_c_f*1.0e-3+T_c_s*1.0e-3)/2
	if 4 == 'inf':
		n_func = abs(n)+ 4
	else:
		n_func = ((4-4)/math.pi)*math.atan(n)+(4+4)/2
	residual = y-(G_func*T_c_func*(1.0/n_func))*(1.0-(np.array(x)/T_c_func)**n_func)
	return(residual)

result = opt.leastsq(fit_func,para0,args = (T_bath_list_K,P_j_W))	#パラメタ列の初期値から始めて、(Ib_1,Vout_1)に対して,
																	#residualが最小になるようなa,bを返す
G_fit = result[0][0]
T_c_fit = result[0][1]
n_fit = result[0][2]

#print(G_fit)
#print(n_fit)
#print(T_c_fit)

if G_f == 'inf':
	a = abs(G_fit)+ G_s
	print( 'G_fit = ',a,'[W/K]')
else:
	a = ((G_f-G_s)/math.pi)*math.atan(G_fit)+(G_f+G_s)/2
	print( 'G_fit = ',a,'[W/K]')
if T_c_f == 'inf':
	b = abs(T_c_fit)+ T_c_s
	print( 'T_c_fit = ',b,'[mK]')
else:
	b = ((T_c_f-T_c_s)/math.pi)*math.atan(T_c_fit)+(T_c_f+T_c_s)/2
	print( 'T_c_fit = ',b,'[mK]')
if 4 == 'inf':
	c = abs(n_fit)+ 4
	print( 'n_fit = ',c,'[-]')
else:
	c = ((4-4)/math.pi)*math.atan(n_fit)+(4+4)/2
	print( 'n_fit = ',c,'[-]')

x1 = np.arange(np.array(T_bath_list)[0]-2.0,np.array(T_bath_list)[T_n-1]+2.0,0.00001)
f1 = (a*(b*1.0e-3)*(1.0/c))*(1.0-(np.array(x1)*1.0e-3/(b*1.0e-3))**c)
#-----------------------------------------------------------------------------------------------------------------------------------------
def fit_func(para,x,y):               #fit関数の定義。paraはパラメタ列。
	G = para[0]
	T_c = para[1]
	n = para[2]
	if G_f == 'inf':
		G_func = abs(G)+ G_s
	else:
		G_func = ((G_f-G_s)/math.pi)*math.atan(G)+(G_f+G_s)/2
	if T_c_f == 'inf':
		T_c_func = abs(T_c)+ T_c_s*1.0e-3
	else:
		T_c_func = ((T_c_f*1.0e-3-T_c_s*1.0e-3)/math.pi)*math.atan(T_c)+(T_c_f*1.0e-3+T_c_s*1.0e-3)/2
	if 5 == 'inf':
		n_func = abs(n)+ 5
	else:
		n_func = ((5-5)/math.pi)*math.atan(n)+(5+5)/2
	residual = y-(G_func*T_c_func*(1.0/n_func))*(1.0-(np.array(x)/T_c_func)**n_func)
	return(residual)

result = opt.leastsq(fit_func,para0,args = (T_bath_list_K,P_j_W))	#パラメタ列の初期値から始めて、(Ib_1,Vout_1)に対して,
																	#residualが最小になるようなa,bを返す
G_fit = result[0][0]
T_c_fit = result[0][1]
n_fit = result[0][2]

#print(G_fit)
#print(n_fit)
#print(T_c_fit)

if G_f == 'inf':
	a = abs(G_fit)+ G_s
	print( 'G_fit = ',a,'[W/K]')
else:
	a = ((G_f-G_s)/math.pi)*math.atan(G_fit)+(G_f+G_s)/2
	print( 'G_fit = ',a,'[W/K]')
if T_c_f == 'inf':
	b = abs(T_c_fit)+ T_c_s
	print( 'T_c_fit = ',b,'[mK]')
else:
	b = ((T_c_f-T_c_s)/math.pi)*math.atan(T_c_fit)+(T_c_f+T_c_s)/2
	print( 'T_c_fit = ',b,'[mK]')
if 5 == 'inf':
	c = abs(n_fit)+ 5
	print( 'n_fit = ',c,'[-]')
else:
	c = ((5-5)/math.pi)*math.atan(n_fit)+(5+5)/2
	print( 'n_fit = ',c,'[-]')

x2 = np.arange(np.array(T_bath_list)[0]-2.0,np.array(T_bath_list)[T_n-1]+2.0,0.00001)
f2 = (a*(b*1.0e-3)*(1.0/c))*(1.0-(np.array(x2)*1.0e-3/(b*1.0e-3))**c)
#------------------------------------------------------------------------------------------------------------------------------------

P_j_nW = np.array(P_j_W)*1.0e+9

a_round = round(a,9)
b_round = round(b,2)
c_round = round(c,3)

G0 = a/(b*1.0e-3)**(c-1.0)
print( 'G0 = ',G0,'[W/K]')

G0_round = round(G0,8)

y = np.arange(80.0,150.0,0.0001)
G_T = G0*(y*1.0e-3)**(c-1.0)


#----------------ファイルに出力--------------------------------

folder = os.path.exists('./output_n=3,4,5(50%)_test')
if not folder:
	os.mkdir('output_n=3,4,5(50%)_test')
os.chdir('./output_n=3,4,5(50%)_test')

#plt.title('G = '+str(a_round)+' [W/K]'+'\nT$_c$ = '+str(b_round)+' [mK]'+' , n = '+str(c_round)+' [-]')
plt.plot(T_bath_list,P_j_nW,'o',color ='red',label = 'Experiment',markersize = 7)
plt.plot(x0,f0*1.0e+9,color ='black',label = 'n = 3',linestyle ='-',linewidth = 2)
plt.plot(x1,f1*1.0e+9,color ='blue',label ='n = 4',linestyle ='--',linewidth = 2)
plt.plot(x2,f2*1.0e+9,color ='green',label ='n = 5',linestyle =':',linewidth = 3)
plt.xlabel(r"$T_{\rm{bath}}$ [mK]",fontsize = 20)
plt.ylabel(r"$P_{\rm{J}}$ [nW]",fontsize = 20)
plt.ylim(ymin=0)
plt.grid(True)
plt.tick_params(labelsize=16)
plt.subplots_adjust(left=0.2,bottom=0.15)
plt.legend(loc ='best',fancybox = True,shadow = True,numpoints = 1)
plt.savefig("P-T.pdf",format = 'pdf')
plt.savefig("P-T.png",format = 'png')
plt.show()
"""
plt.title('G(T) = '+str(G0_round)+' T$^n$$^-$$^1$'+' , n = '+str(c_round)+' [-]')
plt.plot(b,a*1.0e+9,'o',color ='red',label = 'G at T$_c$',markersize = 7)
plt.plot(y,G_T*1.0e+9,color ='black',linestyle ='-',linewidth = 2)
plt.xlabel('Temperature [mK]',fontsize = 20)
plt.ylabel('G [nW/K]',fontsize = 20)
plt.ylim(ymin=0)
plt.grid(True)
plt.tick_params(labelsize=16)
plt.subplots_adjust(left=0.2,bottom=0.15)
plt.legend(loc ='best',fancybox = True,shadow = True,numpoints = 1)
plt.savefig("G-T.pdf",format = 'pdf')
plt.savefig("G-T.png",format = 'png')
plt.show()

f = open('data.txt','w')
f.write('List of Hot Bath = '+str(T_bath_list)+' [mK]\n')
f.write('The Number of Hot Bath = '+str(T_n)+' [-]\n')
f.write('Operating Point = '+str(PoR_n)+' %\n')
f.write('Resistance Value at Operating Point = '+str(R_set)+' [ohm]\n')
f.write('G_fit = '+str(a)+' [W/K]\n')
f.write('T_c_fit = '+str(b)+'[mK]\n')
f.write('n_fit = '+str(c)+' [-]\n')
f.write('G0_fit = '+str(G0)+' [W/K]')
f.flush()
f.close()

f = open('P-T.txt','w')
f.write('#Tbath[mK]'+'\t'+'P[nW]'+'\n')
for i in range(len(T_bath_list)):
	f.write(str(T_bath_list[i])+'\t'+ str(P_j_nW[i])+'\n')
f.flush()
f.close()

os.chdir('..')
"""
#---------------------------------------------------------------------

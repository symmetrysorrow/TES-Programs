import math
import glob
from sys import argv
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import getpara as gp
from natsort import natsorted
import scipy.optimize as opt
import sys

path = "/Volumes/Extreme Pro/matsumi/data/20230616_IV/room1-ch3"
R_n = 0.023
R_SH = 3.9e-3
PoR_n = 30.0
eta = 87
remove = [160,200,210]
G_s = 0.0					#fittng時のパラメータの取り得る範囲。G[W/K],T_c[mK],n[-]
G_f = 'inf'					#例えば、パラメタGを0[W/K]以上で動かしたい場合は、G_s = 0.0,G_f = 'inf'。
T_c_s = 180				#同様に、パラメタT_cを150~160[mK]で動かしたい場合は、T_c_s = 150.0,T_c_f = 160.0。
T_c_f = 200				#nについても同様。
n_s = 3.0					#注意として、全てのパラメタで下限は0.0である。
n_f = 6.0		#またT_cについては、相転移温度±5.0[mK]くらいが良い。

p0 = [1.0e-8,0.190,3]		#fitting時のパラメタ列の初期値。[G,T_c,n]であるが、T_cとnにおいては、0で割ることは許されていない。


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
	if n_f == 'inf':
		n_func = abs(n)+ n_s
	else:
		n_func = ((n_f-n_s)/math.pi)*math.atan(n)+(n_f+n_s)/2
	residual = y-(G_func*T_c_func*(1.0/n_func))*(1.0-(np.array(x)/T_c_func)**n_func)
	return(residual)


def main():
	temps = natsorted(glob.glob(f'{path}/calibration/IV_*mK.txt'))

	R_set = R_n*PoR_n/100
	P_j_list = []
	T_bath_list = []

	for i in temps:
		data = np.loadtxt(i)
		I_bias = data[0]
		V_out = data[1]
		I_tes = eta * V_out
		I_sh = I_bias - I_tes
		V_tes = I_sh * R_SH
		R_tes = V_tes[1:] / I_tes[1:]
		R_tes =  np.append(0.0,R_tes)
		T_bath = int(gp.num(os.path.basename(i))[0])

		if T_bath in remove:
			continue
		r = np.argmax(R_tes >= R_set) # the most close R_tes to bias point 
	
		P_j = eta*V_out[r]*R_SH*(I_bias[r]-eta*V_out[r]) #[pW]
		P_j_list.append(P_j)
		T_bath_list.append(T_bath)
	
		

	T_bath_list = np.array(T_bath_list)*1e-3 #[mK]
	P_j_list = np.array(P_j_list)*1e-12 #[W]
	
	result = opt.leastsq(fit_func,p0,args = (T_bath_list,P_j_list))
	G_fit = result[0][0]
	T_c_fit = result[0][1]
	n_fit = result[0][2]

	if G_f == 'inf':
		a = abs(G_fit)+ G_s
		print ('G_fit = ',a,'[W/K]')
	else:
		a = ((G_f-G_s)/math.pi)*math.atan(G_fit)+(G_f+G_s)/2
		print ('G_fit = ',a,'[W/K]')
	if T_c_f == 'inf':
		b = abs(T_c_fit)+ T_c_s
		print('T_c_fit = ',b,'[mK]')
	else:
		b = ((T_c_f-T_c_s)/math.pi)*math.atan(T_c_fit)+(T_c_f+T_c_s)/2
		print('T_c_fit = ',b,'[mK]')
	if n_f == 'inf':
		c = abs(n_fit)+ n_s
		print('n_fit = ',c,'[-]')
	else:
		c = ((n_f-n_s)/math.pi)*math.atan(n_fit)+(n_f+n_s)/2
		print('n_fit = ',c,'[-]')

	x = np.arange(0.15,0.21,0.0001)
	f = (a*(b*1.0e-3)*(1.0/c))*(1.0-(np.array(x)/(b*1.0e-3))**c)

	P_j_nW = np.array(P_j_list)*1.0e+9

	G0 = a/(b*1.0e-3)**(c-1.0)

	y = np.arange(160.0,210.0,0.0001)
	G_T = G0*(y*1.0e-3)**(c-1.0)

	#plt.title('G = '+str(a_round)+' [W/K]'+'\nT$_c$ = '+str(b_round)+' [mK]'+' , n = '+str(c_round)+' [-]',fontsize=16)
	plt.plot(T_bath_list*1e3,P_j_nW,'o',color ='lightsalmon',label = 'Experiment',markersize = 6)
	plt.plot(x*1e3,f*1.0e+9,color ='black',label = 'Fitting',linestyle ='-',linewidth = 1)
	plt.xlabel(r"$T_{\rm{bath}}$ [mK]",fontsize = 16)
	plt.ylabel(r"$P_{\rm{J}}$ [nW]",fontsize = 16)
	plt.grid()
	plt.legend(fontsize=12)
	plt.tight_layout()
	plt.savefig(f"{path}/P-T.pdf",format = 'pdf')
	plt.show()

	#plt.title('G(T) = '+str(G0_round)+' T$^n$$^-$$^1$'+' , n = '+str(c_round)+' [-]',fontsize=16)
	plt.plot(b,a*1.0e+9,'o',color ='lightsalmon',label = 'G at T$_c$',markersize = 6)
	plt.plot(y,G_T*1.0e+9,color ='black',linestyle ='-',linewidth = 1)
	plt.xlabel('Temperature [mK]',fontsize = 16)
	plt.ylabel('G [nW/K]',fontsize = 16)
	plt.grid()
	plt.legend(fontsize=12)
	plt.tight_layout()
	plt.savefig(f"{path}/G-T.pdf",format = 'pdf')
	plt.show()

	
main()
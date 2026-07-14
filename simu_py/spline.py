import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

def spline_model(x, a, b):
    return a * x + b

def FindPosition(Fit_para,ratios):
    print("A")
    B_right = Fit_para[:, 1]
    B_left = Fit_para[:, 0]

    result=[]

    

def Reso(Data_path):
    fit_para_path=f"f:/hata/662_142_136/ratios.txt"
    fit_para=np.loadtxt(fit_para_path, delimiter=',')

    ratio_base = fit_para[:, 1]
    posi_base = fit_para[:, 0]

    sorted_indices = np.argsort(ratio_base)
    ratio_base = ratio_base[sorted_indices]
    posi_base = posi_base[sorted_indices]

    spline = UnivariateSpline(ratio_base, posi_base, k=3, s=0)

    value = [(float(i)) / float(len(posi_base) + 1) for i in range(len(posi_base))]

    x_new = np.array([2.5,2])
    y_new = spline(x_new)

    x_fine = np.linspace(min(ratio_base), max(ratio_base), 500)  # スムーズな線を描画するためのx値
    y_fine = spline(x_fine)

    plt.scatter(ratio_base, posi_base, label="Data", c=value, cmap="hsv")
    plt.plot(x_fine, y_fine, label="Spline Fit", color="gray", linestyle='--')
    #plt.scatter(x_new, y_new, label=f"Interpolated y", color="green")
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Spline Interpolation")
    plt.show()

Reso("A")

    

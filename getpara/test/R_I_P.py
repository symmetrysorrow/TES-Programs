import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

R_sh = 3.9e-3
file = sorted(glob.glob("IV_140mK.txt"))

if not os.path.exists("./output"):
    os.mkdir("./output")
else:
    shutil.rmtree("./output")
    os.mkdir("./output")

for i in file:
    with open(i,"r") as f:
        data = np.loadtxt(f,skiprows=1)
        I_b = data[:,0]
        V_out = data[:,1]
        R_tes = data[:,2]
        eta = data[:,3]
        P_j = P_j = eta*V_out*R_sh*(I_b-eta*V_out)

        path = i.replace("IV","IP").replace("txt","png")

        plt.title("I_b vs P_j")
        plt.plot(I_b,P_j,'o',color="red",markersize=7)
        plt.xlabel("I_b(uA)",fontsize=10)
        plt.ylabel("P_j(J)",fontsize=10)
        plt.grid()
        plt.savefig("./output/"+path)
        plt.clf()

        f.close()
#-----------------------------------------------------------------------
for i in file:
    with open(i,"r") as f:
        data = np.loadtxt(f,skiprows=1)
        I_b = data[:,0]
        V_out = data[:,1]
        R_tes = data[:,2]
        eta = data[:,3]
        P_j = P_j = eta*V_out*R_sh*(I_b-eta*V_out)

        path = i.replace("IV","IR").replace("txt","png")

        plt.title("R_tes vs I_b")
        plt.plot(R_tes,I_b,'o',color="red",markersize=7)
        plt.xlabel("R_tes",fontsize=10)
        plt.ylabel("I_b",fontsize=10)
        plt.grid()
        plt.savefig("./output/"+path)
        plt.clf()

        f.close()


for i in file:
    with open(i,"r") as f:
        data = np.loadtxt(f,skiprows=1)
        I_b = data[:,0]
        V_out = data[:,1]
        R_tes = data[:,2]
        eta = data[:,3]
        P_j = P_j = eta*V_out*R_sh*(I_b-eta*V_out)

        path = i.replace("IV","RP").replace("txt","png")

        plt.title("R_tes vs P_j")
        plt.plot(R_tes,P_j,'o',color="red",markersize=7)
        plt.xlabel("R_tes",fontsize=10)
        plt.ylabel("P_j",fontsize=10)
        plt.grid()
        plt.savefig("./output/"+path)
        plt.clf()

        f.close()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import getpara as gp

num = input("num: ")
path = f"/Volumes/Untitled/scope_{num}.csv"

data = np.loadtxt(path,skiprows=2,delimiter=',').T

print(data)

ch1 = data[1]
ch2 = data[2]
print(ch1)
base1,ch1 = gp.baseline(ch1,1000,100,100)
base2,ch2 = gp.baseline(ch2,1000,100,100)

plt.plot(ch1)
plt.plot(ch2)
plt.show()
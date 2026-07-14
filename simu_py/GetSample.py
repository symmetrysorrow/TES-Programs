import numpy as np

path=input("dat path:")
data = np.loadtxt(path, comments="#", dtype=np.float64)
print(f"読み込んだサンプル数: {data.size}")
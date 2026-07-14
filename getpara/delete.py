import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import libs.getpara as gp
import sys
import shutil
import pprint
import json
import glob
import re


path = "E:/tsuruta/20230616/room1-ch2-3_180mK_570uA_100kHz_g10/CH0_pulse/output/2/energy_resolution/delete.txt"

numbers = []
with open(path,"r") as f:
    for i in f.readlines():
        num = int(i.split(":")[1])
        numbers.append(num)
print(numbers)

set = gp.loadJson()
ch = set["Config"]["channel"]
os.chdir(set["Config"]["path"])

output = f'CH{set["Config"]["channel"]}_pulse/output/{set["Config"]["output"]}'
df = pd.read_csv((f'{output}/output.csv'),index_col=0)

for num in numbers:
    print(num)
    index = f"CH{ch}_pulse/rawdata\CH{ch}_{num}.dat"
    df.at[index,"quality"] = 0

df.to_csv(f'{output}/output.csv')
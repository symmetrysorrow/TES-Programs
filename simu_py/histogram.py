import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from InquirerPy import inquirer

folder=input("target folder: ")

channels=[0,1]
columns = ['height', 'rise', 'ST_Height']

ch=inquirer.select(
    message="Channel:",
    choices=channels
).execute()
target = inquirer.select(
    message="target:",
    choices=columns
).execute()

df=pd.read_csv(folder+f"/output_TES{ch}.csv",header=0)
selectedId=np.loadtxt(f"{folder}/selected_ids.txt")
df = df[df['id'].isin(selectedId)]

data=df[target].values
plt.hist(data,bins=30)
plt.show()
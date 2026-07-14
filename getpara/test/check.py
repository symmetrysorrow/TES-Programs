
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import shutil
from natsort import natsorted
import json
import pprint
import sys
import re

def main():
    # Get Setting
    set = gp.loadJson()
    
    path = set["Config"]["path"]
    ch = set["Config"]["channel"]
    os.chdir(path)
    trig_ch = np.loadtxt('channel.txt')
    output = f'CH{set["Config"]["channel"]}_pulse/output/{set["Config"]["output"]}'
    df = pd.read_csv((f'{output}/output.csv'),index_col=0)
    path_data = natsorted(glob.glob(f'CH{ch}_pulse/rawdata/CH{ch}_*.dat'))
    cnt=0
    for index,row in df.iterrows():
        num =  re.findall(r'\d+', os.path.basename(index))[1]
        df.at[index,"trigger"] = trig_ch[int(num)-1]
        print(num)

    df.to_csv(os.path.join(output,"output.csv"))
        
if __name__ == '__main__':
    main()
    
    


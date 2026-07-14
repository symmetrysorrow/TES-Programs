import numpy as np
import getpara as gp
import sys

def main():
    parent_path = '/Volumes/Extreme Pro/matsumi/data/20230703_boomeran/180mK_550uA_Difftrig5E-5_100kHz_gain10/CH0_pulse/output/run11/alpha/selected_index.txt'
    child_path = '/Volumes/Extreme Pro/matsumi/data/20230703_boomeran/180mK_550uA_Difftrig5E-5_100kHz_gain10/CH0_pulse/output/run11/alpha/selected_indexのコピー.txt'
    parent = gp.loadIndex(parent_path)
    child = gp.loadIndex(child_path)
    
    remain = []
    for i in parent:
        if i in child:
            continue

        remain.append(int(i))
    print(remain)
    np.savetxt('/Volumes/Extreme Pro/matsumi/data/20230703_boomeran/180mK_550uA_Difftrig5E-5_100kHz_gain10/CH0_pulse/output/run11/alpha/selected_index_remain.txt',remain,fmt="%s")

main()
    
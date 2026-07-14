import os
import glob
import numpy as np
from natsort import natsorted
import libs.getpara as gp





def main():
    set = gp.loadJson()
    path = input('Path: ')
    ch = input('Ch: ')
   
    os.chdir(path) 

    files = natsorted(glob.glob("*.dat"))
    print(files)
    num = 0
    
    for i in files:
        os.rename(i,f"CH{ch}_{num}.dat")
        num += 1
        print(i)

    print(files)
    

if __name__ == "__main__":
    main()






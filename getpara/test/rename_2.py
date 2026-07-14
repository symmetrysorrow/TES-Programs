import os
import glob
import numpy as np
from natsort import natsorted
import libs.getpara as gp





def main():
    set = gp.loadJson()
    path = input('Path: ')
    ch = input('Ch: ')
    data = input('pulse[0], noise[1]: ')
    if data == '0':
        os.chdir(f'{path}/CH{ch}_pulse/rawdata') 
    if data == '1':
        os.chdir(f'{path}/CH{ch}_noize/rawdata') 

    files = natsorted(glob.glob("*.dat"))

    start = int(input("start: "))
    number = np.arange(start,start + len(files),1)

    # ファイル名を変更時、ディレクトリ内に同名のファイルを作らないように仮名に変更
    num = 0
    for i in files:
        os.rename(i,f"{number[num]}.dat")
        num += 1
        print(i)
    
    # 正しいファイルをに変更
    files = natsorted(glob.glob("*.dat")) 
    num = 0
    for i in files:
        os.rename(i,f"CH{ch}_{number[num]}.dat")
        num += 1
        print(i)

    print(files)

if __name__ == "__main__":
    main()






import os
import glob
import natsort
import getpara as gp
import json

def main():

    with open("setting.json") as f:
        set = json.load(f)
    set_config = set["Config"]
    path = set_config["path"]
    ch = set_config["channel"]
    os.chdir(f"{path}/CH{ch}_pulse/rawdata")
    rawdatas = natsort.natsorted(glob.glob(f"CH{ch}_*.dat")) 
    file_num = len(rawdatas)
    print(rawdatas[-1])
    
    print(f"There are {file_num} files!\nbut......")
    num = 1
    count = 0
    for fle in rawdatas:
        if fle == f"CH{ch}_{num}.dat":
            num += 1
            continue
        else:
            print(fle)
            num = int(gp.num(fle)[1])
            count += 1

        num += 1
    print(f"{count} files don't exist!")



if __name__ == '__main__':
    main()
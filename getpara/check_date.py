
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
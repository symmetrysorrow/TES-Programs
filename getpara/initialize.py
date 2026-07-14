import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import shutil
from natsort import natsorted
import libs.getpara as gp
import libs.fft_spectrum as sp
import json
from tkinter import filedialog
import sys
import re

#---Initialize setting------------------------------------------------
main_para = {"main": {
        "base_x": 1000,
        "base_w": 500,
        "peak_max": 1000,
        "peak_x": 3,
        "peak_w": 10,
        "rise_high" : 0.90,
        "rise_low" : 0.10,
        "decay_high" : 0.90,
        "decay_low" :0.10,
        "area_x":1000,
        "area_w":5000,
        "fit_func": "monoExp",
        "fit_x": 5000,
        "fit_w": 50000,
        "fit_p0": [
            0.1,
            1e-05
        ],
        "mv_w": 0,
        "cutoff": 0.0
    }
}


main2_para = {"main2": {
        "base_x": 1000,
        "base_w": 500,
        "peak_max": 4000,
        "peak_x": 100,
        "peak_w": 500,
        "rise_high" : 0.90,
        "rise_low" : 0.10,
        "decay_high" : 0.90,
        "decay_low" :0.10,
        "area_x":1000,
        "area_w":5000,
        "fit_func": "monoExp",
        "fit_x": 5000,
        "fit_w": 50000,
        "fit_p0": [
            0.1,
            1e-05
        ],
        "mv_w": 0,
        "cutoff": 0.0
    }
}



# Json形式へ変更

def main():

    path = filedialog.askdirectory()
    output = input("output name:")

    post_ch = glob.glob(f'{path}/CH*_*')
    print(f"{path}/setting.txt")
    try:
        setting = np.loadtxt(f"{path}/setting.txt",skiprows = 10)
        setting_json = {
            "Config":{
                "path" : path,
                "channel":0,
                "rate" : int(setting[2]),
                "samples" : int(setting[4]),
                "presamples" : int(setting[5]),
                "threshold" : setting[6],
                "output" : output,
                "type": 'binary'
            }}
    except:
        print("Not exist setting.txt !")
        setting_json = {
            "Config":{
                "path" : path,
                "channel":0,
                "rate" : int(250000),
                "samples" : int(50000),
                "presamples" : int(5000),
                "threshold" : 0.03,
                "output" : output,
                "type": 'text'

            }
        }
        

    setting_json.update(main_para)
    if len(post_ch) > 2:
        setting_json.update(main2_para)


    #----- Select Data Parameter ----------------------------------------------
    para = {'select':{
	'samples-=':setting_json['Config']['samples'],
	'height->':setting_json['Config']['threshold'],
	'rSquared->':0,
	'error-!':0
	}}
    #--------------------------------------------------------------------------

    #----- Graugh Parameter ---------------------------------------------------
    xlim = setting_json["Config"]["samples"]/setting_json["Config"]["rate"]
    graugh_para = {"graugh": {
        "xlim->":0.0,
        "xlim-<":xlim,
        "log":False,
        "diff":False
        }}
    #--------------------------------------------------------------------------

    setting_json.update(para)
    setting_json.update(graugh_para)
    jsn = json.dumps(setting_json,indent=4)

    # current directry
    with open("setting.json", 'w') as file:
        file.write(jsn)

    for i in post_ch:
        os.makedirs(f'{i}/output/{output}',exist_ok=True)
        with open(f'{i}/output/{output}/setting.json', 'w') as file:
            file.write(jsn)

    RED = '\033[33m'
    END = '\033[0m'
    print(RED + "Hello GetPara!" + END)


if __name__ == "__main__":
    main()
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from natsort import natsorted
import  pandas as pd



def imshow(path):
    img = mpimg.imread(path)
    plt.imshow(img)
    plt.show()


os.chdir('/Users/matsumi/teion/data/20220902/CH1_150mK_710uA_trig0.1_gain10_10kHz')

df = pd.read_csv('CH1/output/output.csv',index_col=0)

path = natsorted(glob.glob('CH1/output/fig/CH1_*.png'))
#path = 'CH1/output/fig/CH1_1.png'


for i in path:
    print(os.path.basename(i))
    imshow(i)
    input('pulse -> 1 ,noise -> 2')
    
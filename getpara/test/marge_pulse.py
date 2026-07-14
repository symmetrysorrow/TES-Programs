import matplotlib.pyplot as plt
import getpara as gp
import numpy as np


def main():
    data1 = gp.loadbi(
        "/Volumes/Extreme 2TB/Matsumi/data/20231121/boomeran_room2-ch2-100mK_410uA_0.06V_g10_100kHz/CH0_pulse/rawdata/CH0_11.dat",
        "binary",
    )
    data2 = gp.loadbi(
        "/Volumes/Extreme 2TB/Matsumi/data/20231121/boomeran_room2-ch2-100mK_410uA_0.06V_g10_100kHz/CH0_pulse/rawdata/CH0_12.dat",
        "binary",
    )
    data3 = gp.loadbi(
        "/Volumes/Extreme 2TB/Matsumi/data/20231121/boomeran_room2-ch2-100mK_410uA_0.06V_g10_100kHz/CH0_pulse/rawdata/CH0_13.dat",
        "binary",
    )

    data_marge = np.concatenate([data1, data2, data3])

    plt.plot(data_marge,'o',markersize = 0.5)
    plt.show()


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

posis=[30,60,90,120,150]



INPUT_FILE = Path("PoST_Center_time.txt")
OUTPUT_FILE = Path("PoST_center_time_plot.png")


def load_post_center(path: Path) -> np.ndarray:
    """Load COMSOL text table while ignoring comment/header lines."""

    data = np.loadtxt(path, comments="%", encoding="utf-8")
    if data.ndim != 2 or data.shape[1] < 5:
        raise ValueError(
            "Expected at least 5 numeric columns: time, temperatures, and currents."
        )
    return data


def main() -> None:

    Times=[]
    Currents_Tes1=[]
    Currents_Tes2=[]

    for posi in posis:
        data = load_post_center(f"PoST_{posi}.txt")

        time_ms = data[:, 0]
        Times.append(time_ms)

        Temp_abs=data[:, 1]
        Temp_TES=data[:, 2]

        Current_TEs1=data[:, 3]
        Current_TEs2=data[:, 4]

        Currents_Tes1.append(Current_TEs1)
        Currents_Tes2.append(Current_TEs2)


    for posi in posis:
        plt.plot(Times[posis.index(posi)], Currents_Tes1[posis.index(posi)], label=f"{posi}")
    plt.xlabel("Time (ms)")
    plt.ylabel("Current (A)")
    plt.title("Current vs Time for Different Positions")
    plt.grid()
    plt.legend(title="Position")
    plt.savefig("TES1.png", dpi=300,transparent=True)
    plt.close()

    for posi in posis:
        plt.plot(Times[posis.index(posi)], Currents_Tes2[posis.index(posi)], label=f"{posi}")
    plt.xlabel("Time (ms)")
    plt.ylabel("Current (A)")
    plt.title("Current vs Time for Different Positions")
    plt.grid()
    plt.legend(title="Position")
    plt.savefig("TES2.png", dpi=300,transparent=True)
    plt.close()
        

    

if __name__ == "__main__":
    main()

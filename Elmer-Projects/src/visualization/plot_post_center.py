# -*- coding: utf-8 -*-

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_FILE = REPO_ROOT / "docs" / "Single-Pixel.txt"
OUTPUT_DIR = REPO_ROOT / "artifacts" / "plots"


def load_post_center(path: Path) -> np.ndarray:
    """Load COMSOL text table while ignoring comment/header lines."""

    data = np.loadtxt(path, comments="%", encoding="utf-8")
    if data.ndim != 2 or data.shape[1] < 5:
        raise ValueError(
            "Expected at least 5 numeric columns: time, temperatures, and currents."
        )
    return data


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_post_center(INPUT_FILE)

    time_ms = data[:, 0]

    Temp_abs=data[:, 1]
    Temp_Stycast=data[:, 2]
    Temp_TES=data[:, 3]


    Current_TES=data[:, 4]

    Resistance_TES=data[:, 5]

    plt.plot(time_ms, Temp_abs, label="Absorber")
    plt.plot(time_ms, Temp_TES, label="TES")
    plt.grid()
    plt.xlabel("Time (ms)")
    plt.ylabel("Temperature (K)")
    plt.title("Temperature vs Time")
    plt.legend()
    plt.savefig(OUTPUT_DIR / "Temperature_vs_Time.png")
    plt.close()

    plt.plot(time_ms, Current_TES, label="TES Current")
    plt.grid()
    plt.xlabel("Time (ms)")
    plt.ylabel("Current (A)")
    plt.title("TES Current vs Time")
    plt.savefig(OUTPUT_DIR / "TES_Current_vs_Time.png")
    plt.close()
        

    

if __name__ == "__main__":
    main()

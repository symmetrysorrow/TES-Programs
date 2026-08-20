# Steve Smith PhD p42~

import os
import numpy as np
import matplotlib.pyplot as plt
import plt_config
from numpy import linalg as LA
import math
import cmath
import sys
import matplotlib.cm as cm
import pandas as pd
import json
import time
import getpara as gp
import tqdm


k_b = 1.381e-23
ptfn_Flink = 0.5


def create_output_directry(*path):
    if path != None:
        os.makedirs(f"./output/{path[0]}", exist_ok=True)
        filenumber = path[0]
    else:
        for i in range(1000):
            folder = os.path.exists(f"./output/{str(i+1)}")
            if not folder:
                os.makedirs(f"./output/{str(i+1)}", exist_ok=True)
                filenumber = i + 1
                break
            else:
                filenumber = 0
                continue
    return str(filenumber)


def main():
    with open("h:/tagawa/sim/input.json", "r") as f:
        para = json.load(f)

    output_number = create_output_directry(para["output"])
    print(output_number)

    n_abs = para["n_abs"]

    # Total absorber heat capacity is divided among the thermal nodes.
    C_abs = para["C_abs"] / n_abs
    C_tes = para["C_tes"]

    # G_abs-abs in input.json is the end-to-end Pb conductance.
    # For N absorber nodes there are N-1 links in series.
    G_abs_abs = float(para["G_abs-abs"]) * (n_abs - 1)

    G_abs_tes = para["G_abs-tes"]
    G_tes_bath = para["G_tes-bath"]

    R = para["R"]
    R_l = para["R_l"]
    T_c = para["T_c"]
    T_bath = para["T_bath"]
    a = para["alpha"]
    b = para["beta"]
    L = para["L"]
    n = para["n"]
    E = para["E"]
    length = para["length"]
    rate = int(para["rate"])
    samples = int(para["samples"])

    output = output_number
    para["output"] = output

    time = np.linspace(0, samples / rate, samples)
    frequency = np.arange(0, rate, rate / samples)

    I = np.sqrt(
        (
            G_tes_bath
            * T_c
            * (1 - (T_bath / T_c) ** n)
        )
        / (n * R)
    )

    t_el = L / (R_l + R * (1 + b))
    L_I = (a * I**2 * R) / (G_tes_bath * T_c)
    t_I = C_tes / ((1 - L_I) * G_tes_bath)

    # ============================================================
    # Noise amplitudes
    # ============================================================

    ptfn_tes_bath = np.sqrt(
        4 * k_b * T_c**2 * G_tes_bath * ptfn_Flink
    )

    ptfn_abs_tes = np.sqrt(
        4 * k_b * T_c**2 * G_abs_tes * ptfn_Flink
    )

    ptfn_abs_abs = np.sqrt(
        4 * k_b * T_c**2 * G_abs_abs * ptfn_Flink
    )

    # Johnson noise
    enj = np.sqrt(
        4 * k_b * T_c * R * (1 + 2 * b + b**2)
    )

    enj_R = np.sqrt(
        4 * k_b * T_bath * R_l
    )

    # ============================================================
    # Noise source matrix
    #
    # Each ROW before transpose corresponds to one independent
    # physical noise source.
    #
    # For a thermal link, the same fluctuation is applied with
    # opposite signs to the two nodes joined by the link.
    # ============================================================

    def matrix_N(n_abs):
        X = np.zeros(
            (n_abs + 7, n_abs + 4),
            dtype=np.complex128,
        )

        for i in range(n_abs + 7):

            # TES1 Johnson
            if i == 0:
                X[i, 0] = -enj / L
                X[i, 1] = I * enj / C_tes

            # Load1 Johnson
            elif i == 1:
                X[i, 0] = enj_R / L

            # TES1-bath TFN
            elif i == 2:
                X[i, 1] = ptfn_tes_bath / C_tes

            # TES1-absorber TFN
            elif i == 3:
                X[i, 1] = +ptfn_abs_tes / C_tes
                X[i, 2] = -ptfn_abs_tes / C_abs

            # TES2-absorber TFN
            elif i == n_abs + 3:
                X[i, n_abs + 1] = -ptfn_abs_tes / C_abs
                X[i, n_abs + 2] = +ptfn_abs_tes / C_tes

            # TES2-bath TFN
            elif i == n_abs + 4:
                X[i, n_abs + 2] = ptfn_tes_bath / C_tes

            # Load2 Johnson
            elif i == n_abs + 5:
                X[i, n_abs + 3] = enj_R / L

            # TES2 Johnson
            elif i == n_abs + 6:
                X[i, n_abs + 3] = -enj / L
                X[i, n_abs + 2] = I * enj / C_tes

            # Internal absorber-absorber TFN
            else:
                left_node = i - 2
                right_node = i - 1

                X[i, left_node] = +ptfn_abs_abs / C_abs
                X[i, right_node] = -ptfn_abs_abs / C_abs

        return X

    # ============================================================
    # Electrothermal matrix
    # ============================================================

    def matrix_M(n_abs, omega):
        X = np.zeros(
            (n_abs + 4, n_abs + 4),
            dtype=np.complex128,
        )

        for i in range(n_abs + 4):

            # TES1 electrical
            if i == 0:
                X[i, 0] = 1 / t_el + 1j * omega
                X[i, 1] = L_I * G_tes_bath / (I * L)

            # TES1 thermal
            elif i == 1:
                X[i, 0] = -I * R * (2 + b) / C_tes
                X[i, 1] = (
                    1 / t_I
                    + G_abs_tes / C_tes
                    + 1j * omega
                )
                X[i, 2] = -G_abs_tes / C_tes

            # First absorber node
            elif i == 2:
                X[i, 1] = -G_abs_tes / C_abs
                X[i, 2] = (
                    G_abs_tes / C_abs
                    + G_abs_abs / C_abs
                    + 1j * omega
                )
                X[i, 3] = -G_abs_abs / C_abs

            # Last absorber node
            elif i == n_abs + 1:
                X[i, n_abs] = -G_abs_abs / C_abs
                X[i, n_abs + 1] = (
                    (G_abs_tes + G_abs_abs) / C_abs
                    + 1j * omega
                )
                X[i, n_abs + 2] = -G_abs_tes / C_abs

            # TES2 thermal
            elif i == n_abs + 2:
                X[i, n_abs + 1] = -G_abs_tes / C_tes
                X[i, n_abs + 2] = (
                    1 / t_I
                    + G_abs_tes / C_tes
                    + 1j * omega
                )
                X[i, n_abs + 3] = -I * R * (2 + b) / C_tes

            # TES2 electrical
            elif i == n_abs + 3:
                X[i, n_abs + 2] = (
                    L_I * G_tes_bath / (I * L)
                )
                X[i, n_abs + 3] = (
                    1 / t_el + 1j * omega
                )

            # Interior absorber node
            else:
                X[i, i - 1] = -G_abs_abs / C_abs

                X[i, i] = (
                    2 * G_abs_abs / C_abs
                    + 1j * omega
                )

                X[i, i + 1] = -G_abs_abs / C_abs

        return X

    # ============================================================
    # Transfer functions
    # ============================================================

    # After transpose:
    # shape = (state, physical_noise_source)
    N = matrix_N(n_abs).T

    omega = frequency * 2 * math.pi

    transfer_tes1 = []
    transfer_tes2 = []

    for omg in tqdm.tqdm(omega):

        M = matrix_M(n_abs, omg)

        # H[state, source]
        H = np.linalg.solve(M, N)

        # Keep COMPLEX transfer functions.
        transfer_tes1.append(H[0, :])
        transfer_tes2.append(H[n_abs + 3, :])

    # shape = (source, frequency)
    transfer_tes1 = np.asarray(transfer_tes1).T
    transfer_tes2 = np.asarray(transfer_tes2).T

    # ============================================================
    # Individual source ASD
    # ============================================================

    noise_tes1 = np.abs(transfer_tes1)
    noise_tes2 = np.abs(transfer_tes2)

    # ============================================================
    # IMPORTANT:
    #
    # Independent physical noise sources add in PSD:
    #
    # ASD_total = sqrt(sum(ASD_i^2))
    #
    # NOT sum(ASD_i).
    # ============================================================

    noise_total_tes1 = np.sqrt(
        np.sum(noise_tes1**2, axis=0)
    )

    noise_total_tes2 = np.sqrt(
        np.sum(noise_tes2**2, axis=0)
    )

    # ============================================================
    # Internal absorber TFN only
    #
    # source index:
    #   0              TES1 Johnson
    #   1              Load1 Johnson
    #   2              TES1-bath
    #   3              TES1-absorber
    #   4...n_abs+2    absorber-absorber
    # ============================================================

    absorber_slice = slice(4, n_abs + 3)

    noise_absorber_tes1 = np.sqrt(
        np.sum(
            noise_tes1[absorber_slice] ** 2,
            axis=0,
        )
    )

    noise_absorber_tes2 = np.sqrt(
        np.sum(
            noise_tes2[absorber_slice] ** 2,
            axis=0,
        )
    )

    # ============================================================
    # TES1-TES2 cross PSD
    # ============================================================

    cross_psd = np.sum(
        transfer_tes1
        * np.conjugate(transfer_tes2),
        axis=0,
    )

    psd_tes1 = noise_total_tes1**2
    psd_tes2 = noise_total_tes2**2

    # Sum channel
    noise_sum = np.sqrt(
        np.maximum(
            psd_tes1
            + psd_tes2
            + 2.0 * np.real(cross_psd),
            0.0,
        )
    )

    # Difference channel
    noise_diff = np.sqrt(
        np.maximum(
            psd_tes1
            + psd_tes2
            - 2.0 * np.real(cross_psd),
            0.0,
        )
    )

    # ============================================================
    # Save
    # ============================================================

    np.savetxt(
        f"./output/{output_number}/"
        "noise_spectral_total_TES1.dat",
        noise_total_tes1,
    )

    np.savetxt(
        f"./output/{output_number}/"
        "noise_spectral_total_TES2.dat",
        noise_total_tes2,
    )

    np.savetxt(
        f"./output/{output_number}/"
        "noise_spectral_absorber_TES1.dat",
        noise_absorber_tes1,
    )

    np.savetxt(
        f"./output/{output_number}/"
        "noise_spectral_absorber_TES2.dat",
        noise_absorber_tes2,
    )

    np.savetxt(
        f"./output/{output_number}/"
        "noise_cross_psd_TES1_TES2.dat",
        np.column_stack(
            (
                frequency,
                np.real(cross_psd),
                np.imag(cross_psd),
            )
        ),
        header=(
            "frequency_Hz "
            "cross_PSD_real_A2_per_Hz "
            "cross_PSD_imag_A2_per_Hz"
        ),
    )

    np.savetxt(
        f"./output/{output_number}/"
        "noise_spectral_sum_diff.dat",
        np.column_stack(
            (
                frequency,
                noise_sum,
                noise_diff,
            )
        ),
        header=(
            "frequency_Hz "
            "ASD_TES1_plus_TES2 "
            "ASD_TES1_minus_TES2"
        ),
    )

    # Complex transfer functions.
    np.savez_compressed(
        f"./output/{output_number}/"
        "noise_transfer_distributed.npz",
        frequency=frequency,
        transfer_tes1=transfer_tes1,
        transfer_tes2=transfer_tes2,
    )

    jsn = json.dumps(para, indent=4)

    with open(
        f"./output/{output_number}/input.json",
        "w",
    ) as file:
        file.write(jsn)

    # ============================================================
    # Plot TES1
    # ============================================================

    plt.figure(figsize=(8, 8))

    plt.plot(
        frequency,
        noise_tes1[0],
        linewidth=2,
        label="Johnson Noise (TES1)",
    )

    plt.plot(
        frequency,
        noise_tes1[n_abs + 6],
        linewidth=2,
        label="Johnson Noise (TES2)",
    )

    plt.plot(
        frequency,
        noise_tes1[1],
        linewidth=2,
        label="Johnson Noise (Load1)",
    )

    plt.plot(
        frequency,
        noise_tes1[n_abs + 5],
        linewidth=2,
        label="Johnson Noise (Load2)",
    )

    plt.plot(
        frequency,
        noise_tes1[2],
        linewidth=2,
        label="TFN (TES1-Bath)",
    )

    plt.plot(
        frequency,
        noise_tes1[n_abs + 4],
        linewidth=2,
        label="TFN (TES2-Bath)",
    )

    plt.plot(
        frequency,
        noise_tes1[3],
        linewidth=2,
        label="TFN (TES1-Absorber)",
    )

    plt.plot(
        frequency,
        noise_tes1[n_abs + 3],
        linewidth=2,
        label="TFN (TES2-Absorber)",
    )

    plt.plot(
        frequency,
        noise_absorber_tes1,
        linewidth=2,
        label="TFN RSS (Absorber-Absorber)",
    )

    plt.plot(
        frequency,
        noise_total_tes1,
        linewidth=3,
        label="Total Noise (TES1)",
    )

    plt.xlabel("Frequency [Hz]", fontsize=20)
    plt.ylabel("Noise Spectral Density [A/rtHz]", fontsize=20)

    plt.ylim(1e-13, 1e-9)
    plt.xlim(1e-1, 1e5)

    plt.loglog()
    plt.grid()

    plt.legend(
        loc="best",
        fancybox=True,
        fontsize=10,
        ncol=2,
    )

    plt.tight_layout()

    plt.savefig(
        f"./output/{output_number}/"
        "noise_spectrum_corrected.png",
        dpi=350,
    )

    plt.show()
    plt.cla()


if __name__ == "__main__":
    main()
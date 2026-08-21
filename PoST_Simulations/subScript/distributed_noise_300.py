"""Evaluate the 300-node electrothermal noise model near the mid-frequency knee.

This is a fast diagnostic version of ``noise_distributed.py``.  Rather than
solving every source column separately, it solves M.T x = e_CH0 and obtains
the complete CH0 source row as x.T @ N.  It is therefore practical to test a
300-node absorber on a logarithmic frequency grid.
"""

import json
from pathlib import Path

import numpy as np
from scipy import sparse


INPUT_PATH = Path(r"H:\hata2025\new\input.json")
OUTPUT_PATH = Path(r"D:\desktop\distributed_noise_300_midband.csv")


def main():
    para = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
    n_abs = 300
    c_abs = float(para["C_abs"]) / n_abs
    c_tes = float(para["C_tes"])
    g_link = float(para["G_abs-abs"]) * (n_abs - 1)
    g_abs_tes = float(para["G_abs-tes"])
    g_tes_bath = float(para["G_tes-bath"])
    resistance = float(para["R"])
    r_load = float(para["R_l"])
    t_c = float(para["T_c"])
    t_bath = float(para["T_bath"])
    alpha = float(para["alpha"])
    beta = float(para["beta"])
    inductance = float(para["L"])
    exponent = float(para["n"])
    k_b = 1.381e-23
    flink = 0.5

    current = np.sqrt(
        g_tes_bath * t_c * (1 - (t_bath / t_c) ** exponent) / (exponent * resistance)
    )
    tau_el = inductance / (r_load + resistance * (1 + beta))
    loop_gain = alpha * current**2 * resistance / (g_tes_bath * t_c)
    tau_i = c_tes / ((1 - loop_gain) * g_tes_bath)
    tfn_bath = np.sqrt(4 * k_b * t_c**2 * g_tes_bath * flink)
    tfn_abs_tes = np.sqrt(4 * k_b * t_c**2 * g_abs_tes * flink)
    tfn_link = np.sqrt(4 * k_b * t_c**2 * g_link * flink)
    johnson_tes = np.sqrt(4 * k_b * t_c * resistance * (1 + beta) ** 2)
    johnson_load = np.sqrt(4 * k_b * t_bath * r_load)

    states = n_abs + 4
    sources = n_abs + 7
    # State order: I1,T1,absorber[0..N-1],T2,I2.
    n_matrix = np.zeros((states, sources), dtype=complex)
    n_matrix[0, 0] = -johnson_tes / inductance
    n_matrix[1, 0] = current * johnson_tes / c_tes
    n_matrix[0, 1] = johnson_load / inductance
    n_matrix[1, 2] = tfn_bath / c_tes
    n_matrix[1, 3] = tfn_abs_tes / c_tes
    n_matrix[2, 3] = -tfn_abs_tes / c_abs
    # Internal absorber links: sources 4..N+2.
    for link in range(n_abs - 1):
        source = 4 + link
        left = 2 + link
        n_matrix[left, source] = tfn_link / c_abs
        n_matrix[left + 1, source] = -tfn_link / c_abs
    tes2_abs_source = n_abs + 3
    tes2_bath_source = n_abs + 4
    load2_source = n_abs + 5
    tes2_johnson_source = n_abs + 6
    n_matrix[2 + n_abs - 1, tes2_abs_source] = -tfn_abs_tes / c_abs
    n_matrix[2 + n_abs, tes2_abs_source] = tfn_abs_tes / c_tes
    n_matrix[2 + n_abs, tes2_bath_source] = tfn_bath / c_tes
    n_matrix[3 + n_abs, load2_source] = johnson_load / inductance
    n_matrix[3 + n_abs, tes2_johnson_source] = -johnson_tes / inductance
    n_matrix[2 + n_abs, tes2_johnson_source] = current * johnson_tes / c_tes

    # Exact comparison frequencies plus a log grid.  251 points keeps the
    # 300-state solve practical while resolving a broad 5--20 kHz feature.
    fixed = np.array([1e3, 3e3, 5e3, 7e3, 1e4, 1.5e4, 2e4, 3e4])
    frequency = np.unique(np.r_[np.geomspace(500, 50_000, 243), fixed])
    transfer = np.empty((len(frequency), sources), dtype=complex)
    select_ch0 = np.zeros(states)
    select_ch0[0] = 1.0
    for index, freq in enumerate(frequency):
        omega = 2 * np.pi * freq
        matrix = sparse.lil_matrix((states, states), dtype=complex)
        matrix[0, 0] = 1 / tau_el + 1j * omega
        matrix[0, 1] = loop_gain * g_tes_bath / (current * inductance)
        matrix[1, 0] = -current * resistance * (2 + beta) / c_tes
        matrix[1, 1] = 1 / tau_i + g_abs_tes / c_tes + 1j * omega
        matrix[1, 2] = -g_abs_tes / c_tes
        for node in range(n_abs):
            row = 2 + node
            if node == 0:
                matrix[row, 1] = -g_abs_tes / c_abs
                matrix[row, row] = (g_abs_tes + g_link) / c_abs + 1j * omega
                matrix[row, row + 1] = -g_link / c_abs
            elif node == n_abs - 1:
                matrix[row, row - 1] = -g_link / c_abs
                matrix[row, row] = (g_link + g_abs_tes) / c_abs + 1j * omega
                matrix[row, row + 1] = -g_abs_tes / c_abs
            else:
                matrix[row, row - 1] = -g_link / c_abs
                matrix[row, row] = 2 * g_link / c_abs + 1j * omega
                matrix[row, row + 1] = -g_link / c_abs
        tes2_temp = 2 + n_abs
        tes2_current = 3 + n_abs
        matrix[tes2_temp, tes2_temp - 1] = -g_abs_tes / c_tes
        matrix[tes2_temp, tes2_temp] = 1 / tau_i + g_abs_tes / c_tes + 1j * omega
        matrix[tes2_temp, tes2_current] = -current * resistance * (2 + beta) / c_tes
        matrix[tes2_current, tes2_temp] = loop_gain * g_tes_bath / (current * inductance)
        matrix[tes2_current, tes2_current] = 1 / tau_el + 1j * omega
        left_solution = sparse.linalg.spsolve(matrix.T.tocsc(), select_ch0)
        transfer[index] = left_solution @ n_matrix

    components = np.abs(transfer)
    total = np.sqrt(np.sum(components**2, axis=1))
    internal = np.sqrt(np.sum(components[:, 4 : n_abs + 3] ** 2, axis=1))
    bath = components[:, 2]
    load = components[:, 1]
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        OUTPUT_PATH,
        np.column_stack((frequency, total, bath, internal, load)),
        delimiter=",",
        header="frequency_hz,total_asd,bath_tfn_asd,internal_absorber_tfn_asd,load_johnson_asd",
        comments="",
    )
    print(f"n_abs={n_abs}, states={states}, sources={sources}")
    for target in fixed:
        i = np.abs(frequency - target).argmin()
        print(
            f"{frequency[i]:8.0f} Hz total={total[i]:.4e} "
            f"bath_fraction={bath[i]**2/total[i]**2:.3f} "
            f"internal_fraction={internal[i]**2/total[i]**2:.3f}"
        )
    print(f"saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()

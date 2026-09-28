"""Plot how the distributed absorber noise converges with node count.

All physical parameters are held fixed from ``input.json``.  The absorber
heat capacity is divided evenly across nodes, and the internal-link
conductance is scaled so the end-to-end absorber conductance stays fixed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib.tes_noise_model import F_LINK, noise_components  # noqa: E402

DEFAULT_INPUT = ROOT / "input.json"
DEFAULT_OUTPUT = ROOT / "figures" / "absorber_partition_convergence.png"
NODE_COUNTS = (1, 2, 4, 8, 16, 32, 64, 128, 300)


def total_asd(parameters: dict, node_count: int, frequency_hz: np.ndarray) -> np.ndarray:
    """Return CH0 total ASD for the symmetric two-TES distributed model."""
    c_abs = float(parameters["C_abs"]) / node_count
    c_tes = float(parameters["C_tes"])
    g_link = float(parameters["G_abs-abs"]) * (node_count - 1)
    g_abs_tes = float(parameters["G_abs-tes"])
    g_tes_bath = float(parameters["G_tes-bath"])
    resistance = float(parameters["R"])
    r_load = float(parameters["R_l"])
    t_c = float(parameters["T_c"])
    t_bath = float(parameters["T_bath"])
    alpha = float(parameters["alpha"])
    beta = float(parameters["beta"])
    inductance = float(parameters["L"])
    exponent = float(parameters["n"])
    k_b = 1.381e-23
    flink = F_LINK

    current = np.sqrt(
        g_tes_bath * t_c * (1.0 - (t_bath / t_c) ** exponent)
        / (exponent * resistance)
    )
    tau_el = inductance / (r_load + resistance * (1.0 + beta))
    loop_gain = alpha * current**2 * resistance / (g_tes_bath * t_c)
    tau_i = c_tes / ((1.0 - loop_gain) * g_tes_bath)

    tfn_bath = np.sqrt(4.0 * k_b * t_c**2 * g_tes_bath * flink)
    tfn_abs_tes = np.sqrt(4.0 * k_b * t_c**2 * g_abs_tes * flink)
    tfn_link = np.sqrt(4.0 * k_b * t_c**2 * g_link * flink)
    excess_m = float(parameters.get("excess_johnson_M", 0.0))
    johnson_tes = np.sqrt(
        4.0 * k_b * t_c * resistance * (1.0 + 2.0 * beta) * (1.0 + excess_m**2)
    )
    johnson_load = np.sqrt(4.0 * k_b * t_bath * r_load)

    states = node_count + 4
    sources = node_count + 7
    # States: I1, T1, absorber nodes, T2, I2.
    # Sources: TES/load Johnson, TES-bath TFN, TES-absorber TFN,
    # absorber-link TFNs, then the symmetric channel-2 sources.
    noise = np.zeros((states, sources), dtype=complex)
    noise[0, 0] = -johnson_tes / inductance
    noise[1, 0] = current * johnson_tes / c_tes
    noise[0, 1] = johnson_load / inductance
    noise[1, 2] = tfn_bath / c_tes
    noise[1, 3] = tfn_abs_tes / c_tes
    noise[2, 3] = -tfn_abs_tes / c_abs
    for link in range(node_count - 1):
        source = 4 + link
        left = 2 + link
        noise[left, source] = tfn_link / c_abs
        noise[left + 1, source] = -tfn_link / c_abs
    source_tes2_abs = node_count + 3
    source_tes2_bath = node_count + 4
    source_load2 = node_count + 5
    source_tes2_johnson = node_count + 6
    noise[node_count + 1, source_tes2_abs] = -tfn_abs_tes / c_abs
    noise[node_count + 2, source_tes2_abs] = tfn_abs_tes / c_tes
    noise[node_count + 2, source_tes2_bath] = tfn_bath / c_tes
    noise[node_count + 3, source_load2] = johnson_load / inductance
    noise[node_count + 3, source_tes2_johnson] = -johnson_tes / inductance
    noise[node_count + 2, source_tes2_johnson] = current * johnson_tes / c_tes

    select_ch0 = np.zeros(states)
    select_ch0[0] = 1.0
    result = np.empty(len(frequency_hz))
    for index, frequency in enumerate(frequency_hz):
        omega = 2.0 * np.pi * frequency
        matrix = sparse.lil_matrix((states, states), dtype=complex)
        matrix[0, 0] = 1.0 / tau_el + 1j * omega
        matrix[0, 1] = loop_gain * g_tes_bath / (current * inductance)
        matrix[1, 0] = -current * resistance * (2.0 + beta) / c_tes
        matrix[1, 1] = 1.0 / tau_i + g_abs_tes / c_tes + 1j * omega
        matrix[1, 2] = -g_abs_tes / c_tes
        for node in range(node_count):
            row = 2 + node
            if node == 0:
                matrix[row, 1] = -g_abs_tes / c_abs
                matrix[row, row] = (g_abs_tes + g_link) / c_abs + 1j * omega
                if node_count > 1:
                    matrix[row, row + 1] = -g_link / c_abs
            elif node == node_count - 1:
                matrix[row, row - 1] = -g_link / c_abs
                matrix[row, row] = (g_link + g_abs_tes) / c_abs + 1j * omega
                matrix[row, row + 1] = -g_abs_tes / c_abs
            else:
                matrix[row, row - 1] = -g_link / c_abs
                matrix[row, row] = 2.0 * g_link / c_abs + 1j * omega
                matrix[row, row + 1] = -g_link / c_abs
        tes2_temp = 2 + node_count
        tes2_current = 3 + node_count
        matrix[tes2_temp, tes2_temp - 1] = -g_abs_tes / c_tes
        matrix[tes2_temp, tes2_temp] = 1.0 / tau_i + g_abs_tes / c_tes + 1j * omega
        matrix[tes2_temp, tes2_current] = -current * resistance * (2.0 + beta) / c_tes
        matrix[tes2_current, tes2_temp] = loop_gain * g_tes_bath / (current * inductance)
        matrix[tes2_current, tes2_current] = 1.0 / tau_el + 1j * omega

        transfer = sparse.linalg.spsolve(matrix.T.tocsc(), select_ch0) @ noise
        result[index] = np.sqrt(np.sum(np.abs(transfer) ** 2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-nodes", type=int, default=300)
    args = parser.parse_args()

    parameters = json.loads(args.input.read_text(encoding="utf-8"))
    node_counts = tuple(count for count in NODE_COUNTS if count <= args.max_nodes)
    if not node_counts or node_counts[0] != 1:
        raise ValueError("--max-nodes must be at least 1")
    sample_rate = float(parameters["rate"])
    sample_count = float(parameters["samples"])
    frequency = np.geomspace(sample_rate / sample_count, sample_rate / 2.0, 241)
    spectra = {count: total_asd(parameters, count, frequency) for count in node_counts}
    # The production five-state model represents the absorber as one thermal
    # body with the effective TES-to-absorber conductance.
    reference = noise_components(parameters, frequency)["total_ch0"]

    fig, ax_spectrum = plt.subplots(figsize=(9.2, 5.8))
    palette = plt.get_cmap("viridis")
    colors = palette(np.linspace(0.05, 0.95, len(node_counts)))
    ax_spectrum.loglog(frequency, reference * 1e6, color="black", linestyle="--",
                       linewidth=2.0, label="Lumped absorber (5-state)")
    for count, color in zip(node_counts, colors):
        ax_spectrum.loglog(frequency, spectra[count] * 1e6, color=color,
                           linewidth=2.4 if count == 1 else 1.35,
                           label=f"{count} node" if count == 1 else f"{count} nodes")

    ax_spectrum.set_ylabel(r"Total current ASD ($\mu$A/$\sqrt{\rm Hz}$)")
    ax_spectrum.set_title("Absorber subdivision convergence")
    ax_spectrum.set_xlabel("Frequency (Hz)")
    ax_spectrum.grid(True, which="both", alpha=0.25)
    ax_spectrum.legend(ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)

    print(f"Input: {args.input}")
    print(f"Saved: {args.output}")
    for count in node_counts:
        delta = 100.0 * (spectra[count] / reference - 1.0)
        print(
            f"{count:>3} nodes: max |difference| = {np.max(np.abs(delta)):.6g}%, "
            f"RMS = {np.sqrt(np.mean(delta**2)):.6g}%"
        )


if __name__ == "__main__":
    main()

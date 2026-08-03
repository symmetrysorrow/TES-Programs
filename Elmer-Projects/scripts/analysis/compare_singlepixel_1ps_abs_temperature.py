"""Compare the new COMSOL 1-ps-pulse absorber temperature to Phase32."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/hybrid_prism_diagnostics/singlepixel_1ps_abs_temperature"
PULSE_S = 20.020e-3


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    c = np.loadtxt(ROOT / "reference/SinglePixel_1ps.txt", comments="%", usecols=(0, 1), encoding="utf-8")
    ct, ctemp = c[:, 0] * 1e-3, c[:, 1]
    cb = float(np.mean(ctemp[(ct >= 19.5e-3) & (ct < PULSE_S)]))
    e = np.loadtxt(ROOT / "artifacts/hybrid_prism_diagnostics/phase32_abs_temperature_splitpulse/absorber_temperature_comparison.csv", delimiter=",", skiprows=1)
    et, etemp = e[:, 0], e[:, 1]
    # Phase32 was extracted against the previous COMSOL file; replace only
    # the COMSOL columns here with the new 1-ps COMSOL data.
    cb_delta = ctemp - cb
    eb = float(np.mean(etemp[et < 0]))
    metrics = {
        "comsol_baseline_K": cb, "elmer_baseline_K": eb,
        "comsol_peak_K": float(ctemp.max()), "elmer_peak_K": float(etemp.max()),
        "comsol_peak_time_us": float((ct[np.argmax(ctemp)] - PULSE_S) * 1e6),
        "elmer_peak_time_us": float(et[np.argmax(etemp)]),
        "comsol_temperature_rise_K": float(ctemp.max() - cb),
        "elmer_temperature_rise_K": float(etemp.max() - eb),
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    with (OUT / "absorber_temperature_comparison.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["time_from_pulse_us", "comsol_abs_T_K", "comsol_delta_T_K", "elmer_abs_T_K", "elmer_delta_T_K"])
        for t, temp in zip(et, etemp):
            w.writerow([t, float(np.interp(t, (ct-PULSE_S)*1e6, ctemp)), float(np.interp(t, (ct-PULSE_S)*1e6, cb_delta)), temp, temp-eb])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True, constrained_layout=True)
    ax1.plot((ct-PULSE_S)*1e6, ctemp, "o-", color="#2166ac", label="COMSOL SinglePixel_1ps", markersize=2, linewidth=1.2)
    ax1.plot(et, etemp, "o-", color="#6a3d9a", label="Elmer Phase32 (10 ps pulse / 1 ns after)", markersize=3, linewidth=1.8)
    ax2.plot((ct-PULSE_S)*1e6, cb_delta, "o-", color="#2166ac", label="COMSOL", markersize=2, linewidth=1.2)
    ax2.plot(et, etemp-eb, "o-", color="#6a3d9a", label="Elmer Phase32", markersize=3, linewidth=1.8)
    for ax in (ax1, ax2):
        ax.set_xscale("symlog", linthresh=0.1, linscale=1); ax.set_xlim(-10, 10); ax.axvline(0, color="#555", linestyle="--", linewidth=.9); ax.grid(True, alpha=.25); ax.legend(frameon=False, loc="best")
    ax1.set_ylabel("Absorber temperature [K]"); ax1.set_title("Absorber temperature: COMSOL 1 ps pulse vs Elmer")
    ax2.set_ylabel("Temperature rise [K]"); ax2.set_xlabel("Time from pulse [µs] (symlog)")
    fig.savefig(OUT / "absorber_temperature_singlepixel_1ps_vs_elmer.png", dpi=240); fig.savefig(OUT / "absorber_temperature_singlepixel_1ps_vs_elmer.svg"); plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

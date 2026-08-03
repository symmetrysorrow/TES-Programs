"""Run the absorber-temperature extractor for Phase32."""
from __future__ import annotations

import numpy as np
import scripts.analysis.extract_abs_temperature_comparison as base

base.CASE = "case_p19_pulse_phase32_abs_temperature_splitpulse"
base.OUT = base.ROOT / "artifacts/hybrid_prism_diagnostics/phase32_abs_temperature_splitpulse"


def split_times():
    stages = [(18e-6, 1, 1), (1999.9995e-9, 1, 1), (0.1e-12, 10, 1), (10e-12, 100, 10), (0.1e-12, 10, 1), (1e-9, 100, 10)]
    times, t = [], base.PULSE_S - 20e-6
    for dt, count, interval in stages:
        for j in range(1, count + 1):
            t += dt
            if j % interval == 0:
                times.append(t)
    return np.asarray(times)


base.output_times = split_times

if __name__ == "__main__":
    base.main()

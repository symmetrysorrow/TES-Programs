"""Run the absorber-temperature extractor for the Phase31 fine-step case."""
from __future__ import annotations

import scripts.analysis.extract_abs_temperature_comparison as base

base.CASE = "case_p19_pulse_phase31_abs_temperature_fine1us"
base.OUT = base.ROOT / "artifacts/hybrid_prism_diagnostics/phase31_abs_temperature_fine1us"


def fine_times():
    stages = [(18e-6, 1, 1), (1999.9995e-9, 1, 1), (0.1e-12, 10, 1), (999e-12, 1, 1), (0.1e-12, 10, 1), (1e-9, 100, 1)]
    times, t = [], base.PULSE_S - 20e-6
    for dt, count, interval in stages:
        for j in range(1, count + 1):
            t += dt
            if j % interval == 0:
                times.append(t)
    return __import__("numpy").asarray(times)


base.output_times = fine_times

if __name__ == "__main__":
    base.main()

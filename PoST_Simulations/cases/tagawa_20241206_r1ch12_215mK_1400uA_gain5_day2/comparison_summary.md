# Target comparison record

Status: **blocked — no target simulation is reported**. This is the
**post-analysis** spectrum; the pre-analysis spectrum is not back-calculated.

The experimental values below are recomputed from the target `CH0_noise/rawdata` with the shared production estimator (345 accepted records) on its native 5 Hz grid (`500000 / 100000`). They remain in raw CH0 voltage ASD units because no target voltage-to-current calibration is available. The requested intrinsic-physical-noise comparison keeps TES Johnson, load Johnson, TES-bath TFN, and TES-absorber TFN enabled, while empirical white/readout floors, residual sources/poles, Lorentzian resistance fluctuation, and hanging TES are disabled. However, no independently sourced target operating point is available. Therefore simulation ASD, normalized simulation ASD, and simulation/experiment are `—`; filling them from the generic input or from a residual would invalidate the comparison.

| Frequency | Experimental post-analysis ASD (raw V/√Hz) | Experimental post-analysis / 1 kHz | Simulation post-analysis ASD | Simulation post-analysis / 1 kHz | Sim / exp |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 Hz | 7.3987686e-4 | 46.6210487 | — | — | — |
| 100 Hz | 6.4502958e-5 | 4.06445411 | — | — | — |
| 1 kHz | 1.5870018e-5 | 1.00000000 | — | — | — |
| 3 kHz | 1.3661775e-5 | 0.860854465 | — | — | — |
| 5 kHz | 1.0540628e-5 | 0.664185034 | — | — | — |
| 7 kHz | 7.2576767e-6 | 0.457320017 | — | — | — |
| 10 kHz | 3.4008741e-6 | 0.214295548 | — | — | — |

The values and provenance are machine-readable in [`comparison_summary.json`](comparison_summary.json) and [`provenance.json`](provenance.json). An absolute physical comparison remains disallowed until the target output calibration is independently established; the pre-existing `modelnoise.txt` is not used as a substitute for this recomputation.

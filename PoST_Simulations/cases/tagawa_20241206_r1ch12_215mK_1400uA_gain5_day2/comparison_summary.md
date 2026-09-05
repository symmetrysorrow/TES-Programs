# Target comparison record

Status: **blocked — no target simulation is reported**. The experimental
pre-analysis spectrum is now reconstructed directly from raw CH0 records;
it is not obtained by inverse filtering.

Both spectra use the same 345 accepted records and native 5 Hz grid
(`500000 / 100000`). Post-analysis is mean removal → 10 kHz digital Bessel
`filtfilt` → Hann → power average. Pre-analysis is mean removal → Hann → power
average. Values are raw CH0 voltage ASD units because no target voltage-to-
current calibration or readout transfer is independently established.

| Frequency | Post ASD (raw V/√Hz) | Post / 1 kHz | Pre ASD (raw V/√Hz) | Pre / 1 kHz | Simulation / ratio |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 Hz | 7.3987686e-4 | 46.6210487 | 7.39877399e-4 | 46.1561965 | — / — |
| 100 Hz | 6.4502958e-5 | 4.06445411 | 6.45093835e-5 | 4.02432591 | — / — |
| 1 kHz | 1.5870018e-5 | 1.00000000 | 1.6029860e-5 | 1.00000000 | — / — |
| 3 kHz | 1.3661775e-5 | 0.860854465 | 1.4998462e-5 | 0.935657667 | — / — |
| 5 kHz | 1.0540628e-5 | 0.664185034 | 1.3826595e-5 | 0.862552418 | — / — |
| 7 kHz | 7.2576767e-6 | 0.457320017 | 1.2546959e-5 | 0.782724169 | — / — |
| 10 kHz | 3.4008741e-6 | 0.214295548 | 1.0202737e-5 | 0.636483211 | — / — |

The intrinsic physical-source policy remains unchanged: TES Johnson, load
Johnson, TES-bath TFN, and TES-absorber TFN are enabled; empirical white or
readout floors, residual-derived sources/poles, Lorentzian resistance
fluctuation, and hanging TES are disabled. Simulation columns remain null
because the target physical operating point is not independently admissible.

The machine-readable record is [`comparison_summary.json`](comparison_summary.json).
Linkage and provenance are in [`detector_channel_linkage.json`](detector_channel_linkage.json)
and [`provenance.json`](provenance.json).

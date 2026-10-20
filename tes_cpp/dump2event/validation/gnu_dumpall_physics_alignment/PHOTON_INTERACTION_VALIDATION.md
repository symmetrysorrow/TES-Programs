# Photon interaction validation

PHITS 3.341 internal diagnostics for the existing runs report:

| Case | Histories | Internal score calls | Sum tlv (MeV) |
|---|---:|---:|---:|
| 2 MeV, thick Pb | 20 | 123,789 | 22.093269303402575 |
| 511 keV, thick Pb | 20 | 52,897 | 8.0676441116014352 |
| 511 keV, thin Pb plus vacuum | 20 | 12,440 | 1.8818912517211674 |

The transport runs completed with exit code 0 in the earlier diagnostic campaign, and baseline/diagnostic output hashes matched. Those facts show instrumentation did not perturb those runs, but do not pair internal score calls to raw rows. This turn’s GNU Dumpall parse stops at the first NCOL=4 particle header in each case. Photoelectric, Compton, Rayleigh, pair production, fluorescence and escape have 0 newly raw-mapped deposit additions here. The counts/sums above are PHITS internal values only, not raw-derived validation.

# Target-case continuation: r1ch12, 215 mK, 1400 uA, gain 5

This note records the continuation for commit `26806c78cba15f8a8467c2e465a30be286255887` and the target acquisition
`G:\\tagawa\\20241206\\r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2`.

The exact acquisition metadata is independently available from the target
`PulseConfig.json` and `Setting.txt`: `500 kS/s`, `100,000` samples, `10 kHz`
cutoff, and therefore `df=5 Hz`. The directory name provides the `215mK`,
`1400uA1400uA`, and `gain5` labels, but these are not promoted to `T_c`,
calibrated TES current, or output calibration.

The repository and target-folder audit found no exact target source for `T_c`,
TES `R`, `alpha/beta`, `L`, `n`, the heat capacities/conductances, `n_abs`, or
the CH0 voltage-to-current calibration. A same-day campaign file,
`G:\\tagawa\\20241206\\room1-ch1-iv3\\calibration\\IV_215mK.txt`, does recover
the bath-setpoint meaning and provides a conditional TES-current candidate:
the existing superconducting-slope derivation gives `eta=100.725337 uA/V`
and `I_TES=251.847592 uA` at the 1400 uA bias row. The channel linkage is not
encoded, so this remains provenance-only and is not copied into `input.json`.
Generic `PoST_Simulations/input.json`, `R_SH` constants in general IV scripts,
and unrelated Elmer cases are listed as rejected candidates in the case
`provenance.json`; none is used as target input.

Consequently the target no-added-noise comparison is saved as a transparent
blocked table in
`cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.md`
and `.json`: experimental ASD is listed at 10, 100, 1k, 3k, 5k, 7k, and 10k
Hz, while simulation and ratio fields are null. This is the only valid result
until an independent target operating-point/calibration record is supplied.

The C++ `posi2pulse --dump-linearization <input.json>` debug path and
`tests/test_cpp_python_linearization.py` regression test verify the shared
electrical/TES intrinsic sub-block against Python. The absorber boundary term
is reported separately because the C++ distributed model uses `G_abs-tes`,
whereas the Python five-state reduction uses `G_eff`; that documented model
reduction is not silently treated as a calibration mismatch.

The generic fixture parity result is passing for operating current, `tau_el`,
loop gain, `tau_i`, `dI/dI`, `dI/dT`, `dT/dI`, intrinsic TES thermal diagonal,
boundary term, sign convention, and left/right symmetry. The target parity test
is conditionally enabled by `cpp_parity_ready` and therefore does not run while
the target input contains unresolved values. No target operating-point or
stability eigenvalue is reported.

The target report is explicitly post-analysis: `experimental_post_analysis_asd`
and `experimental_post_analysis_normalized` are based on the 10 kHz digital
Bessel `filtfilt` path. `experimental_pre_analysis_asd` is null because inverse
filtering in the stopband is not safe. The target normalized no-added-noise
comparison is consequently not available yet; the physical source semantics
are preserved in `input.json` (TES Johnson, load Johnson, TES-bath TFN, and
TES-absorber TFN enabled; only empirical/residual/hanging terms disabled).

## Readiness and final judgment

| Capability | Ready | Blocking information |
| --- | --- | --- |
| Operating point | No | `T_c`, target TES `R`, `G_tes-bath`, `n` |
| Python stability | No | operating point plus thermal/electrical parameters |
| Reduced intrinsic noise | No | stability set plus FFT settings (FFT settings are known) |
| C++ TES parity | No for target | same set plus `n_abs`; `E` is not required |
| Normalized post-analysis comparison | No | physical reduced-noise input; calibration is not required |
| Absolute ASD comparison | No | normalized comparison plus independent CH0/readout calibration |

Final classification: **C — existing information is insufficient to construct
the target physical case**. The work has nevertheless reached the partial
provenance stage described by B: acquisition settings, bath-setpoint meaning,
and a conditional same-campaign TES-current candidate were recovered. The
minimum next information is prioritized in `parameter_dependency.json`:
first target-linked `T_c`, `R`, `G_tes-bath`, and thermal-law `n`; next the
remaining electrothermal/circuit values; calibration is only needed for the
absolute comparison. No residual fit, 1 kHz normalization, empirical noise
term, or arbitrary filter is used.

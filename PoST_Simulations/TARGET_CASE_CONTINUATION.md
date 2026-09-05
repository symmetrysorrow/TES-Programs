# Target-case continuation: r1ch12, 215 mK, 1400 uA, gain 5

This note records the continuation for commit `26806c78cba15f8a8467c2e465a30be286255887` and the target acquisition
`G:\\tagawa\\20241206\\r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2`.

The exact acquisition metadata is independently available from the target
`PulseConfig.json` and `Setting.txt`: `500 kS/s`, `100,000` samples, `10 kHz`
cutoff, and therefore `df=5 Hz`. The directory name provides the `215mK`,
`1400uA1400uA`, and `gain5` labels, but these are not promoted to `T_c`,
calibrated TES current, or output calibration.

The repository and target-folder audit found no exact target source for `T_c`,
`T_bath` as a measured TES operating temperature, TES `I/R`, `alpha/beta`,
`L`, `n`, the heat capacities/conductances, `n_abs`, or the CH0
voltage-to-current calibration. Generic `PoST_Simulations/input.json`,
`R_SH` constants in general IV scripts, and unrelated Elmer cases are listed
as rejected candidates in the case `provenance.json`; none is used as target
input.

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

## Final judgment

Primary judgment: **6 — current information is insufficient to identify the
cause**. Ranked next are **1 — noise-source provenance/definition mismatch**
(including the now-corrected estimator path and possible unit/ASD provenance),
**2 — operating-point inconsistency**, and **5 — readout/bias transfer
omission**. No residual fit, 1 kHz normalization, or empirical noise term is
used to choose among them.

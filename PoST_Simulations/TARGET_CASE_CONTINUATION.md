# Target-case continuation: r1ch12, 215 mK, 1400 uA, gain 5

This note records the continuation for commit `95fab67b51f71e97bb47e01e998c1019ff3a1dc9` and the target acquisition
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

The target report contains post-analysis `experimental_post_analysis_asd` and
`experimental_post_analysis_normalized` based on the 10 kHz digital Bessel
`filtfilt` path. `experimental_pre_analysis_asd` is directly reconstructed
from the same 345 accepted raw records using mean removal, Hann, power average,
and the shared one-sided ASD helper. No inverse filtering is used. The target
normalized no-added-noise comparison is consequently not available yet; the physical source semantics
are preserved in `input.json` (TES Johnson, load Johnson, TES-bath TFN, and
TES-absorber TFN enabled; only empirical/residual/hanging terms disabled).

## Readiness and final judgment

| Capability | Ready | Conditionally ready | Blocking information |
| --- | --- | --- | --- |
| Operating point | No | No | `T_c`, target TES `R`, `G_tes-bath`, `n`; `T_bath=0.215 K` is setpoint-only |
| Python stability | No | No | operating point plus thermal/electrical parameters |
| Reduced intrinsic noise | No | No | stability set plus FFT settings (FFT settings are known) |
| C++ TES parity | No for target | No | same set plus `n_abs`; `E` is not required |
| Normalized comparison | No | No | physical reduced-noise input; calibration is not required |
| Absolute comparison | No | No | independent CH0 calibration and readout gain/transfer in addition |

Final classification: **C — existing information is insufficient to construct
the target physical case**. The work has nevertheless reached the partial
provenance stage described by B: acquisition settings, bath-setpoint meaning,
and a conditional same-campaign TES-current candidate were recovered. The
minimum next information is prioritized in `parameter_dependency.json`:
first target-linked `T_c`, `R`, `G_tes-bath`, and thermal-law `n`; next the
remaining electrothermal/circuit values; calibration is only needed for the
absolute comparison. No residual fit, 1 kHz normalization, empirical noise
term, or arbitrary filter is used.

## Audit report

### Newly confirmed facts

- Target `setting.xml` maps the acquisition to `PXI2Slot2/ai0:1`; target CH0/CH1 are paired records.
- Target acquisition is `500000 Hz`, `100000` samples, `10 kHz` cutoff, `df=5 Hz`; 345 CH0 records pass the unchanged production mask.
- Direct pre-analysis is available at all requested frequencies and uses the same accepted mask.
- The C++ dump now reads the actual `make_matrix()` TES sub-blocks.

### Conditional candidates and rejected candidates

The same-day `room1-ch1-iv3` file gives the conditional candidate
`eta=100.725337 uA/V`, `I_TES=251.847592 uA` at the 1400 uA row. Generic
`R_SH=3.8/3.9 mOhm`, the generic input, and all residual-derived parameters
remain rejected.

### Detector/channel linkage result

`likely_but_unproven`: same date and `room1/ch1` naming are suggestive, but
the IV folder has no DAQ channel string, detector ID/serial, or explicit map
to target `CH0`; the candidate is not admissible. See
`cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/detector_channel_linkage.json`.

### Multi-temperature IV availability and parameter provenance

Only `IV_215mK.txt` was found under the same-day campaign. There is no second
temperature, RT sweep, or thermal-law record, so `T_c`, `G_tes-bath`, and `n`
cannot be independently constrained. `R_SH` is also unresolved from campaign
circuit/specification records. `T_bath=0.215 K` is `setpoint_only`, not an
independent thermometer measurement.

### Direct pre-analysis spectrum

The post/pre normalized values at 10, 100, 1000, 3000, 5000, 7000, and
10000 Hz are recorded in `comparison_summary.json` and `comparison_summary.md`.
At 10 kHz the direct pre-analysis normalized value is `0.636483211`, while
the post-analysis value is `0.214295548`; this is the known digital Bessel
analysis difference, not an empirical correction.

### C++ actual matrix parity and hanging fix

The generic fixture tests pass against actual C++ matrix entries. The hanging
TES1 branch was missing `1/t_i + G_abs-tes/C_tes`; it now has the same
`1/t_i + G_abs-tes/C_tes + G_tes-hanging/C_tes` diagonal structure as TES2,
and regression coverage checks both TES blocks and hanging rows. Target
hanging remains disabled.

### Operating point, stability, and intrinsic noise comparison

Not run: the provenance-aware operating-point gate is closed, so stability
and intrinsic physical-noise simulation are not admissible. No noise fitting,
white floor, pole-zero, Lorentzian, hanging parameter, or arbitrary transfer
was introduced.

### Remaining blockers and next best action

The highest-value existing-data action is to locate a target-linked detector
mapping or calibration record that proves the IV linkage and supplies
`T_c`, `R`, `G_tes-bath`, `n`, and the campaign `R_SH`. Until then the correct
final judgment is **C: existing data do not determine the target physics**.

# 5--20 kHz noise-shape investigation handoff

## Objective

Explain and, if possible, correct the mismatch between experimental and
simulated normalized CH0 noise ASD around 5--20 kHz.  Both spectra are
normalized at 1 kHz.  The experimental spectrum is above the simulation near
10 kHz; changing the model to raise 10 kHz has tended to over-predict 30 kHz.

## Authoritative paths

| Purpose | Path |
| --- | --- |
| Comparison script | `D:\Github\Analysis\CompareNoise-double.py` |
| Experimental data root | `G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2` |
| Experimental configuration | `G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2\PulseConfig.json` |
| Experimental CH0 model noise | `G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2\CH0_noise\modelnoise.txt` |
| Experimental raw CH0 noise | `G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2\CH0_noise\rawdata` |
| Active simulation output/input | `H:\hata2025\new` |
| Simulation program | `D:\github\TES-Programs\PoST_Simulations\PoST_Simulation.py` |
| Simulation input | `H:\hata2025\new\input.json` |
| Simulation component output | `H:\hata2025\new\noise.h5` |
| Experimental analysis implementation | `D:\github\TES-Programs\Analyze_Experimental_Data\tes_analysis\operations.py` |

## Current comparison and analysis conditions

- Acquisition: 500 kS/s, 100,000 samples, Nyquist 250 kHz.
- Hardware stage: 100 kHz, fourth-order analog Bessel in the comparison.
- Software analysis: second-order 10 kHz Bessel applied with `filtfilt`, then
  Hann window and FFT.
- 345 CH0 noise records pass the current selection used for direct rebuilds.
- Run the comparison:

```powershell
python D:\Github\Analysis\CompareNoise-double.py --no-show
```

The high-frequency (`>225 kHz`) mismatch was solved before this investigation:
the old experimental `modelnoise.txt` had a high-frequency analysis artifact.
The current regenerated spectrum agrees well there.  Do not revisit that issue
unless `modelnoise.txt` is regenerated with different settings.

## Mid-band observations

After removing the common 10 kHz analysis-filter response, the normalized
input ASD is:

| Frequency | Experiment | Simulation | Experiment / simulation |
| ---: | ---: | ---: | ---: |
| 3 kHz | 0.947 | 0.797 | 1.19 |
| 5 kHz | 0.861 | 0.605 | 1.42 |
| 7 kHz | 0.787 | 0.473 | 1.66 |
| 10 kHz | 0.634 | 0.349 | 1.82 |
| 15 kHz | 0.420 | 0.239 | 1.76 |
| 20 kHz | 0.272 | 0.179 | 1.52 |
| 30 kHz | 0.138 | 0.117 | 1.18 |

Therefore the mismatch is not an analysis Bessel mismatch: both experimental
and simulation paths have the same `filtfilt` magnitude of about 0.333 at
10 kHz.  It is a broad 5--20 kHz shape difference before that filter.

## CH0--CH1 correlation result

Cross spectra from 400 simultaneous CH0/CH1 records show no significant
common-mode term in 1--30 kHz:

- Coherence: about 0.0016--0.0019.
- Correlation coefficient: near zero.
- For 400 averages, this is at/below the random estimator floor.
- At 10.5 kHz, normalized sum/difference ASD: 0.553 / 0.554.

Result file: `D:\desktop\CH0_CH1_noise_common_mode.csv`.

This rules out a *dominant common-mode* readout/bias noise explanation.  It
does not rule out independent electronics noise in each channel.

## What MakeNoise actually calculates

`PoST_Simulation.py:MakeNoise()` is not the spatial PoST thermal simulation.
It implements the explicitly documented reduced five-state linearized
electrothermal model:

```text
I_TES1, T_TES1, T_abs(center), T_TES2, I_TES2
```

For each angular frequency it constructs `M(omega)` and physical noise source
matrix `N`, then computes `H = solve(M, N)`.  The output HDF5 attribute is
`noise_model = five_state_effective_conductance`.

The individual source ASDs in `noise.h5/components_ch0` are consequently
available for diagnosis.  From 1--30 kHz their simulated CH0 power fractions
are approximately:

| Source | Power fraction |
| --- | ---: |
| `phonon_tes1_bath` | 96.6--96.7% |
| `phonon_tes1_absorber_effective` | 1.9% |
| `johnson_load1` | 1.2% |
| `johnson_tes1` | 0.2--0.3% |

The dominant TES1--bath TFN and the non-negligible sources have almost the
same frequency shape because they propagate through the same linearized TES
current response.  Scaling existing source amplitudes cannot form the needed
localized shoulder.

## 300-node distributed thermal result

`subScript/noise_distributed.py` contains a distributed absorber model but is
an old standalone program with hard-coded paths.  A fast compatible 300-node
diagnostic was added at:

```text
D:\github\TES-Programs\PoST_Simulations\subScript\distributed_noise_300.py
```

It uses the same equations as `noise_distributed.py`, evaluates 300 absorber
nodes (304 state variables, 307 physical noise sources), and writes:

```text
D:\desktop\distributed_noise_300_midband.csv
```

Run it with:

```powershell
python D:\github\TES-Programs\PoST_Simulations\subScript\distributed_noise_300.py
```

Result: 300-node / five-state total-ASD ratio is only `1.000001--1.000004`
from 1--50 kHz.  Internal absorber-link TFN has power fraction below
`5.4e-6`.  The absorber internal thermal-mode scale is about 45 MHz with the
current values, so it cannot cause a 10 kHz feature.  The five-state absorber
reduction is therefore not the cause in this parameter set.

## Parameter sweeps already performed

All sweeps regenerated only isolated desktop outputs; `H:\hata2025\new` was
not overwritten.

Reusable sweep script:

```text
D:\github\TES-Programs\PoST_Simulations\subScript\sweep_noise_lr.py
```

The script rebuilds the experimental CH0 spectrum from raw data, runs
`MakeNoise()` in an isolated output directory, applies the hardware and
analysis responses, and scores 1--30 kHz log-ASD residuals.

| Sweep | Result | CSV |
| --- | --- | --- |
| `L` / `R` | Improves total residual but cannot match 10 and 30 kHz simultaneously. Low-temperature best: L=1.0 uH, R=38.08 mOhm, 10k ratio=0.813, 30k ratio=1.345. | `D:\desktop\post_noise_lr_sweep.csv` |
| `alpha` / `beta` (±50%) | No material improvement. | `D:\desktop\post_noise_etf_sweep.csv` |
| `C_tes`, `G_abs-tes`, `G_tes-bath` (0.1--10x individually) | No material improvement. | `D:\desktop\post_noise_thermal_sweep.csv` |
| `T_bath` 15--180 mK | No useful improvement; superseded by high-temperature scan. | `D:\desktop\post_noise_bath_sweep.csv` |
| `T_bath` 195--210 mK at fixed L/R | Raises 10 kHz but over-raises 30 kHz. | `D:\desktop\post_noise_bath_high_sweep.csv` |
| Joint L/R at 195 and 200 mK | Best at 200 mK, L=0.75 uH, R=22.85 mOhm; score=0.07562, but still 10k ratio=0.812 and 30k ratio=1.339. | `D:\desktop\post_noise_bath_lr_sweep.csv` |

`T_bath=195--210 mK` is the physically relevant range (the earlier low-
temperature scan was not physically appropriate).  The joint scan confirms
that using the proper bath range does not resolve the shape mismatch.

## Current interpretation

Within the currently implemented model and tested parameter ranges, existing
`input.json` values cannot simultaneously make the normalized 10 kHz and
30 kHz regions agree.  This does **not** prove a coding error.  It means the
model lacks a degree of freedom with the observed broad, non-monotonic
mid-band shape, or a relevant parameter/physical configuration is not being
represented by the current equations.

Most productive next investigations:

1. Audit the electrical/readout transfer function represented by `matrix_M`:
   especially whether the real bias/readout chain has an additional pole,
   zero, or impedance absent from the model.
2. Derive the residual ASD/PSD after the best physical model and determine
   whether its shape is band-pass-like, resonance-like, or consistent with an
   independent electronics term.
3. If allowed by the experiment, measure the readout transfer function with a
   swept injected signal or an off-TES/reference input.  The CH0--CH1 result
   says a dominant *common-mode* term is absent, but an independent per-channel
   electronics term remains possible.

## Electrical/readout audit and residual inversion (2026-08-20)

A reproducible diagnostic was added at:

```text
D:\github\TES-Programs\PoST_Simulations\subScript\midband_residual.py
```

Run it without touching the active simulation output:

```powershell
python D:\github\TES-Programs\PoST_Simulations\subScript\midband_residual.py `
  --work-dir D:\desktop\post_midband_residual_work `
  --csv D:\desktop\post_midband_residual.csv `
  --plot D:\desktop\post_midband_residual.png
```

It rebuilds the 345-record experimental spectrum, removes the common digital
10 kHz response, regenerates the best physically plausible sweep point
(`T_bath=200 mK`, `L=0.75 uH`, `R=22.8493 mOhm`) in an isolated directory,
and evaluates both additive-PSD and multiplicative-transfer explanations.

### What `matrix_M` does and does not contain

Writing `matrix_M = A + i omega I`, the implemented coefficients and signs
match the usual linearized stiff-voltage-bias TES equations.  CH0 is then
selected directly as state 0.  Consequently, the model contains:

- a frequency-independent series load resistance `R_l`;
- one frequency-independent series inductance `L`;
- the TES electrothermal impedance generated by the five states;
- a fixed current readout (`C = [1, 0, 0, 0, 0]`, no output dynamics).

The 100 kHz fourth-order Bessel and the post-filter white ASD are applied only
after the detector calculation.  There is no representation of a
frequency-dependent Thevenin/Norton bias impedance, shunt capacitance, wiring
resonance, SQUID input impedance, flux-locked-loop transfer, amplifier
voltage/current noise shaped by source impedance, or channel-specific readout
pole/zero.  Thus an independently measured 5--30 kHz injected-signal transfer
function can be compared directly against the correction below.

One amplitude audit item is not capable of explaining the shoulder but should
be resolved: the TES Johnson source uses
`4 k_B T R (1 + beta)^2`.  If the intended convention is the usual
non-equilibrium Johnson correction, its voltage-noise PSD factor is
`1 + 2 beta`, without `beta^2`.  With the current `beta=3.994`, the code's
source PSD is 2.77 times that convention, but TES Johnson noise is only a few
tenths of one percent of total modeled power in this band.

### Stability finding (must be resolved before another unconstrained sweep)

The frequency-domain magnitude calculation does not check whether `A` is a
stable time-domain matrix.  The best 200 mK point has:

- loop gain `L_I = 18.356`;
- electrical rate `1/tau_el = 1.5734e5 s^-1`;
- open-loop TES thermal diagonal rate `-5.2896e6 s^-1`;
- four right-half-plane poles: two near `+3.7380e4 s^-1`
  (`5.949 kHz`) and two near `+5.0949e6 s^-1` (`810.9 kHz`);
- one stable slow absorber pole near `-5.356 s^-1` (`0.852 Hz`).

The same issue exists in the active base input: its right-half-plane pole
scales are about 3.60 kHz and 3.91 MHz.  The C++ distributed pulse matrix uses
the same coefficients and signs, so this is not a disagreement between the
noise and pulse implementations.

For the best point, the local TES electrical/thermal trace changes sign only
around `L = 22.3 nH`, versus the fitted `750 nH`.  Equivalently, holding `L`
fixed requires approximately `C_tes = 2.50e-12 J/K`, 33.6 times the current
value, to reach that boundary.  This does not by itself identify which input
parameter is wrong; it shows that the points compared so far are not valid
stable operating points.  A magnitude response can still look smooth because
the imaginary-axis magnitude of a real left- or right-half-plane pole is the
same.  Future sweeps should reject any trial with `max(real(eig(-A))) >= 0`.

### Multiplicative missing-transfer result

The experimental/five-state ASD ratio from 1--30 kHz is reproduced to an RMS
factor of 1.035 by the minimum useful pole/zero form

```text
H(s) = (1 + s/w_z)^2 / (1 + s/(Q w_p) + (s/w_p)^2)
```

normalized at 1 kHz, with:

- pole-pair natural frequency `f_p = 12.944 kHz`;
- pole quality factor `Q = 0.810`;
- two real zeros at `f_z = 17.626 kHz`;
- asymptotic high-frequency amplitude gain `(f_p/f_z)^2 = 0.539`.

This is a quantitative target response, not proof of a particular circuit.
At 10 kHz it requires a magnitude around 1.27, while at 30 kHz it requires
about 0.74, relative to 1 kHz.  A swept injected-signal measurement should
accept or reject this interpretation cleanly.

### Additive residual-PSD result

This result is a **diagnostic hypothesis, not an approved production noise
model**.  Its parameters were inferred from the same CH0 spectrum used for the
comparison, so agreement after adding it is circular and cannot validate the
physical origin.  Do not enable it in `PoST_Simulation.py` until its amplitude
and transfer shape have been fixed by an independent electronics/reference
measurement and then checked against held-out TES noise data.

Because both spectra were normalized independently at 1 kHz, simply computing
`PSD_exp - PSD_sim` at equal 1 kHz amplitude is not a physical extra source:
the residual is negative in 4526 of 7801 bins from 1--40 kHz.  The relative
absolute scale must be supplied before PSD subtraction.

Using the median experimental/simulation ASD ratio over 28--32 kHz as the
anchor gives a five-state amplitude scale of 0.7206.  The exact non-negative
residual then requires an independent source ASD (in units of experimental
ASD at 1 kHz) of approximately:

| Frequency | Required extra ASD | Extra fraction of total PSD |
| ---: | ---: | ---: |
| 1 kHz | 0.693 | 48% |
| 3 kHz | 0.687 | 53% |
| 5 kHz | 0.655 | 58% |
| 7 kHz | 0.631 | 64% |
| 10 kHz | 0.514 | 66% |
| 15 kHz | 0.325 | 60% |
| 20 kHz | 0.182 | 45% |
| 30 kHz | 0.036 | 7% |

Its shape is therefore **low-pass-like, not band-pass-like**.  A generalized
low-pass ASD

```text
A / sqrt(1 + (f/f_c)^(2 n))
```

with the 1 kHz PSD closure enforced has `A=0.693`, `f_c=11.45 kHz`, and
`n=2.41`; the combined model has an RMS factor of 1.053 over log-spaced
1--30 kHz samples.  In physical terms, the shoulder can be formed by the
cross-over between the existing detector noise and a substantial independent
per-channel source that rolls off around 10--15 kHz.  Absolute current/voltage
calibration is needed to distinguish that mixture from the multiplicative
pole/zero explanation.

### Recommended order now

1. Add a stability rejection to every physical-parameter sweep and establish
   a stable operating point from independently measured `L`, `C_tes`,
   `alpha`, `beta`, and the actual bias impedance.
2. Measure the injected-signal or complex-impedance response from 1--50 kHz.
   Compare it with the `12.94 kHz / Q=0.81 / 17.63 kHz` target above.
3. Obtain an absolute CH0 input-current calibration (including SQUID/FLL
   gain).  Then repeat PSD subtraction without an arbitrary 30 kHz anchor.
4. If no multiplicative feature is measured, test independent electronics
   voltage/current noise with an approximately 11.5 kHz, order-2.4 roll-off.

The diagnostic plot was simplified accordingly.  It now compares only the
final spectra after identical digital analysis in the upper panel, and shows
`simulation / experiment` in the lower panel with a +/-10% band.  The current
five-state model and the residual-derived *provisional* candidate are labeled
separately; the 1--30 kHz region is explicitly identified as a fit/comparison
region rather than independent validation.

## Stable internal-source-only existence test (2026-08-20)

Because independent SQUID/off-TES noise is expected to be small, a subsequent
search explicitly prohibited any added readout ASD and required all eigenvalues
of the time-domain five-state matrix to have negative real parts.

Within the previously tested, near-current parameter envelope, the best stable
fit still had a 1--30 kHz RMS factor of 1.129.  It also sat on many bounds
(`C_tes=10x`, `G_tes-bath=0.5x`, `alpha=0.5x`, `beta=0.5x`, `T_bath=210 mK`),
so the failure is not just a small local adjustment.

A deliberately broad simultaneous search did demonstrate that the five-state
model is mathematically capable of the measured shape without any added noise.
The stable candidate had an RMS factor of 1.019 and used:

| Parameter | Candidate | Relative to active input |
| --- | ---: | ---: |
| `T_bath` | 210 mK | physically allowed upper edge |
| `L` | 0.145 uH | 0.113x |
| `R` | 10.0 mOhm | 0.328x |
| `C_tes` | `1.393e-12 J/K` | 18.7x |
| `G_abs-tes` | `4.428e-8 W/K` | 100x |
| `G_tes-bath` | `2.272e-9 W/K` | 0.1x |
| `alpha` | 861.9 | 3x |
| `beta` | 0.399 | 0.1x |

Several values reached the broad search bounds.  This is therefore an
existence proof, not a credible parameter determination from noise alone.
The mechanism is clear: the current model is about 96.7% TES--bath TFN, while
the broad candidate becomes about 94% TES1--absorber TFN at 1 kHz, 93% at
10 kHz, and 80% at 30 kHz.  The shoulder is generated by changing the internal
TFN balance, chiefly through the much larger
`G_abs-tes / G_tes-bath` ratio, rather than by adding an electronics source.

The same candidate was evaluated with 300 absorber nodes.  Its distributed /
five-state ASD ratio was 0.994--1.012 from 1--30 kHz (RMS deviation 0.0064),
and distributed internal-link TFN contributed only 2.6% at 1 kHz, falling to
0.6% at 30 kHz.  Thus absorber subdivision still does not create the shoulder;
the fitted interface and bath conductances do.

The clear internal-only comparison is saved as:

```text
D:\desktop\post_noise_stable_internal_validation.png
```

The lower panel shows the stable internal-only candidate within +/-10% for
100% of the 250 Hz bins from 1--30 kHz.  This is still the fit band.  The next
validity test is not another noise-only fit, but fixing `L`, `R`, `C_tes`,
`G_abs-tes`, `G_tes-bath`, `alpha`, and `beta` from pulse, DC, or complex-
impedance information and checking whether the internal TFN-dominant regime
survives.

## Worktree note

The repository was already dirty before this investigation.  The diagnostic
scripts added during this work are:

- `PoST_Simulations/subScript/sweep_noise_lr.py`
- `PoST_Simulations/subScript/distributed_noise_300.py`

Do not discard unrelated existing worktree changes.

## Estimator/Jacobian audit and correction (2026-09-05)

### Baseline

Before this change, the repository test suite passed `42` tests.  The
production `tes_analysis.operations.NoiseAnalysis()` already used a Hann
window and power averaging, but `subScript/midband_residual.py` rebuilt the
experimental spectrum by averaging `abs(FFT)` record-by-record.  That was a
real estimator-parity error.  `finite_record_simulation_spectrum()` had the
right power-average intent, but duplicated the normalization and preprocessing
instead of sharing the production implementation.

### Confirmed and corrected errors

The following changes were made without fitting an experimental residual or
adding white/readout noise:

* `Analyze_Experimental_Data/tes_analysis/noise_utils.py` now owns record
  preprocessing, windowed rFFT power, one-sided ASD normalization, and the
  finite-record estimator.
* `NoiseAnalysis()`, `noise_main.py`, `midband_residual.experimental_asd()`,
  and `PoST_Simulation.finite_record_simulation_spectrum()` use the shared
  implementation.  The order is per-record mean removal, digital Bessel
  `filtfilt`, Hann window, power average, square root, then DC/Nyquist/interior
  one-sided factors.  `df = rate/sample` and the Hann power gain are applied
  once.
* The analytic reduced TES Jacobian is now explicit as a time-domain matrix
  `A`, with `M(omega) = -A + i*omega*I`.  A canonical nonlinear five/seven
  state RHS and central finite-difference Jacobian are available for audit.
* A stability gate reports all eigenvalues, maximum real part, unstable-mode
  status, and pole-frequency scales.  `stability_mode` defaults to `warn`;
  `strict` rejects an unstable operating point.  Invalid `T_bath >= T_c`
  operating points are rejected instead of producing NaNs.
* The hanging model's right TES thermal `dT2/dI2` element was found missing
  in the new explicit matrix and restored.  The existing left/right symmetry
  test catches regressions of this kind.

### Jacobian and stability result

For the repository nominal `PoST_Simulations/input.json` (`I=57.1173 uA`,
`T=T_c=142 mK`), the analytic and finite-difference Jacobians agree to a
maximum relative error of `1.78e-9` (maximum absolute error `4.69e-2 s^-1`,
from finite-difference roundoff on entries of order `1e8 s^-1`).  The time
domain eigenvalues are approximately

```text
-16.1647, -5066.72, -5070.76, -8.00792e8, -8.00792e8  [s^-1]
```

so this nominal point is stable; its pole scales are `2.57 Hz`, `806.4 Hz`,
`807.0 Hz`, and `127.45 MHz` (the last pair is electrical).  Earlier broad
search points reported in this file were not stability-valid; the new gate
will warn or reject them before a spectrum is treated as physical.

### Estimator parity result

The synthetic regression uses a sloped known ASD and 512 finite records.  The
experimental-style and simulation-style estimates agree in absolute level
within 4--5% in the passband and in normalized spectral shape within 3% (the
remaining variation is finite-record scatter).  The test explicitly covers
DC, Nyquist, interior `sqrt(2)`, `df`, Hann power correction, mean removal,
`filtfilt`, and power-vs-magnitude averaging.  The full repository suite now
passes `45` tests.

### Recomputed nominal comparison (no fitted parameters)

`noise_model_audit.py` and `midband_residual.py` were run using the repository
nominal input, the experiment's `500 kS/s / 100,000 sample / 10 kHz` analysis
configuration, and `white/readout ASD = 0`.  The table reports the pre-analysis
normalized five-state CH0 ASD divided by the rebuilt experimental ASD; the
common digital Bessel response cancels in this ratio below the diagnostic
cutoff.

The generated CSV/PNG under
`PoST_Simulations/diagnostics/noise_comparison_nominal_20260905/` are local
diagnostics, not committed artifacts: `.gitignore` excludes `*.csv` and the
PoST-specific ignore excludes `*.png`. The versioned target-case summary is
instead under
`PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/`.

| Frequency | Experiment normalized | Simulation normalized | Simulation / experiment |
| ---: | ---: | ---: | ---: |
| 10 Hz | 46.1561 | 1.49635 | 0.0324 |
| 100 Hz | 4.02432 | 1.48496 | 0.3690 |
| 1 kHz | 1.00000 | 1.00000 | 1.0000 |
| 3 kHz | 0.93566 | 0.57709 | 0.6168 |
| 5 kHz | 0.86256 | 0.49675 | 0.5759 |
| 7 kHz | 0.78273 | 0.47017 | 0.6007 |
| 10 kHz | 0.63648 | 0.45394 | 0.7132 |

This is a diagnostic comparison, not a parameter fit.  It shows that fixing
the estimator parity does not by itself remove the broad 1--10 kHz difference.

### Physical-noise and circuit audit

The current source matrix treats each independent source column as an ASD,
then adds independent columns in PSD.  TES Johnson noise remains
`sqrt(4*k_B*T*R*(1+2*beta)*(1+M^2))`; the old `(1+beta)^2` expression was not
reintroduced.  The TES voltage source's electrical and Joule-heating entries
retain their correlated relative sign.  Thermal-link sources enter the two
connected nodes with opposite signs.  `alpha` and `beta` are used as the local
logarithmic derivatives `d ln R/d ln T` and `d ln R/d ln I`, respectively; the
new nonlinear RHS makes this convention explicit.

The modeled electrical/readout chain contains `R_l`, one series `L`, the TES
electrothermal states, and direct TES-current output.  Repository search found
no configured frequency-dependent bias impedance, shunt/stray capacitance,
wiring resonance, SQUID input impedance, FLL transfer, amplifier voltage or
current noise, source-impedance-dependent readout noise, or additional
pole/zero.  No such term was added, because no independent circuit constant
in the repository uniquely determines one.

### Unresolved items and most likely causes

1. The previous `experimental_asd()` magnitude-average bug was confirmed and
   fixed, but the remaining nominal mismatch is much larger than its expected
   estimator bias; it is not sufficient as the sole explanation.
2. The nominal comparison uses repository parameters, while the external
   experiment's independently measured operating-point/circuit parameter set
   is not versioned here.  The stability-validity and parameter-definition
   parity must be established from DC/pulse/impedance information.
3. The model has no independently specified readout/bias transfer function.
   Such a transfer may be a physical cause, but adding a pole/zero or fitting
   it to this spectrum is intentionally left unresolved pending an external
   injected-signal or circuit measurement.

The most likely root causes at present are therefore: (a) the former
experimental estimator inconsistency, now corrected; (b) operating-point or
parameter-definition mismatch, especially stability and the actual bias
impedance; and (c) an unmodeled but independently measurable readout/bias
chain.  No empirical noise source or fitted transfer was promoted to the
production model.

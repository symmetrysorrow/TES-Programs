# Pre-analysis readout biquad falsification diagnostic

## Purpose

The 10 kHz digital analysis transfer has been independently validated from the
raw CH0 records, and replacing the stored `modelnoise.txt` with a fresh
canonical reconstruction does not remove the characteristic 5--15 kHz model
deficit plus 40--100 kHz model excess.

This diagnostic therefore moves one stage upstream. It compares the frozen
detector model directly with the fresh **pre-analysis** experimental ASD and
tests whether a small, physically structured pre-ADC readout/electrical
transfer can supply the missing shape.

The production detector parameters, production 100 kHz analog Bessel response,
alias bookkeeping, and post-filter white term are frozen.

## Experimental target

The target is reconstructed directly from the exact
`accepted_record_indices` stored in the target-case
`comparison_summary.json`.

For every accepted CH0 raw record the target estimator is:

```text
mean removal
-> Hann window
-> rFFT power
-> power average
-> one-sided ASD
```

No 10 kHz digital Bessel is applied. The target and model are normalized at
1 kHz, so unresolved absolute voltage/current calibration does not enter the
shape diagnostic.

## Model placement

The frozen model path is:

```text
intrinsic detector noise
-> production 100 kHz / order-4 / mag-normalized analog Bessel
-> diagnostic pre-ADC biquad section(s)
-> ADC first-alias fold
-> production post-filter white term
-> normalized pre-analysis ASD
```

The diagnostic biquad is evaluated separately at the main frequency and at the
first-alias frequency before the two PSD contributions are folded.

The production `post_filter_white` term retains its existing semantics and is
added after the hardware/alias stage. It is therefore deliberately not shaped
by the diagnostic pre-ADC section.

## Biquad parameterization

Each section is a DC-normalized analog second-order pole/zero pair:

```text
[(s/wz)^2 + s/(Qz*wz) + 1]
--------------------------------
[(s/wp)^2 + s/(Qp*wp) + 1]
```

with four parameters:

- pole center `f_p`;
- pole Q `Q_p`;
- zero center `f_z`;
- zero Q `Q_z`.

Both pole and zero Q are constrained above 0.5, so each section represents
complex-conjugate left-half-plane roots. With positive center frequencies this
is stable and minimum phase.

If `f_p = f_z` and `Q_p = Q_z`, the section is exactly unity. The DC gain is
always one, so there is no free arbitrary amplitude parameter.

Default search box:

```text
center frequency: 1 kHz -- 300 kHz
Q:                0.55 -- 20
```

Only three model families are tested:

- no added section;
- one biquad, four free parameters;
- two biquads, eight free parameters.

This intentionally avoids returning to an unconstrained generic high-order
pole/zero correction.

## Optimization and outputs

Each non-unity family is fit with differential evolution followed by a bounded
least-squares polish. The detector model is evaluated only once to construct
the frozen intrinsic context; the optimizer changes only the diagnostic
readout sections.

The JSON reports:

- baseline pre-analysis shape score and bands;
- one- and two-biquad fitted parameters;
- score ratio to baseline;
- 1--5, 5--15, 15--40, 40--100, and 100--200 kHz residuals;
- whether each family moves the mid/high residuals in the required direction;
- whether it improves both mid and high absolute residuals;
- whether it also avoids worsening the 100--200 kHz tail;
- parameter positions relative to the search bounds;
- a configurable 1 dB RMS / 3 dB max residual screen;
- representative end-to-end model corrections at 1, 5, 10, 20, 40, 70,
  100, 150, and 200 kHz.

The screen tolerance is diagnostic only and is not a circuit prior.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/preanalysis_readout_biquad_diagnostic.py --summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json --comparison-summary PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/preanalysis_readout_biquad_diagnostic.json
```

If the experimental drive is mounted elsewhere, supply `--experiment-path`.

The center/Q search box and numerical effort can also be changed explicitly
with `--center-min-hz`, `--center-max-hz`, `--q-min`, `--q-max`, and
`--de-maxiter`.

## Interpretation

A strong one-biquad improvement would be especially informative: it would show
that a single structured pre-ADC second-order feature has enough leverage to
produce the observed shoulder/dip pattern without changing detector physics.

If one section fails but two sections succeed, the result supports a more
structured readout/electrical-transfer explanation but is less identifying and
requires independent hardware evidence before assigning circuit meaning to the
fitted centers or Q values.

If neither one nor two biquads jointly improve the 5--15 and 40--100 kHz
residuals, the simple single-path pre-ADC resonance/lead-lag hypothesis is
weakened substantially.

A successful magnitude fit does **not** identify the actual circuit or phase
response. The fitted centers and Q values remain diagnostic until checked
against hardware schematics, measured transfer functions, or readout
configuration data.

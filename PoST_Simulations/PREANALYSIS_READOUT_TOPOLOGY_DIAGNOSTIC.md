# Pre-analysis readout topology comparison diagnostic

## Purpose

The previous pre-analysis biquad diagnostic showed that a structured transfer
placed before ADC alias folding has very large leverage on the residual shape:

- one complex biquad reduced the shape score by roughly an order of magnitude;
- two complex biquads reduced it by roughly two orders of magnitude;
- several fitted Q values landed at the imposed lower bound, Q = 0.55.

That boundary behavior makes the next question topological rather than merely
parametric: does the data actually prefer a resonant complex-conjugate section,
or is it asking for an overdamped / real-pole network more like an RC/RL,
anti-alias, amplifier, or other readout transfer?

This diagnostic answers that question while keeping detector physics frozen.

## Shared placement and target

Every tested transfer is placed at the same location:

```text
frozen intrinsic detector ASD
-> production 100 kHz order-4 mag-normalized analog Bessel
-> diagnostic pre-ADC transfer topology
-> first ADC alias fold
-> production post-filter white term
-> normalized pre-analysis ASD
```

The experimental target is reconstructed from the exact accepted CH0 raw
record indices stored in `comparison_summary.json`:

```text
mean removal
-> Hann window
-> rFFT power
-> power average
-> one-sided ASD
```

No 10 kHz digital analysis filter is applied to either target or model.

All diagnostic transfer families have DC gain exactly one. There is no free
overall amplitude parameter.

## Matched four-parameter core comparison

Three topologies have exactly four free parameters each.

### 1. `complex_biquad_4p`

A second-order pole pair and second-order zero pair:

```text
[(s/wz)^2 + s/(Qz*wz) + 1]
--------------------------------
[(s/wp)^2 + s/(Qp*wp) + 1]
```

with:

- `f_p`;
- `Q_p`;
- `f_z`;
- `Q_z`.

For this family Q is restricted to 0.55--20, so both pole and zero roots remain
complex-conjugate.

This is the direct matched-complexity reference to the previous one-biquad
diagnostic.

### 2. `general_second_order_4p`

The transfer equation and four parameters are identical, but Q is allowed over
0.10--20.

For a positive-Q second-order factor:

- Q > 0.5 gives a complex-conjugate pair;
- Q = 0.5 is critically damped;
- Q < 0.5 gives two negative real roots.

When Q <= 0.5, the output JSON also reports the equivalent two real-root corner
frequencies. This makes an overdamped fit directly interpretable as a real-pole
or real-zero network rather than a resonance.

### 3. `real_2p2z_4p`

Two independent real first-order poles and two independent real first-order
zeros:

```text
(1 + s/wz1)(1 + s/wz2)
-----------------------
(1 + s/wp1)(1 + s/wp2)
```

The four corner frequencies are the only free parameters.

This is a direct stable/minimum-phase RC/RL-like reference with the same number
of free parameters as both second-order families.

## Secondary six-parameter family

`hybrid_general_plus_leadlag_6p` combines:

- one general second-order pole/zero section;
- one additional real first-order lead/lag pole/zero pair.

It is included only as a complexity-ladder check. It is **not** considered a
matched-DOF winner against the four-parameter core families.

## Search box

Defaults:

```text
all center/corner frequencies: 1 kHz -- 300 kHz
complex biquad Q:             0.55 -- 20
general second-order Q:       0.10 -- 20
```

Each non-unity family is optimized with differential evolution followed by a
bounded least-squares polish.

## Output

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/preanalysis_readout_topology_diagnostic.json
```

Important fields include:

- `best_matched_4p_family`;
- `best_overall_family`;
- every fitted parameter and named boundary position;
- root regime for each general second-order pole/zero factor;
- equivalent real-root frequencies whenever Q <= 0.5;
- score ratio to the frozen baseline;
- band residuals in 1--5, 5--15, 15--40, 40--100, and 100--200 kHz;
- `topology_evidence.general_second_order_prefers_Q_below_0p5`;
- `topology_evidence.general_score_ratio_to_complex_4p`;
- `topology_evidence.real_2p2z_score_ratio_to_complex_4p`;
- `topology_evidence.best_4p_is_nonresonant_or_overdamped`;
- representative correction values from 1 to 200 kHz.

The same diagnostic screen as the previous readout test is retained:

```text
overall RMS residual <= 1 dB
overall max |residual| <= 3 dB
```

This is a diagnostic tolerance, not a physical prior.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/preanalysis_readout_topology_diagnostic.py --summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json --comparison-summary PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json
```

If the raw-data drive is mounted elsewhere, add `--experiment-path`.

## Interpretation

### General second order wins and Q < 0.5

This is strong evidence that the previous Q = 0.55 boundary was not merely a
numerical accident. The preferred transfer is overdamped and can be represented
by real roots. The next step should be to compare those inferred real corner
frequencies with known RC/RL, anti-alias, amplifier, wiring, SQUID/readout, or
other hardware constants.

### `real_2p2z_4p` wins

This is even cleaner evidence against a required resonance. A four-corner real
network explains the shape better than a same-parameter-count complex biquad.
The next diagnostic should move from generic real poles/zeros to circuit
families constrained by the actual schematic.

### Complex biquad remains best

Then the resonance-like interpretation survives the matched-complexity test.
The fitted center frequencies and Q values should be checked against hardware
features before any production-model change.

### Six-parameter hybrid improves substantially beyond every four-parameter model

Then one simple four-parameter topology is insufficient, but the required
complexity is still lower than the previous two-biquad eight-parameter model.
The fitted second-order and lead/lag scales should be compared with the
~10--20 kHz and ~80--150 kHz features seen in earlier diagnostics.

## Guardrail

A successful topology fit demonstrates shape leverage at the stated pre-ADC
location. It does not identify a unique physical circuit, phase response, or
component value. Production code should not adopt a diagnostic transfer until
the fitted topology is supported by schematic information, direct transfer
measurements, or cross-dataset consistency.

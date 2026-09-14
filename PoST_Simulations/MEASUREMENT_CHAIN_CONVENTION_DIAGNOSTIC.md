# Measurement-chain convention diagnostic

## Purpose

The detector-side parameter screens did not produce a stable way to raise the
5--15 kHz deficit while lowering the 40--100 kHz excess. A subsequent
low-order residual-transfer fit showed strong leverage, with fitted structure
clustered near the known 10 kHz digital analysis stage and the 80--100 kHz
analog hardware region.

This diagnostic therefore freezes every detector parameter and changes only
the known measurement-chain conventions. It asks whether reasonable
cutoff/order/pass/normalization alternatives can reproduce the residual shape
without introducing a new detector state or noise source.

It is a falsification screen, not a production refit.

## Frozen detector physics

The script reads `best_case_parameters` from the optimizer `summary.json`.
TES operating point, heat capacities, thermal conductances, alpha, beta,
inductance, effective series resistance, all physical noise-source amplitudes,
and the post-filter white term are held fixed.

The intrinsic TES ASD is computed once. Every grid point reuses that same
intrinsic spectrum and changes only the analog/digital measurement-chain
response.

## Production baseline

The production chain is explicitly included as one grid point:

- analog Bessel cutoff: 100 kHz;
- analog order: 4;
- analog normalization: `mag`;
- digital analysis Bessel cutoff: 10 kHz;
- digital order: 2;
- digital passes: 2, representing the forward/backward analysis filter.

The diagnostic independently reconstructs this baseline and compares it with
`Opt_noise.deterministic_simulated_spectrum`. The JSON reports the maximum
relative difference and a boolean requiring agreement to 1e-10.

## Grid

Analog stage:

- cutoff: 80, 100, 120 kHz;
- order: 2, 4, 6;
- normalization: `phase`, `mag`, `delay`.

Digital analysis stage:

- cutoff: 8, 10, 12 kHz;
- order: 1, 2, 3, 4;
- passes: 1 or 2.

This produces 648 deterministic cases.

The non-production variants are convention diagnostics. Their inclusion does
not imply that the actual acquisition used those settings.

## Run

From the repository root in PowerShell:

```powershell
python PoST_Simulations/subScript/measurement_chain_convention_diagnostic.py --summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/measurement_chain_convention_diagnostic.json
```

## Decision fields

For every chain configuration the JSON records the production shape score and
the same five broad-band residual diagnostics used elsewhere in the noise
analysis.

It also reports whether any measurement-chain-only case:

- raises the 5--15 kHz deficit while lowering the 40--100 kHz excess;
- reduces the absolute mean residual in both bands;
- does both without worsening 100--200 kHz.

The output includes the global best row, best correct-direction row, best
strict row, and the best row within each analog-normalization and
digital-pass convention.

## Interpretation guardrail

A positive result is evidence that the known measurement/analysis chain has
enough leverage to explain part of the residual. It is not permission to
change the production convention until the actual acquisition and analysis
implementation is verified independently.

A negative result would substantially weaken the hypothesis that the
remaining discrepancy is merely a cutoff/order/pass/normalization
interpretation error in the existing single-path Bessel chain.

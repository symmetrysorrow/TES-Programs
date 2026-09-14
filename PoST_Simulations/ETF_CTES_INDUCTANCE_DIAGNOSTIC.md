# ETF / C_tes / inductance boundary-stress diagnostic

## Purpose

The IV-power-preserving thermal-path scan found no stable combination of
`T_c`, `n`, derived `G_tes-bath`, and `G_tes-stycast` that raised the
5--15 kHz deficit while lowering the 40--100 kHz excess. The production
best-fit also places `alpha` at its upper bound, `C_tes` at its lower
bound, and `L` near its upper bound.

`subScript/etf_ctes_inductance_diagnostic.py` therefore tests whether moving
those dynamical parameters back into the allowed interior can create the
required spectral shape without changing the DC/thermal operating point.

## What is varied

- `alpha`: from a configurable fraction of the best-fit value up to the
  best-fit value. The default lower fraction is 0.05. The scan never extends
  above the fitted upper-bound solution.
- `C_tes`: the full production physical range,
  `0.25--5 times C_TES_MATERIAL_J_PER_K`.
- `L`: the full production inductance range.

The exact best-fit value is inserted on every axis.

All thermal operating parameters, TES resistance, `beta`, effective series
resistance, noise-source amplitudes, and transfer filters remain fixed.
The output records TES-current and Joule-power ratios as guardrails.

## Run

From the repository root in PowerShell:

```powershell
python PoST_Simulations/subScript/etf_ctes_inductance_diagnostic.py --summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/etf_ctes_inductance_diagnostic.json
```

Resolution can be changed with `--alpha-points`, `--ctes-points`, and
`--l-points`. The inward alpha extent can be changed with
`--alpha-min-fraction`.

## Decision fields

The JSON reports whether any stable trial:

- moves the 5--15 kHz deficit upward and 40--100 kHz excess downward;
- improves the absolute mean residual in both bands;
- does both without worsening 100--200 kHz.

It also reports the best global shape-score row and the best correct-direction
row.

A positive result shows that the present reduced topology has dynamical
leverage when ETF strength, TES thermal inertia, and electrical inertia are
changed together. It is not a physical measurement of `alpha`, `C_tes`, or
`L`; those parameters still require pulse or complex-impedance constraints.

A negative result, after the earlier source, RC, thermal-node, thermal-path,
and IV-compatible thermal-law screens, is evidence that the reduced
electrothermal topology or independently assumed transfer model needs to be
revisited rather than adding more unconstrained fit freedom.

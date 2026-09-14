# IV-power-preserving thermal-path diagnostic

## Why this diagnostic exists

The previous thermal-path partition scan found no stable point where changing
only `G_tes-stycast` and `G_stycast-abs` simultaneously improved the
5--15 kHz deficit and the 40--100 kHz excess. The next question is whether the
dominant direct TES--bath branch and electrothermal feedback can supply the
missing spectral leverage without violating the measured DC operating point.

## What is varied

`subScript/iv_power_thermal_path_diagnostic.py` scans:

- `T_c` across the same target-case range recorded in `summary.json`;
- thermal exponent `n` across the production bounds;
- `G_tes-stycast` across the production physical bounds.

For every `T_c, n` pair, `G_tes-bath` is recomputed with
`Opt_noise.g_tes_bath_from_joule_power()` using the same-campaign IV
`P_J`. Therefore the trial TES current and Joule power remain tied to the IV
operating point rather than becoming new noise-fit degrees of freedom.

`T_bath`, TES resistance, `alpha`, `beta`, `L`,
`R_series_eff`, `G_stycast-abs`, heat capacities, noise-source amplitudes,
and measurement filters remain fixed at the best-fit values.

## Run

From the repository root in PowerShell:

```powershell
python PoST_Simulations/subScript/iv_power_thermal_path_diagnostic.py --summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json
```

The default output is:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/iv_power_thermal_path_diagnostic.json
```

Default resolution is 7 points on each axis before insertion of the exact
best-fit value. It can be changed with `--tc-points`, `--n-points`, and
`--g-points`.

## Required checks

The output explicitly records:

- maximum IV Joule-power closure error;
- number of stable points;
- best global shape score;
- whether 5--15 kHz and 40--100 kHz move in the required opposite directions;
- whether both absolute band residuals improve;
- whether that improvement also avoids worsening 100--200 kHz.

A positive screen only establishes leverage inside the reduced model. It does
not independently determine `T_c`, `n`, `G_tes-bath`, or
`G_tes-stycast`; promotion still requires pulse, DC, or complex-impedance
constraints.

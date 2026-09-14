# Thermal-path partition diagnostic

## Purpose

The constrained noise fit leaves a 5--15 kHz model deficit together with a
40--100 kHz model excess. Existing source-class rescaling, the TES hanging-body
screen, and the TES--Stycast series thermalization-node screen do not solve
that pair without a high-frequency penalty.

The next test therefore adds no source, state, or production fit parameter.
It asks whether the existing internal thermal path has the required spectral
leverage when its two fitted conductances are moved while every other best-fit
parameter is frozen.

## Scan

`subScript/thermal_path_partition_diagnostic.py` scans:

- `G_tes-stycast`
- `G_stycast-abs`

over the same physical bounds already defined in `Opt_noise.py`. The exact
best-fit values are inserted into both logarithmic grids.

`G_tes-bath`, `T_c`, `T_bath`, `alpha`, `beta`, `L`,
`R_series_eff`, heat capacities, Johnson terms, and filter settings remain
fixed. The diagnostic records current and Joule-power ratios as a guardrail
against silently moving the reduced-model operating point.

## Run

Run `Opt_noise.py` first, then:

```text
python subScript/thermal_path_partition_diagnostic.py ^
  --summary .noise_optimization_work_rsh_sweep/summary.json
```

The default output is:

```text
.noise_optimization_work_rsh_sweep/thermal_path_partition_diagnostic.json
```

Use `--grid-points N` to change the logarithmic resolution on each axis.

## Decision fields

The output reports three progressively stronger tests:

1. `can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction`
2. `can_improve_5_15k_and_40_100k_absolute_error`
3. `can_improve_5_15k_and_40_100k_without_worsening_100_200k`

A negative result means that simple repartitioning of the existing
TES--Stycast--absorber conductance path cannot repair the alternating residual.
A positive result is only a leverage screen; it is not a physical conductance
measurement and still requires independent pulse, DC, or complex-impedance
constraints before promotion into the production model.

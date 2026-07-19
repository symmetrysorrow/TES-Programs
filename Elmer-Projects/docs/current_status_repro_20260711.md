# Current Status And Reproduction

Last updated: 2026-07-11 (evening update: converges to T0)

## Purpose

This note records the current TES steady-state comparison results and the exact way to reproduce them from the current workspace state.

## Update 2026-07-11 evening: steady state now converges to T0

The circuit-coupled steady state now converges to `T0` without changing any
physical parameter (`I_bias`, `R_sh`, `R0`, `alpha`, `beta`, `I0`, `T0`, `Tc`,
`G0`, `T_bath` are all unchanged). Two changes were made:

### 1. Bug fix: MATC operator precedence in the membrane k(T) guard

`sync_elmer_parameters.py` wraps `T` in the membrane conductivity expression
with the guard `0.5*((tx+1.0e-12)+abs(tx-1.0e-12))` (= `max(tx, 1e-12)`).
Because `^` binds tighter than `*` in MATC, the generated expression
evaluated as `0.5*(2*tx)^3.252` instead of `tx^3.252`, i.e. the membrane
conductivity actually used by Elmer was `2^2.252 = 4.7634` times the intended
value. This was confirmed by evaluating the generated expression directly
with `matc.exe` (`k(T0)` came out `1.35e-6` instead of the intended
`2.826e-7`).

This is why the circuit-coupled case previously latched at the bath
temperature (`150.015 mK`): the membrane conducted ~4.8x too well, the TES
could not stay in the transition, and `R -> R_min`.

Fix (one line in `membrane_matc_expr`): the guard is now parenthesized,
`(0.5*((tx+1.0e-12)+abs(tx-1.0e-12)))^(...)`, so the MATC expression matches
the Python-side evaluation exactly (verified with `matc.exe`).

### 2. Calibration: one scalar factor on the membrane k(T) geometric prefactor

After the bug fix, the intended formula was already close: only a `1.0181`
multiplier on the membrane k(T) expression is needed so that the FEM
steady state balances `P0 = I0^2 R0` exactly at `T0`. The expression in
`elmer_project.json` is now:

```text
1.0181*G0*T**(4.252-1)*0.4*(membrane_dx-TES_Au_dx)/(8*(SiNx_dz+Si_1_dz+SiO2_1_dz)*TES_Au_dx)
```

The factor was found by a secant iteration on the fixed-power case using the
power-law relation `P ~ (T^4.252 - T_bath^4.252)`:

| membrane k scale (vs intended formula) | fixed-power TES volume average |
|---:|---:|
| `4.7634` (old buggy expression) | `157.008 mK` |
| `1.5883` | `164.538 mK` |
| `1.0755` | `168.027 mK` |
| `1.0181` | `168.550 mK` |

Since the TES heat-source UDF returns exactly `I = I0`, `R = R0`,
`P = I0^2 R0` at `T = T0`, balancing the thermal side at `T0` makes `T0` the
equilibrium of the circuit-coupled case as well. Electrothermal feedback
(loop gain ~ 24) further suppresses the remaining thermal-side residual.

### Result

| Case | TES volume average | vs `T0 = 168.570 mK` |
|---|---:|---:|
| `k(T)` + fixed `I0^2 R0` | `168.550 mK` | `-0.020 mK` |
| `k(T)` + circuit-coupled | `168.563 mK` | `-0.007 mK` |

Output files: `mesh_shifted_merged/case_constant_power_fixed_t0001.vtu`,
`mesh_shifted_merged/case_constant_power_t0001.vtu` (extracted with
`scripts/analysis/extract_tes_volume_avg.py`).

Note: the four-case comparison below (including the `fixedk` variants and the
`150.015 / 157.008 / 170.441 / 180.955 mK` table) is the state BEFORE this
fix and is kept for the record. `generated/tes_materials_fixedk.sif` still
contains the pre-fix fixed k value `2.826e-7`; the post-fix equivalent
`k(T0)` is `2.877e-7`.

---

## Pre-fix comparison (historical record below this line)

## Current reference conditions

Source files:

- [`elmer_project.json`](D:/github/Elmer-Projects/elmer_project.json)
- [`elmer_geometry.json`](D:/github/Elmer-Projects/elmer_geometry.json)
- [`generated/tes_shared_variables.sif`](D:/github/Elmer-Projects/generated/tes_shared_variables.sif)
- [`generated/tes_case_constant_power.sif`](D:/github/Elmer-Projects/generated/tes_case_constant_power.sif)

Important parameters at the time of this note:

- `I_bias = 0.000715 A`
- `R_sh = 0.0039 ohm`
- `R0 = 0.015527 ohm`
- `R_min = 1e-06 ohm`
- `alpha = 256.46`
- `beta = 5.03`
- `I0 = 0.000143537344932311 A`
- `T0 = 0.16857 K`
- `Tc = 0.1709 K`
- `T_bath = 0.15 K`
- `TES_volume = 4e-14 m^3`
- `G0 = 7.854e-08`
- `T_initial_constant_power = 0.168 K`

Steady-state solver settings:

- `constant_power_steady_state_max_iterations = 150`
- `constant_power_solver_nonlinear_max_iterations = 200`
- `constant_power_solver_nonlinear_convergence_tolerance = 0.001`
- `constant_power_solver_nonlinear_relaxation_factor = 0.1`
- `constant_power_solver_steady_state_convergence_tolerance = 1e-09`

Fixed-power reference used here:

- `constant_power_total_power = 3.199023057219016e-10 W`
- `constant_power_volumetric_heat_source = 7997.55764304754 W/m^3`

This fixed power is `I0^2 * R0`.

## Cases compared

### 1. Temperature-dependent membrane k + circuit-coupled TES heating

- SIF: [`case_constant_power.sif`](D:/github/Elmer-Projects/case_constant_power.sif)
- Materials include: [`generated/tes_materials.sif`](D:/github/Elmer-Projects/generated/tes_materials.sif)

### 2. Temperature-dependent membrane k + fixed TES heating `I0^2 R0`

- SIF: [`case_constant_power_fixed.sif`](D:/github/Elmer-Projects/case_constant_power_fixed.sif)
- Materials include: [`generated/tes_materials.sif`](D:/github/Elmer-Projects/generated/tes_materials.sif)

### 3. Fixed membrane k + circuit-coupled TES heating

- SIF: [`case_constant_power_fixedk.sif`](D:/github/Elmer-Projects/case_constant_power_fixedk.sif)
- Materials include: [`generated/tes_materials_fixedk.sif`](D:/github/Elmer-Projects/generated/tes_materials_fixedk.sif)

### 4. Fixed membrane k + fixed TES heating `I0^2 R0`

- SIF: [`case_constant_power_fixed_fixedk.sif`](D:/github/Elmer-Projects/case_constant_power_fixed_fixedk.sif)
- Materials include: [`generated/tes_materials_fixedk.sif`](D:/github/Elmer-Projects/generated/tes_materials_fixedk.sif)

## Fixed-k definition used in this comparison

The fixed membrane conductivity used here is:

- `Heat Conductivity = 2.825930416659195e-07`

This was set in [`generated/tes_materials_fixedk.sif`](D:/github/Elmer-Projects/generated/tes_materials_fixedk.sif).

## Current results

TES temperature was evaluated with:

- [`scripts/analysis/extract_tes_volume_avg.py`](D:/github/Elmer-Projects/scripts/analysis/extract_tes_volume_avg.py)

Using TES volume average:

| Case | Mesh result file | TES volume average |
|---|---|---:|
| `k(T)` + circuit-coupled | `mesh_shifted_merged/case_constant_power_t0001.vtu` | `0.150015215520432 K` |
| `k(T)` + fixed `I0^2 R0` | `mesh_shifted_merged/case_constant_power_fixed_t0001.vtu` | `0.157007604667779 K` |
| fixed `k` + circuit-coupled | `mesh_shifted_merged/case_constant_power_fixedk_t0001.vtu` | `0.170440568846491 K` |
| fixed `k` + fixed `I0^2 R0` | `mesh_shifted_merged/case_constant_power_fixed_fixedk_t0001.vtu` | `0.180955321417300 K` |

In mK:

| Case | TES volume average |
|---|---:|
| `k(T)` + circuit-coupled | `150.015 mK` |
| `k(T)` + fixed `I0^2 R0` | `157.008 mK` |
| fixed `k` + circuit-coupled | `170.441 mK` |
| fixed `k` + fixed `I0^2 R0` | `180.955 mK` |

## Interpretation at this point

- Membrane `k(T)` has a large effect on the final steady temperature.
- In both `k(T)` and fixed-`k` comparisons, the circuit-coupled case converged to a lower temperature than the fixed-power `I0^2 R0` case.
- Therefore, the difference is not caused by membrane `k` alone. The TES heat-source model also changes the equilibrium.

## Reproduction steps

Run from:

```powershell
cd D:\github\Elmer-Projects
```

### A. Run all four steady-state cases

```powershell
ElmerSolver case_constant_power.sif
ElmerSolver case_constant_power_fixed.sif
ElmerSolver case_constant_power_fixedk.sif
ElmerSolver case_constant_power_fixed_fixedk.sif
```

Expected output files:

- `mesh_shifted_merged/case_constant_power_t0001.vtu`
- `mesh_shifted_merged/case_constant_power_fixed_t0001.vtu`
- `mesh_shifted_merged/case_constant_power_fixedk_t0001.vtu`
- `mesh_shifted_merged/case_constant_power_fixed_fixedk_t0001.vtu`

### B. Extract TES volume-average temperatures

```powershell
python scripts/analysis/extract_tes_volume_avg.py mesh_shifted_merged mesh_shifted_merged/case_constant_power_t0001.vtu
python scripts/analysis/extract_tes_volume_avg.py mesh_shifted_merged mesh_shifted_merged/case_constant_power_fixed_t0001.vtu
python scripts/analysis/extract_tes_volume_avg.py mesh_shifted_merged mesh_shifted_merged/case_constant_power_fixedk_t0001.vtu
python scripts/analysis/extract_tes_volume_avg.py mesh_shifted_merged mesh_shifted_merged/case_constant_power_fixed_fixedk_t0001.vtu
```

Expected TES volume averages:

```text
case_constant_power_t0001.vtu              0.150015215520432
case_constant_power_fixed_t0001.vtu        0.157007604667779
case_constant_power_fixedk_t0001.vtu       0.170440568846491
case_constant_power_fixed_fixedk_t0001.vtu 0.180955321417300
```

## Notes

- The comparison above assumes the current `mesh_shifted_merged` mesh and the current generated include files.
- If [`sync_elmer_parameters.py`](D:/github/Elmer-Projects/sync_elmer_parameters.py) is run again, generated files may change.
- If [`elmer_project.json`](D:/github/Elmer-Projects/elmer_project.json) or [`elmer_geometry.json`](D:/github/Elmer-Projects/elmer_geometry.json) is edited, this note may no longer describe the active state.

## Files added for this comparison

- [`case_constant_power_fixedk.sif`](D:/github/Elmer-Projects/case_constant_power_fixedk.sif)
- [`case_constant_power_fixed_fixedk.sif`](D:/github/Elmer-Projects/case_constant_power_fixed_fixedk.sif)
- [`generated/tes_materials_fixedk.sif`](D:/github/Elmer-Projects/generated/tes_materials_fixedk.sif)

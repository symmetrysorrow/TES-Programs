# Phase24 absolute operating-point offset diagnostic

Generated: `2026-09-16T10:18:17.600264+00:00`

Baseline window: `-2..0 us` relative to the event; common post-event window is reported per comparison.

## Baselines

| series | T [K] | I [uA] | R [ohm] | P [W] |
|---|---:|---:|---:|---:|
| HYPRE 5e-7 | 0.168448309447 | 147.80143942 | 0.0149725153889 | 3.27078188213e-10 |
| same Phase24 Direct MUMPS | 0.168549546881 | 144.268506298 | 0.0154416808299 | 3.21393818706e-10 |
| CPU/MUMPS Phase23 | 0.168459765942 | 147.375245042 | 0.0150112483744 | 3.26036275378e-10 |

## Pairwise decomposition

| comparison | quantity | baseline offset | max raw difference | max baseline-corrected difference |
|---|---|---:|---:|---:|
| hypre_vs_mumps | T [K] | -0.000101237434 | 0.000108812747 | 7.57531362e-06 |
| hypre_vs_mumps | I [uA] | 3.53293312 | 3.44622496 | 0.142527967 |
| hypre_vs_mumps | R [ohm] | -0.000469165441 | 0.000709752426 | 0.000240586985 |
| hypre_vs_mumps | P [W] | 5.68436951e-12 | 5.9161009e-12 | 5.7025344e-12 |
| hypre_vs_cpu_reference | T [K] | -1.14564949e-05 | 2.47622663e-05 | 3.62187612e-05 |
| hypre_vs_cpu_reference | I [uA] | 0.426194378 | 0.68234266 | 1.10853704 |
| hypre_vs_cpu_reference | R [ohm] | -3.87329855e-05 | 0.000262186828 | 0.000252408214 |
| hypre_vs_cpu_reference | P [W] | 1.04191283e-12 | 5.27446795e-12 | 6.31638079e-12 |
| mumps_vs_cpu_reference | T [K] | 8.97809388e-05 | 0.000121590474 | 3.18095354e-05 |
| mumps_vs_cpu_reference | I [uA] | -3.10673874 | 4.07274782 | 0.966009071 |
| mumps_vs_cpu_reference | R [ohm] | 0.000430432455 | 0.000656235075 | 0.000225802619 |
| mumps_vs_cpu_reference | P [W] | -4.64245667e-12 | 5.27352222e-12 | 6.31065545e-13 |

## Current checkpoints

| comparison | 0.1 us corrected [uA] | 0.5 us corrected [uA] | 0.9 us corrected [uA] |
|---|---:|---:|---:|
| hypre_vs_mumps | -0.0899181804 | -0.104204376 | -0.142400227 |
| hypre_vs_cpu_reference | -0.319033447 | -0.557508027 | -1.10704105 |
| mumps_vs_cpu_reference | -0.229115266 | -0.45330365 | -0.964640822 |

Interpretation: a large current baseline offset with small corrected T/R/P differences indicates readout or circuit operating-point extraction. A large corrected temperature or resistance difference indicates a state/solution mismatch. A large corrected current difference only, with matching state variables, indicates current extraction or circuit bookkeeping.

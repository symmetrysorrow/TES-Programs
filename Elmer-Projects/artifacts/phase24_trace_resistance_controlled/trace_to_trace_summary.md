# Phase24 trace-to-trace resistance decomposition

Trace temperatures are area-weighted named-boundary averages from the captured CPU/MUMPS solutions. R is the differential d(Delta T_trace)/dP from the 0.95/1.00/1.05 P bracket.

The Membrane-to-bath endpoint is the active second interval. The Si1/SiNx and Si2/SiO2_2 rows are parallel downstream branches; they are reported for localization and must not be summed as a series chain.

The direct TES/Membrane trace audit is: historical 4274/6702 faces at mean
edge 11.8/12.6 um, refined parent 308/308 faces at 43.4 um, and the usable
Membrane/substrate control 308/310 faces. Thus the control leaves the direct
TES/Membrane trace unchanged. The refined parent and the controls have an
essentially zero temperature jump on this conforming trace; the historical
nonconforming trace retains a finite jump. That jump is a coupling diagnostic,
not an additional material contact resistance.

- `TES_to_Membrane`: historical `1.559330710e+07 K/W`, refined mortar `8.665480596e-07 K/W`, ratio `0.000000`
- `TES_to_Stycast`: historical `-1.027839965e+03 K/W`, refined mortar `-1.096745981e+03 K/W`, ratio `1.067040`
- `Stycast_to_substrate`: historical `2.487114248e-01 K/W`, refined mortar `2.746749377e-01 K/W`, ratio `1.104392`
- `Membrane_to_bath_endpoint`: historical `3.764650845e+07 K/W`, refined mortar `4.926328285e+07 K/W`, ratio `1.308575`
- `Membrane_SiNx_to_Membrane_Si1`: historical `7.123035672e+04 K/W`, refined mortar `1.384436479e+07 K/W`, ratio `194.360458`
- `Membrane_Si1_to_SiO2_1`: historical `8.095475840e+05 K/W`, refined mortar `3.367011813e+07 K/W`, ratio `41.591277`
- `SiO2_1_to_Si_1`: historical `1.000868633e+06 K/W`, refined mortar `9.422211551e+05 K/W`, ratio `0.941403`
- `Si_1_to_SiNx`: historical `-7.653436518e+00 K/W`, refined mortar `-6.661238989e+00 K/W`, ratio `0.870359`
- `SiO2_1_to_Si_2`: historical `3.676544527e+07 K/W`, refined mortar `9.421838589e+05 K/W`, ratio `0.025627`
- `Si_2_to_SiO2_2`: historical `5.402318055e+01 K/W`, refined mortar `5.402319528e+01 K/W`, ratio `1.000000`
- `SiO2_2_to_bath`: historical `2.312155049e+02 K/W`, refined mortar `-2.442967593e+02 K/W`, ratio `-1.056576`

## One-at-a-time controlled tests

| case | G_eff (W/K) | historical ratio | status |
|---|---:|---:|---|
| refined mortar parent | 2.029899315e-08 | 1.080714 | reference |
| TES/Membrane trace refinement | 2.685719122e-08 | 1.429871 | rejected: Stycast top trace also changed (4579 -> 4070 faces) |
| Membrane/substrate trace refinement | 1.959713948e-08 | 1.043348 | accepted one-at-a-time diagnostic; direct TES/Membrane trace retained |

The accepted Membrane/substrate test moves the conductance 46.2% of the way
from the refined parent toward historical. It therefore identifies the
downstream trace neighborhood as a real contributor, but it is not yet
historical-equivalent within the 2--3% criterion.

The controlled one-at-a-time runs are in trace_to_trace_resistance.csv. Their strict validity still depends on the mesh audit: the TES/Membrane variant changed the Stycast top trace face count, so it is not accepted as a clean one-factor result; the Membrane/substrate variant retained the direct TES/Membrane trace and is the usable controlled test.

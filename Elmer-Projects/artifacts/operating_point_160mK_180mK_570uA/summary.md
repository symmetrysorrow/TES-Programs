# TES operating point: 160 mK bath, 180 mK reference, 570 uA reference current

## Requested operating condition

- Bath temperature: `T_bath = 160 mK`
- TES reference temperature: `T_0 = 180 mK`
- TES reference current: `I_0 = 570 uA`
- Existing electrical parameters retained: `R_0 = 15.527 mOhm`, `R_sh = 3.9 mOhm`, `alpha = 256.46`, `beta = 5.03`
- Bias current was set to `I_bias = 2.83933076923077 mA`, so that the shunted steady circuit gives exactly `I_TES = I_0` at `T = T_0` and `R = R_0`.

At the requested reference point, the Joule power is

`P_0 = I_0^2 R_0 = 5.0447223 nW`.

Heat capacities were left unchanged because they do not change the mathematical steady-state equilibrium; they affect transient response and numerical relaxation. The parameter varied here is the effective membrane thermal-conductivity prefactor `G0`.

## Elmer steady-state results

| Case | Membrane scale | Initial T | TES volume-average T | Error from 180 mK | I_TES | R_TES | P_TES | Convergence |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Existing membrane model | 1.000 | 180 mK | 189.023835 mK | +9.023835 mK | 73.332 uA | 147.104 mOhm | 0.791059 nW | `6.47e-10` |
| First thermal-balance estimate | 11.847 | 180 mK | 181.052571 mK | +1.052571 mK | 441.375 uA | 21.188 mOhm | 4.127749 nW | `7.78e-10` |
| Tuned membrane model | 15.402 | 180 mK | 180.058588 mK | +0.058588 mK | 562.445 uA | 15.788 mOhm | 4.994430 nW | `9.25e-10` |
| Tuned model, low initial T | 15.402 | 170 mK | 180.058586 mK | +0.058586 mK | 562.445 uA | 15.788 mOhm | 4.994432 nW | `8.30e-8` |

The 170 mK and 180 mK initial-temperature solutions differ by only `2.74 nK` in TES volume-average temperature. Within the tested basin, the tuned solution is therefore independent of the initial condition.

## Tuned condition

The condition that reproduces the requested operating point on the current geometry and current `R(T,I)` model is

- `G0 = 1.20965674443295e-4 W/K`
- scale relative to the current main project: `15.4017920095868`
- effective membrane conductivity at 180 mK: `k_eff(180 mK) = 5.38752640943619e-4 W/(m K)`
- static conductance implied by `P_0/(T_0-T_bath)`: `2.52236115e-7 W/K`

`G0` and `k_eff` are effective parameters for the combined membrane path represented by the model, not a claim about a single homogeneous material.

## Is this a valid convergence point?

Numerically, yes:

1. The strict run reached an Elmer steady-state relative change of `9.25e-10`.
2. The TES temperature span across the final mesh is only `136.081 uK` (`179.944809–180.080890 mK`), so the lumped TES-temperature approximation is consistent at this operating point.
3. The final circuit state is close to the requested reference state: `T = 180.058588 mK`, `I = 562.445 uA`, `R = 15.787950 mOhm`.
4. A local heat-balance derivative check gives `d(P_J-P_bath)/dT = -1.159e-6 W/K`, which is negative. The point is locally stable against small temperature perturbations.
5. Repeating from `170 mK` returns to the same TES temperature to nanokelvin-level agreement.

## Physical-model caveats

- The main project currently has `T_c = 170.9 mK`, below the requested `T_0 = 180 mK`. The steady `TESHeatSource` UDF uses the linearized `R(T,I)` reference parameters `T_0`, `I_0`, `R_0`, `alpha`, and `beta`; it does not use `T_c`. The computed point is therefore mathematically consistent with the implemented UDF, but `T_c` and the transition model should be updated before treating it as a fully physical TES design.
- The required membrane prefactor is about 15.4 times the current project value. This is an effective calibration for the present geometry/contact model. Before adopting it, compare the implied total conductance with measured IV/load-curve data or a literature-based layer-by-layer boundary-conductance model.
- Only the steady state was calibrated. Pulse shape, recovery time, and stability under the transient/inductive circuit still require a restarted transient calculation using this steady field.

## Reproduction files

- `elmer_project_op160_180_570_base.json`
- `elmer_project_op160_180_570_k11p847.json`
- `elmer_project_op160_180_570_k15p402.json`
- `elmer_project_op160_180_570_k15p402_init170.json`
- Detailed numerical values: `artifacts/operating_point_160mK_180mK_570uA/metrics.csv`

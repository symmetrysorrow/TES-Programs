# Historical versus current circuit semantics

Both SIFs select `TES Inner Circuit Update = True` and publish the heat source as
`TES Parallel Power / TES Volume`.  The common active law is:

`R = R0 * (1 + alpha*(T-T0)/T0 - beta) + R0*beta*abs(I)/I0`, clipped only at `Rmin`;
`I*(Rsh + R) = Ibias*Rsh` in steady state; `P = I^2 R`.

The current native implementation samples an element-equal nodal TES temperature and deliberately uses the completed prior assembly sweep (`CircuitSweepTemperature`) on the next call.  It owns power relaxation and checkpoint restore, then writes `TES Parallel Power`.  The body-force UDF only divides that scalar by `TES Volume`.

The historical waveform artifact is a transient launched with `--reuse-known-serial-steady`; its archived summary proves waveform agreement but does **not** prove that its 143.777852 uA seed was a full nonlinear steady fixed point.  Therefore it cannot be used as evidence for a distinct exact 143-uA solution until the historical state/SIF/binary is recovered and its residual is measured.

No difference in the affine TES resistance law, nominal circuit constants, or volumetric normalization was found.  The confirmed semantic difference to audit further is the coupling authority/timing (external historical provenance is incomplete versus current native HeatSolve hook), not a different circuit equation.

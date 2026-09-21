# Windows build/runtime audit and native outer capture

## Scope

No production physics, material, geometry, boundary-condition, timestep, or
solver-configuration source was changed.  The only source edit is the opt-in
HeatSolve diagnostic capture hook; it saves an already assembled outer system
immediately before and after `DefaultSolve` and restores the caller's ordinary
`Linear System Save Prefix` afterwards.

## Build failure

- The original full verbose CMake/Ninja command and its empty stderr are saved
  in `cmake_build_verbose.log` and `cmake_build_exit.txt`.
- The checked-in `test_heat.f90` is zero bytes, so it is not a valid minimal
  compiler test.  A real `hello.f90` fails in the inherited PowerShell
  environment with the same silent exit 1.
- Direct `f951.exe --version` exits `-1073741515` (`STATUS_DLL_NOT_FOUND`).
  `f951_ldd_stdout.txt` shows that `C:\Program Files\Elmer 26.1-Release\bin`
  supplied incompatible `libgcc_s_seh-1.dll` / `libwinpthread-1.dll`, while
  UCRT GCC dependencies (`libisl-23.dll`, `libmpc-3.dll`, `libmpfr-6.dll`) were
  unresolved.  The `gfortran` driver itself can print a version even when its
  `f951` child cannot load.
- With `MSYSTEM=UCRT64` and `/ucrt64/bin:/usr/bin` ahead of inherited PATH,
  `gfortran -cpp -E hello.f90` succeeds (`msys_hello_*`).  The HeatSolve target
  then builds with exact verbose commands in `cmake_build_ucrt64_stdout.log`.
- Cache/compiler/include/flag/PATH evidence is saved as
  `cmake_cache_selected.txt`, `environment.txt`, `command_resolution.txt`, and
  `ninja_relevant_rules.txt`.

## Runtime isolation

The first clean-UCRT run still segfaulted because inherited `ELMER_HOME` was
`C:\Program Files\Elmer 26.1-Release`.  GDB proves it loaded the old
`HeatSolve.dll` from that tree next to the fresh `libelmersolver.dll`; see
`capture_off_gdb_stdout.log`.  This is an ABI mismatch, independent of capture
or model physics.  The launcher used for the successful run sets:

```sh
MSYSTEM=UCRT64
ELMER_HOME=/d/Github/TES-Programs/tools/elmer-hypre/install-outer-capture-v2
ELMER_LIB=$ELMER_HOME/share/elmersolver/lib
PATH=$ELMER_HOME/bin:/ucrt64/bin:/usr/bin:/c/Program Files/CMake/bin:/c/Windows/System32
```

`capture_off_cleanruntime_stdout.log` reaches nonlinear iteration 2 normally
with capture disabled.  The successful capture uses the same environment and
is logged in `capture_on_v2_stdout.log`.

## Verified outer systems

`run_v2/ts0001_nl0001`, `nl0002`, and `nl0003` each contain the exact outer
`outer_before_{a,b,sol,sizes}.dat` and
`outer_after_{a,b,sol,sizes}.dat` files.  `metadata.json` was generated for
each directory by `phase24_outer_capture.py materialize`; all have 96,769
rows and 1,187,615 matrix records.

The run is the pulse-off, BDF1, one-step MUMPS restart diagnostic and retains
its original SIF physics.  It was intentionally stopped after capture nl3;
iteration 4 writes only the SIF's ordinary `fallback_*` files and cannot
overwrite a capture directory.

## Direct checks after capture

- NL2 independent SciPy LU residual: `6.359666272857715e-11` relative.
- NL2 direct minus saved outer-before maximum: `1.2359227999133027 mK`.
- NL2 direct minus saved outer-after maximum: `1.246956043716524 mK`.
- NL2 native outer-after minus outer-before maximum: `0.026100000000001122 mK`.
- NL1 -> NL2 crossed sensitivity is operator-dominated: max `0.0315806 mK`
  operator-only versus `0.00001188 mK` RHS-only.
- NL2 -> NL3 crossed sensitivity is RHS-dominated: max `0.0101499 mK`
  RHS-only versus `0.0004561 mK` operator-only.

Raw JSON is retained in `direct_nl2_nl3_stdout.json`,
`direct_nl2_native_after_stdout.json`, and `crossed_sensitivity_stdout.json`.

## Native residual / post-solve audit

`native_outer_system_audit.json` performs no new solver run.  It reads the
captured before/after files, independently solves each before system, and uses
the exact `TESInnerCircuitUpdate` aggregation: each TES (body 101) element is
given equal weight and each of its nodes has weight `1 / node_count`.  The
resulting weights cover 409 of 96,769 DOFs and sum exactly to one.

For nl1/nl2/nl3 respectively, the native outer-after residual relative to the
saved RHS is `2.0632e-2`, `1.9362e-2`, and `1.5446e-2`.  The residual is wholly
primal: this captured outer system has 96,769 nonzero-diagonal primal rows and
zero algebraic constraint rows.  The corresponding direct-vs-native null
residual `A (x_native - x_direct)` has the same relative magnitudes, so the
difference is not a near-nullspace/ill-conditioning ambiguity.

Raw SHA-256 is identical for before/after A and b at each nonlinear iteration.
Thus DefaultSolve did not change the saved assembled system; the nonzero
residual and solution discrepancy arise in the native solve/solution path.
The matched TES averages are:

| iteration | x_before (K) | x_direct (K) | x_after (K) | after - before (mK) |
| --- | ---: | ---: | ---: | ---: |
| nl1 | 0.168569126 | 0.166569454 | 0.167800483 | -0.768643 |
| nl2 | 0.167800483 | 0.166592884 | 0.167810210 | +0.009727 |
| nl3 | 0.167810210 | 0.166603131 | 0.167814144 | +0.003934 |

The `-0.768643 mK` observed jump is therefore the native post-DefaultSolve
field transfer/update for nl1; the same A/b's direct solution calls for a
larger `-1.999672 mK` TES-average correction.  Further investigation should
stay on the DefaultSolve/MUMPS solution-to-variable path (including nonlinear
relaxation and permutation/copy handling), rather than assembly or a
near-nullspace explanation.

## Mortar comparison status

This newly captured run is the mortar-on reference.  No new conformal
no-mortar physics run was launched as part of this Windows build/capture
isolation.  The existing physically conformal shared-node no-mortar short
campaign remains available at `artifacts/phase24_gate4_5_nomortar/short` and
reports a normal MUMPS completion; it is not an identical-mesh replacement for
the present 96,769-row mortar reference.  A fresh matched conformal control is
therefore still required before making a mortar-causality claim.

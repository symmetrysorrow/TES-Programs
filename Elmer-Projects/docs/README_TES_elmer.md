# TES shunt thermal-equilibrium prototype

This directory now contains a first Elmer prototype for the mesh in `gmsh/project.msh`.

## Reproducible Layout

The repository should now be interpreted with three roles:

- Source of truth:
  - `elmer_project.json`
- Generated artifacts:
  - `generated/*.sif`
  - `case_*.sif`
- Frozen reproducible runs:
  - `runs/<run_name>/`
  - `runs/<run_name>.zip`

The main rule is:

- edit `elmer_project.json` (numeric fields are self-healing: `sync_elmer_parameters.py`
  and `generate_project_geometry.py` always recompute them from the paired
  `*_expr`/`expression` fields via `scripts/support/reconcile_project.py` and write
  the result back, so editing an expression is enough)
- do not hand-edit `generated/*.sif` inside a reproducible workflow
- freeze confirmed states into `runs/`

To freeze the current confirmed state into a reproducible bundle:

```powershell
python freeze_repro_run.py current_reference
```

This creates:

- `runs/current_reference/`
- `runs/current_reference.zip`

Each frozen run contains copied inputs, copied generated files, copied mesh, copied results, and a `manifest.json`.

## Mesh IDs

The active Gmsh mesh is `gmsh/project.msh`; the original centered mesh is kept
as `gmsh/project_centered.msh`. The active Elmer mesh is in
`mesh_shifted_merged/`. The important physical IDs are:

- `100`: `abs`
- `101`: `TES`
- `102`: `Stycast`
- `103`: `Membrane_SiNx`
- `104`: `SiO2_1`
- `105`: `Si_1`
- `106`: `SiNx`
- `107`: `Si_2`
- `108`: `SiO2_2`
- `109`: `Membrane_Si1`

The TES volume computed from the Gmsh mesh is approximately `4.0e-14 m^3`.

## Cases

- `case_constant_power.sif`: steady heat solve using the steady-state TES/shunt circuit to compute TES Joule heating from temperature.
- `case_tes_shunt_internal.sif`: first internal Elmer/MATC version of the TES + shunt model.
- `case_tes_shunt_transient.sif`: transient Elmer case with the TES-branch
  inductance and the temperature/current-dependent TES resistance.
- `case_tes_pulse_1332kev.sif`: 1332 keV absorber-center pulse response,
  restarted from the converged electrothermal state.

The internal prototype currently uses the steady-state limit of the TES/shunt circuit:

```text
I_TES = I_bias * R_sh / (R_sh + R_TES)
P_TES = I_TES^2 * R_TES
R_TES(T,I) = max(R_min, R0 * (1 + alpha * (T_TES - T0) / T0 + beta * (I_TES - I0) / I0))
```

The current-dependent resistance is solved in `tes_heat_source.f90`. It evaluates the positive root of the TES/shunt circuit equation and applies the `R_min` branch from the JSON resistance law.
The TES-branch inductor from `tes_test2.json` is `L1 = 1.23e-8 H`. The
transient user function `tes_transient_heat_source.f90` advances the branch
current inside Elmer using backward Euler:

```text
L_TES * (I_n - I_(n-1)) / dt
  + I_n * (R_sh + R_TES(T_(n-1), I_n)) = I_bias * R_sh
```

The resulting `I_n^2 R_TES` is applied as the TES volumetric heat source. The
coupling is implicit within each timestep: the circuit is re-solved on every
nonlinear heat iteration (detected via `GetNonlinIter()`) using the running
TES nodal-average temperature, and the Joule power fixed point is converged
with a trust-region-capped Aitken (secant) relaxation. This is required for
stability: with the old once-per-step staggered update the explicit
electrothermal feedback loop has an amplification `|dP/dT|*dt / (C + G*dt)`
of order 10 at `dt = 1 ms`, which produced a growing period-2 oscillation
saturating at the transition edges. The transient cases therefore need
`Nonlinear System Max Iterations` of ~15 (steps converge in 2-4 iterations
once the relaxation factor has adapted). The CSV series file name can be set
per case with `TES Series File = String "..."` in the `Constants` block; one
row is written per converged timestep (the final timestep is not committed).

The supplied transient case is a short smoke test: `dt = 1.0e-7 s` and 10
steps. Increase `n_steps` near the top of the SIF to follow the thermal response
for longer times. The tested 10-step run completed successfully; TES current
rose from `23.0495 uA` at step 1 to `201.088 uA` at step 10.

The 1332 keV case deposits `2.1341e-13 J` at the Pb absorber center using a
3D Gaussian with `sigma = 50 um` during the first timestep.
`tes_pulse_series.csv` contains time, TES temperature, TES current, resistance,
and Joule power. The 100-ms simulation uses five timestep stages:

```text
0-20 us       dt = 0.5 us
20-200 us     dt = 2 us
0.2-2 ms      dt = 20 us
2-10 ms       dt = 0.1 ms
10-100 ms     dt = 0.5 ms
```

The completed 480-step run reached `0.1 s`. The TES temperature peaked at
`168.648618 mK` at `0.36 ms`, while the TES current reached its minimum of
`142.801170 uA` at the same time. At `100 ms`, the temperature was
`168.461491 mK` and the current was `149.403052 uA`, effectively returning to
the initial electrothermal state.

The parameters and material constants are taken from `D:\Github\Thermal-and-Electoric-Sim\tes_test2.json`. The `abs` body is assigned to the `Pb` material.

The membrane proxy has been replaced by two explicit center-island bodies:
`Membrane_Si1` (`body 109`) and `Membrane_SiNx` (`body 103`). They reuse the
bulk `Si` and `SiNx` material properties so the suspended island now carries
its own layer-resolved thermal resistance instead of a single effective
membrane conductivity.

The SiO2, Si, and SiNx substrate mother boxes are shifted by `-1 mm` in Y.
The membrane, TES, Stycast, Pb absorber, and all substrate subtraction boxes
(the cavity) remain centered at `Y = 0`. The source Gmsh mesh contains
separately meshed bodies. `mesh_shifted_merged/` merges coincident nodes so the
substrate layers and the explicit Si/SiNx membrane-island bodies form one
connected component. The remaining TES-Membrane_SiNx, TES-Stycast, and
Stycast-Pb interfaces are coupled with Elmer mortar boundary conditions.

Stycast is a Z-axis cylinder centered at `(0, 0)`, with radius `249 um` and
height `20 um`. Because the upstream geometry package supports only boxes and
wedges directly, `generate_project_geometry.py` adds cylinder support while
preserving the existing physical IDs.

## Mesh Modes

There are currently two distinct Elmer workflows in this repository.

### 1. Merged mesh

This is the stable baseline workflow used for the main steady-state runs:

```powershell
python generate_project_geometry.py
ElmerGrid 14 2 gmsh\project.msh -merge 1e-10 -out mesh_shifted_merged
ElmerSolver case_constant_power.sif
```

This merges coincident nodes and therefore collapses some nominal mortar
interfaces into conforming ones.

### 2. Unmerged mesh with minimal mortar

This is the current experimental workflow for keeping the mesh pieces separate
while only coupling the Z-stack interfaces with mortar:

```powershell
python generate_project_geometry.py
ElmerGrid 14 2 gmsh\project.msh -out mesh_shifted_unmerged
ElmerSolver case_constant_power_unmerged_minimal_mortar.sif
```

The current reference case is:

- `case_constant_power_unmerged_minimal_mortar.sif`

This case intentionally omits the cavity sidewall mortar pairs. During
debugging, those sidewall pairs were found to destabilize the solve much more
strongly than the Z-stack-only configuration.

At the moment, the unmerged minimal-mortar case should be treated as an
investigation scaffold rather than a validated production case.

## Run

```powershell
python generate_project_geometry.py
ElmerGrid 14 2 gmsh\project.msh -merge 1e-10 -out mesh_shifted_merged
elmerf90 tes_heat_source.f90 -o tes_heat_source_t0.dll
elmerf90 tes_transient_heat_source.f90 -o tes_transient_heat_source_t0.dll
ElmerSolver case_constant_power.sif
ElmerSolver case_tes_shunt_internal.sif
ElmerSolver case_tes_shunt_transient.sif
ElmerSolver case_tes_pulse_1332kev.sif
```

Open the resulting `.vtu` files with ParaView.

To print TES-only temperature summaries:

```powershell
python scripts\analysis\summarize_tes_temperature.py
```

## Temperature Reporting

TES temperature should be reported with the following priority:

- Official representative temperature: TES volume average
- Secondary check value: TES nodal average
- Do not use the old EP-only `tes` element-name average as a representative value

Why:

- The volume average is weighted by TES element volume, so it is the most natural
  representative temperature for thermal energy balance and for coupling back
  into the TES circuit model.
- The nodal average is still useful as a quick check, but it weights each node
  equally and therefore depends slightly on the mesh discretization.
- In the current merged mesh, the TES is nearly isothermal, so the nodal and
  volume averages are almost identical.

The analysis scripts now follow this convention:

- `scripts/analysis/extract_tes_volume_avg.py`: prints both TES nodal average
  and TES volume average from the VTU result
- `scripts/analysis/extract_tes_avg_ep.py`: reads temperature from the EP file,
  but uses `mesh.elements` body ID `101` to identify TES nodes safely
- `scripts/analysis/summarize_tes_temperature.py`: treats TES volume average as
  the reference temperature for `T - T0`, `R_TES`, `I_TES`, and `P_TES`

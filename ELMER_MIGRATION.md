# Elmer code migration

Elmer simulation code, projects and the PoST fitting workflow have moved to
the sibling repository `../Elmer-TES-Runner`.

Run from `D:/Github/Elmer-TES-Runner`:

```powershell
./Run_DualTESFast.ps1 -Action Prepare
./Run_DualTESFast.ps1 -Action Run
```

The entry point is `Fitting_LSE_Elmer.py`, with the sequence in
`config/dualtes_fast_steps.json` and instructions in
`docs/DualTESFast_fitting.md`. The original `PoST_Simulations/Fitting_LSE.py`
remains here; the Elmer fitter uses the same fixed experimental data location.

The retired `Elmer-Projects` directory, including its local generated meshes,
results and uncommitted snapshots, is preserved outside Git under
`D:/Github/Elmer-TES-Runner/archives/tes-programs-elmer-projects-20261014`.
The destination's existing newer source and project changes were preserved.
Historical provenance records referring to Elmer-Projects retain their original
names; active diagnostic configuration paths point to Elmer-TES-Runner.

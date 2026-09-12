# Archived PoST noise solution: R_SH = 3.85 mOhm

This directory archives the accepted 2024-12-06, 215 mK, 1400 uA target-case
noise solution.  The fixed TES resistance was derived from the same-campaign
IV data using the 3.85 mOhm shunt branch.

Archived solution parameters:

- `R_SH = 3.85 mOhm`
- `R_TES = 17.5518326 mOhm`
- `alpha = 102.0962702`
- `beta = 2.2638623`
- `L = 2.1065141 uH`
- `T_bath = 215 mK`
- `hardware_bessel_cutoff = 179527.6505 Hz`
- `post_filter_white_fraction = 0.0079696321`

## Reproduce generation and comparison

From the repository root:

```powershell
python PoST_Simulations\cases\tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2\noise_solution_rsh_3p85\run_noise_comparison.py
```

The default measured file is the target CH0 `modelnoise.txt` on the original
G: drive.  On another machine, pass its location explicitly:

```powershell
python PoST_Simulations\cases\tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2\noise_solution_rsh_3p85\run_noise_comparison.py `
  --measured-noise G:\\path\\to\\modelnoise.txt
```

The script copies the archived input into a temporary run directory, invokes
`PoST_Simulations/PoST_Simulation.py --noise-only`, copies the generated
`noise_total-bessel100k.dat` and its spectrum PNG into this archive, and
creates the full positive-frequency comparison plot.  The temporary
`noise.h5` is removed afterwards.

The optimizer used to obtain this fixed solution is
`PoST_Simulations/Opt_noise.py`; the simulator and its shared model helpers
remain in their normal repository locations.

The fit score uses 1 kHz--200 kHz.  The plot displays the full positive
frequency range; DC is omitted because both axes are logarithmic.

The archived standalone comparison gives a 1 kHz--200 kHz mean-squared
log10-residual score of 0.009343 and a median model/measured ratio of 0.9862.

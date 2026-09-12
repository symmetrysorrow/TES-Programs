# Archived PoST noise solution v2: R_SH = 3.85 mOhm

This directory archives the latest near-perfect 2024-12-06, 215 mK, 1400 uA
target-case noise solution.  The fixed TES resistance was derived from the
same-campaign IV data using the 3.85 mOhm shunt branch.

Archived solution parameters:

- `R_SH = 3.85 mOhm`
- `R_TES = 17.5518326 mOhm`
- `alpha = 18.2440912`
- `beta = 0`
- `L = 2.8116853 uH`
- `T_bath = 215 mK`
- `T_c = 241.2532358 mK`
- `hardware_bessel_cutoff = 100000 Hz`
- `hardware_bessel_order = 4`
- `post_filter_white_fraction = 0.0177971816`

## Reproduce generation and comparison

From the repository root:

```powershell
python PoST_Simulations\cases\tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2\noise_solution_rsh_3p85_v2\run_noise_comparison.py
```

The default measured file is the target CH0 `modelnoise.txt` on the original
G: drive.  On another machine, pass its location explicitly:

```powershell
python PoST_Simulations\cases\tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2\noise_solution_rsh_3p85_v2\run_noise_comparison.py `
  --measured-noise G:\\path\\to\\modelnoise.txt
```

The script copies the archived input into a temporary run directory, invokes
`PoST_Simulations/PoST_Simulation.py --noise-only`, copies the generated
`noise_total-bessel100k.dat` and its spectrum PNG into this archive, and
creates the full positive-frequency comparison plot.  The temporary
`noise.h5` is removed afterwards.

The optimizer used to obtain this fixed solution is
`PoST_Simulations/Opt_noise.py`; the simulator and shared model helpers remain
in their normal repository locations.

The archived comparison uses the 1 kHz--200 kHz band.  Its validation score is
0.001303 and the model/measured ratio median is 1.0158, with a 5--95%
interval of 0.9256--1.0931.  The plot displays the full positive frequency
range; DC is omitted because both axes are logarithmic.

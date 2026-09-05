# Noise-blind target-like sensitivity

This is a deterministic sensitivity study around `median_of_pulse_consistent_scenarios`. It is not a target parameter estimate and uses no experimental spectrum.
A separate generic 0.5x–2x reference sensitivity is stored in `reference_sensitivity_summary.json` and is not used for the physical envelope.

| parameter | 10-100 Hz | 100-1000 Hz | 1-3 kHz | 3-10 kHz |
|---|---:|---:|---:|---:|
| T_c | 1.042 | 0.5142 | 0.2879 | 0.2954 |
| T_bath | 1.041 | 0.5138 | 0.2879 | 0.2952 |
| R | 0.08794 | 0.04341 | 0.04673 | 0.1064 |
| R_l | 0.08795 | 0.04341 | 0.04673 | 0.1064 |
| alpha | 0.1868 | 0.09206 | 0.01143 | 0.2614 |
| beta | 0.1226 | 0.06046 | 0.02486 | 0.05399 |
| L | 1.554e-05 | 7.707e-06 | 8.589e-06 | 2.017e-05 |
| n | 0.0513 | 0.02532 | 0.01434 | 0.01442 |
| C_tes | 0.3524 | 0.1746 | 0.07938 | 0.1301 |
| C_abs | 0.001308 | 3.952e-05 | 2.757e-05 | 6.212e-05 |
| G_tes-bath | 0.3119 | 0.1547 | 0.0737 | 0.1073 |
| G_abs-tes | 0.04079 | 0.01942 | 0.00553 | 0.02566 |
| G_abs-abs | 0.001099 | 0.0005234 | 0.000149 | 0.0006916 |

Interpretation: values are finite-difference shape sensitivities around the generic nominal reference. They do not establish target provenance.

# Phase24 explicit mortar saddle sensitivity

- System: `96955 = 96769 + 186`; TES DOFs `409`
- Schur singular range: `4.75157176e-11` .. `8.85227990e-09`
- Schur 2-norm condition number: `1.86302140e+02`
- NumPy numerical rank: `186/186`
- Schur symmetry relative difference: `1.04613790e-13`
- Actual direct-minus-restart TES correction: `-0.76865053 mK`

## Weak-mode reconstruction

- top 1: L2 fraction `0.00142303`, TES signed fraction `-1.07525e-06`, TES reconstructed `-0.00000083 mK`
- top 3: L2 fraction `0.00540788`, TES signed fraction `3.6619e-07`, TES reconstructed `0.00000028 mK`
- top 5: L2 fraction `0.0058012`, TES signed fraction `-5.30238e-08`, TES reconstructed `-0.00000004 mK`
- top 10: L2 fraction `0.00645712`, TES signed fraction `-4.4012e-06`, TES reconstructed `-0.00000338 mK`

## Backward error

- restart augmented primal max eta: `4.07204542e-03`; constraint max eta: `5.14791142e-05`
- direct primal max eta: `1.98673762e-08`; constraint max eta: `1.80362475e-07`
- native primal max eta: `1.27723465e-08`; constraint max eta: `4.76435964e-08`

Interpretation is deferred until the Schur spectrum, weak-mode TES projection, and row classifications are considered together.

## Decision

- Root-cause assessment: the weakest explicit-mortar Schur modes do not explain the 0.76865053 mK direct-minus-restart correction.
- First production mitigation: constraint elimination, preserving the constraint physics but avoiding the explicit saddle path.
- HYPRE scaling/tolerance: secondary numerical tuning; the 5e-8 Gate3 sweep changed the discrepancy by only 0.931 microK.
- Conformal/shared-node no-mortar: long-term structural alternative requiring mesh/interface validation.

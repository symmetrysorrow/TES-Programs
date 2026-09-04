# HYPRE Phase19 Schur diagnostic artifacts

These JSON files are the checkout-visible snapshots used by the Phase19
documentation.  They are independent one-rank CPU diagnostics on the saved
explicit mortar matrix:

- `exact_schur.json`: SuperLU K factorization and dense Schur correctness oracle.
- `diag_schur.json`: `diag(K)^-1` Schur approximation.
- `ilu_proxy.json`: sparse ILU K-action proxy.
- `mgr_best.json`: the targeted MGR best result copied from the handoff status.

The JSON values are not production acceptance claims.  The exact dense Schur
is an oracle, and the MGR result reached its iteration limit above the gate.
The matrix-free Elmer action validator is
`scripts/analysis/validate_matrix_free_schur.py`; it requires the four
`*_vN.dat`/`*_yN.dat` files emitted by an Elmer diagnostic run.
`matrix_free_validation_current.json` is the current four-vector SciPy
cross-backend sensitivity snapshot: vector parity passes, while its strict
different-backend action gate does not (`max action relative error =
5.0716313e-7`). This is not the implementation correctness gate.
`superlu_parity_cpu.json` separates `self_consistency`, `b_vs_bt_transpose`,
`block_equivalence`, and `monolithic_comparison`.  It records K backward
residuals, componentwise backward errors, stage parity, raw/canonical SHA-256
fingerprints, Frobenius differences, and Schur action perturbations.  The
diagnostic emits the corresponding `*_K|B|Bt|D.triplets` and
`*_vN|btN|kuN|bkuN|dvN|yN.dat` files under
`results/case_p19_hypre_block_schur_diag_cpu_time5us/`.  The raw block
fingerprints mismatch for K/B/Bt, but the separate classification is
`NUMERICALLY CLOSE`; only D is raw-exact.
The emitted-stage matvec self-check passes. `same_binary_parity.json` uses
`BlockSchurSuperLUOracle`, linked against Elmer's exact
`block_schur_superlu_solve`; same-binary K and Schur parity pass for 4/4
vectors. Its `cross_backend_comparison` preserves the SciPy diagnostic: the
Schur gate is 2/4 fail and the composite gate is 3/4 fail. The finite
stale-tail result is in `short_bt_regression.json`. Lower/full CPU baseline
execution is separately gated; GPU, MPI, and transient promotion remain out
of scope.
The lower baseline attempt is recorded as `INCOMPLETE` in
`lower_cpu_one_step.json`; no partial iterate is accepted and full CPU was
not started. `full_cpu_one_step.json` records that gate explicitly.

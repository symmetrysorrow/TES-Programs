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
`matrix_free_validation_current.json` is the current four-vector result and
is an explicit failure snapshot: vector parity passes, but the strict action
gate does not (`max action relative error = 5.0716313e-7`).
`superlu_parity_cpu.json` separates `self_consistency`, `b_vs_bt_transpose`,
`block_equivalence`, and `monolithic_comparison`.  It records K backward
residuals, componentwise backward errors, stage parity, raw/canonical SHA-256
fingerprints, Frobenius differences, and Schur action perturbations.  The
diagnostic emits the corresponding `*_K|B|Bt|D.triplets` and
`*_vN|btN|kuN|bkuN|dvN|yN.dat` files under
`results/case_p19_hypre_block_schur_diag_cpu_time5us/`.  The raw block
fingerprints mismatch for K/B/Bt, but the separate classification is
`NUMERICALLY CLOSE`; only D is raw-exact.
The emitted-stage matvec self-check passes, while the SciPy actual-block K
oracle remains strict-fail because the system SuperLU and SciPy SuperLU
backends differ.
Lower/full CPU, GPU, MPI, and transient promotion tests were intentionally not
run after this failed gate.

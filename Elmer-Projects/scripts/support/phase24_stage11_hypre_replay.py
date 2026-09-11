"""Replay an exact Phase24 SaveLinearSystem A/b dump through HYPRE 2.33.

This is a diagnostic-only ctypes harness.  It uses one Microsoft MPI rank and
the same FlexGMRES + BoomerAMG settings as the native Phase24 CPU path.
"""
from __future__ import annotations

import argparse
import ctypes as C
import hashlib
from pathlib import Path

import numpy as np


MPI_COMM_WORLD = 0x44000000
HYPRE_PARCSR = 5555
# HYPRE_config.h has HYPRE_BIGINT and HYPRE_MIXEDINT disabled for this build;
# HYPRE_BigInt therefore resolves to C int, matching the native wrapper's
# int* index buffers.
BIG = C.c_int
INT = C.c_int
REAL = C.c_double
PTR = C.c_void_p


def fn(lib: C.CDLL, name: str, args: list[object], result: object = INT):
    f = getattr(lib, name)
    f.argtypes = args
    f.restype = result
    return f


def check(code: int, name: str) -> None:
    if code != 0:
        raise RuntimeError(f"{name} returned HYPRE/MPI status {code}")


def load_dump(a_path: Path, b_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    a = np.loadtxt(a_path, dtype=np.float64)
    b = np.loadtxt(b_path, dtype=np.float64)
    if a.ndim != 2 or a.shape[1] != 3 or b.ndim != 2 or b.shape[1] != 2:
        raise ValueError("unexpected SaveLinearSystem dump shape")
    rows = a[:, 0].astype(np.int64) - 1
    cols = a[:, 1].astype(np.int64) - 1
    vals = np.ascontiguousarray(a[:, 2], dtype=np.float64)
    rhs_rows = b[:, 0].astype(np.int64) - 1
    rhs = np.zeros(87534, dtype=np.float64)
    rhs[rhs_rows] = b[:, 1]
    return np.ascontiguousarray(rows), np.ascontiguousarray(cols), vals, rhs


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("a", type=Path)
    p.add_argument("b", type=Path)
    p.add_argument("--max-iter", type=int, default=2000)
    p.add_argument("--csv", type=Path, default=Path("phase24_stage11_hypre_residual.csv"))
    p.add_argument("--scaling", choices=("none", "symmetric-jacobi", "row"), default="none")
    p.add_argument("--guess-npy", type=Path)
    p.add_argument("--solution-npy-out", type=Path)
    args = p.parse_args()
    n = 87534
    rows, cols, vals, rhs = load_dump(args.a, args.b)
    original_vals = vals.copy()
    original_rhs = rhs.copy()
    scaling = np.ones(n, dtype=np.float64)
    if args.scaling == "symmetric-jacobi":
        diagonal = np.zeros(n, dtype=np.float64)
        diagonal[rows[rows == cols]] = vals[rows == cols]
        if np.any(diagonal <= 0.0):
            raise ValueError("symmetric Jacobi scaling requires positive diagonal")
        scaling = 1.0 / np.sqrt(diagonal)
        vals = np.ascontiguousarray(vals * scaling[rows] * scaling[cols])
        rhs = np.ascontiguousarray(rhs * scaling)
    elif args.scaling == "row":
        row_norm = np.zeros(n, dtype=np.float64)
        np.add.at(row_norm, rows, np.abs(vals))
        if np.any(row_norm <= 0.0):
            raise ValueError("row scaling requires nonzero rows")
        scaling = 1.0 / row_norm
        vals = np.ascontiguousarray(vals * scaling[rows])
        rhs = np.ascontiguousarray(rhs * scaling)
    if rows.size == 0 or rows.min() < 0 or rows.max() >= n:
        raise ValueError("matrix row range is invalid")
    if cols.min() < 0 or cols.max() >= n:
        raise ValueError("matrix column range is invalid")
    if np.any(rows[1:] < rows[:-1]):
        raise ValueError("matrix dump is not row ordered")

    mpi = C.CDLL(r"C:\Windows\System32\msmpi.dll")
    hypre = C.CDLL(r"C:\msys64\ucrt64\bin\libHYPRE.dll")
    mpi.MPI_Init.argtypes = [PTR, PTR]
    mpi.MPI_Init.restype = INT
    mpi.MPI_Finalize.argtypes = []
    mpi.MPI_Finalize.restype = INT
    mpi.MPI_Comm_size = fn(mpi, "MPI_Comm_size", [INT, C.POINTER(INT)])

    mpi_init = mpi.MPI_Init(None, None)
    check(mpi_init, "MPI_Init")
    try:
        hypre_initialize = fn(hypre, "HYPRE_Initialize", [])
        hypre_finalize = fn(hypre, "HYPRE_Finalize", [])
        hypre_clear_errors = fn(hypre, "HYPRE_ClearAllErrors", [])
        check(hypre_initialize(), "HYPRE_Initialize")
        size = INT()
        check(mpi.MPI_Comm_size(MPI_COMM_WORLD, C.byref(size)), "MPI_Comm_size")
        if size.value != 1:
            raise RuntimeError("replay requires exactly one MPI rank")

        matrix_create = fn(hypre, "HYPRE_IJMatrixCreate", [INT, BIG, BIG, BIG, BIG, C.POINTER(PTR)])
        matrix_type = fn(hypre, "HYPRE_IJMatrixSetObjectType", [PTR, INT])
        matrix_init = fn(hypre, "HYPRE_IJMatrixInitialize", [PTR])
        matrix_set = fn(hypre, "HYPRE_IJMatrixSetValues", [PTR, INT, C.POINTER(INT), C.POINTER(BIG), C.POINTER(BIG), C.POINTER(REAL)])
        matrix_assemble = fn(hypre, "HYPRE_IJMatrixAssemble", [PTR])
        matrix_object = fn(hypre, "HYPRE_IJMatrixGetObject", [PTR, C.POINTER(PTR)])
        matrix_destroy = fn(hypre, "HYPRE_IJMatrixDestroy", [PTR])
        vector_create = fn(hypre, "HYPRE_IJVectorCreate", [INT, BIG, BIG, C.POINTER(PTR)])
        vector_type = fn(hypre, "HYPRE_IJVectorSetObjectType", [PTR, INT])
        vector_init = fn(hypre, "HYPRE_IJVectorInitialize", [PTR])
        vector_set = fn(hypre, "HYPRE_IJVectorSetValues", [PTR, INT, C.POINTER(BIG), C.POINTER(REAL)])
        vector_assemble = fn(hypre, "HYPRE_IJVectorAssemble", [PTR])
        vector_object = fn(hypre, "HYPRE_IJVectorGetObject", [PTR, C.POINTER(PTR)])
        vector_get = fn(hypre, "HYPRE_IJVectorGetValues", [PTR, INT, C.POINTER(BIG), C.POINTER(REAL)])
        vector_destroy = fn(hypre, "HYPRE_IJVectorDestroy", [PTR])

        matrix = PTR()
        check(matrix_create(MPI_COMM_WORLD, 0, n - 1, 0, n - 1, C.byref(matrix)), "HYPRE_IJMatrixCreate")
        check(matrix_type(matrix, HYPRE_PARCSR), "HYPRE_IJMatrixSetObjectType")
        check(matrix_init(matrix), "HYPRE_IJMatrixInitialize")
        start = 0
        while start < rows.size:
            stop = start + 1
            while stop < rows.size and rows[stop] == rows[start]:
                stop += 1
            count = INT(stop - start)
            row = BIG(int(rows[start]))
            row_cols = np.ascontiguousarray(cols[start:stop], dtype=np.int64)
            row_vals = np.ascontiguousarray(vals[start:stop], dtype=np.float64)
            check(matrix_set(matrix, 1, C.byref(count), C.byref(row), row_cols.ctypes.data_as(C.POINTER(BIG)), row_vals.ctypes.data_as(C.POINTER(REAL))), "HYPRE_IJMatrixSetValues")
            start = stop
        check(matrix_assemble(matrix), "HYPRE_IJMatrixAssemble")
        A = PTR()
        check(matrix_object(matrix, C.byref(A)), "HYPRE_IJMatrixGetObject")

        indices = np.arange(n, dtype=np.int64)
        rhs = np.ascontiguousarray(rhs, dtype=np.float64)
        x = np.zeros(n, dtype=np.float64) if args.guess_npy is None else np.ascontiguousarray(np.load(args.guess_npy), dtype=np.float64)
        if x.shape != (n,):
            raise ValueError("initial guess must contain exactly 87534 values")
        vector_b = PTR()
        vector_x = PTR()
        for name, vector in (("b", C.byref(vector_b)), ("x", C.byref(vector_x))):
            check(vector_create(MPI_COMM_WORLD, 0, n - 1, vector), f"HYPRE_IJVectorCreate({name})")
        check(vector_type(vector_b, HYPRE_PARCSR), "HYPRE_IJVectorSetObjectType(b)")
        check(vector_type(vector_x, HYPRE_PARCSR), "HYPRE_IJVectorSetObjectType(x)")
        check(vector_init(vector_b), "HYPRE_IJVectorInitialize(b)")
        check(vector_init(vector_x), "HYPRE_IJVectorInitialize(x)")
        check(vector_set(vector_b, n, indices.ctypes.data_as(C.POINTER(BIG)), rhs.ctypes.data_as(C.POINTER(REAL))), "HYPRE_IJVectorSetValues(b)")
        check(vector_set(vector_x, n, indices.ctypes.data_as(C.POINTER(BIG)), x.ctypes.data_as(C.POINTER(REAL))), "HYPRE_IJVectorSetValues(x)")
        check(vector_assemble(vector_b), "HYPRE_IJVectorAssemble(b)")
        check(vector_assemble(vector_x), "HYPRE_IJVectorAssemble(x)")
        b = PTR()
        xv = PTR()
        check(vector_object(vector_b, C.byref(b)), "HYPRE_IJVectorGetObject(b)")
        check(vector_object(vector_x, C.byref(xv)), "HYPRE_IJVectorGetObject(x)")

        flex_create = fn(hypre, "HYPRE_ParCSRFlexGMRESCreate", [INT, C.POINTER(PTR)])
        flex_destroy = fn(hypre, "HYPRE_ParCSRFlexGMRESDestroy", [PTR])
        flex_setup = fn(hypre, "HYPRE_ParCSRFlexGMRESSetup", [PTR, PTR, PTR, PTR])
        flex_solve = fn(hypre, "HYPRE_ParCSRFlexGMRESSolve", [PTR, PTR, PTR, PTR])
        flex_set = {
            key: fn(hypre, symbol, [PTR, INT if key in ("max_iter", "kdim", "logging", "print") else REAL])
            for key, symbol in {
                "max_iter": "HYPRE_FlexGMRESSetMaxIter",
                "kdim": "HYPRE_FlexGMRESSetKDim",
                "tol": "HYPRE_FlexGMRESSetTol",
                "absolute": "HYPRE_FlexGMRESSetAbsoluteTol",
                "logging": "HYPRE_FlexGMRESSetLogging",
                "print": "HYPRE_FlexGMRESSetPrintLevel",
            }.items()
        }
        flex_get_iter = fn(hypre, "HYPRE_ParCSRFlexGMRESGetNumIterations", [PTR, C.POINTER(INT)])
        flex_get_res = fn(hypre, "HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm", [PTR, C.POINTER(REAL)])
        flex_get_res_generic = fn(hypre, "HYPRE_FlexGMRESGetFinalRelativeResidualNorm", [PTR, C.POINTER(REAL)])
        flex_precond = fn(hypre, "HYPRE_FlexGMRESSetPrecond", [PTR, PTR, PTR, PTR])
        amg_create = fn(hypre, "HYPRE_BoomerAMGCreate", [C.POINTER(PTR)])
        amg_destroy = fn(hypre, "HYPRE_BoomerAMGDestroy", [PTR])
        amg_set = {
            key: fn(hypre, symbol, [PTR, INT if key not in ("tol", "strong") else REAL])
            for key, symbol in {
                "tol": "HYPRE_BoomerAMGSetTol",
                "max_iter": "HYPRE_BoomerAMGSetMaxIter",
                "relax": "HYPRE_BoomerAMGSetRelaxType",
                "coarsen": "HYPRE_BoomerAMGSetCoarsenType",
                "sweeps": "HYPRE_BoomerAMGSetNumSweeps",
                "levels": "HYPRE_BoomerAMGSetMaxLevels",
                "interp": "HYPRE_BoomerAMGSetInterpType",
                "smooth": "HYPRE_BoomerAMGSetSmoothType",
                "cycle": "HYPRE_BoomerAMGSetCycleType",
                "functions": "HYPRE_BoomerAMGSetNumFunctions",
                "strong": "HYPRE_BoomerAMGSetStrongThreshold",
            }.items()
        }
        amg_solve = C.cast(getattr(hypre, "HYPRE_BoomerAMGSolve"), PTR)
        amg_setup = C.cast(getattr(hypre, "HYPRE_BoomerAMGSetup"), PTR)
        solver = PTR()
        precond = PTR()
        check(flex_create(MPI_COMM_WORLD, C.byref(solver)), "HYPRE_ParCSRFlexGMRESCreate")
        check(amg_create(C.byref(precond)), "HYPRE_BoomerAMGCreate")
        check(flex_set["max_iter"](solver, args.max_iter), "HYPRE_FlexGMRESSetMaxIter")
        check(flex_set["kdim"](solver, 100), "HYPRE_FlexGMRESSetKDim")
        check(flex_set["tol"](solver, 5.0e-7), "HYPRE_FlexGMRESSetTol")
        check(flex_set["absolute"](solver, 0.0), "HYPRE_FlexGMRESSetAbsoluteTol")
        check(flex_set["logging"](solver, 1), "HYPRE_FlexGMRESSetLogging")
        check(flex_set["print"](solver, 3), "HYPRE_FlexGMRESSetPrintLevel")
        for key, value in (("tol", 0.0), ("max_iter", 1), ("relax", 18), ("coarsen", 8), ("sweeps", 1), ("levels", 25), ("interp", 6), ("smooth", 0), ("cycle", 1), ("functions", 1), ("strong", 0.25)):
            check(amg_set[key](precond, value), f"HYPRE_BoomerAMGSet{key}")
        check(flex_precond(solver, amg_solve, amg_setup, precond), "HYPRE_FlexGMRESSetPrecond")
        check(flex_setup(solver, A, b, xv), "HYPRE_ParCSRFlexGMRESSetup")
        status = flex_solve(solver, A, b, xv)
        iterations = INT()
        residual = REAL()
        generic_residual = REAL()
        generic_residual_before_clear = REAL()
        generic_residual_before_status = flex_get_res_generic(solver, C.byref(generic_residual_before_clear))
        hypre_clear_errors()
        iter_status = flex_get_iter(solver, C.byref(iterations))
        residual_status = flex_get_res(solver, C.byref(residual))
        generic_residual_status = flex_get_res_generic(solver, C.byref(generic_residual))
        if generic_residual_before_clear.value > 0.0:
            residual.value = generic_residual_before_clear.value
        elif generic_residual_status == 0:
            residual.value = generic_residual.value
        check(vector_get(vector_x, n, indices.ctypes.data_as(C.POINTER(BIG)), x.ctypes.data_as(C.POINTER(REAL))), "HYPRE_IJVectorGetValues(x)")
        if args.scaling == "symmetric-jacobi":
            x_original = scaling * x
        else:
            x_original = x
        original_residual = np.zeros(n, dtype=np.float64)
        np.add.at(original_residual, rows, original_vals * x_original[cols])
        original_residual -= original_rhs
        original_relative_residual = float(np.linalg.norm(original_residual) / np.linalg.norm(original_rhs))
        args.csv.write_text("iteration,relative_residual\n0,1.0\n" f"{iterations.value},{residual.value:.17g}\n", encoding="utf-8")
        solution_hash = hashlib.sha256(x.tobytes()).hexdigest()
        if args.solution_npy_out:
            np.save(args.solution_npy_out, np.ascontiguousarray(x))
        print(f"replay dimensions={n} nnz={rows.size} scaling={args.scaling} rhs_l2={np.linalg.norm(rhs):.17g}")
        print(f"replay status={status} result={'PASS' if status == 0 and residual.value <= 5.0e-7 else 'FAIL'} iterations={iterations.value} final_relative_residual={residual.value:.17g} original_relative_residual={original_relative_residual:.17g} tolerance=5e-7 getter_statuses={iter_status},{residual_status},{generic_residual_before_status},{generic_residual_status}")
        print(f"replay solution_sha256={solution_hash} residual_csv={args.csv}")
        flex_destroy(solver)
        amg_destroy(precond)
        vector_destroy(vector_x); vector_destroy(vector_b); matrix_destroy(matrix)
        result = 0 if status == 0 and residual.value <= 5.0e-7 else 1
        check(hypre_finalize(), "HYPRE_Finalize")
        return result
    finally:
        mpi.MPI_Finalize()


if __name__ == "__main__":
    raise SystemExit(main())

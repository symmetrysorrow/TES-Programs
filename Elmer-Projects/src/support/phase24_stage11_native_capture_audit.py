"""Audit and convert a Phase24 native HYPRE capture.

The native wrapper writes ``PHASE24_CSR_V1`` as a portable, little-endian
single-rank CSR file.  This script converts that exact A/b/x capture to
NumPy artifacts and computes independent residuals without referring to the
canonical SaveLinearSystem dump.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import struct
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix


HEADER = struct.Struct("<16siiqii")
MAGIC = b"PHASE24_CSR_V1"


def read_csr(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[int, ...]]:
    raw = path.read_bytes()
    magic, nrows, ncols, nnz, row_lower, row_upper = HEADER.unpack_from(raw)
    if magic.rstrip(b"\0") != MAGIC:
        raise ValueError(f"unexpected CSR magic: {magic!r}")
    offset = HEADER.size
    row_ids = np.frombuffer(raw, dtype="<i4", count=nrows, offset=offset).copy()
    offset += 4 * nrows
    indptr = np.frombuffer(raw, dtype="<i8", count=nrows + 1, offset=offset).copy()
    offset += 8 * (nrows + 1)
    indices = np.frombuffer(raw, dtype="<i4", count=nnz, offset=offset).copy()
    offset += 4 * nnz
    data = np.frombuffer(raw, dtype="<f8", count=nnz, offset=offset).copy()
    if offset + 8 * nnz != len(raw):
        raise ValueError("CSR file has trailing or truncated data")
    if indptr[-1] != nnz or np.unique(row_ids).size != nrows:
        raise ValueError("CSR row metadata is inconsistent")
    return row_ids, indptr, indices, data, (row_lower, row_upper)


def read_vector(path: Path, n: int) -> np.ndarray:
    records = np.loadtxt(path, dtype=np.float64)
    if records.ndim != 2 or records.shape[1] != 2:
        raise ValueError(f"unexpected vector shape in {path}: {records.shape}")
    vector = np.full(n, np.nan, dtype=np.float64)
    vector[records[:, 0].astype(np.int64)] = records[:, 1]
    if np.isnan(vector).any() or np.isinf(vector).any():
        raise ValueError(f"vector is not a complete finite capture: {path}")
    return vector


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("capture_dir", type=Path)
    args = parser.parse_args()
    d = args.capture_dir
    a_path = d / "native_exact_A_rank0000.csrbin"
    row_ids, indptr, indices, data, row_range = read_csr(a_path)
    n = row_ids.size
    b = read_vector(d / "native_hypre_b_rank0000.dat", n)
    x_before = read_vector(d / "native_hypre_x_before_rank0000.dat", n)
    x_after = read_vector(d / "native_hypre_x_after_rank0000.dat", n)
    matrix = csr_matrix((data, indices, indptr), shape=(n, n))

    def relative_residual(x: np.ndarray) -> float:
        return float(np.linalg.norm(b - matrix @ x) / np.linalg.norm(b))

    np.save(d / "native_exact_indices.npy", row_ids)
    np.save(d / "native_exact_b.npy", b)
    np.save(d / "native_exact_x_before.npy", x_before)
    np.save(d / "native_exact_x_after.npy", x_after)
    np.savez_compressed(
        d / "native_exact_A.npz",
        indptr=indptr,
        indices=indices,
        data=data,
        shape=np.array([n, n], dtype=np.int64),
        row_ids=row_ids,
    )
    metadata = d / "native_hypre_metadata_rank0000.txt"
    metadata_text = metadata.read_text(encoding="utf-8")
    metadata_values: dict[str, str] = {}
    for line in metadata_text.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            metadata_values[key] = value
    result = {
        "capture_dir": str(d),
        "A_native": {
            "artifact": str(d / "native_exact_A.npz"),
            "binary_artifact": str(a_path),
            "rows": int(n),
            "nnz": int(data.size),
            "row_range": list(map(int, row_range)),
            "sha256_binary": hashlib.sha256(a_path.read_bytes()).hexdigest().upper(),
        },
        "b_native": {
            "artifact": str(d / "native_exact_b.npy"),
            "l2": float(np.linalg.norm(b)),
            "max_abs": float(np.max(np.abs(b))),
        },
        "x_before": {
            "artifact": str(d / "native_exact_x_before.npy"),
            "l2": float(np.linalg.norm(x_before)),
            "min": float(np.min(x_before)),
            "max": float(np.max(x_before)),
            "r_before": relative_residual(x_before),
        },
        "x_after": {
            "artifact": str(d / "native_exact_x_after.npy"),
            "l2": float(np.linalg.norm(x_after)),
            "min": float(np.min(x_after)),
            "max": float(np.max(x_after)),
            "r_after": relative_residual(x_after),
        },
        "native_metadata": metadata_values,
    }
    (d / "native_exact_audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

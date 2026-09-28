"""Export the h8 matrices as raw CSR (int32 n, int32 nnz, int32 indptr, int32 indices, float64 data)."""
import numpy as np, scipy.sparse as sp

def write(path, M):
    M = sp.csr_matrix(M)
    M.sort_indices()
    with open(path, "wb") as f:
        np.array([M.shape[0], M.nnz], dtype=np.int32).tofile(f)
        M.indptr.astype(np.int32).tofile(f)
        M.indices.astype(np.int32).tofile(f)
        M.data.astype(np.float64).tofile(f)
    print(path, M.shape[0], M.nnz)

A = sp.load_npz("/root/linsys/A.npz").tocsr()
A = ((A + A.T) * 0.5).tocsr()
Al = sp.load_npz("/root/linsys/A_lumped.npz").tocsr()
write("/root/linsys/A.csr", A)
write("/root/linsys/Al.csr", Al)
rng = np.random.default_rng(3)
b = A @ rng.standard_normal(A.shape[0])
b.astype(np.float64).tofile("/root/linsys/b.bin")
print("b", b.size)

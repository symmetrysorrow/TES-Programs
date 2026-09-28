/* PCG + BoomerAMG probe on exported CSR matrices (CPU or HIP HYPRE).
 *
 *   hypre_probe A.csr Al.csr b.bin  [key=value ...]
 *   keys: gpu=0/1 relax=8 coarsen=8 interp=6 agg=0 aggi=4 pmax=4 tol=3e-11
 *         strong=0.25 reps=3 lap=0 (lap=N: use an N^3 7-point Laplacian instead)
 *
 * Rows are split contiguously over MPI ranks.  BoomerAMG is set up and
 * applied on Al (preconditioning matrix), PCG runs on A, as in Elmer.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "_hypre_utilities.h"

static HYPRE_ParCSRMatrix g_tilde;
static HYPRE_Int amg_solve(HYPRE_Solver p, HYPRE_ParCSRMatrix A, HYPRE_ParVector b, HYPRE_ParVector x)
{ (void)A; return HYPRE_BoomerAMGSolve(p, g_tilde, b, x); }
static HYPRE_Int noop(HYPRE_Solver p, HYPRE_ParCSRMatrix A, HYPRE_ParVector b, HYPRE_ParVector x)
{ (void)p; (void)A; (void)b; (void)x; return 0; }

typedef struct { int n, nnz; int *ip, *ix; double *v; } csr_t;

static void read_csr(const char *path, csr_t *m)
{
  FILE *f = fopen(path, "rb");
  if (!f) { perror(path); MPI_Abort(MPI_COMM_WORLD, 1); }
  int h[2]; fread(h, sizeof(int), 2, f);
  m->n = h[0]; m->nnz = h[1];
  m->ip = malloc(sizeof(int) * (m->n + 1)); m->ix = malloc(sizeof(int) * m->nnz);
  m->v = malloc(sizeof(double) * m->nnz);
  fread(m->ip, sizeof(int), m->n + 1, f); fread(m->ix, sizeof(int), m->nnz, f);
  fread(m->v, sizeof(double), m->nnz, f); fclose(f);
}

static void laplace(int N, csr_t *m)
{
  int n = N * N * N, k = 0;
  m->n = n; m->nnz = 7 * n;
  m->ip = malloc(sizeof(int) * (n + 1)); m->ix = malloc(sizeof(int) * m->nnz); m->v = malloc(sizeof(double) * m->nnz);
  for (int r = 0; r < n; r++) {
    int i = r % N, j = (r / N) % N, l = r / (N * N);
    m->ip[r] = k;
    int nb[6][3] = {{i-1,j,l},{i+1,j,l},{i,j-1,l},{i,j+1,l},{i,j,l-1},{i,j,l+1}};
    m->ix[k] = r; m->v[k++] = 6.0;
    for (int q = 0; q < 6; q++) {
      int a = nb[q][0], b = nb[q][1], c = nb[q][2];
      if (a < 0 || b < 0 || c < 0 || a >= N || b >= N || c >= N) continue;
      m->ix[k] = a + N * (b + N * c); m->v[k++] = -1.0;
    }
  }
  m->ip[n] = k; m->nnz = k;
}

static HYPRE_IJMatrix build_ij(csr_t *m, int lo, int hi, int gpu)
{
  HYPRE_IJMatrix A;
  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, lo, hi, lo, hi, &A);
  HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize_v2(A, HYPRE_MEMORY_HOST);
  for (int r = lo; r <= hi; r++) {
    int nc = m->ip[r + 1] - m->ip[r];
    HYPRE_IJMatrixSetValues(A, 1, &nc, &r, &m->ix[m->ip[r]], &m->v[m->ip[r]]);
  }
  HYPRE_IJMatrixAssemble(A);
  if (gpu) HYPRE_IJMatrixMigrate(A, HYPRE_MEMORY_DEVICE);
  return A;
}

static HYPRE_IJVector build_vec(double *vals, int lo, int hi, int gpu)
{
  HYPRE_IJVector v;
  int n = hi - lo + 1;
  int *idx = malloc(sizeof(int) * n);
  for (int i = 0; i < n; i++) idx[i] = lo + i;
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, lo, hi, &v);
  HYPRE_IJVectorSetObjectType(v, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize_v2(v, HYPRE_MEMORY_HOST);
  HYPRE_IJVectorSetValues(v, n, idx, vals + lo);
  HYPRE_IJVectorAssemble(v);
  if (gpu) HYPRE_IJVectorMigrate(v, HYPRE_MEMORY_DEVICE);
  free(idx);
  return v;
}

static double kv(int argc, char **argv, const char *key, double def)
{
  size_t L = strlen(key);
  for (int i = 4; i < argc; i++)
    if (!strncmp(argv[i], key, L) && argv[i][L] == '=') return atof(argv[i] + L + 1);
  return def;
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);
  int gpu = (int)kv(argc, argv, "gpu", 0), relax = (int)kv(argc, argv, "relax", 8);
  int coarsen = (int)kv(argc, argv, "coarsen", 8), interp = (int)kv(argc, argv, "interp", 6);
  int agg = (int)kv(argc, argv, "agg", 0), aggi = (int)kv(argc, argv, "aggi", 4);
  int pmax = (int)kv(argc, argv, "pmax", 4), reps = (int)kv(argc, argv, "reps", 3);
  int lap = (int)kv(argc, argv, "lap", 0);
  double tol = kv(argc, argv, "tol", 3e-11), strong = kv(argc, argv, "strong", 0.25);

  HYPRE_Initialize();
  if (gpu) {
    HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
    HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
    HYPRE_SetSpGemmUseVendor(0);
    HYPRE_SetUseGpuRand(1);
  }

  csr_t A, Al; double *b;
  double t0 = MPI_Wtime();
  if (lap > 0) { laplace(lap, &A); Al = A; }
  else { read_csr(argv[1], &A); read_csr(argv[2], &Al); }
  b = malloc(sizeof(double) * A.n);
  if (lap > 0) { for (int i = 0; i < A.n; i++) b[i] = 1.0; }
  else { FILE *f = fopen(argv[3], "rb"); fread(b, sizeof(double), A.n, f); fclose(f); }
  int lo = (int)((long)A.n * rank / size), hi = (int)((long)A.n * (rank + 1) / size) - 1;
  double *zero = calloc(A.n, sizeof(double));

  HYPRE_IJMatrix ijA = build_ij(&A, lo, hi, gpu);
  HYPRE_IJMatrix ijT = (lap > 0) ? ijA : build_ij(&Al, lo, hi, gpu);
  HYPRE_IJVector ijb = build_vec(b, lo, hi, gpu), ijx = build_vec(zero, lo, hi, gpu);
  HYPRE_ParCSRMatrix pA, pT; HYPRE_ParVector pb, px;
  HYPRE_IJMatrixGetObject(ijA, (void **)&pA); HYPRE_IJMatrixGetObject(ijT, (void **)&pT);
  HYPRE_IJVectorGetObject(ijb, (void **)&pb); HYPRE_IJVectorGetObject(ijx, (void **)&px);
  MPI_Barrier(MPI_COMM_WORLD);
  double t_build = MPI_Wtime() - t0;
  if (rank == 0) printf("n=%d nnzA=%d nnzT=%d ranks=%d gpu=%d build=%.2fs\n", A.n, A.nnz, Al.nnz, size, gpu, t_build);
  fflush(stdout);

  HYPRE_Solver amg, pcg;
  HYPRE_BoomerAMGCreate(&amg);
  HYPRE_BoomerAMGSetPrintLevel(amg, 1);
  HYPRE_BoomerAMGSetMaxIter(amg, 1); HYPRE_BoomerAMGSetTol(amg, 0.0);
  HYPRE_BoomerAMGSetRelaxType(amg, relax); HYPRE_BoomerAMGSetNumSweeps(amg, 1);
  HYPRE_BoomerAMGSetCoarsenType(amg, coarsen); HYPRE_BoomerAMGSetInterpType(amg, interp);
  HYPRE_BoomerAMGSetStrongThreshold(amg, strong); HYPRE_BoomerAMGSetPMaxElmts(amg, pmax);
  HYPRE_BoomerAMGSetAggNumLevels(amg, agg); HYPRE_BoomerAMGSetAggInterpType(amg, aggi);
  if (gpu) HYPRE_BoomerAMGSetKeepTranspose(amg, 1);

  for (int rep = 0; rep < reps; rep++) {
    HYPRE_ParVectorSetConstantValues(px, 0.0);
    MPI_Barrier(MPI_COMM_WORLD);
    double s0 = MPI_Wtime();
    HYPRE_BoomerAMGSetup(amg, pT, pb, px);
    g_tilde = pT;
    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &pcg);
    HYPRE_PCGSetTol(pcg, tol); HYPRE_PCGSetMaxIter(pcg, 2000); HYPRE_PCGSetTwoNorm(pcg, 0);
    HYPRE_PCGSetPrecond(pcg, (HYPRE_PtrToSolverFcn)amg_solve, (HYPRE_PtrToSolverFcn)noop, amg);
    HYPRE_ParCSRPCGSetup(pcg, pA, pb, px);
    MPI_Barrier(MPI_COMM_WORLD);
    double s1 = MPI_Wtime();
    HYPRE_ParCSRPCGSolve(pcg, pA, pb, px);
    MPI_Barrier(MPI_COMM_WORLD);
    double s2 = MPI_Wtime();
    int its; double res;
    HYPRE_PCGGetNumIterations(pcg, &its); HYPRE_PCGGetFinalRelativeResidualNorm(pcg, &res);
    if (rank == 0)
      printf("rep=%d setup=%.3fs solve=%.3fs its=%d per_it=%.2fms rel=%.2e\n",
             rep, s1 - s0, s2 - s1, its, 1e3 * (s2 - s1) / (its > 0 ? its : 1), res);
    fflush(stdout);
    HYPRE_ParCSRPCGDestroy(pcg);
  }
  HYPRE_BoomerAMGDestroy(amg);
  HYPRE_Finalize();
  MPI_Finalize();
  return 0;
}

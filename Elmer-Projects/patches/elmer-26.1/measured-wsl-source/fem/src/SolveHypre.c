/*****************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! ******************************************************************************
! *
! *  Elmer interface for Hypre - High performance Preconditioners
! *
! *  For more information on Hypre see
! *  https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Thomas Zwinger, Jonas Thies, Peter Råback, Mika Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2000
! *
! *****************************************************************************/


#include "../config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif

/*
 * Optional process-memory telemetry for the long-running Phase24 campaign.
 * It is disabled unless PHASE24_MEMORY_TRACE is set to 1/true, so normal
 * production runs do not incur a process query or extra output.  The public
 * entry point is also a no-op on non-Windows builds to keep the Fortran
 * ListMatrix module portable.  Stages 1-5 bracket List_ToCRSMatrix, 10-12
 * bracket initial HYPRE construction/setup, 20-21 bracket a persistent
 * refresh, 22-26 split its IJ update and Krylov setup phases, and 30-31
 * bracket HYPRE destruction.
 */
static int phase24_memory_trace_enabled = -1;

/* Optional native Krylov trace for a focused convergence diagnostic.  This is
 * deliberately environment-controlled so production runs retain their
 * normal output level and cost. */
static int phase24_hypre_residual_trace_enabled = -1;

/* A refreshed ParCSR matrix can require a fresh outer Krylov setup even when
 * the BoomerAMG hierarchy is intentionally retained.  Keep this diagnostic
 * opt-in because the normal persistent path avoids a known HYPRE workspace
 * growth issue. */
static int phase24_force_krylov_setup = -1;
static int phase24_force_amg_setup = -1;

static int phase24_force_krylov_setup_is_enabled(void)
{
  const char *value;
  if (phase24_force_krylov_setup >= 0)
    return phase24_force_krylov_setup;
  value = getenv("PHASE24_HYPRE_FORCE_KRYLOV_SETUP");
  phase24_force_krylov_setup =
    (value != NULL &&
     (value[0] == '1' || value[0] == 'y' || value[0] == 'Y' ||
      value[0] == 't' || value[0] == 'T')) ? 1 : 0;
  return phase24_force_krylov_setup;
}

static int phase24_force_amg_setup_is_enabled(void)
{
  const char *value;
  if (phase24_force_amg_setup >= 0)
    return phase24_force_amg_setup;
  value = getenv("PHASE24_HYPRE_FORCE_AMG_SETUP");
  phase24_force_amg_setup =
    (value != NULL &&
     (value[0] == '1' || value[0] == 'y' || value[0] == 'Y' ||
      value[0] == 't' || value[0] == 'T')) ? 1 : 0;
  return phase24_force_amg_setup;
}

static int phase24_hypre_residual_trace_is_enabled(void)
{
  const char *value;
  if (phase24_hypre_residual_trace_enabled >= 0)
    return phase24_hypre_residual_trace_enabled;
  value = getenv("PHASE24_HYPRE_RESIDUAL_TRACE");
  phase24_hypre_residual_trace_enabled =
    (value != NULL &&
     (value[0] == '1' || value[0] == 'y' || value[0] == 'Y' ||
      value[0] == 't' || value[0] == 'T')) ? 1 : 0;
  return phase24_hypre_residual_trace_enabled;
}

static int phase24_memory_trace_is_enabled(void)
{
  const char *value;
  if (phase24_memory_trace_enabled >= 0)
    return phase24_memory_trace_enabled;
  value = getenv("PHASE24_MEMORY_TRACE");
  phase24_memory_trace_enabled =
    (value != NULL &&
     (value[0] == '1' || value[0] == 'y' || value[0] == 'Y' ||
      value[0] == 't' || value[0] == 'T')) ? 1 : 0;
  return phase24_memory_trace_enabled;
}

static void phase24_memory_trace_event(int stage, int nrows, int nnz,
                                       int epoch, const void *container)
{
  if (!phase24_memory_trace_is_enabled())
    return;

#ifdef _WIN32
  {
    PROCESS_MEMORY_COUNTERS_EX counters;
    memset(&counters, 0, sizeof(counters));
    if (GetProcessMemoryInfo(GetCurrentProcess(),
        (PROCESS_MEMORY_COUNTERS *)&counters, sizeof(counters))) {
      fprintf(stderr,
        "PHASE24_MEMORY pid=%lu stage=%d rows=%d nnz=%d epoch=%d container=%p "
        "private_bytes=%llu working_set=%llu peak_working_set=%llu "
        "pagefile_bytes=%llu\n",
        (unsigned long)GetCurrentProcessId(), stage, nrows, nnz, epoch,
        container,
        (unsigned long long)counters.PrivateUsage,
        (unsigned long long)counters.WorkingSetSize,
        (unsigned long long)counters.PeakWorkingSetSize,
        (unsigned long long)counters.PagefileUsage);
    } else {
      fprintf(stderr,
        "PHASE24_MEMORY pid=%lu stage=%d rows=%d nnz=%d epoch=%d container=%p "
        "query_failed=1 error=%lu\n",
        (unsigned long)GetCurrentProcessId(), stage, nrows, nnz, epoch,
        container,
        (unsigned long)GetLastError());
    }
    fflush(stderr);
  }
#else
  (void)stage;
  (void)nrows;
  (void)nnz;
  (void)epoch;
  (void)container;
#endif
}

/* BIND(C) entry point used by ListMatrix.F90. */
void phase24_memory_trace(int *stage, int *nrows, int *nnz, int *epoch)
{
  if (stage == NULL || nrows == NULL || nnz == NULL || epoch == NULL)
    return;
  phase24_memory_trace_event(*stage, *nrows, *nnz, *epoch, NULL);
}

#ifdef HAVE_HYPRE
#include <math.h>
#include <stdint.h>
#include "HYPRE_config.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "_hypre_parcsr_mv.h"

/* Opt-in (PHASE24_HYPRE_MIGRATE=1): with the HOST memory location the HIP
   build executes on the CPU.  Migrate the assembled ParCSR matrix and the
   ParVectors to device memory so setup and Krylov run on the GPU; vectors
   are moved back to the host around every host-side IJ access. */
static int phase24_migrate_enabled(void)
{
  static int enabled = -1;
  if (enabled < 0) {
    const char *e = getenv("PHASE24_HYPRE_MIGRATE");
    enabled = (e != NULL && *e != '\0' && *e != '0');
  }
  return enabled;
}
#include "HYPRE_parcsr_ls.h"

#define realtime_ FC_FUNC_(realtime,REALTIME)
typedef struct ElmerHypreContainer ElmerHypreContainer;

/* BoomerAMG built on a separate preconditioning matrix: the V-cycle must use
   that matrix on every level, whatever matrix the Krylov method passes. */
static HYPRE_ParCSRMatrix phase24_tilde_matrix = NULL;

static HYPRE_Int phase24_tilde_amg_solve(HYPRE_Solver precond, HYPRE_ParCSRMatrix A,
                                         HYPRE_ParVector b, HYPRE_ParVector x)
{
  (void)A;
  return HYPRE_BoomerAMGSolve(precond, phase24_tilde_matrix, b, x);
}

static HYPRE_Int phase24_noop_setup(HYPRE_Solver precond, HYPRE_ParCSRMatrix A,
                                    HYPRE_ParVector b, HYPRE_ParVector x)
{
  (void)precond; (void)A; (void)b; (void)x;
  return 0;
}

/* Zero the owned part of an initialized IJ vector before AddToValues. */
static void phase24_zero_ij_vector(HYPRE_IJVector v)
{
  HYPRE_ParVector pv;
  HYPRE_IJVectorGetObject(v, (void **) &pv);
  HYPRE_ParVectorSetConstantValues(pv, 0.0);
}

typedef struct {
  ElmerHypreContainer *container;
  HYPRE_Solver boomer;
} Phase24BoomerAMGContext;

struct ElmerHypreContainer {

  int ilower, iupper;

  HYPRE_IJMatrix A;
  HYPRE_IJMatrix Atilde;
  HYPRE_ParCSRMatrix parcsr_A;

  HYPRE_IJVector x,b;
  
  int hypre_method;
  int hypre_use_gpu;
  HYPRE_Solver solver, precond;
  HYPRE_Solver mgr_fsolver;
  int *mgr_marker;

  /* Reusable host-side scratch buffers for the persistent IJ path. */
  int *index_buf;
  int index_capacity;
  double *vec_buf;
  int vec_capacity;
  double *rhs_buf;
  int *rhs_index_buf;
  int rhs_capacity;
  int rhs_count;
  int structure_epoch;
  int matrix_epoch;
  int preconditioner_epoch;
  /* The numeric setup epoch is distinct from the matrix epoch.  A solve is
     allowed to enter the Krylov method directly only when these match. */
  int setup_epoch;
  int vectors_initialized;
  int phase24_keep_alive;

  /* Native-call counters for the persistent Phase24 path.  These are kept at
     the wrapper boundary so lifecycle labels cannot hide an actual setup. */
  long long ij_matrix_initialize_calls;
  long long ij_matrix_set_values_calls;
  long long ij_matrix_assemble_calls;
  long long ij_matrix_get_object_calls;
  long long ij_vector_initialize_calls;
  long long ij_vector_set_values_calls;
  long long ij_vector_assemble_calls;
  long long ij_vector_get_object_calls;
  long long matrix_refresh_calls;
  long long krylov_setup_calls;
  long long amg_setup_calls;
  long long solve_calls;
  long long case_a_direct_solve_calls;
  long long case_b_reused_solve_calls;
  long long case_b_fresh_solve_calls;
  long long preconditioner_age_forced_rebuild_calls;
  long long degradation_trigger_rebuild_calls;
  long long recovery_trigger_rebuild_calls;
  long long recovery_setup_calls;
  int last_krylov_iterations;
  int setup_since_last_solve;
  int last_solved_matrix_epoch;
  int recovery_solve_pending;
  int preconditioner_setup_requested;
  int gpu_migration_failed;
  Phase24BoomerAMGContext *phase24_boomer_context;

  /* Explicit device-residency diagnostics for the Phase24 GPU path. */
  long long gpu_matrix_h2d_refresh_calls;
  long long gpu_rhs_h2d_calls;
  long long gpu_solution_h2d_calls;
  long long gpu_solution_d2h_calls;
  long long gpu_full_object_migration_calls;
  long long gpu_device_allocation_calls;

  /* AMS specific stuff */
  HYPRE_IJMatrix G, Pi;

};

static const char *Phase24HypreBackendName(void)
{
#if defined(HYPRE_USING_CUDA)
  return "CUDA";
#elif defined(HYPRE_USING_HIP)
  return "HIP";
#else
  return "CPU";
#endif
}

/* Diagnostic-only native vector provenance.  The selector is deliberately
   narrow and deterministic for the current 31 us failing state, while the
   files retain global HYPRE indices so the capture is not tied to Elmer's
   SaveLinearSystem numbering. */
static int phase24_capture_attempted = 0;
static int phase24_capture_active = 0;
static int phase24_capture_selector_logged = 0;
static char phase24_capture_dir[512];

typedef struct {
  int nan_count;
  int inf_count;
  int unchanged_count;
  double min_value;
  double max_value;
  double l2_norm;
} Phase24CaptureStats;

static Phase24CaptureStats phase24_x_before_stats;
static Phase24CaptureStats phase24_b_stats;
static Phase24CaptureStats phase24_x_after_ij_stats;
static Phase24CaptureStats phase24_x_after_par_stats;
static Phase24CaptureStats phase24_x_after_selected_stats;
static HYPRE_Int phase24_error_before = 0;
static HYPRE_Int phase24_global_error_before = 0;
static HYPRE_Int phase24_conv_before = 0;
static HYPRE_Int phase24_clear_conv_status = 0;
static HYPRE_Int phase24_error_after_clear = 0;
static HYPRE_Int phase24_global_error_after_clear = 0;
static HYPRE_Int phase24_conv_after_clear = 0;
static HYPRE_Int phase24_ij_after_status = -1;
static HYPRE_Int phase24_par_after_status = -1;
static HYPRE_Int phase24_matrix_capture_status = -1;
static int64_t phase24_matrix_capture_nnz = 0;

#define PHASE24_CAPTURE_SENTINEL (-9.87654321e307)

static Phase24CaptureStats Phase24CaptureMeasure(int count,
                                                 const double *values)
{
  Phase24CaptureStats stats;
  int i;
  stats.nan_count = 0;
  stats.inf_count = 0;
  stats.unchanged_count = 0;
  stats.min_value = HUGE_VAL;
  stats.max_value = -HUGE_VAL;
  stats.l2_norm = 0.0;
  for (i = 0; i < count; i++) {
    double value = values[i];
    if (value == PHASE24_CAPTURE_SENTINEL) stats.unchanged_count++;
    if (isnan(value)) {
      stats.nan_count++;
      continue;
    }
    if (isinf(value)) {
      stats.inf_count++;
      continue;
    }
    if (value < stats.min_value) stats.min_value = value;
    if (value > stats.max_value) stats.max_value = value;
    stats.l2_norm += value * value;
  }
  stats.l2_norm = sqrt(stats.l2_norm);
  if (stats.min_value == HUGE_VAL) stats.min_value = 0.0;
  if (stats.max_value == -HUGE_VAL) stats.max_value = 0.0;
  return stats;
}

static void Phase24CaptureLogStats(const char *label, int status, int count,
                                   const Phase24CaptureStats *stats)
{
  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_VECTOR label=%s status=%d count=%d "
          "nan=%d inf=%d unchanged=%d min=%.17e max=%.17e l2=%.17e\n",
          label, status, count, stats->nan_count, stats->inf_count,
          stats->unchanged_count, stats->min_value, stats->max_value,
          stats->l2_norm);
  fflush(stderr);
}

static void Phase24CaptureLog(const char *message)
{
  fprintf(stderr, "%s\n", message);
  fflush(stderr);
}

static int Phase24CaptureWriteIndexed(const char *suffix, int count,
                                      const int *indices,
                                      const double *values, int rank)
{
  char path[768];
  FILE *fp;
  int i;

  snprintf(path, sizeof(path), "%s/native_hypre_%s_rank%04d.dat",
           phase24_capture_dir, suffix, rank);
  fprintf(stderr, "NATIVE_XVEC_CAPTURE_OPEN label=%s path=%s count=%d\n",
          suffix, path, count);
  fflush(stderr);
  fp = fopen(path, "w");
  if (fp == NULL) {
    fprintf(stderr, "NATIVE_XVEC_CAPTURE_OPEN_FAILURE label=%s path=%s\n",
            suffix, path);
    fflush(stderr);
    return 1;
  }
  for (i = 0; i < count; i++)
    fprintf(fp, "%d %.17e\n", indices[i], values[i]);
  fclose(fp);
  fprintf(stderr, "NATIVE_XVEC_CAPTURE_DUMP_COMPLETED label=%s path=%s\n",
          suffix, path);
  fflush(stderr);
  return 0;
}

static int Phase24CaptureWriteIndices(int count, const int *indices, int rank)
{
  char path[768];
  FILE *fp;
  int i;

  snprintf(path, sizeof(path), "%s/native_hypre_indices_rank%04d.dat",
           phase24_capture_dir, rank);
  fprintf(stderr, "NATIVE_XVEC_CAPTURE_OPEN label=indices path=%s count=%d\n",
          path, count);
  fflush(stderr);
  fp = fopen(path, "w");
  if (fp == NULL) {
    fprintf(stderr, "NATIVE_XVEC_CAPTURE_OPEN_FAILURE label=indices path=%s\n",
            path);
    fflush(stderr);
    return 1;
  }
  for (i = 0; i < count; i++) fprintf(fp, "%d\n", indices[i]);
  fclose(fp);
  fprintf(stderr, "NATIVE_XVEC_CAPTURE_DUMP_COMPLETED label=indices path=%s\n",
          path);
  fflush(stderr);
  return 0;
}

/* Save the exact IJ matrix presented to the ParCSR solve.  The file is a
   small, self-describing binary CSR container: 16-byte magic, int32 rows /
   cols, int64 nnz, int32 row range, then global row ids, uint64 row pointer,
   int32 global columns, and float64 values.  MPI-1 is the bounded target;
   retaining global ids keeps the format unambiguous if this is later split
   into one file per rank. */
static int Phase24CaptureWriteMatrix(ElmerHypreContainer *Container,
                                     int count, int *row_ids, int rank)
{
  char path[768];
  FILE *fp;
  HYPRE_Int *row_counts;
  HYPRE_BigInt *cols;
  HYPRE_Complex *matrix_values;
  int64_t *row_ptr;
  int64_t nnz = 0;
  int i;
  HYPRE_Int status;
  int32_t nrows = (int32_t)count;
  int32_t ncols = (int32_t)(Container->iupper - Container->ilower + 1);
  int32_t row_lower = (int32_t)Container->ilower;
  int32_t row_upper = (int32_t)Container->iupper;
  char magic[16] = "PHASE24_CSR_V1";

  row_counts = (HYPRE_Int *)calloc((size_t)count, sizeof(HYPRE_Int));
  row_ptr = (int64_t *)calloc((size_t)count + 1, sizeof(int64_t));
  if (row_counts == NULL || row_ptr == NULL) {
    free(row_counts);
    free(row_ptr);
    phase24_matrix_capture_status = -1;
    return 1;
  }
  status = HYPRE_IJMatrixGetRowCounts(Container->A, count,
                                      (HYPRE_BigInt *)row_ids, row_counts);
  if (status != 0) {
    fprintf(stderr,
            "NATIVE_XVEC_CAPTURE_MATRIX_ROWCOUNTS_FAILURE status=%d\n",
            status);
    fflush(stderr);
    free(row_counts);
    free(row_ptr);
    phase24_matrix_capture_status = status;
    return status;
  }
  for (i = 0; i < count; i++) {
    nnz += (int64_t)row_counts[i];
    row_ptr[i + 1] = nnz;
  }
  cols = (HYPRE_BigInt *)malloc((size_t)nnz * sizeof(HYPRE_BigInt));
  matrix_values = (HYPRE_Complex *)malloc((size_t)nnz * sizeof(HYPRE_Complex));
  if (cols == NULL || matrix_values == NULL) {
    free(row_counts);
    free(row_ptr);
    free(cols);
    free(matrix_values);
    phase24_matrix_capture_status = -1;
    return 1;
  }
  status = HYPRE_IJMatrixGetValues(Container->A, count, row_counts,
                                   (HYPRE_BigInt *)row_ids, cols,
                                   matrix_values);
  if (status != 0) {
    fprintf(stderr, "NATIVE_XVEC_CAPTURE_MATRIX_VALUES_FAILURE status=%d\n",
            status);
    fflush(stderr);
    free(row_counts);
    free(row_ptr);
    free(cols);
    free(matrix_values);
    phase24_matrix_capture_status = status;
    return status;
  }

  snprintf(path, sizeof(path), "%s/native_exact_A_rank%04d.csrbin",
           phase24_capture_dir, rank);
  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_OPEN label=A_native path=%s rows=%d nnz=%lld\n",
          path, count, (long long)nnz);
  fflush(stderr);
  fp = fopen(path, "wb");
  if (fp == NULL) {
    fprintf(stderr, "NATIVE_XVEC_CAPTURE_OPEN_FAILURE label=A_native path=%s\n",
            path);
    fflush(stderr);
    free(row_counts);
    free(row_ptr);
    free(cols);
    free(matrix_values);
    phase24_matrix_capture_status = 1;
    return 1;
  }
  fwrite(magic, sizeof(magic), 1, fp);
  fwrite(&nrows, sizeof(nrows), 1, fp);
  fwrite(&ncols, sizeof(ncols), 1, fp);
  fwrite(&nnz, sizeof(nnz), 1, fp);
  fwrite(&row_lower, sizeof(row_lower), 1, fp);
  fwrite(&row_upper, sizeof(row_upper), 1, fp);
  fwrite(row_ids, sizeof(int32_t), (size_t)count, fp);
  fwrite(row_ptr, sizeof(int64_t), (size_t)count + 1, fp);
  fwrite(cols, sizeof(int32_t), (size_t)nnz, fp);
  fwrite(matrix_values, sizeof(double), (size_t)nnz, fp);
  if (fclose(fp) != 0) {
    phase24_matrix_capture_status = 1;
  } else {
    phase24_matrix_capture_status = 0;
  }
  phase24_matrix_capture_nnz = nnz;
  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_DUMP_COMPLETED label=A_native path=%s "
          "status=%d nnz=%lld\n",
          path, phase24_matrix_capture_status, (long long)nnz);
  fflush(stderr);
  free(row_counts);
  free(row_ptr);
  free(cols);
  free(matrix_values);
  return phase24_matrix_capture_status;
}

/* Capture the matrix object that is actually passed to the Krylov routine.
   Container->A is normally the same object, but the persistent lifecycle can
   temporarily retain a different ParCSR object after a refresh.  Reading
   through the public ParCSR row API makes that distinction explicit. */
static int Phase24CaptureWriteSolveMatrix(HYPRE_ParCSRMatrix matrix,
                                          ElmerHypreContainer *Container,
                                          int count, int *row_ids, int rank)
{
  char path[768];
  char magic[16] = "PHASE24_CSR_V1";
  FILE *fp;
  int64_t *row_ptr;
  HYPRE_BigInt *cols = NULL;
  HYPRE_Complex *matrix_values = NULL;
  int64_t capacity = 0;
  int64_t nnz = 0;
  int i;
  int32_t nrows = (int32_t)count;
  int32_t ncols = (int32_t)(Container->iupper - Container->ilower + 1);
  int32_t row_lower = (int32_t)Container->ilower;
  int32_t row_upper = (int32_t)Container->iupper;
  HYPRE_Int status;

  row_ptr = (int64_t *)calloc((size_t)count + 1, sizeof(int64_t));
  if (row_ptr == NULL) {
    phase24_matrix_capture_status = 1;
    return 1;
  }
  for (i = 0; i < count; i++) {
    HYPRE_BigInt row = (HYPRE_BigInt)row_ids[i];
    HYPRE_Int row_size = 0;
    HYPRE_BigInt *row_cols = NULL;
    HYPRE_Complex *row_values = NULL;
    status = HYPRE_ParCSRMatrixGetRow(matrix, row, &row_size,
                                      &row_cols, &row_values);
    if (status != 0 || row_size < 0) {
      fprintf(stderr,
              "NATIVE_XVEC_CAPTURE_SOLVE_MATRIX_ROW_FAILURE row=%d status=%d size=%d\n",
              row_ids[i], status, row_size);
      fflush(stderr);
      free(row_ptr);
      free(cols);
      free(matrix_values);
      phase24_matrix_capture_status = status != 0 ? status : 1;
      return phase24_matrix_capture_status;
    }
    if (nnz + row_size > capacity) {
      int64_t next_capacity = capacity > 0 ? capacity : 1024;
      while (next_capacity < nnz + row_size) next_capacity *= 2;
      cols = (HYPRE_BigInt *)realloc(cols,
                                     (size_t)next_capacity * sizeof(HYPRE_BigInt));
      matrix_values = (HYPRE_Complex *)realloc(
          matrix_values, (size_t)next_capacity * sizeof(HYPRE_Complex));
      if (cols == NULL || matrix_values == NULL) {
        HYPRE_ParCSRMatrixRestoreRow(matrix, row, &row_size, &row_cols,
                                     &row_values);
        free(row_ptr);
        free(cols);
        free(matrix_values);
        phase24_matrix_capture_status = 1;
        return 1;
      }
      capacity = next_capacity;
    }
    if (row_size > 0) {
      memcpy(cols + nnz, row_cols, (size_t)row_size * sizeof(HYPRE_BigInt));
      memcpy(matrix_values + nnz, row_values,
             (size_t)row_size * sizeof(HYPRE_Complex));
    }
    nnz += row_size;
    row_ptr[i + 1] = nnz;
    status = HYPRE_ParCSRMatrixRestoreRow(matrix, row, &row_size, &row_cols,
                                          &row_values);
    if (status != 0) {
      fprintf(stderr,
              "NATIVE_XVEC_CAPTURE_SOLVE_MATRIX_ROW_RESTORE_FAILURE row=%d status=%d\n",
              row_ids[i], status);
      fflush(stderr);
      free(row_ptr);
      free(cols);
      free(matrix_values);
      phase24_matrix_capture_status = status;
      return status;
    }
  }

  snprintf(path, sizeof(path), "%s/native_exact_A_rank%04d.csrbin",
           phase24_capture_dir, rank);
  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_OPEN label=A_native_solve_object path=%s rows=%d nnz=%lld\n",
          path, count, (long long)nnz);
  fflush(stderr);
  fp = fopen(path, "wb");
  if (fp == NULL) {
    fprintf(stderr,
            "NATIVE_XVEC_CAPTURE_OPEN_FAILURE label=A_native_solve_object path=%s\n",
            path);
    fflush(stderr);
    free(row_ptr);
    free(cols);
    free(matrix_values);
    phase24_matrix_capture_status = 1;
    return 1;
  }
  fwrite(magic, sizeof(magic), 1, fp);
  fwrite(&nrows, sizeof(nrows), 1, fp);
  fwrite(&ncols, sizeof(ncols), 1, fp);
  fwrite(&nnz, sizeof(nnz), 1, fp);
  fwrite(&row_lower, sizeof(row_lower), 1, fp);
  fwrite(&row_upper, sizeof(row_upper), 1, fp);
  fwrite(row_ids, sizeof(int32_t), (size_t)count, fp);
  fwrite(row_ptr, sizeof(int64_t), (size_t)count + 1, fp);
  fwrite(cols, sizeof(int32_t), (size_t)nnz, fp);
  fwrite(matrix_values, sizeof(double), (size_t)nnz, fp);
  phase24_matrix_capture_status = fclose(fp) == 0 ? 0 : 1;
  phase24_matrix_capture_nnz = nnz;
  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_DUMP_COMPLETED label=A_native_solve_object path=%s status=%d nnz=%lld\n",
          path, phase24_matrix_capture_status, (long long)nnz);
  fflush(stderr);
  free(row_ptr);
  free(cols);
  free(matrix_values);
  return phase24_matrix_capture_status;
}

static void Phase24CaptureBefore(ElmerHypreContainer *Container,
                                 int local_size, int *globaldofs,
                                 HYPRE_ParCSRMatrix parcsr_A,
                                 int rank, int rounds, double tol,
                                 int hypre_sol)
{
  const char *requested_dir = getenv("PHASE24_HYPRE_CAPTURE_DIR");
  int x_status, b_status;

  phase24_capture_attempted = 1;
  phase24_capture_active = 1;
  if (requested_dir != NULL && requested_dir[0] != '\0')
    snprintf(phase24_capture_dir, sizeof(phase24_capture_dir), "%s",
             requested_dir);
  else
    snprintf(phase24_capture_dir, sizeof(phase24_capture_dir),
             "results/case_phase24_adaptive_output_smoke");

  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_ENTER solver=%p method=%d hypre_sol=%d "
          "global_size=%d local_size=%d local_range=%d:%d matrix_epoch=%d "
          "preconditioner_epoch=%d rounds=%d tol=%.17e requested_dir=%s\n",
          (void *)Container->solver, Container->hypre_method, hypre_sol,
          Container->iupper - Container->ilower + 1, local_size,
          Container->ilower, Container->iupper, Container->matrix_epoch,
          Container->preconditioner_epoch, rounds, tol, phase24_capture_dir);
  fflush(stderr);

  Phase24CaptureWriteIndices(local_size, globaldofs, rank);
  phase24_matrix_capture_status = Phase24CaptureWriteSolveMatrix(
      parcsr_A, Container, local_size, globaldofs, rank);
  for (int i = 0; i < local_size; i++) {
    Container->vec_buf[i] = PHASE24_CAPTURE_SENTINEL;
    Container->rhs_buf[i] = PHASE24_CAPTURE_SENTINEL;
  }
  x_status = HYPRE_IJVectorGetValues(Container->x, local_size,
                                     globaldofs, Container->vec_buf);
  b_status = HYPRE_IJVectorGetValues(Container->b, local_size,
                                     globaldofs, Container->rhs_buf);
  phase24_x_before_stats = Phase24CaptureMeasure(local_size,
                                                  Container->vec_buf);
  phase24_b_stats = Phase24CaptureMeasure(local_size, Container->rhs_buf);
  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_GET_BEFORE x_status=%d b_status=%d\n",
          x_status, b_status);
  fflush(stderr);
  Phase24CaptureLogStats("x_before", x_status, local_size,
                         &phase24_x_before_stats);
  Phase24CaptureLogStats("b", b_status, local_size, &phase24_b_stats);
  if (x_status == 0)
    Phase24CaptureWriteIndexed("x_before", local_size, globaldofs,
                               Container->vec_buf, rank);
  if (b_status == 0)
    Phase24CaptureWriteIndexed("b", local_size, globaldofs,
                               Container->rhs_buf, rank);
}

static int Phase24CaptureAfter(ElmerHypreContainer *Container,
                               int count, int *indices, double *values,
                               HYPRE_ParVector par_x, MPI_Comm comm, int rank,
                               int solve_status, int iterations,
                               double final_residual)
{
  char path[768];
  FILE *fp;
  int i;
  const double sentinel = PHASE24_CAPTURE_SENTINEL;

  if (!phase24_capture_active)
    return HYPRE_IJVectorGetValues(Container->x, count, indices, values);

  phase24_error_before = HYPRE_GetError();
  phase24_global_error_before = HYPRE_GetGlobalError(comm);
  phase24_conv_before = HYPRE_CheckError(phase24_error_before,
                                         HYPRE_ERROR_CONV);
  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_ERROR_BEFORE local=%d global=%d conv=%d "
          "solve_status=%d iterations=%d final_relative_residual=%.17e\n",
          phase24_error_before, phase24_global_error_before,
          phase24_conv_before, solve_status, iterations, final_residual);
  fflush(stderr);

  phase24_clear_conv_status = HYPRE_ClearError(HYPRE_ERROR_CONV);
  phase24_error_after_clear = HYPRE_GetError();
  phase24_global_error_after_clear = HYPRE_GetGlobalError(comm);
  phase24_conv_after_clear = HYPRE_CheckError(phase24_error_after_clear,
                                              HYPRE_ERROR_CONV);
  fprintf(stderr,
          "NATIVE_XVEC_CAPTURE_ERROR_CLEAR_CONV status=%d local_after=%d "
          "global_after=%d conv_after=%d\n",
          phase24_clear_conv_status, phase24_error_after_clear,
          phase24_global_error_after_clear, phase24_conv_after_clear);
  fflush(stderr);

  for (i = 0; i < count; i++) values[i] = sentinel;
  phase24_ij_after_status = HYPRE_IJVectorGetValues(Container->x, count,
                                                     indices, values);
  phase24_x_after_ij_stats = Phase24CaptureMeasure(count, values);
  Phase24CaptureLogStats("x_after_ij", phase24_ij_after_status, count,
                         &phase24_x_after_ij_stats);

  if (phase24_ij_after_status == 0)
    memcpy(Container->rhs_buf, values, (size_t)count * sizeof(double));
  for (i = 0; i < count; i++) values[i] = sentinel;
  phase24_par_after_status = HYPRE_ParVectorGetValues(par_x, count, indices,
                                                       values);
  phase24_x_after_par_stats = Phase24CaptureMeasure(count, values);
  Phase24CaptureLogStats("x_after_parvector", phase24_par_after_status,
                         count, &phase24_x_after_par_stats);

  if (phase24_ij_after_status == 0 && phase24_par_after_status == 0) {
    double max_abs_diff = 0.0;
    double diff_norm = 0.0;
    for (i = 0; i < count; i++) {
      double diff = values[i] - Container->rhs_buf[i];
      double abs_diff = fabs(diff);
      if (abs_diff > max_abs_diff) max_abs_diff = abs_diff;
      diff_norm += diff * diff;
    }
    fprintf(stderr,
            "NATIVE_XVEC_CAPTURE_COMPARE ij_vs_parvector_max_abs=%.17e "
            "ij_vs_parvector_l2=%.17e\n",
            max_abs_diff, sqrt(diff_norm));
    fflush(stderr);
  }

  /* Prefer the IJ values when both APIs succeed; retain the ParVector values
     as the fallback because the latter bypasses IJ's error state. */
  if (phase24_ij_after_status == 0) {
    memcpy(values, Container->rhs_buf, (size_t)count * sizeof(double));
    phase24_x_after_selected_stats = phase24_x_after_ij_stats;
  } else if (phase24_par_after_status == 0) {
    phase24_x_after_selected_stats = phase24_x_after_par_stats;
  } else {
    phase24_x_after_selected_stats = Phase24CaptureMeasure(count, values);
  }
  {
    int selected_status = phase24_ij_after_status == 0 ? 0 :
                          phase24_par_after_status;
    fprintf(stderr,
            "NATIVE_XVEC_CAPTURE_GET_AFTER ij_status=%d par_status=%d "
            "selected_status=%d solve_status=%d iterations=%d "
            "final_relative_residual=%.17e\n",
            phase24_ij_after_status, phase24_par_after_status, selected_status,
            solve_status, iterations, final_residual);
    fflush(stderr);
    if (selected_status == 0)
      Phase24CaptureWriteIndexed("x_after", count, indices, values, rank);

    snprintf(path, sizeof(path), "%s/native_hypre_metadata_rank%04d.txt",
             phase24_capture_dir, rank);
    fprintf(stderr, "NATIVE_XVEC_CAPTURE_OPEN label=metadata path=%s\n", path);
    fflush(stderr);
    fp = fopen(path, "w");
    if (fp == NULL) {
      fprintf(stderr, "NATIVE_XVEC_CAPTURE_OPEN_FAILURE label=metadata path=%s\n",
              path);
      fflush(stderr);
    } else {
      fprintf(fp, "global_size=%d\n", Container->iupper - Container->ilower + 1);
      fprintf(fp, "local_size=%d\n", count);
      fprintf(fp, "local_range=%d %d\n", Container->ilower, Container->iupper);
      fprintf(fp, "matrix_epoch=%d\n", Container->matrix_epoch);
      fprintf(fp, "preconditioner_epoch=%d\n", Container->preconditioner_epoch);
      fprintf(fp, "matrix_capture_status=%d\n", phase24_matrix_capture_status);
      fprintf(fp, "matrix_capture_nnz=%lld\n",
              (long long)phase24_matrix_capture_nnz);
      fprintf(fp, "solve_status=%d\n", solve_status);
      fprintf(fp, "error_before=%d\n", phase24_error_before);
      fprintf(fp, "global_error_before=%d\n", phase24_global_error_before);
      fprintf(fp, "conv_before=%d\n", phase24_conv_before);
      fprintf(fp, "clear_conv_status=%d\n", phase24_clear_conv_status);
      fprintf(fp, "error_after_clear=%d\n", phase24_error_after_clear);
      fprintf(fp, "global_error_after_clear=%d\n", phase24_global_error_after_clear);
      fprintf(fp, "conv_after_clear=%d\n", phase24_conv_after_clear);
      fprintf(fp, "ij_after_status=%d\n", phase24_ij_after_status);
      fprintf(fp, "parvector_after_status=%d\n", phase24_par_after_status);
      fprintf(fp, "selected_after_status=%d\n", selected_status);
      fprintf(fp, "iterations=%d\n", iterations);
      fprintf(fp, "final_relative_residual=%.17e\n", final_residual);
      fprintf(fp, "solver_method=%d\n", Container->hypre_method);
      fprintf(fp, "x_before_nan=%d\n", phase24_x_before_stats.nan_count);
      fprintf(fp, "x_before_inf=%d\n", phase24_x_before_stats.inf_count);
      fprintf(fp, "x_before_unchanged=%d\n", phase24_x_before_stats.unchanged_count);
      fprintf(fp, "b_nan=%d\n", phase24_b_stats.nan_count);
      fprintf(fp, "b_inf=%d\n", phase24_b_stats.inf_count);
      fprintf(fp, "b_unchanged=%d\n", phase24_b_stats.unchanged_count);
      fprintf(fp, "x_after_ij_nan=%d\n", phase24_x_after_ij_stats.nan_count);
      fprintf(fp, "x_after_ij_inf=%d\n", phase24_x_after_ij_stats.inf_count);
      fprintf(fp, "x_after_ij_unchanged=%d\n", phase24_x_after_ij_stats.unchanged_count);
      fprintf(fp, "x_after_par_nan=%d\n", phase24_x_after_par_stats.nan_count);
      fprintf(fp, "x_after_par_inf=%d\n", phase24_x_after_par_stats.inf_count);
      fprintf(fp, "x_after_par_unchanged=%d\n", phase24_x_after_par_stats.unchanged_count);
      fprintf(fp, "x_after_selected_nan=%d\n", phase24_x_after_selected_stats.nan_count);
      fprintf(fp, "x_after_selected_inf=%d\n", phase24_x_after_selected_stats.inf_count);
      fprintf(fp, "x_after_selected_unchanged=%d\n", phase24_x_after_selected_stats.unchanged_count);
      fclose(fp);
      fprintf(stderr, "NATIVE_XVEC_CAPTURE_DUMP_COMPLETED label=metadata path=%s\n",
              path);
      fflush(stderr);
    }
    phase24_capture_active = 0;
    return selected_status;
  }
}

static HYPRE_Int Phase24BoomerAMGSetup(HYPRE_Solver context_solver,
                                       HYPRE_Matrix matrix,
                                       HYPRE_Vector b,
                                       HYPRE_Vector x)
{
  Phase24BoomerAMGContext *context = (Phase24BoomerAMGContext *)context_solver;
  ElmerHypreContainer *container;
  if (context == NULL || context->container == NULL || context->boomer == NULL)
    return -1;
  container = context->container;
  if (!container->preconditioner_setup_requested &&
      !phase24_force_amg_setup_is_enabled())
    return 0;
  {
    HYPRE_Int status = HYPRE_BoomerAMGSetup(
      context->boomer, (HYPRE_ParCSRMatrix)matrix,
      (HYPRE_ParVector)b, (HYPRE_ParVector)x);
    if (status == 0) {
      container->amg_setup_calls++;
      container->preconditioner_setup_requested = 0;
    }
    return status;
  }
}

static HYPRE_Int Phase24BoomerAMGSolve(HYPRE_Solver context_solver,
                                       HYPRE_Matrix matrix,
                                       HYPRE_Vector b,
                                       HYPRE_Vector x)
{
  Phase24BoomerAMGContext *context = (Phase24BoomerAMGContext *)context_solver;
  if (context == NULL || context->boomer == NULL)
    return -1;
  return HYPRE_BoomerAMGSolve(
    context->boomer, (HYPRE_ParCSRMatrix)matrix,
    (HYPRE_ParVector)b, (HYPRE_ParVector)x);
}

static int RefreshHypreSetup(ElmerHypreContainer *Container,
                             HYPRE_ParCSRMatrix parcsr_A,
                             HYPRE_ParVector par_b,
                             HYPRE_ParVector par_x)
{
  int hypre_sol = Container->hypre_method / 100;
  int hypre_pre = Container->hypre_method % 100;
  int status;

  /* For the Phase24 reuse path the IJ/ParCSR matrix object remains the same
     object and only its values have been refreshed.  Re-running the
     FlexGMRES setup here allocates a new Krylov workspace on every refresh;
     the Windows telemetry showed a persistent ~67 MiB increase per call.
     Keep the existing solver workspace for reuse, while retaining the real
     setup for initial construction, requested AMG rebuilds, and recovery. */
  if (hypre_sol == 9 && hypre_pre == 1 &&
      !Container->preconditioner_setup_requested &&
      Container->setup_epoch >= 0 &&
      !phase24_force_krylov_setup_is_enabled() &&
      !phase24_force_amg_setup_is_enabled())
    return 0;

  Container->krylov_setup_calls++;
  /* FlexGMRES invokes the registered BoomerAMG setup callback internally.
     Count that native callback at the same wrapper boundary as the outer
     setup; this is the actual setup event seen by the persistent path. */
  if ((hypre_pre == 1 || hypre_sol == 1) &&
      !(hypre_sol == 9 && hypre_pre == 1)) Container->amg_setup_calls++;
  /* HYPRE returns its process-global error bit from subsequent API calls.
     Clear a prior convergence bit so this status belongs to this setup. */
  HYPRE_ClearAllErrors();

  switch (hypre_sol) {
  case 1:
    status = HYPRE_BoomerAMGSetup(Container->solver, parcsr_A, par_b, par_x);
    break;
  case 2:
    status = HYPRE_AMSSetup(Container->solver, parcsr_A, par_b, par_x);
    break;
  case 6:
    status = HYPRE_ParCSRPCGSetup(Container->solver, parcsr_A, par_b, par_x);
    break;
  case 7:
    status = HYPRE_ParCSRBiCGSTABSetup(Container->solver, parcsr_A, par_b, par_x);
    break;
  case 8:
    status = HYPRE_ParCSRGMRESSetup(Container->solver, parcsr_A, par_b, par_x);
    break;
  case 9:
    status = HYPRE_ParCSRFlexGMRESSetup(Container->solver, parcsr_A, par_b, par_x);
    break;
  case 10:
    status = HYPRE_ParCSRLGMRESSetup(Container->solver, parcsr_A, par_b, par_x);
    break;
  case 11:
    status = HYPRE_ParCSRCOGMRESSetup(Container->solver, parcsr_A, par_b, par_x);
    break;
  default:
    status = -1;
    break;
  }
  return status;
}


/* The interface calls Hypre step-wise separating phases for 
   procedure of setup, solve and cleanup.
   The original monolithic one that was not used for years has been removed. 

   Phase24 SolveHYPRE3 updates the persistent current operator in place and
   can retain the BoomerAMG hierarchy through the setup callback wrapper. */
   
/* initialization for a new matrix.
   - convert matrix
   - setup solver and preconditioner
   - return a pointer 'Container' which the calling fortran
   program should not alter but pass back into subsequent
   SolveHYPRE2, ~3 and ~4 calls.

 This function has an additional feature compared to the SolveHYPRE call above,
 namely to use a block diagonal approximation of A for the preconditioner setup.
 This mimics the behavior of the BILUn preconditioners in Elmer, although any   
 preconditioner (like ParaSails or BoomerAMG) can still be used in combination  
 with block diagonal approximation. 
 BILU=0 or 1 - use A. 
 BILU=k - assume k equations and use block diagonal A with k blocks.
*/
void STDCALLBULL FC_FUNC(solvehypre1,SOLVEHYPRE1)
 (
  int *nrows,int *rows, int *cols, double *vals, int *precflag, double *precvals, 
  int *globaldofs, int *owner, int *ILUn, int *BILU, int *hypre_method,
  int *hypre_intpara, double *hypre_dppara,
  int *Rounds, double *TOL, int *verbosityPtr, int *structureEpochPtr,
  int *matrixEpochPtr, int *phase24KeepPtr, int** ContainerPtr,
  int *fcomm
 )
{
   int i, j, k;
   int myid, num_procs;
   int N, n, csize=128;

   int ilower, iupper;
   int local_size, extra;
   int hypre_sol, hypre_pre, hypre_use_gpu;
   MPI_Comm comm=MPI_Comm_f2c(*fcomm);
   ElmerHypreContainer* Container;

   HYPRE_IJMatrix A, Atilde;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;

   HYPRE_Solver solver, precond;

   int num_iterations;
   int setup_status = 0;
   double final_res_norm;
  
   double st, realtime_();
   
   int verbosity = *verbosityPtr, myverb;
   HYPRE_Int init_status;

   phase24_memory_trace_event(10, *nrows, rows[*nrows] - rows[0],
                              matrixEpochPtr ? *matrixEpochPtr : -1, NULL);

   /* HYPRE 3 requires an explicit initialization before any other call.
      Initialization is process-global; do not repeat it for every rebuilt
      container during a long transient run. */
   if (!HYPRE_Initialized()) {
     init_status = HYPRE_Initialize();
     if (init_status != 0) {
       fprintf(stderr, "SolveHypre: HYPRE_Initialize failed with status %d\n",
               (int) init_status);
       return;
     }
   }

   /* which process number am I? */
   MPI_Comm_rank(comm, &myid);

   if( myid == 0 ) 
     myverb = verbosity;
   else
     myverb = 0;
   
   if (myverb > 8) fprintf(stdout,"SolveHypre: Performing HYPRE Setup\n");
   
   st  = realtime_();
   
   /* How many rows do I have? */
   local_size = *nrows;
   hypre_sol = *hypre_method / 100;
   hypre_pre = *hypre_method % 100;
   hypre_use_gpu = hypre_intpara[18];

#if !defined(HYPRE_USING_GPU)
   if (hypre_use_gpu) {
     if (myverb > 0) fprintf(stderr,
       "SolveHypre: GPU requested but this HYPRE library has no GPU backend; using CPU fallback\n");
     hypre_use_gpu = 0;
   }
#endif

   /* Keep IJ insertion in host memory.  The persistent GPU branch below
      migrates the assembled ParCSR objects once and stages only changed
      values at the Phase24 boundary. */
   /* HYPRE's CUDA/HIP IJ interface accepts host-side FEM insertion buffers
      and stages them for device kernels.  Keep that supported boundary mode;
      explicit IJ migration is not safe with the current non-host-memory
      HYPRE build (it faults before the first solve). */
   HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
   HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
   if (hypre_use_gpu) {
     /* Explicitly create the CUDA/HIP stream and bind HYPRE to the process's
        already-selected device before any IJ allocation. */
     HYPRE_DeviceInitialize();
     if (myverb > 4) HYPRE_PrintDeviceInfo();
   }
   if (hypre_use_gpu && myverb > 0) fprintf(stdout,
     "SolveHypre: Phase24 HYPRE backend=%s, persistent device execution enabled\n",
     Phase24HypreBackendName());

   if(hypre_sol == 2 || hypre_pre == 2)  {
     if (*ContainerPtr == NULL) {
       fprintf( stdout, "Hypre pointer should not be zero at start for AMS\n");
     }   
     Container = (ElmerHypreContainer*)(*ContainerPtr);
   }
   else {   
      Container = (ElmerHypreContainer*)malloc(sizeof(ElmerHypreContainer));   
      *ContainerPtr=(int*)(Container);
      memset(Container, 0, sizeof(ElmerHypreContainer));
   }
   if (Container->phase24_boomer_context == NULL) {
     Container->phase24_boomer_context =
       (Phase24BoomerAMGContext *)malloc(sizeof(Phase24BoomerAMGContext));
     if (Container->phase24_boomer_context == NULL) {
       fprintf(stderr, "SolveHypre: unable to allocate Phase24 AMG context\n");
       return;
     }
     Container->phase24_boomer_context->container = Container;
     Container->phase24_boomer_context->boomer = NULL;
   }
   /* Make the method available to the native setup counters before the first
      solver setup call.  The epoch is committed below only after setup has
      returned successfully. */
   Container->hypre_method = *hypre_method;
   Container->setup_epoch = -1;
   Container->preconditioner_epoch = -1;
   Container->last_solved_matrix_epoch = -1;
   ilower =  1000000000;
   iupper = -1;
   for( i=0; i<local_size; i++ ) {
     if ( owner[i] ) {
       if ( iupper < globaldofs[i] ) iupper = globaldofs[i];
       if ( ilower > globaldofs[i] ) ilower = globaldofs[i];
     }
   }


   /* if the partition doesn't own any of the dofs, apply null range (with valid indices) */
   if ( iupper == -1 ) { ilower = 1; iupper = 0; }

   /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);
   Container->ij_matrix_initialize_calls++;

   /* Now go through my local rows and set the matrix entries.
      Note that here we are setting one row at a time, though
      one could set all the rows together (see the User's Manual).
   */
   {
      int nnz,irow,i,j,k,*rcols;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++) {
	nnz = rows[i+1]-rows[i];
	if ( nnz>csize ) {
          csize = nnz+csize;
          rcols = (int *)realloc( rcols, csize*sizeof(int) );
	}
	irow=globaldofs[i];
        for( k=0,j=rows[i]; j<rows[i+1]; j++,k++) {
	  rcols[k] = globaldofs[cols[j-1]-1];
	}
        HYPRE_IJMatrixAddToValues(A, 1, &nnz, &irow, rcols, &vals[rows[i]-1]);
      }
      free( rcols );
   }

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);
   Container->ij_matrix_assemble_calls++;
   phase24_memory_trace_event(11, local_size, rows[local_size] - rows[0],
                              matrixEpochPtr ? *matrixEpochPtr : -1,
                              Container);

   if (!*precflag && *BILU <= 1) {
     /* Standard version - use A as preconditioner */
     Atilde = A;
   } else if ( *precflag ) {
     /* We have another matrix that is used as preconditioning matrix! */
     int nnz,irow,jcol,i,j,k,*rcols;
     double *dbuf;

     HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &Atilde);
     HYPRE_IJMatrixSetObjectType(Atilde, HYPRE_PARCSR);
     HYPRE_IJMatrixInitialize(Atilde);
     {
        int nnz,irow,i,j,k,*rcols;

        rcols = (int *)malloc( csize*sizeof(int) );
        for (i = 0; i < local_size; i++) {
          nnz = rows[i+1]-rows[i];
          if ( nnz>csize ) {
            csize = nnz+csize;
            rcols = (int *)realloc( rcols, csize*sizeof(int) );
          }
          irow=globaldofs[i];
          for( k=0,j=rows[i]; j<rows[i+1]; j++,k++) {
             rcols[k] = globaldofs[cols[j-1]-1];
          }
          HYPRE_IJMatrixAddToValues(Atilde, 1, &nnz, &irow, rcols, &precvals[rows[i]-1]);
        }
        free( rcols );
     }
     /* Assemble after setting the coefficients */
     HYPRE_IJMatrixAssemble(Atilde);     
   } else {
     /* We only take the block diagonal values of the original matrix for our preconditioner */
     int nnz,irow,jcol,i,j,k,*rcols;
     double *dbuf;
     if (myverb > 6) fprintf(stdout,"HYPRE: using BILU(%d) approximation for preconditioner\n",*BILU);
     
     HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &Atilde);
     HYPRE_IJMatrixSetObjectType(Atilde, HYPRE_PARCSR);
     HYPRE_IJMatrixInitialize(Atilde);
     
     rcols = (int *)malloc( csize*sizeof(int) );
     dbuf = (double *)malloc( csize*sizeof(double) );
     for (i = 0; i < local_size; i++) {
       irow=globaldofs[i];
       nnz = 0;
       for (j=rows[i];j<rows[i+1];j++) {
         jcol = globaldofs[cols[j-1]-1];
         /*TODO - is the block ordering preserved in the linear numbering?
	   Here we assume it is.
	 */
         if ((irow%*BILU)==(jcol%*BILU)) {
	   rcols[nnz] = jcol;
           dbuf[nnz] = vals[j-1];
           nnz++;
	 }
       }
       HYPRE_IJMatrixAddToValues(Atilde, 1, &nnz, &irow, rcols, dbuf);
     }
     free( rcols );
     free( dbuf );
     /* Assemble after setting the coefficients */
     HYPRE_IJMatrixAssemble(Atilde);     
   }

   /* Get the parcsr matrix object to use */
   /* note: this is only used for setup,  */
   /* so we put in the possibly approxima-*/
   /* ted matrix Atilde                   */
   HYPRE_IJMatrixGetObject(Atilde, (void**) &parcsr_A);
   Container->ij_matrix_get_object_calls++;
   
   HYPRE_IJVectorCreate(comm, ilower, iupper,&b);
   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(b);
   HYPRE_IJVectorAssemble(b);
   Container->ij_vector_assemble_calls++;
   HYPRE_IJVectorGetObject(b, (void **) &par_b);
   Container->ij_vector_get_object_calls++;

   HYPRE_IJVectorCreate(comm, ilower, iupper,&x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);
   HYPRE_IJVectorAssemble(x);
   Container->ij_vector_assemble_calls++;
   HYPRE_IJVectorGetObject(x, (void **) &par_x);
   Container->ij_vector_get_object_calls++;

   if (hypre_use_gpu && phase24_migrate_enabled()) {
     hypre_ParCSRMatrixMigrate(parcsr_A, HYPRE_MEMORY_DEVICE);
     hypre_ParVectorMigrate(par_b, HYPRE_MEMORY_DEVICE);
     hypre_ParVectorMigrate(par_x, HYPRE_MEMORY_DEVICE);
     if (myverb > 0) fprintf(stdout, "SolveHypre: Phase24 ParCSR/ParVector migrated to device\n");
   }

   if (hypre_use_gpu) HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);

   if (hypre_use_gpu) {
     /* The persistent solver, preconditioner, IJ handles, and ParCSR object
        are created once.  HYPRE owns the device-side work buffers behind this
        supported host-insertion/device-execution boundary. */
     Container->gpu_device_allocation_calls += (Atilde == A) ? 3 : 4;
   }

   if(myverb > 12) {     
     fprintf(stdout,"SolveHypre: NonZero Hypre parameters\n");
     for(i=0;i<=12;i++) {
       if(hypre_intpara[i]) fprintf(stdout,"  intpara %d: %d\n",i+1,hypre_intpara[i]);
     }
     for(i=0;i<=5;i++) {
       if(fabs(hypre_dppara[i]) > 1.0e-20) fprintf(stdout,"   dppara %d: %lg\n",i+1,hypre_dppara[i]);
     }
   }
   
    /* This is copy-pasted from SParIterSolver for convenience.
    !---------------------------------------
    !              No  Prec
    ! none         0    x   -
    ! BoomerAMG    1    x   x 
    ! AMS          2    x   x
    ! ILU          3    x   x
    ! Parasails    4    x   -
    ! FSAI         5    x   -
    ! PCG          6    -   x
    ! BiCGStab     7    -   x
    ! GMRes        8    -   x
    ! FlexGMRes    9    -   x
    ! LGMRes       10   -   x
    ! COGMRes      11   -   x
    !---------------------------------------  */  

   
   /* Create preconditioner for Krylov methods. */
   /* Some methods may act as preconditioners and solvers and are also initialized here
      and later changes the pointer to solvers. */

   if ( hypre_pre == 1 || hypre_sol == 1 )  {
     if( myverb > 8 ) 
       fprintf( stdout,"SolveHypre: Creating Hypre BoomerAMG\n");
     if( myverb > 10 ) {
       fprintf( stdout,"RelaxType = %d\n",hypre_intpara[0]); 
       fprintf( stdout,"CoarsenType = %d\n",hypre_intpara[1]); 
       fprintf( stdout,"NumSweeps = %d\n",hypre_intpara[2]); 
       fprintf( stdout,"MaxLevels = %d\n",hypre_intpara[3]); 
       fprintf( stdout,"Interpolation Type = %d\n",hypre_intpara[4]); 
       fprintf( stdout,"Smooth Type = %d\n",hypre_intpara[5]);
       fprintf( stdout,"Cycle Type = %d\n",hypre_intpara[6]);
       fprintf( stdout,"DOFs = %d\n",hypre_intpara[7]);
       fprintf( stdout,"StrongThreshold = %g\n",hypre_dppara[0]);
     }
     HYPRE_BoomerAMGCreate(&precond);
     Container->phase24_boomer_context->boomer = precond;
     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_BoomerAMGSetNumFunctions(precond, hypre_intpara[7]); /* No. of PDE's */

     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_BoomerAMGSetPrintLevel(precond, i); /* print amg solution info */

     HYPRE_BoomerAMGSetNumSweeps(precond, hypre_intpara[2]); /* fixed for preconditioner to 1 */
     HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
     HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
     HYPRE_BoomerAMGSetRelaxType(precond, hypre_intpara[0]);   /* G-S/Jacobi hybrid relaxation */
     HYPRE_BoomerAMGSetCoarsenType(precond, hypre_intpara[1]);  /* coarsening type */
     
     HYPRE_BoomerAMGSetMaxLevels(precond, hypre_intpara[3]); /* levels of coarsening */
     HYPRE_BoomerAMGSetInterpType(precond, hypre_intpara[4]);  /* interpolation type */
     HYPRE_BoomerAMGSetSmoothType(precond, hypre_intpara[5]);  /* smoother type */
     HYPRE_BoomerAMGSetCycleType(precond, hypre_intpara[6]);  /* coarsening type */
     /* threshold for strong coupling (default 0.25 recommended for 2D Laplace, 0.5-0.6 
	for 3D Laplace, 0.9 for elasticity) */
     HYPRE_BoomerAMGSetStrongThreshold(precond, hypre_dppara[0]);
     /* Optional operator-complexity controls (negative = HYPRE default). */
     if (hypre_intpara[14] >= 0) HYPRE_BoomerAMGSetPMaxElmts(precond, hypre_intpara[14]);
     if (hypre_dppara[6] >= 0.0) HYPRE_BoomerAMGSetTruncFactor(precond, hypre_dppara[6]);
     if (hypre_intpara[15] >= 0) HYPRE_BoomerAMGSetAggNumLevels(precond, hypre_intpara[15]);
     if (hypre_intpara[16] >= 0) HYPRE_BoomerAMGSetAggInterpType(precond, hypre_intpara[16]);
     if (hypre_intpara[17] >= 0) HYPRE_BoomerAMGSetPrintLevel(precond, hypre_intpara[17]);

     if( myverb > 10) fprintf(stdout,"SolveHypre: Created BoomerAMG preconditioner!\n");

   } else if ( hypre_pre == 2 || hypre_sol == 2 ) {
     precond = Container->precond;
     if(precond) {
       if(myverb > 10) fprintf( stdout,"SolveHypre: Using previously defined AMS!\n");
     } else {
       fprintf( stdout,"SolveHypre: Pointer to AMS preconditioner is NULL!\n");       
       exit(EXIT_FAILURE);
     }
   
   } else if ( hypre_pre == 3 || hypre_sol == 3 ) {
     int ilu_type, max_iter, reordering, print_level;
     int max_nnz_row, schur_max_iter,tri_solve, ljac_iters,ujac_iters; 
     double tol, threshold;
       
     if (myverb > 8) fprintf( stdout,"SolveHypre: using ILU%d as preconditioner\n",*ILUn); 
          
     /* (Required) Create ILU solver */
     HYPRE_ILUCreate(&precond);

     /* (Recommended) General solver options */
     HYPRE_ILUSetType(precond, ilu_type=*ILUn); /* 0, 1, 10, 11, 20, 21, 30, 31, 40, 41, 50 */
     HYPRE_ILUSetMaxIter(precond, max_iter=1);
#if 0 
     HYPRE_ILUSetTol(precond, tol);
     HYPRE_ILUSetLocalReordering(precond, reordering); /* 0: none, 1: RCM */
     HYPRE_ILUSetPrintLevel(precond, print_level);

     /* (Optional) Function calls for ILUK variants */
     HYPRE_ILUSetLevelOfFill(precond, fill);

     /* (Optional) Function calls for ILUT variants */
     HYPRE_ILUSetMaxNnzPerRow(precond, max_nnz_row);
     HYPRE_ILUSetDropThreshold(precond, threshold);

     /* (Optional) Function calls for GMRES-ILU or NSH-ILU */
     HYPRE_ILUSetNSHDropThreshold(precond, threshold);
     HYPRE_ILUSetSchurMaxIter(precond, schur_max_iter);

     /* (Optional) Function calls for iterative ILU variants */
     HYPRE_ILUSetTriSolve(precond, tri_solve);
     HYPRE_ILUSetLowerJacobiIters(precond, ljac_iters);
     HYPRE_ILUSetUpperJacobiIters(precond, ujac_iters);
     
     /* (Exclusively required) Function calls for using ILU as standalone solver */
     HYPRE_ILUSetup(precond, parcsr_M, b, x);
     HYPRE_ILUSolve(precond, parcsr_A, b, x);

     /* (Exclusively required) Function calls for using ILU as preconditioner to GMRES */
#endif

   }

   else if ( hypre_pre == 4 ) {
     if (myverb > 8) fprintf( stdout,"SolveHypre: using ParaSails as preconditioner\n"); 

     /* Now set up the ParaSails preconditioner and specify any parameters */
     HYPRE_ParaSailsCreate(comm, &precond);

     /* Set some parameters (See Reference Manual for more parameters) */
     /* threshold = dppara[0]; maxlevels= intpara[1] */
     HYPRE_ParaSailsSetParams(precond, hypre_dppara[0], hypre_intpara[1]);
     /* filter = dppara[1] */
     HYPRE_ParaSailsSetFilter(precond, hypre_dppara[1]);
     /* symmetry = intpara[0] */
     HYPRE_ParaSailsSetSym(precond, hypre_intpara[0]);

     i = 3*(verbosity >= 6 );
     HYPRE_ParaSailsSetLogging(precond, i);

   } else if ( hypre_pre == 5 ) {
     int max_steps=5,max_step_size=3;
     double kap_tolerance=1.0e-3;
#if 0
     if (myverb > 8) fprintf( stdout,"SolveHypre: using SPAI as preconditioner\n"); 
     HYPRE_FSAICreate(&precond);
     HYPRE_FSAISetMaxSteps(precond, max_steps);
     HYPRE_FSAISetMaxStepSize(precond, max_step_size);
     HYPRE_FSAISetKapTolerance(precond, kap_tolerance);     
     i = 3*(verbosity >= 6 );
     HYPRE_SPAISetLogging(precond, i);
#else
     fprintf( stdout,"Hypre preconditioning method FSAI not compiled with!\n");
     exit(EXIT_FAILURE);
#endif  
   } else if ( hypre_pre == 12 ) {
     int mgr_cstart = hypre_intpara[12];
     int mgr_ccount = hypre_intpara[13];
     HYPRE_Int mgr_num_cpoints[1] = {1};
     HYPRE_Int mgr_cpoint_tags[1] = {1};
     HYPRE_Int *mgr_tag_ptr[1] = {mgr_cpoint_tags};

     if (mgr_cstart < 1 || mgr_ccount < 1 ||
         mgr_cstart + mgr_ccount - 1 > local_size) {
       fprintf(stderr,
          "SolveHypre: invalid MGR C-point range start=%d count=%d local_rows=%d\n",
          mgr_cstart, mgr_ccount, local_size);
       exit(EXIT_FAILURE);
     }
     if (myverb > 8) fprintf(stdout,
          "SolveHypre: creating MGR preconditioner (C rows %d..%d)\n",
          mgr_cstart, mgr_cstart + mgr_ccount - 1);

     Container->mgr_marker = (int *)malloc((size_t)local_size * sizeof(int));
     if (!Container->mgr_marker) {
       fprintf(stderr, "SolveHypre: MGR marker allocation failed\n");
       exit(EXIT_FAILURE);
     }
     for (i = 0; i < local_size; i++) {
       Container->mgr_marker[i] =
          (i + 1 >= mgr_cstart && i + 1 < mgr_cstart + mgr_ccount) ? 1 : 0;
     }

     HYPRE_MGRCreate(&precond);
     /* Tags 0/1 identify F/C points; block_size=2 avoids MGR's scalar
        shortcut while point_marker_array supplies the actual local split. */
     HYPRE_MGRSetCpointsByPointMarkerArray(precond, 2, 1,
          mgr_num_cpoints, mgr_tag_ptr, Container->mgr_marker);
     HYPRE_MGRSetNonCpointsToFpoints(precond, 1);
     /* Use a dedicated BoomerAMG object for the F block.  MGR only consumes
        HYPRE_MGRSetFSolver when the F-relaxation type is AMG (2); setting a
        Jacobi relax type here silently bypasses the user solver. */
     HYPRE_BoomerAMGCreate(&Container->mgr_fsolver);
     HYPRE_BoomerAMGSetTol(Container->mgr_fsolver, 1.0e-4);
     HYPRE_BoomerAMGSetMaxIter(Container->mgr_fsolver, 2);
     HYPRE_BoomerAMGSetCoarsenType(Container->mgr_fsolver, 8);
     HYPRE_BoomerAMGSetInterpType(Container->mgr_fsolver, 6);
     HYPRE_BoomerAMGSetRelaxType(Container->mgr_fsolver, 18);
     HYPRE_BoomerAMGSetNumSweeps(Container->mgr_fsolver, 1);
     HYPRE_MGRSetFSolver(precond,
          (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve,
          (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup,
          Container->mgr_fsolver);
     HYPRE_MGRSetRelaxType(precond, 2);
     HYPRE_MGRSetNumRelaxSweeps(precond, 2);
     /* Match the GPU-capable BoomerAMG interpolation used by no-mortar. */
     /* Approximate-inverse restriction is the first Schur-quality variant;
        injection is retained in the history as the baseline. */
     HYPRE_MGRSetRestrictType(precond, 3);
     HYPRE_MGRSetNumRestrictSweeps(precond, 2);
     HYPRE_MGRSetInterpType(precond, 2);
     HYPRE_MGRSetNumInterpSweeps(precond, 2);
     HYPRE_MGRSetMaxIter(precond, 1);
     HYPRE_MGRSetTol(precond, 0.0);
     i = (verbosity >= 6);
     if (verbosity >= 10) i = 3;
     HYPRE_MGRSetPrintLevel(precond, i);

   } else if( hypre_pre != 0) {
     fprintf( stdout,"Hypre preconditioning method %d not implemented\n",hypre_pre);
     exit(EXIT_FAILURE);
   }

   /* Create Hypre solver */

   if ( hypre_sol == 1 ) { /* boomer AMG */     
     solver = precond;
     precond = NULL;
     /* Now setup - note that the input vectors are ignored so we can pass in NULLs */
     if (myverb > 10 ) fprintf(stdout,"SolveHypre: construct BoomerAMG solver");

     /* When BoomerAMG is used as a solver these override the definitions related to AMS as a preconditioner. */
     HYPRE_BoomerAMGSetTol(solver, *TOL);      
     HYPRE_BoomerAMGSetMaxIter(solver, *Rounds); 

     Container->amg_setup_calls++;
     setup_status = HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
     

   } else if ( hypre_sol == 2 ) { /* AMS */
     solver = precond;
     precond = NULL;
     
     /* Now setup - note that the input vectors are ignored so we can pass in NULLs */
     if (myverb > 10 ) fprintf(stdout,"SolveHypre: construct AMS solver");

     /* When AMS is used as a solver these override the definitions related to AMS as a preconditioner. */
     HYPRE_AMSSetMaxIter(solver,*Rounds);
     HYPRE_AMSSetTol(solver, *TOL);

     HYPRE_AMSSetup(solver, parcsr_A, par_b, par_x);


   } else if ( hypre_sol == 6) { /* PCG */
     /* Create solver */
     HYPRE_ParCSRPCGCreate(comm, &solver);
     
     HYPRE_PCGSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_PCGSetTol(solver, *TOL);        /* conv. tolerance */

     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_PCGSetAbsoluteTol(solver, hypre_dppara[5]);       

     HYPRE_PCGSetTwoNorm(solver,hypre_intpara[10]);  /* use the two norm as the stopping criteria */
     //HYPRE_PCGSetFlex(solver,hypre_intpara[11]);     /* use flexible method for robustness */
     
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_PCGSetPrintLevel(solver, i);   /* print solve info */

     i = (verbosity >= 6);
     HYPRE_PCGSetLogging(solver, i);      /* needed to get run info later */
     
     /* Set the PCG preconditioner */
     if( hypre_pre == 1 && Atilde != A && *precflag ) {
       /* Separate preconditioning matrix: Krylov on the true A, BoomerAMG
          built and applied on Atilde only (upstream passed Atilde to the
          Krylov solve as well, i.e. solved the wrong system). */
       HYPRE_ParCSRMatrix parcsr_true;
       double t_amg = realtime_();
       HYPRE_BoomerAMGSetup(precond, parcsr_A, par_b, par_x);
       phase24_tilde_matrix = parcsr_A;
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) phase24_tilde_amg_solve,
                           (HYPRE_PtrToSolverFcn) phase24_noop_setup, precond);
       HYPRE_IJMatrixGetObject(A, (void**) &parcsr_true);
       parcsr_A = parcsr_true;
       if (getenv("PHASE24_HYPRE_REPORT") != NULL)
         fprintf(stdout, "SolveHypre: Phase24 BoomerAMG setup on preconditioning matrix %.3f s\n",
                 realtime_() - t_amg);
     } else if( hypre_pre == 1){
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre == 3) {
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if ( hypre_pre == 4) { 
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
       
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up PCG linear system");
     {
       double t_setup = realtime_();
       HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
       if (getenv("PHASE24_HYPRE_REPORT") != NULL)
         fprintf(stdout, "SolveHypre: Phase24 PCG+AMG setup %.3f s (total since entry %.3f s)\n",
                 realtime_() - t_setup, realtime_() - st);
     }


   } else if ( hypre_sol == 7) { /* BiGSTAB methods */
     /* Create solver */
     HYPRE_ParCSRBiCGSTABCreate(comm, &solver);

     HYPRE_BiCGSTABSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_BiCGSTABSetTol(solver, *TOL);       /* conv. tolerance */

     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_BiCGSTABSetAbsoluteTol(solver, hypre_dppara[5]);       
     
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_BiCGSTABSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity >= 6);
     HYPRE_BiCGSTABSetLogging(solver, i);      /* needed to get run info later */
 
     /* Set the BiCGStabl preconditioner */
     if(hypre_pre == 1) {
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre == 3) {
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if (hypre_pre == 4) { 
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up BiCGStab linear system");
     HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x);     

     
   } else if ( hypre_sol == 8) { /* GMRES */
     /* Create solver */
     HYPRE_ParCSRGMRESCreate(comm, &solver);

     HYPRE_GMRESSetMaxIter(solver, *Rounds); /* max GMRES iterations */
     HYPRE_GMRESSetTol(solver, *TOL);        /* GMRES conv. tolerance */

     HYPRE_GMRESSetKDim(solver, hypre_intpara[9]);
     HYPRE_GMRESSetAbsoluteTol(solver, hypre_dppara[5]);       
     
     
     /* Set some parameters (See Reference Manual for more parameters) */
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_GMRESSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity > 6);
     HYPRE_GMRESSetLogging(solver, i);      /* needed to get run info later */

     if( hypre_pre == 1){
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre  == 3) {
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if (hypre_pre == 4) { 
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
       
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up GMRes linear system");
     HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x);   

     
   } else if ( hypre_sol == 9) { /* FlexGMRes */
     /* Create solver */
     HYPRE_ParCSRFlexGMRESCreate(comm, &solver);

     HYPRE_FlexGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_FlexGMRESSetTol(solver, *TOL);       /* conv. tolerance */

     HYPRE_FlexGMRESSetKDim(solver,hypre_intpara[9]);
     HYPRE_FlexGMRESSetAbsoluteTol(solver, hypre_dppara[5]);       
     
     /* Set some parameters (See Reference Manual for more parameters) */
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_FlexGMRESSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity >= 6);
     HYPRE_FlexGMRESSetLogging(solver, i);      /* needed to get run info later */
     if (phase24_hypre_residual_trace_is_enabled()) {
       /* HYPRE print level 3 emits the native Krylov residual progression;
          logging keeps the solver's iteration counters available as well. */
       HYPRE_FlexGMRESSetPrintLevel(solver, 3);
       HYPRE_FlexGMRESSetLogging(solver, 1);
       fprintf(stderr,
               "PHASE24_HYPRE_RESIDUAL_TRACE_ENABLED method=%d max_iterations=%d tolerance=%.17e\n",
               *hypre_method, *Rounds, *TOL);
       fflush(stderr);
     }

     /* Set the FlexGMRES preconditioner */
     if( hypre_pre == 1) {
       HYPRE_FlexGMRESSetPrecond(solver,
                                 Phase24BoomerAMGSolve,
                                 Phase24BoomerAMGSetup,
                                 (HYPRE_Solver)Container->phase24_boomer_context);
     } else if(hypre_pre == 2 ) {
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
				 (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre == 3) {
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if ( hypre_pre == 4) { 
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     } else if ( hypre_pre == 12) {
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_MGRSolve,
                                 (HYPRE_PtrToSolverFcn) HYPRE_MGRSetup, precond);
     }
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up FlexGMRes linear system");
     Container->krylov_setup_calls++;
     if (hypre_pre == 1) Container->preconditioner_setup_requested = 1;
     setup_status = HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);

     
   } else if ( hypre_sol == 10) { /* LGMRes */
     /* Create solver */
     HYPRE_ParCSRLGMRESCreate(comm, &solver);

     HYPRE_LGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_LGMRESSetTol(solver, *TOL);       /* conv. tolerance */

     HYPRE_LGMRESSetKDim(solver,hypre_intpara[9]);     
     HYPRE_LGMRESSetAbsoluteTol(solver, hypre_dppara[5]);
     
     HYPRE_LGMRESSetAugDim(solver, hypre_intpara[10]);       
     
     /* Set some parameters (See Reference Manual for more parameters) */
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_LGMRESSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity >= 6);
     HYPRE_LGMRESSetLogging(solver, i);      /* needed to get run info later */

     /* Set the LGMRES preconditioner */
     if( hypre_pre == 1){
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
				 (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     }
     else if ( hypre_pre  == 3) {
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if ( hypre_pre == 4) { 
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
       
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up LGMRes linear system");
     HYPRE_ParCSRLGMRESSetup(solver, parcsr_A, par_b, par_x);

     
   } else if ( hypre_sol == 11) { /* COGMRES */
     /* Create solver */
     HYPRE_ParCSRCOGMRESCreate(comm, &solver);

     /* Set some parameters (See Reference Manual for more parameters) */
     //HYPRE_COGMRESSetStopCrit(solver, 0);     /* use the two norm as the stopping criteria */
     HYPRE_COGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_COGMRESSetTol(solver, *TOL);       /* conv. tolerance */

     HYPRE_COGMRESSetKDim(solver,hypre_intpara[9]);
     HYPRE_COGMRESSetAbsoluteTol(solver, hypre_dppara[5]);       

     HYPRE_COGMRESSetUnroll(solver, hypre_intpara[10]);       
     HYPRE_COGMRESSetCGS(solver, hypre_intpara[11]);       


     
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_COGMRESSetPrintLevel(solver, i);   /* print solve info */
     HYPRE_COGMRESSetLogging(solver,verbosity>=6);      /* needed to get run info later */

     /* Set the COGMRES preconditioner */
     if(hypre_pre == 1) {
       HYPRE_COGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_COGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre == 3) {
       HYPRE_COGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if (hypre_pre == 4) { 
       HYPRE_COGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
       
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up COGMRes linear system");
     HYPRE_ParCSRCOGMRESSetup(solver, parcsr_A, par_b, par_x);

     
   } else {
     fprintf( stdout,"SolveHypre: Hypre solver method %d not implemented!\n",hypre_sol);
     exit(EXIT_FAILURE);
   }

   Container->ilower = ilower;
   Container->iupper = iupper;     
   Container->hypre_method = *hypre_method;
   Container->hypre_use_gpu = hypre_use_gpu;
   Container->A = A;
   Container->b = b;
   Container->x = x;
   Container->Atilde = Atilde;
   Container->parcsr_A = parcsr_A;
   Container->structure_epoch = structureEpochPtr ? *structureEpochPtr : -1;
   Container->matrix_epoch = matrixEpochPtr ? *matrixEpochPtr : -1;
   /* SolveHYPRE1 has completed the initial native setup at this point. */
   Container->setup_epoch = (setup_status == 0) ? Container->matrix_epoch : -1;
   Container->preconditioner_epoch = Container->setup_epoch;
   Container->setup_since_last_solve = (setup_status == 0);
   if (setup_status != 0 && myverb > 0)
     fprintf(stderr, "SolveHypre: initial native setup failed with status %d\n",
             setup_status);
   Container->vectors_initialized = 0;
   Container->phase24_keep_alive = phase24KeepPtr ? *phase24KeepPtr : 0;
   Container->solver = solver;
   Container->precond = precond;

   phase24_memory_trace_event(12, local_size, rows[local_size] - rows[0],
                              Container->matrix_epoch, Container);
   
   if (myverb > 5) fprintf( stdout, "SolveHypre: setup time (method %d): %g\n", 
			    Container->hypre_method, realtime_()-st );
   
} /* SolveHypre1 - matrix conversion and solver setup */


/* Update the stopping tolerance of a previously constructed solver */
void STDCALLBULL FC_FUNC(updatehypre,UPDATEHYPRE)
  ( double *TOL,  int *hypre_method,  int** ContainerPtr, int *verbosityPtr, int *fcomm )
{
  HYPRE_Solver solver;
  ElmerHypreContainer *Container;
  int hypre_sol;
  
  Container = (ElmerHypreContainer*)(*ContainerPtr);
  solver = Container->solver;

  hypre_sol = *hypre_method / 100;

  if ( hypre_sol == 1 ) { /* boomer AMG */
    HYPRE_BoomerAMGSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 2 ) { /* AMS */
    HYPRE_AMSSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 6) { /* CG */
    HYPRE_PCGSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 7) { /* BiGSTAB method */
    HYPRE_BiCGSTABSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 8) { /* GMRES */
    HYPRE_GMRESSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 9) { /* FlexGMRes */
    HYPRE_FlexGMRESSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 10) { /* LGMRes */
    HYPRE_LGMRESSetTol(solver, *TOL);
  }
  if ( hypre_sol == 11) { /* COGMRES method */
    HYPRE_COGMRESSetTol(solver, *TOL);
  }
  else {
    fprintf( stdout,"SolveHypre: Hypre solver method not implemented\n");
    exit(EXIT_FAILURE);
  }
}


/*////////////////////////////////////////////////////////////////////////////////////////////////*/

/* Return the epochs associated with the currently persistent HYPRE setup. */
void STDCALLBULL FC_FUNC(solvehypre6,SOLVEHYPRE6)
  ( int** ContainerPtr, int *structureEpochPtr, int *matrixEpochPtr,
    int *preconditionerEpochPtr ) {
  ElmerHypreContainer *Container = (ElmerHypreContainer*)(*ContainerPtr);
  if (Container == NULL) {
    *structureEpochPtr = -1;
    *matrixEpochPtr = -1;
    *preconditionerEpochPtr = -1;
    return;
  }
  *structureEpochPtr = Container->structure_epoch;
  *matrixEpochPtr = Container->matrix_epoch;
  *preconditionerEpochPtr = Container->preconditioner_epoch;
}

/* Tell Elmer's matrix owner whether this container is intentionally retained
   across its wrapper-matrix teardown. */
void STDCALLBULL FC_FUNC(solvehypre5,SOLVEHYPRE5)
  ( int** ContainerPtr, int *keepAlivePtr ) {
  ElmerHypreContainer *Container = (ElmerHypreContainer*)(*ContainerPtr);
  *keepAlivePtr = (Container != NULL) ? Container->phase24_keep_alive : 0;
}

/* Exceptional recovery only: invoke a fresh numeric setup on the persistent
   fixed-pattern IJ/ParCSR objects after a failed refreshed solve. */
int STDCALLBULL FC_FUNC(solvehypre8,SOLVEHYPRE8)
  ( int** ContainerPtr, int *verbosityPtr ) {
  ElmerHypreContainer *Container = (ElmerHypreContainer*)(*ContainerPtr);
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_ParVector par_b, par_x;
  if (Container == NULL) return 1;
  /* A failed Krylov cycle can leave the assembled IJ RHS state unusable for
     the next native setup.  Reassemble the same persistent vector from the
     last RHS values; this does not recreate the IJ object. */
  if (Container->rhs_count > 0 && Container->rhs_buf != NULL) {
#if defined(HYPRE_USING_GPU)
    if (Container->hypre_use_gpu) {
      HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
      HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
    }
#endif
    HYPRE_IJVectorInitialize(Container->b);
    phase24_zero_ij_vector(Container->b);
    HYPRE_IJVectorAddToValues(Container->b, Container->rhs_count,
                              Container->rhs_index_buf, Container->rhs_buf);
    Container->ij_vector_set_values_calls++;
    HYPRE_IJVectorAssemble(Container->b);
    Container->ij_vector_assemble_calls++;
    if (Container->hypre_use_gpu) {
      Container->gpu_rhs_h2d_calls++;
      HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
    }
  }
  parcsr_A = Container->parcsr_A;
  HYPRE_IJVectorGetObject(Container->b, (void **) &par_b);
  HYPRE_IJVectorGetObject(Container->x, (void **) &par_x);
  {
    if (Container->hypre_method == 901)
      Container->preconditioner_setup_requested = 1;
    int status = RefreshHypreSetup(Container, parcsr_A, par_b, par_x);
    /* A failed FlexGMRES solve can leave its internal cycle state dirty.  A
       second native setup on the same persistent IJ/ParCSR objects is a
       cheap, allocation-preserving recovery and is still a real setup event. */
    if (status != 0) {
      int retry_status = RefreshHypreSetup(Container, parcsr_A, par_b, par_x);
      if (retry_status == 0) status = 0;
    }
    if (status == 0) {
      Container->setup_epoch = Container->matrix_epoch;
      Container->preconditioner_epoch = Container->matrix_epoch;
      Container->setup_since_last_solve = 1;
      Container->recovery_setup_calls++;
      Container->recovery_trigger_rebuild_calls++;
      Container->recovery_solve_pending = 1;
    } else {
      Container->setup_epoch = -1;
      Container->preconditioner_epoch = -1;
      fprintf(stderr, "SolveHypre: recovery native setup failed with status %d\n", status);
    }
    return status;
  }
}

/* Restore the Elmer iterate before an exceptional retry; normal solves keep
   the persistent x vector as their warm start. */
void STDCALLBULL FC_FUNC(solvehypre9,SOLVEHYPRE9)
  ( int *nrows, int *globaldofs, double *xvec, int** ContainerPtr ) {
  ElmerHypreContainer *Container = (ElmerHypreContainer*)(*ContainerPtr);
  int i;
  if (Container == NULL) return;
  if (Container->index_capacity < *nrows) {
    Container->index_buf = (int *)realloc(Container->index_buf,
                                          (size_t)(*nrows) * sizeof(int));
    Container->index_capacity = *nrows;
  }
  for (i = 0; i < *nrows; i++) Container->index_buf[i] = globaldofs[i];
  HYPRE_IJVectorSetValues(Container->x, *nrows, Container->index_buf, xvec);
  Container->ij_vector_set_values_calls++;
  HYPRE_IJVectorAssemble(Container->x);
  Container->ij_vector_assemble_calls++;
}

/* Return native counters for Phase24 diagnostics.  The array is deliberately
   plain C integers so the Fortran caller can log and assert call deltas
   without reaching into HYPRE's private solver objects.
   [0] matrix refreshes, [1] Krylov setups, [2] Boomer setups,
   [3] solves, [4] direct Case-A solves, [5] recovery setups,
   [6] IJ matrix Initialize, [7] IJ matrix SetValues,
   [8] IJ matrix Assemble, [9] IJ matrix GetObject,
   [10] IJ vector SetValues, [11] IJ vector Assemble,
   [12] IJ vector GetObject, [13] last Krylov iterations,
   [14] Case-B reused AMG solves, [15] Case-B fresh AMG solves,
   [16] age-triggered rebuilds, [17] degradation-triggered rebuilds,
   [18] recovery-triggered rebuilds, [19] current matrix epoch,
   [20] current preconditioner epoch, [21] GPU backend active,
   [22] matrix host-to-device refreshes, [23] RHS host-to-device transfers,
   [24] solution host-to-device transfers, [25] solution device-to-host syncs,
   [26] full object migrations, [27] device allocation events. */
void STDCALLBULL FC_FUNC(solvehypre10,SOLVEHYPRE10)
  ( int** ContainerPtr, int *counters ) {
  ElmerHypreContainer *Container = (ElmerHypreContainer*)(*ContainerPtr);
  int i;
  for (i = 0; i < 28; i++) counters[i] = 0;
  if (Container == NULL) return;
  counters[0] = (int)Container->matrix_refresh_calls;
  counters[1] = (int)Container->krylov_setup_calls;
  counters[2] = (int)Container->amg_setup_calls;
  counters[3] = (int)Container->solve_calls;
  counters[4] = (int)Container->case_a_direct_solve_calls;
  counters[5] = (int)Container->recovery_setup_calls;
  counters[6] = (int)Container->ij_matrix_initialize_calls;
  counters[7] = (int)Container->ij_matrix_set_values_calls;
  counters[8] = (int)Container->ij_matrix_assemble_calls;
  counters[9] = (int)Container->ij_matrix_get_object_calls;
  counters[10] = (int)Container->ij_vector_set_values_calls;
  counters[11] = (int)Container->ij_vector_assemble_calls;
  counters[12] = (int)Container->ij_vector_get_object_calls;
  counters[13] = Container->last_krylov_iterations;
  counters[14] = (int)Container->case_b_reused_solve_calls;
  counters[15] = (int)Container->case_b_fresh_solve_calls;
  counters[16] = (int)Container->preconditioner_age_forced_rebuild_calls;
  counters[17] = (int)Container->degradation_trigger_rebuild_calls;
  counters[18] = (int)Container->recovery_trigger_rebuild_calls;
  counters[19] = Container->matrix_epoch;
  counters[20] = Container->preconditioner_epoch;
  counters[21] = Container->hypre_use_gpu ? 1 : 0;
  counters[22] = (int)Container->gpu_matrix_h2d_refresh_calls;
  counters[23] = (int)Container->gpu_rhs_h2d_calls;
  counters[24] = (int)Container->gpu_solution_h2d_calls;
  counters[25] = (int)Container->gpu_solution_d2h_calls;
  counters[26] = (int)Container->gpu_full_object_migration_calls;
  counters[27] = (int)Container->gpu_device_allocation_calls;
}

/*////////////////////////////////////////////////////////////////////////////////////////////////*/

/* solve a linear system with previously constructed solver and preconditioner */
int STDCALLBULL FC_FUNC(solvehypre2,SOLVEHYPRE2)
 (
  int *nrows, int *globaldofs, int *owner,  double *xvec,
  double *rhsvec, int *Rounds, double *TOL,
  int *verbosityPtr, int** ContainerPtr, int *fcomm
 )
{

   int i, j, k;
   int myid, num_procs;
   int N, n;

   int ilower, iupper;
   int local_size, extra;

   int print_solution, print_system;

   double st, realtime_();
   
   HYPRE_Solver solver, precond;

   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_ParVector par_b;
   HYPRE_ParVector par_x;

   MPI_Comm comm=MPI_Comm_f2c(*fcomm);
   
   ElmerHypreContainer *Container;

   int verbosity = *verbosityPtr, myverb;
   int hypre_sol, hypre_pre;

   /* Keep failure telemetry defined even when a backend returns before its
      optional GetNumIterations/GetFinalRelativeResidualNorm calls. */
   int num_iterations = -1;
   int solve_status = 0;
   int solve_kind = 0;
   double final_res_norm = -1.0;
   
   Container = (ElmerHypreContainer*)(*ContainerPtr);

   /* which process number am I? */
   MPI_Comm_rank(comm, &myid); 

   if(0) myid = 0;
   
   if( myid == 0 ) 
     myverb = verbosity;
   else
     myverb = 0;

   if (myverb > 10) fprintf(stdout,"SolveHypre: solving the linear system\n");
   if (Container==NULL) {
     if(myverb) fprintf( stdout, "ID no. %i: pointer passed into SolveHypre2 is NULL, not solving",myid);
     return 1;
   }
   st = realtime_();

   if (Container->matrix_epoch < 0 || Container->preconditioner_epoch < 0 ||
       Container->preconditioner_epoch > Container->matrix_epoch) {
     if (myverb > 0)
       fprintf(stderr,
         "SolveHypre: refusing solve with invalid epochs preconditioner=%d matrix=%d\n",
         Container->preconditioner_epoch, Container->matrix_epoch);
     return 1;
   }

   parcsr_A = Container->parcsr_A;
   solver = Container->solver;
   precond = Container->precond;

   ilower = Container->ilower;
   iupper = Container->iupper;
   local_size = *nrows;

   hypre_sol = Container->hypre_method / 100;
   hypre_pre = Container->hypre_method % 100;


   /* Reuse the IJ vectors owned by the persistent container.  RHS values are
      overwritten and assembled for every solve; x is initialized from Elmer
      once, then left intact as the warm-start vector for subsequent solves. */
   if (Container->index_capacity < local_size) {
     Container->index_buf = (int *)realloc(Container->index_buf,
                                           (size_t)local_size * sizeof(int));
     Container->index_capacity = local_size;
   }
   if (Container->vec_capacity < local_size) {
     Container->vec_buf = (double *)realloc(Container->vec_buf,
                                            (size_t)local_size * sizeof(double));
     Container->vec_capacity = local_size;
   }
   if (Container->rhs_capacity < local_size) {
     Container->rhs_buf = (double *)realloc(Container->rhs_buf,
                                             (size_t)local_size * sizeof(double));
     Container->rhs_index_buf = (int *)realloc(Container->rhs_index_buf,
                                                (size_t)local_size * sizeof(int));
     Container->rhs_capacity = local_size;
   }
   for( i=0; i<local_size; i++ ) Container->index_buf[i] = globaldofs[i];
   memcpy(Container->rhs_index_buf, globaldofs, (size_t)local_size * sizeof(int));
   memcpy(Container->rhs_buf, rhsvec, (size_t)local_size * sizeof(double));
   Container->rhs_count = local_size;

#if defined(HYPRE_USING_GPU)
   if (Container->hypre_use_gpu) {
     HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
     HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
   }
#endif

   /* HYPRE 2.33's IJ interface requires Initialize() before SetValues()/
      Assemble() to modify an already-assembled vector, even when only
      replacing existing entries -- this call was missing here (every
      persistent-path solve refreshes the RHS, so this ran unconditionally
      on every one of the ~324+ solves observed before the Stage 11C-2 40ns
      window-length check crashed after 705 matrix refreshes; suspected to
      leave the IJ vector's internal assembly bookkeeping in an
      undefined/stale state that eventually corrupts a later HYPRE AMG
      setup rather than failing immediately). */
   if (Container->hypre_use_gpu && phase24_migrate_enabled()) {
     HYPRE_ParVector pv;
     HYPRE_IJVectorGetObject(Container->b, (void **) &pv);
     hypre_ParVectorMigrate((hypre_ParVector *) pv, HYPRE_MEMORY_HOST);
     HYPRE_IJVectorGetObject(Container->x, (void **) &pv);
     hypre_ParVectorMigrate((hypre_ParVector *) pv, HYPRE_MEMORY_HOST);
   }
   {
     HYPRE_Int ierr_init_b = HYPRE_IJVectorInitialize(Container->b);
     Container->ij_vector_initialize_calls++;
     if (ierr_init_b != 0 && myverb)
       fprintf(stderr, "SolveHypre: HYPRE_IJVectorInitialize(b) returned %d before RHS refresh\n",
               (int)ierr_init_b);
   }
   /* Elmer hands over the unsummed local RHS: interface rows carry only this
      rank's element contributions.  As in upstream Elmer, zero b and add every
      local row so HYPRE sums the partial values on the owner at Assemble.
      (SetValues made the last writer win and broke MPI runs.) */
   phase24_zero_ij_vector(Container->b);
   HYPRE_IJVectorAddToValues(Container->b, local_size, Container->rhs_index_buf, rhsvec);
   Container->ij_vector_set_values_calls++;
   int n_owned = 0;
   int *owned_idx = (int *)malloc((size_t)(local_size > 0 ? local_size : 1) * sizeof(int));
   double *owned_val = (double *)malloc((size_t)(local_size > 0 ? local_size : 1) * sizeof(double));
   for (i = 0; i < local_size; i++) {
     if (owner[i]) {
       owned_idx[n_owned] = globaldofs[i];
       n_owned++;
     }
   }
   HYPRE_IJVectorAssemble(Container->b);
   Container->ij_vector_assemble_calls++;

   if (!Container->vectors_initialized) {
     for (i = 0, k = 0; i < local_size; i++)
       if (owner[i]) owned_val[k++] = xvec[i];
     HYPRE_IJVectorSetValues(Container->x, n_owned, owned_idx, owned_val);
     Container->ij_vector_set_values_calls++;
     HYPRE_IJVectorAssemble(Container->x);
     Container->ij_vector_assemble_calls++;
     Container->vectors_initialized = 1;
   }
   free(owned_idx);
   free(owned_val);
   if (Container->hypre_use_gpu) HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
   HYPRE_IJVectorGetObject(Container->b, (void **) &par_b);
   Container->ij_vector_get_object_calls++;
   HYPRE_IJVectorGetObject(Container->x, (void **) &par_x);
   Container->ij_vector_get_object_calls++;
   if (Container->hypre_use_gpu && phase24_migrate_enabled()) {
     hypre_ParVectorMigrate((hypre_ParVector *) par_b, HYPRE_MEMORY_DEVICE);
     hypre_ParVectorMigrate((hypre_ParVector *) par_x, HYPRE_MEMORY_DEVICE);
   }

   if (!phase24_capture_attempted) {
     const char *capture_disabled = getenv("PHASE24_DISABLE_NATIVE_CAPTURE");
     if (capture_disabled != NULL &&
         (strcmp(capture_disabled, "1") == 0 ||
          strcmp(capture_disabled, "true") == 0 ||
          strcmp(capture_disabled, "TRUE") == 0)) {
       /* The recovery diagnostic needs lifecycle output only.  Keep the
          epoch-44 capture available by default for its dedicated replay, but
          allow lightweight runs to disable all native x/b/A capture before
          any capture-side HYPRE vector reads or files are created. */
       phase24_capture_attempted = 1;
       fprintf(stderr,
               "NATIVE_XVEC_CAPTURE_DISABLED env=PHASE24_DISABLE_NATIVE_CAPTURE\n");
       fflush(stderr);
     }
   }

   if (!phase24_capture_attempted) {
     /* Capture the first full-size 2000-iteration candidate.  This avoids
        replaying the whole 31 us campaign merely to reach the later
        epoch-44 failure; the metadata records the actual epoch and solve
        status for provenance. */
     if (local_size == 87534 && *Rounds == 2000 &&
         Container->matrix_epoch == 44) {
       Phase24CaptureBefore(Container, local_size, globaldofs, parcsr_A, myid,
                            *Rounds, *TOL, hypre_sol);
     } else if (!phase24_capture_selector_logged) {
       fprintf(stderr,
               "NATIVE_XVEC_CAPTURE_CANDIDATE_SEEN_BUT_SELECTOR_MISMATCH "
               "global_size=%d local_size=%d matrix_epoch=%d rounds=%d\n",
               Container->iupper - Container->ilower + 1, local_size,
               Container->matrix_epoch, *Rounds);
       fflush(stderr);
       phase24_capture_selector_logged = 1;
     }
   }

   Container->solve_calls++;
   if (Container->recovery_solve_pending) {
     solve_kind = 3;
     Container->recovery_solve_pending = 0;
   } else if (Container->last_solved_matrix_epoch == Container->matrix_epoch) {
     solve_kind = 1;
     Container->case_a_direct_solve_calls++;
   } else if (Container->last_solved_matrix_epoch >= 0 &&
              Container->preconditioner_epoch < Container->matrix_epoch) {
     solve_kind = 2;
     Container->case_b_reused_solve_calls++;
   } else if (Container->last_solved_matrix_epoch >= 0) {
     solve_kind = 4;
     Container->case_b_fresh_solve_calls++;
   }
   Container->last_solved_matrix_epoch = Container->matrix_epoch;
   Container->setup_since_last_solve = 0;
   HYPRE_ClearAllErrors();
   
   /* Now setup and solve! */
   if( hypre_sol == 1) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with BoomerAMG (method %d)\n",Container->hypre_method);
     solve_status = HYPRE_BoomerAMGSolve(Container->solver, parcsr_A, par_b, par_x);     
     if (myverb > 5 ) {
       HYPRE_BoomerAMGGetNumIterations(Container->solver, &num_iterations);
       HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }
   
   else if( hypre_sol == 2 ) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with AMS (method %d)\n",Container->hypre_method);
     solve_status = HYPRE_AMSSolve(Container->solver, parcsr_A, par_b, par_x);
     if (myverb > 5 ) {
       HYPRE_AMSGetNumIterations(Container->solver, &num_iterations);
       HYPRE_AMSGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }
   
   else if ( hypre_sol == 6) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with PCG (method %d)\n",Container->hypre_method);
//     HYPRE_ParCSRPCGSetMaxIter(solver, *Rounds); /* max iterations */
     solve_status = HYPRE_ParCSRPCGSolve(Container->solver, parcsr_A, par_b, par_x);

     if (myverb > 5 || getenv("PHASE24_HYPRE_REPORT") != NULL) {
       HYPRE_PCGGetNumIterations(Container->solver, &num_iterations);
       HYPRE_PCGGetFinalRelativeResidualNorm(Container->solver, &final_res_norm);
       if (getenv("PHASE24_HYPRE_REPORT") != NULL)
         fprintf(stdout, "SolveHypre: Phase24 PCG iterations=%d final_rel_res=%.3e\n",
                 (int)num_iterations, (double)final_res_norm);
     }
   }

   if ( hypre_sol == 7) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with BiCGStab (method %d)\n",Container->hypre_method);
     //     HYPRE_ParCSRBiCGSTABSetMaxIter(solver, *Rounds); /* max iterations */
     solve_status = HYPRE_ParCSRBiCGSTABSolve(Container->solver, parcsr_A, par_b, par_x);

     if (myverb > 5 ) {
       HYPRE_BiCGSTABGetNumIterations(Container->solver, &num_iterations);
       HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }

   else if ( hypre_sol == 8) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with GMRes (method %d)\n",Container->hypre_method);
//     HYPRE_GMRESSetMaxIter(solver, *Rounds); /* max GMRES iterations */
     solve_status = HYPRE_ParCSRGMRESSolve(Container->solver, parcsr_A, par_b, par_x);

     if (myverb > 5 ) {
       HYPRE_GMRESGetNumIterations(Container->solver, &num_iterations);
       HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }

   else if ( hypre_sol == 9) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with FlexGMRes (method %d)\n",Container->hypre_method);
//     HYPRE_ParCSRFlexGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     solve_status = HYPRE_ParCSRFlexGMRESSolve(Container->solver, parcsr_A, par_b, par_x);
     
     HYPRE_FlexGMRESGetNumIterations(Container->solver, &num_iterations);
     HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
   }

   else if ( hypre_sol == 10) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with LGMRes (method %d)\n",Container->hypre_method);     
//     HYPRE_ParCSRLGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     solve_status = HYPRE_ParCSRLGMRESSolve(Container->solver, parcsr_A, par_b, par_x);

     if (myverb > 5 ) {
       HYPRE_LGMRESGetNumIterations(Container->solver, &num_iterations);
       HYPRE_LGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }
   else if ( hypre_sol == 11) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with COGMRes (method %d)\n",Container->hypre_method);     
//     HYPRE_ParCSRLGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     solve_status = HYPRE_ParCSRCOGMRESSolve(Container->solver, parcsr_A, par_b, par_x);
          
     if (myverb > 5 ) {
       HYPRE_COGMRESGetNumIterations(Container->solver, &num_iterations);
       HYPRE_COGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }
   
   for( k=0,i=0; i<local_size; i++ )
     if ( owner[i] ) Container->index_buf[k++] = globaldofs[i];

   /* Return the current HYPRE iterate even on a convergence failure.  The
      native diagnostic path uses this to distinguish the Elmer warm start
      from the failed Krylov iterate before the recovery branch restores x. */
   if (Container->hypre_use_gpu) Container->gpu_solution_d2h_calls++;
   {
     int x_after_status;

     /* The diagnostic path must inspect and selectively clear HYPRE's
        convergence bit before attempting either vector getter.  The normal
        path retains the historical IJ getter behavior. */
     if (Container->hypre_use_gpu && phase24_migrate_enabled()) {
       HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
       hypre_ParVectorMigrate((hypre_ParVector *) par_x, HYPRE_MEMORY_HOST);
       hypre_ParVectorMigrate((hypre_ParVector *) par_b, HYPRE_MEMORY_HOST);
     }
     if (phase24_capture_active) {
       x_after_status = Phase24CaptureAfter(
           Container, k, Container->index_buf, Container->vec_buf, par_x,
           comm, myid, solve_status, num_iterations, final_res_norm);
     } else {
       x_after_status = HYPRE_IJVectorGetValues(Container->x, k,
                                                 Container->index_buf,
                                                 Container->vec_buf);
     }
     (void)x_after_status;
   
     for( i=0,k=0; i<local_size; i++ )
       if ( owner[i] ) xvec[i] = Container->vec_buf[k++];

     Container->last_krylov_iterations = num_iterations;
     if (solve_status != 0) {
       /* Do not let the Fortran wrapper collapse every native failure into
          the same fatal message.  In particular, a strict tolerance can set
          HYPRE's convergence/error state while still leaving a useful
          iteration count and residual available from the solver object. */
       fprintf(stderr,
               "PHASE24_HYPRE_SOLVE_FAILURE phase=backend_status method=%d "
               "solve_status=%d iterations=%d final_relative_residual=%.17e "
               "requested_tolerance=%.17e max_iterations=%d matrix_epoch=%d "
               "preconditioner_epoch=%d case=%d hypre_error=%d "
               "hypre_global_error=%d\n",
               Container->hypre_method, solve_status, num_iterations,
               final_res_norm, *TOL, *Rounds, Container->matrix_epoch,
               Container->preconditioner_epoch, solve_kind,
               HYPRE_GetError(), HYPRE_GetGlobalError(comm));
       fflush(stderr);
       return solve_status;
     }

     /* Never return a silently inaccurate Krylov iterate.  Elmer's generic
        option is advisory in some nonlinear paths, so enforce the HYPRE
        convergence contract at this boundary. */
     if (hypre_sol == 9 && num_iterations >= *Rounds && final_res_norm > *TOL) {
       fprintf(stderr, "SolveHypre: FlexGMRES failed to converge (%d iterations, relative residual %g; tolerance %g)\n",
               num_iterations, final_res_norm, *TOL);
       fprintf(stderr,
               "PHASE24_HYPRE_SOLVE_FAILURE phase=flexgmres_contract "
               "method=%d solve_status=1 iterations=%d "
               "final_relative_residual=%.17e requested_tolerance=%.17e "
               "max_iterations=%d matrix_epoch=%d preconditioner_epoch=%d "
               "case=%d hypre_error=%d hypre_global_error=%d\n",
               Container->hypre_method, num_iterations, final_res_norm,
               *TOL, *Rounds, Container->matrix_epoch,
               Container->preconditioner_epoch, solve_kind,
               HYPRE_GetError(), HYPRE_GetGlobalError(comm));
       fflush(stderr);
       return 1;
     }

     if(myverb > 4) fprintf(stdout,"SolveHypre: Required iterations %d (method %d) to norm %lg\n",
                            num_iterations,Container->hypre_method, final_res_norm);
     if (myverb > 4) fprintf( stdout, "SolveHypre: Solution time (method %d): %g\n",
                              Container->hypre_method, realtime_()-st );
     return 0;
   }
}


/* Update matrix coefficients in an existing IJ matrix and run the supported
   numeric Krylov setup while keeping IJ objects and solver allocations alive.
   For Phase24 B2, the registered AMG setup callback returns success without
   rebuilding the existing hierarchy; B1 and recovery request the real setup. */
/* Refresh the values of an assembled IJ matrix in place from the unsummed
   local CSR rows (zero the owned values, then add every local row). */
static void phase24_refresh_ij_values(HYPRE_IJMatrix M, int local_size, int *rows,
                                      int *cols, double *vals, int *globaldofs,
                                      ElmerHypreContainer *Container)
{
  int i, j, k, nnz, irow, max_nnz = 0;
  for (i = 0; i < local_size; i++) {
    nnz = rows[i+1] - rows[i];
    if (nnz > max_nnz) max_nnz = nnz;
  }
  if (Container->index_capacity < max_nnz) {
    Container->index_buf = (int *)realloc(Container->index_buf, (size_t)max_nnz * sizeof(int));
    Container->index_capacity = max_nnz;
  }
  HYPRE_IJMatrixInitialize(M);
  HYPRE_IJMatrixSetConstantValues(M, 0.0);
  for (i = 0; i < local_size; i++) {
    nnz = rows[i+1] - rows[i];
    if (nnz <= 0) continue;
    irow = globaldofs[i];
    for (k = 0, j = rows[i]; j < rows[i+1]; j++, k++)
      Container->index_buf[k] = globaldofs[cols[j-1]-1];
    HYPRE_IJMatrixAddToValues(M, 1, &nnz, &irow, Container->index_buf, &vals[rows[i]-1]);
  }
  HYPRE_IJMatrixAssemble(M);
}

/* Phase24 reuse with a separate preconditioning matrix: the Krylov matrix A
   is refreshed on every call; the preconditioning matrix and its BoomerAMG
   hierarchy only when *rebuildPtr is set.  The IJ graphs, PCG object and
   vectors are kept (fixed sparsity pattern). */
void STDCALLBULL FC_FUNC(solvehypre11,SOLVEHYPRE11)
 (
  int *nrows, int *rows, int *cols, double *vals, double *precvals,
  int *globaldofs, int *rebuildPtr, int** ContainerPtr, int *fcomm
 )
{
   ElmerHypreContainer *Container = (ElmerHypreContainer*)(*ContainerPtr);
   HYPRE_ParCSRMatrix parcsr_A, parcsr_T;
   HYPRE_ParVector par_b, par_x;
   int local_size = *nrows;
   double st, realtime_();

   if (Container == NULL) return;
   st = realtime_();
   /* Start from Elmer's current iterate as a fresh container would: the
      linear system is rescaled for every solve, so HYPRE's previous x is in
      a different scaling and is a poor warm start. */
   Container->vectors_initialized = 0;
   phase24_refresh_ij_values(Container->A, local_size, rows, cols, vals, globaldofs, Container);
   HYPRE_IJMatrixGetObject(Container->A, (void **) &parcsr_A);
   Container->parcsr_A = parcsr_A;
   Container->matrix_refresh_calls++;
   if (*rebuildPtr && Container->Atilde != Container->A) {
     double sa;
     phase24_refresh_ij_values(Container->Atilde, local_size, rows, cols, precvals, globaldofs, Container);
     HYPRE_IJMatrixGetObject(Container->Atilde, (void **) &parcsr_T);
     HYPRE_IJVectorGetObject(Container->b, (void **) &par_b);
     HYPRE_IJVectorGetObject(Container->x, (void **) &par_x);
     sa = realtime_();
     HYPRE_BoomerAMGSetup(Container->precond, parcsr_T, par_b, par_x);
     phase24_tilde_matrix = parcsr_T;
     if (getenv("PHASE24_HYPRE_REPORT") != NULL)
       fprintf(stdout, "SolveHypre: Phase24 reuse: preconditioner rebuilt (AMG setup %.3f s, total %.3f s)\n",
               realtime_() - sa, realtime_() - st);
   } else if (getenv("PHASE24_HYPRE_REPORT") != NULL) {
     fprintf(stdout, "SolveHypre: Phase24 reuse: matrix refreshed, AMG kept (%.3f s)\n",
             realtime_() - st);
   }
}

void STDCALLBULL FC_FUNC(solvehypre3,SOLVEHYPRE3)
 (
  int *nrows, int *rows, int *cols, double *vals, int *globaldofs,
  int *verbosityPtr, int *structureEpochPtr, int *matrixEpochPtr,
  int *refreshSetupPtr, int *setupReasonPtr,
  int** ContainerPtr, int *fcomm
 )
{
   int i, j, k, nnz, irow;
   MPI_Comm comm=MPI_Comm_f2c(*fcomm);
   ElmerHypreContainer *Container;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_ParVector par_b, par_x;
   int verbosity = *verbosityPtr, myverb;
   int myid;
   int local_size = *nrows;
   int max_nnz = 0;
   int refresh_setup = refreshSetupPtr ? *refreshSetupPtr : 1;
   int setup_reason = setupReasonPtr ? *setupReasonPtr : 0;

   Container = (ElmerHypreContainer*)(*ContainerPtr);
   if (Container == NULL) return;
   phase24_memory_trace_event(20, local_size, rows[local_size] - rows[0],
                              Container->matrix_epoch, Container);
   MPI_Comm_rank(comm, &myid);
   myverb = (myid == 0) ? verbosity : 0;

   for (i = 0; i < local_size; i++) {
     nnz = rows[i+1] - rows[i];
     if (nnz > max_nnz) max_nnz = nnz;
   }
   if (Container->index_capacity < max_nnz) {
     Container->index_buf = (int *)realloc(Container->index_buf,
                                           (size_t)max_nnz * sizeof(int));
     Container->index_capacity = max_nnz;
   }

#if defined(HYPRE_USING_GPU)
   if (Container->hypre_use_gpu) {
     HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
     HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
   }
#endif

   /* Re-open the assembled IJ matrix for coefficient replacement. */
   /* CORRECTED (Stage 11C-2 40ns window-length check crash investigation):
      the previous version of this comment argued that re-Initializing the
      IJ handle here would invalidate the live ParCSR object used by the
      solver/preconditioner, and omitted the call on that basis. Per the
      HYPRE 2.33 IJ interface specification, HYPRE_IJMatrixInitialize() IS
      required before SetValues()/Assemble() on an already-assembled
      matrix, even when only replacing existing fixed-pattern values -- it
      resets the matrix's internal assembly bookkeeping (not the underlying
      ParCSR data this call intentionally keeps live) for the new round of
      SetValues calls. Running this refresh path repeatedly without it
      (every accepted matrix-coefficient-change epoch; the crash was
      reproduced after 705 refreshes in a 40ns run that never happens in
      the ~50-trial bounded regressions) is suspected to leave that
      bookkeeping stale until a subsequent HYPRE AMG setup call reads it
      and crashes. HYPRE_IJMatrixGetObject() below still returns the same
      live parcsr_A either way -- Initialize() does not tear that down. */
   {
     HYPRE_Int ierr_init_A = HYPRE_IJMatrixInitialize(Container->A);
     Container->ij_matrix_initialize_calls++;
     if (ierr_init_A != 0 && myverb)
       fprintf(stderr, "SolveHypre: HYPRE_IJMatrixInitialize(A) returned %d before matrix refresh\n",
               (int)ierr_init_A);
   }
   phase24_memory_trace_event(22, local_size, rows[local_size] - rows[0],
                              Container->matrix_epoch, Container);
   /* Local rows are unsummed element contributions (interface rows are shared
      between ranks), so zero the owned values and add, as SolveHYPRE1 does. */
   HYPRE_IJMatrixSetConstantValues(Container->A, 0.0);
   for (i = 0; i < local_size; i++) {
     nnz = rows[i+1] - rows[i];
     if (nnz <= 0) continue;
     irow = globaldofs[i];
     for (k = 0, j = rows[i]; j < rows[i+1]; j++, k++)
       Container->index_buf[k] = globaldofs[cols[j-1]-1];
     HYPRE_IJMatrixAddToValues(Container->A, 1, &nnz, &irow,
                               Container->index_buf, &vals[rows[i]-1]);
     Container->ij_matrix_set_values_calls++;
   }
   HYPRE_IJMatrixAssemble(Container->A);
   Container->ij_matrix_assemble_calls++;
   phase24_memory_trace_event(23, local_size, rows[local_size] - rows[0],
                              Container->matrix_epoch, Container);
   if (Container->hypre_use_gpu) {
     Container->gpu_matrix_h2d_refresh_calls++;
     HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
   }
   HYPRE_IJMatrixGetObject(Container->A, (void **) &parcsr_A);
   Container->ij_matrix_get_object_calls++;
   Container->parcsr_A = parcsr_A;
   HYPRE_IJVectorGetObject(Container->b, (void **) &par_b);
   Container->ij_vector_get_object_calls++;
   HYPRE_IJVectorGetObject(Container->x, (void **) &par_x);
   Container->ij_vector_get_object_calls++;
   phase24_memory_trace_event(24, local_size, rows[local_size] - rows[0],
                              Container->matrix_epoch, Container);

   Container->structure_epoch = structureEpochPtr ? *structureEpochPtr : Container->structure_epoch;
   Container->matrix_epoch = matrixEpochPtr ? *matrixEpochPtr : Container->matrix_epoch;
   Container->matrix_refresh_calls++;
   if (Container->hypre_method == 901)
     Container->preconditioner_setup_requested = refresh_setup ? 1 : 0;
   phase24_memory_trace_event(25, local_size, rows[local_size] - rows[0],
                              Container->matrix_epoch, Container);
   {
     int refresh_status = RefreshHypreSetup(Container, parcsr_A, par_b, par_x);
     phase24_memory_trace_event(26, local_size, rows[local_size] - rows[0],
                                Container->matrix_epoch, Container);
     if (refresh_status != 0) {
       Container->setup_epoch = -1;
       Container->preconditioner_epoch = -1;
       Container->setup_since_last_solve = 0;
       if (myverb) fprintf(stderr,
         "SolveHypre: persistent numeric setup is unavailable for method %d\n",
         Container->hypre_method);
       return;
     }
   }
   phase24_memory_trace_event(21, local_size, rows[local_size] - rows[0],
                              Container->matrix_epoch, Container);
   if (refresh_setup) {
     if (setup_reason == 1) Container->preconditioner_age_forced_rebuild_calls++;
     if (setup_reason == 2) Container->degradation_trigger_rebuild_calls++;
     Container->setup_epoch = Container->matrix_epoch;
     Container->preconditioner_epoch = Container->matrix_epoch;
     Container->setup_since_last_solve = 1;
   } else {
     /* Keep the current ParCSR operator and run the supported Krylov setup;
        only the Phase24 AMG callback is suppressed for B2. */
     Container->setup_since_last_solve = 0;
   }
   if (myverb > 8) fprintf(stdout,
       "SolveHypre: refreshed matrix values in persistent IJ objects (%s preconditioner)\n",
       refresh_setup ? "fresh" : "reused");
}

/* destroy HYPRE data structure stored in a fortran environment */
void STDCALLBULL FC_FUNC(solvehypre4,SOLVEHYPRE4)
  ( int** ContainerPtr, int *verbosityPtr ) {
  
   ElmerHypreContainer* Container = (ElmerHypreContainer*)(*ContainerPtr);
   int verbosity = *verbosityPtr, myverb;

   int hypre_sol, hypre_pre, myid;

   if (Container==0) return;

   phase24_memory_trace_event(30, -1, -1, Container->matrix_epoch, Container);

   myverb = verbosity;
   
   hypre_sol = Container->hypre_method / 100;
   hypre_pre = Container->hypre_method % 100;

   
   if(myverb > 10 ) fprintf(stdout,"SolveHypre: Destroying Hypre solver structures!\n");

   /* Destroy the outer Krylov/solver object before the preconditioner.
      FlexGMRES retains the preconditioner callbacks and data pointer, so
      destroying BoomerAMG first can leave a dangling callback target.  Keep
      all solver-dependent objects alive until the outer solver is gone. */
   if ( hypre_sol == 1) { /* boomer AMG */
     if(myverb > 10) fprintf(stdout,"SolveHypre: Destroying BoomerAMG solver\n");
     HYPRE_BoomerAMGDestroy(Container->solver);
   }
   else if ( hypre_sol == 2 ) {
     if(myverb > 10) fprintf(stdout,"SolveHypre: Destroying AMS solver\n");
     HYPRE_AMSDestroy(Container->solver);
     if (Container->G) {
       // This leads to core dump...
       // HYPRE_IJMatrixDestroy(Container->G);
     }
   }
   else if ( hypre_sol == 6 ) {
     HYPRE_ParCSRPCGDestroy(Container->solver);
   }
   else if ( hypre_sol == 7) {
     HYPRE_ParCSRBiCGSTABDestroy(Container->solver);
   }
   else if ( hypre_sol == 8 ) {
     HYPRE_ParCSRGMRESDestroy(Container->solver);
   }
   else if ( hypre_sol == 9 ) {
     HYPRE_ParCSRFlexGMRESDestroy(Container->solver);
   }
   else if ( hypre_sol == 10 ) {
     HYPRE_ParCSRLGMRESDestroy(Container->solver);
   }
   else if ( hypre_sol == 11 ) {
     HYPRE_ParCSRCOGMRESDestroy(Container->solver);
   }

   /* Destroy Hypre preconditioner after its owning solver. */
   if ( hypre_pre == 1 ) {
     if(myverb > 10) fprintf(stdout,"SolveHypre: Destroying BoomerAMG preconditioner\n");
     HYPRE_BoomerAMGDestroy(Container->precond);
     if (Container->phase24_boomer_context)
       Container->phase24_boomer_context->boomer = NULL;
   }
   else if ( hypre_pre == 2 ) {
     if(myverb > 10) fprintf(stdout,"SolveHypre: Destroying AMS preconditioner\n");
     HYPRE_AMSDestroy(Container->precond);
     if (Container->G) {
       // This leads to core dump...
       // HYPRE_IJMatrixDestroy(Container->G);
     }
   }
   else if ( hypre_pre == 3 ) {
     HYPRE_ILUDestroy(Container->precond);
   } else if ( hypre_pre == 4 ) {
     HYPRE_ParaSailsDestroy(Container->precond);
   } else if ( hypre_pre == 12 ) {
     HYPRE_MGRDestroy(Container->precond);
     if (Container->mgr_fsolver) HYPRE_BoomerAMGDestroy(Container->mgr_fsolver);
   }

   /* IJ vectors must outlive the solver and preconditioner because their
      ParVector objects are passed into setup/solve callbacks. */
   if(Container->x) {
     HYPRE_IJVectorDestroy(Container->x);
     Container->x = NULL;
   }
   if(Container->b) {
     HYPRE_IJVectorDestroy(Container->b);
     Container->b = NULL;
   }

   if (Container->Atilde != Container->A) {
     HYPRE_IJMatrixDestroy(Container->Atilde);
   }

   if (Container->A) {
     HYPRE_IJMatrixDestroy(Container->A);
   }
   phase24_memory_trace_event(31, -1, -1, Container->matrix_epoch, Container);
   if (Container->mgr_marker) free(Container->mgr_marker);
   if (Container->index_buf) free(Container->index_buf);
   if (Container->vec_buf) free(Container->vec_buf);
   if (Container->rhs_buf) free(Container->rhs_buf);
   if (Container->rhs_index_buf) free(Container->rhs_index_buf);
   if (Container->phase24_boomer_context) free(Container->phase24_boomer_context);

   free(Container);
   *ContainerPtr = NULL;
}

/* HYPRE is process-global.  Containers are rebuilt many times during a
   transient run, so finalization must not be tied to SolveHypre4; it belongs
   after the last container has been released and before MPI_Finalize. */
void STDCALLBULL FC_FUNC(finalizehypre,FINALIZEHYPRE)(void) {
   if (HYPRE_Initialized()) {
     HYPRE_Int status = HYPRE_Finalize();
     if (status != 0)
       fprintf(stderr, "SolveHypre: HYPRE_Finalize failed with status %d\n",
               (int) status);
   }
}



void STDCALLBULL FC_FUNC(createhypreams,CREATEHYPREAMS)
 (
  int *nrows,int *rows, int *cols, double *vals, int *nnodes,
  int *grows, int *gcols, double *gvals,
  int *pirows, int *picols, double *pivals, 
  int *perm, int *invperm, int *globaldofs, int *owner,  int *globalnodes, 
  int *nodeowner, double *xvec,
  double *rhsvec, int *pe, int *ILUn, int *Rounds, double *TOL,
  double *xx_d, double *yy_d, double *zz_d, 
  int *hypre_method, int *hypre_intpara, double *hypre_dppara,
  int *verbosityPtr, int** ContainerPtr, int *fcomm
 )
{
   int i, j, k, *rcols;
   int myid, num_procs;
   int N, n;

   int ilower, iupper, nlower, nupper;
   int local_size, local_nodes, extra;

   int solver_id;
   int print_solution, print_system;

   double  *txvec, st, realtime_();
   ElmerHypreContainer *Container;
   HYPRE_ParCSRMatrix parcsr_A,parcsr_G, parcsr_Pi;
   HYPRE_IJMatrix A,G, Pi;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;
   HYPRE_IJVector xx,yy,zz;
   HYPRE_ParVector par_xx,par_yy,par_zz;

   HYPRE_Solver solver, precond;
   int verbosity = *verbosityPtr, myverb;
   MPI_Comm comm=MPI_Comm_f2c(*fcomm);

   
   Container = (ElmerHypreContainer*)malloc(sizeof(ElmerHypreContainer));   
   *ContainerPtr=(int*)(Container);
   
   st  = realtime_();
   MPI_Comm_rank(comm, &myid);

   if(myid==0)
     myverb = verbosity;
   else
     myverb = 0;
   
   if(myverb>6) fprintf(stdout,"SolveHypre: Setting up AMS preconditioner / solver\n");
   
   /* How many rows do I have? */
   local_size = *nrows;
   local_nodes = *nnodes;

   ilower=1000000000;
   iupper=0;
   for( i=0; i<local_size; i++ )
     {
       if ( owner[i] ) {
	 if ( iupper < globaldofs[i] ) iupper=globaldofs[i];
	 if ( ilower > globaldofs[i] ) ilower=globaldofs[i];
       }
     }

#if 0
   nlower=1000000000;
   nupper=0;
   for( i=0; i<local_nodes; i++ )
     {
       if ( nodeowner[i] ) {
	 if ( nupper < globalnodes[i] ) nupper=globalnodes[i];
	 if ( nlower > globalnodes[i] ) nlower=globalnodes[i];
       }
     }
#else
   nlower=1000000000;
   nupper=0;
   for( i=0; i<local_nodes; i++ )
       if  (nodeowner[i])  {
         k = 3*globalnodes[i]+2;
         if ( nupper < k ) nupper = k;
         k = 3*globalnodes[i];
         if ( nlower > k ) nlower = k;
       }
#endif
//   fprintf( stderr, "%d %d %d %d\n", ilower, iupper, nlower, nupper );

   HYPRE_IJMatrixCreate(comm, ilower, iupper, nlower, nupper, &G);
   HYPRE_IJMatrixSetObjectType(G, HYPRE_PARCSR);
   HYPRE_IJMatrixInitialize(G);
   
   {
      int nnz,irow,i,j,k,l,p,q,*rcols,csize=32;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++)
      {
         if( !owner[i] ) continue;
         nnz = grows[i+1] - grows[i];
         if ( nnz>csize ) {
           rcols = (int *)realloc( rcols, nnz*sizeof(int) );
           csize = nnz;
         }
         irow = globaldofs[i];
         for( k=0,j=grows[i]; j<grows[i+1]; j++,k++)
         {
           l = gcols[j-1]-1;
           p = l % 3;
           q = l / 3;
           rcols[k] = 3*globalnodes[q]+p;
         }
         HYPRE_IJMatrixAddToValues(G, 1, &nnz, &irow, rcols, &gvals[grows[i]-1]);
      }
      free( rcols );
   }
   
   HYPRE_IJMatrixAssemble(G);
   HYPRE_IJMatrixGetObject(G, (void**) &parcsr_G);


   nlower=1000000000;
   nupper=0;
   for( i=0; i<local_nodes; i++ )
     {
       if  (nodeowner[i])  {
         k = 3*globalnodes[i]+2;
         if ( nupper < k ) nupper = k;
         k = 3*globalnodes[i];
         if ( nlower > k ) nlower = k;
       }
     }
//   fprintf( stderr, "%d %d %d %d\n", ilower, iupper, nlower, nupper );
   HYPRE_IJMatrixCreate(comm, ilower, iupper, nlower, nupper, &Pi);
   HYPRE_IJMatrixSetObjectType(Pi, HYPRE_PARCSR);
   HYPRE_IJMatrixInitialize(Pi);
   
   {
      int nnz,irow,i,j,k,l,p,q,*rcols,csize=32;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++)
      {
         if( !owner[i] ) continue;
         nnz =  pirows[i+1] - pirows[i];
         if ( nnz>csize ) {
           rcols = (int *)realloc( rcols, nnz*sizeof(int) );
           csize = nnz;
         }
         irow = globaldofs[i];
         for( k=0,j=pirows[i]; j<pirows[i+1]; j++,k++)
         {
           l = picols[j-1]-1;
           p = l % 3;
           q = l / 3;
           rcols[k] = 3*globalnodes[q]+p;
         }
         HYPRE_IJMatrixAddToValues(Pi, 1, &nnz, &irow, rcols, &pivals[pirows[i]-1]);
      }
      free( rcols );
   }
   
   HYPRE_IJMatrixAssemble(Pi);
   HYPRE_IJMatrixGetObject(Pi, (void**) &parcsr_Pi);

   
#if 0
   for( k=0,i=0; i<local_nodes; i++ ) rcols[k++] = globalnodes[i];

   HYPRE_IJVectorCreate(comm, nlower, nupper,&xx);
   HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(xx);
   HYPRE_IJVectorSetValues(xx, local_nodes, rcols,xx_d);
   HYPRE_IJVectorAssemble(xx);
   HYPRE_IJVectorGetObject(xx, (void **) &par_xx);

   HYPRE_IJVectorCreate(comm, nlower, nupper,&yy);
   HYPRE_IJVectorSetObjectType(yy, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(yy);
   HYPRE_IJVectorSetValues(yy, local_nodes, rcols, yy_d);
   HYPRE_IJVectorAssemble(yy);
   HYPRE_IJVectorGetObject(yy, (void **) &par_yy);

   HYPRE_IJVectorCreate(comm, nlower, nupper,&zz);
   HYPRE_IJVectorSetObjectType(zz, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(zz);
   HYPRE_IJVectorSetValues(zz, local_nodes, rcols, zz_d);
   HYPRE_IJVectorAssemble(zz);
   HYPRE_IJVectorGetObject(zz, (void **) &par_zz);
#else
   HYPRE_IJVectorCreate(comm, ilower, iupper,&xx);
   HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(xx);
   HYPRE_IJVectorSetValues(xx, local_size, rcols, xx_d);
   HYPRE_IJVectorAssemble(xx);
   HYPRE_IJVectorGetObject(xx, (void **) &par_xx);

   HYPRE_IJVectorCreate(comm, ilower, iupper,&yy);
   HYPRE_IJVectorSetObjectType(yy, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(yy);
   HYPRE_IJVectorSetValues(yy, local_size, rcols, yy_d);
   HYPRE_IJVectorAssemble(yy);
   HYPRE_IJVectorGetObject(yy, (void **) &par_yy);

   HYPRE_IJVectorCreate(comm, ilower, iupper,&zz);
   HYPRE_IJVectorSetObjectType(zz, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(zz);
   HYPRE_IJVectorSetValues(zz, local_size, rcols, zz_d);
   HYPRE_IJVectorAssemble(zz);
   HYPRE_IJVectorGetObject(zz, (void **) &par_zz);
#endif

   HYPRE_AMSCreate(&precond); 
   HYPRE_AMSSetDiscreteGradient(precond,parcsr_G);
   HYPRE_AMSSetInterpolations(precond, parcsr_Pi, NULL, NULL, NULL);
//   HYPRE_AMSSetEdgeConstantVectors(precond,par_xx,par_yy,par_zz);
//   HYPRE_AMSSetCoordinateVectors(precond,par_xx,par_yy,par_zz);

   // AMS Parameters
   HYPRE_AMSSetMaxIter(precond,hypre_intpara[0]);
   HYPRE_AMSSetTol(precond,hypre_dppara[0]);
   
   HYPRE_AMSSetCycleType(precond, hypre_intpara[1]); // 1-14
   HYPRE_AMSSetSmoothingOptions(precond, hypre_intpara[2], hypre_intpara[3],
				hypre_dppara[1], hypre_dppara[2]);
   HYPRE_AMSSetAlphaAMGOptions(precond, 10, 1, 3, hypre_dppara[3], 0, 0);
   HYPRE_AMSSetBetaAMGOptions(precond, 10, 1, 3, hypre_dppara[4], 0, 0);

   if(hypre_intpara[4]) 
     HYPRE_AMSSetBetaPoissonMatrix(precond,NULL);
   
   i = (verbosity >= 6);
   if(verbosity >= 10) i=3;
   HYPRE_AMSSetPrintLevel(precond, i); 
   
   Container->precond = precond;
   Container->G = G; 
   Container->Pi = Pi; 
   
   if(myverb>5) fprintf( stdout, "SolveHypre: AMS preconditioner setup time: %g\n", realtime_()-st );
}

/* Bounded production-matrix symmetry audit.  It consumes finalized Elmer CSR
 * values, after boundary/Dirichlet elimination and before HYPRE conversion.
 * This is diagnostic evidence for optional PCG only; FlexGMRES remains the
 * safe default when SPD has not been independently established. */
void STDCALLBULL FC_FUNC(phase24csr_symmetry_metric,PHASE24CSR_SYMMETRY_METRIC)
 (int *nrows, int *rows, int *cols, double *vals, int *globaldofs, double *metric)
{
  int n=*nrows, i, p, k, cap=1, slot, row_local, col_local;
  int *keys=NULL, *values=NULL;
  double norm2=0.0, diff2=0.0;
  const int empty=-2147483647;
  *metric=0.0;
  if (n<=0) return;
  while (cap<2*n) cap<<=1;
  keys=(int*)malloc((size_t)cap*sizeof(int));
  values=(int*)malloc((size_t)cap*sizeof(int));
  if (!keys || !values) { free(keys); free(values); return; }
  for (i=0;i<cap;++i) { keys[i]=empty; values[i]=-1; }
  for (i=0;i<n;++i) {
    unsigned int h=(unsigned int)globaldofs[i]*2654435761u;
    slot=(int)(h&(unsigned int)(cap-1));
    while (keys[slot]!=empty && keys[slot]!=globaldofs[i]) slot=(slot+1)&(cap-1);
    keys[slot]=globaldofs[i]; values[slot]=i;
  }
  for (i=0;i<n;++i) {
    int start=rows[i]-1, end=rows[i+1]-1;
    for (p=start;p<end;++p) {
      double a=vals[p], at=0.0;
      norm2+=a*a;
      col_local=cols[p]-1;
      if (col_local<0 || col_local>=n) { diff2+=a*a; continue; }
      {
        unsigned int h=(unsigned int)globaldofs[col_local]*2654435761u;
        slot=(int)(h&(unsigned int)(cap-1));
        while (keys[slot]!=empty && keys[slot]!=globaldofs[col_local]) slot=(slot+1)&(cap-1);
        row_local=(keys[slot]==globaldofs[col_local]) ? values[slot] : -1;
      }
      if (row_local>=0) {
        int tstart=rows[row_local]-1, tend=rows[row_local+1]-1;
        for (k=tstart;k<tend;++k) {
          int tc=cols[k]-1;
          if (tc>=0 && tc<n && globaldofs[tc]==globaldofs[i]) { at=vals[k]; break; }
        }
      }
      diff2+=(a-at)*(a-at);
    }
  }
  if (norm2>1.0e-300) *metric=sqrt(diff2/norm2);
  free(keys); free(values);
}

#endif

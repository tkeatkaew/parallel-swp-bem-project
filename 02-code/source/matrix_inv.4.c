/*----------------------------------------------------------------------------------*/
/*--------------------------------- matrix_inv.c ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on matrix structures */
/*----------------------------------------------------------------------------------*/
#include <immintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*----------------------------------------------------------------------------------*/
#include "matrix_types.h"

#include <time.h>
#include <omp.h>
#include <x86intrin.h>

#include "cblas.h"
#include "lapacke.h"
#include "sys/time.h"
#include <sys/resource.h>
#include "openblas_config.h"

// OpenBLAS runtime query functions
extern int openblas_get_num_threads(void);
extern int openblas_get_parallel(void);

extern void dgemm_(char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);

/*----------------------------------------------------------------------------------*/
/* Global timing accumulators for matrix operations */
/*----------------------------------------------------------------------------------*/
typedef struct {
    double total_dgemm_time;
    double total_mat_inv_time;
    long total_dgemm_calls;
    long total_mat_inv_calls;
    long long total_flops;
    size_t peak_memory_bytes;
    size_t current_allocated_bytes;
} matrix_perf_stats;

static matrix_perf_stats g_perf_stats = {0};

/*----------------------------------------------------------------------------------*/
/* Memory tracking */
/*----------------------------------------------------------------------------------*/
void update_memory_stats(size_t bytes_allocated) {
    g_perf_stats.current_allocated_bytes += bytes_allocated;
    if (g_perf_stats.current_allocated_bytes > g_perf_stats.peak_memory_bytes) {
        g_perf_stats.peak_memory_bytes = g_perf_stats.current_allocated_bytes;
    }
}

void free_memory_stats(size_t bytes_freed) {
    g_perf_stats.current_allocated_bytes -= bytes_freed;
}

/*----------------------------------------------------------------------------------*/
/* Get current memory usage */
/*----------------------------------------------------------------------------------*/
void get_memory_usage_kb(long *vmrss_kb, long *vmsize_kb) {
    FILE* file = fopen("/proc/self/status", "r");
    char line[128];
    
    *vmrss_kb = 0;
    *vmsize_kb = 0;
    
    if (file) {
        while (fgets(line, 128, file) != NULL) {
            if (strncmp(line, "VmRSS:", 6) == 0) {
                sscanf(line + 6, "%ld", vmrss_kb);
            }
            if (strncmp(line, "VmSize:", 7) == 0) {
                sscanf(line + 7, "%ld", vmsize_kb);
            }
        }
        fclose(file);
    }
}

/*----------------------------------------------------------------------------------*/
/* Update statistics from multiply_matrix() calls */
/*----------------------------------------------------------------------------------*/
void update_multiply_matrix_stats(double duration, long long flops) {
    g_perf_stats.total_dgemm_time += duration;
    g_perf_stats.total_dgemm_calls++;
    g_perf_stats.total_flops += flops;
}

/*----------------------------------------------------------------------------------*/
/* Update statistics from invert_this_matrix() calls */
/*----------------------------------------------------------------------------------*/
void update_matrix_inversion_stats(double duration) {
    g_perf_stats.total_mat_inv_time += duration;
    g_perf_stats.total_mat_inv_calls++;
}

/*----------------------------------------------------------------------------------*/
/* Print performance summary */
/*----------------------------------------------------------------------------------*/
void print_matrix_performance_summary() {
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    
    long vmrss, vmsize;
    get_memory_usage_kb(&vmrss, &vmsize);
    
    printf("\n");
    printf("================================================================================\n");
    printf("                    MATRIX OPERATIONS PERFORMANCE SUMMARY\n");
    printf("================================================================================\n");
    printf("\n");
    printf("TIMING STATISTICS:\n");
    printf("------------------\n");
    printf("  Total DGEMM (matrix multiply) time:  %.6f seconds\n", g_perf_stats.total_dgemm_time);
    printf("  Total Matrix Inversion time:         %.6f seconds\n", g_perf_stats.total_mat_inv_time);
    printf("  Total computation time:              %.6f seconds\n", 
           g_perf_stats.total_dgemm_time + g_perf_stats.total_mat_inv_time);
    printf("\n");
    
    printf("OPERATION COUNTS:\n");
    printf("-----------------\n");
    printf("  DGEMM calls:                         %ld\n", g_perf_stats.total_dgemm_calls);
    printf("  Matrix Inversion calls:              %ld\n", g_perf_stats.total_mat_inv_calls);
    printf("\n");
    
    if (g_perf_stats.total_dgemm_time > 0) {
        printf("PERFORMANCE METRICS:\n");
        printf("--------------------\n");
        printf("  DGEMM GFLOPS:                        %.2f\n", 
               (double)g_perf_stats.total_flops / g_perf_stats.total_dgemm_time / 1.0e9);
        printf("  Average DGEMM time per call:         %.6f seconds\n", 
               g_perf_stats.total_dgemm_time / g_perf_stats.total_dgemm_calls);
    }
    
    if (g_perf_stats.total_mat_inv_time > 0 && g_perf_stats.total_mat_inv_calls > 0) {
        printf("  Average Matrix Inv time per call:    %.6f seconds\n", 
               g_perf_stats.total_mat_inv_time / g_perf_stats.total_mat_inv_calls);
    }
    printf("\n");
    
    printf("MEMORY USAGE:\n");
    printf("-------------\n");
    printf("  Peak allocated (tracked):            %.2f MB\n", 
           g_perf_stats.peak_memory_bytes / (1024.0 * 1024.0));
    printf("  VmRSS (resident set):                %.2f MB\n", vmrss / 1024.0);
    printf("  VmSize (virtual memory):             %.2f MB\n", vmsize / 1024.0);
    printf("  Max RSS (rusage):                    %.2f MB\n", r_usage.ru_maxrss / 1024.0);
    printf("\n");
    
    printf("OPENMP/OPENBLAS CONFIGURATION:\n");
    printf("------------------------------\n");
    printf("  OMP_NUM_THREADS:                     %s\n", getenv("OMP_NUM_THREADS") ?: "not set");
    printf("  OPENBLAS_NUM_THREADS:                %s\n", getenv("OPENBLAS_NUM_THREADS") ?: "not set");
    printf("  OMP_PROC_BIND:                       %s\n", getenv("OMP_PROC_BIND") ?: "not set");
    printf("  OMP_PLACES:                          %s\n", getenv("OMP_PLACES") ?: "not set");
    printf("\n");
    printf("================================================================================\n");
    printf("\n");
}

/*----------------------------------------------------------------------------------*/
/* get element from matrix */
/*----------------------------------------------------------------------------------*/
double get_matrix_element1(x, i, j, rows)
double *x;
int i, j, rows;
{
    double value;

    value = x[j * rows + i];

    return (value);
}

void dgemm(double *a, double *b, double *x, int a_row_num, int a_col_num, int b_row_num, int b_col_num)
{
    double value, value_new;
    int i;

    for (int row_a = 0; row_a < a_row_num; row_a++)
    {
        for (int column_b = 0; column_b < b_col_num; column_b++)
        {
            value_new = 0.0;
            for (i = 0; i < a_col_num; i++)
            {
                value = a[i * a_row_num + row_a] * b[column_b * b_row_num + i];
                value_new = value_new + value;
            }
            x[column_b * a_row_num + row_a] = value_new;
        }
    }
}

void dgemm2(double *a, double *b, double *c, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                c[i * n + j] = a[i * n + k] * b[k * n + j] + c[i * n + j];
            }
        }
    }
}

void dgemm_tpl(double *restrict a, double *restrict b, double *restrict c, int n)
{
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < n; j++)
            {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
}

static inline void gemm_tlp_simd(const double *a, const double *b, double *c, int n)
{
    for (int i = 0; i < n; ++i) {
        int j = 0;
        for (; j <= n - 4; j += 4) {
            __m256d acc = _mm256_loadu_pd(c + i*(size_t)n + j);
            for (int k = 0; k < n; ++k) {
                __m256d m1 = _mm256_set1_pd(a[i*(size_t)n + k]);
                __m256d m2 = _mm256_loadu_pd(b + k*(size_t)n + j);
                acc = _mm256_fmadd_pd(m1, m2, acc);
            }
            _mm256_storeu_pd(c + i*(size_t)n + j, acc);
        }
        for (; j < n; ++j) {
            double acc = c[i*(size_t)n + j];
            for (int k = 0; k < n; ++k)
                acc += a[i*(size_t)n + k] * b[k*(size_t)n + j];
            c[i*(size_t)n + j] = acc;
        }
    }
}

/*----------------------------------------------------------------------------------*/
/* Enhanced mat_mul with timing and memory tracking */
/*----------------------------------------------------------------------------------*/
void mat_mul(double *A, double *B, double *X, int m, int k, int n)
{
    int sizeofx = m * n;
    char ta = 'N';
    char tb = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    struct timeval start, finish;
    double duration;
    long vmrss, vmsize;

    // Track memory before operation
    get_memory_usage_kb(&vmrss, &vmsize);
    
    printf("=== DGEMM Operation ===\n");
    printf("Matrix dimensions: m=%d, n=%d, k=%d\n", m, n, k);
    printf("Memory before: VmRSS=%.2f MB, VmSize=%.2f MB\n", vmrss/1024.0, vmsize/1024.0);

    gettimeofday(&start, NULL);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
                m, n, k, alpha, A, k, B, n, beta, X, n);

    gettimeofday(&finish, NULL);

    duration = ((double)(finish.tv_sec - start.tv_sec) * 1000000 + (double)(finish.tv_usec - start.tv_usec)) / 1000000;
    
    // Calculate FLOPS
    long long flops = 2LL * m * n * k;
    double gflops = (double)flops / duration / 1.0e9;

    // Update global statistics
    g_perf_stats.total_dgemm_time += duration;
    g_perf_stats.total_dgemm_calls++;
    g_perf_stats.total_flops += flops;

    // Track memory after operation
    get_memory_usage_kb(&vmrss, &vmsize);
    
    printf("DGEMM completed in %.6f seconds (%.2f GFLOPS)\n", duration, gflops);
    printf("Memory after: VmRSS=%.2f MB, VmSize=%.2f MB\n\n", vmrss/1024.0, vmsize/1024.0);
}

/*----------------------------------------------------------------------------------*/
/* Enhanced mat_inv with timing, memory tracking, and thread diagnostics */
/*----------------------------------------------------------------------------------*/
lapack_int mat_inv(double *A, unsigned n)
{
    int ipiv[n + 1];
    lapack_int ret;
    struct timeval start, finish, start_total;
    double duration, dgetrf_time, dgetri_time;
    long vmrss_before, vmsize_before, vmrss_after, vmsize_after;
    int actual_threads = 1;  // Default to 1 if query fails

    get_memory_usage_kb(&vmrss_before, &vmsize_before);
    
    if (A == NULL) {
        fprintf(stderr, "mat_inv ERROR: A is NULL (n=%u)\n", n);
        return -4;
    }

    // ========================================================================
    // DIAGNOSTIC OUTPUT - Shows actual threading configuration
    // ========================================================================
    
    printf("=== LAPACK Matrix Inversion (DIAGNOSTIC) ===\n");
    printf("Matrix size: %u x %u\n", n, n);
    printf("\nEnvironment Variables:\n");
    printf("  OMP_NUM_THREADS:      %s\n", getenv("OMP_NUM_THREADS") ?: "not set");
    printf("  OPENBLAS_NUM_THREADS: %s\n", getenv("OPENBLAS_NUM_THREADS") ?: "not set");
    printf("  OMP_PROC_BIND:        %s\n", getenv("OMP_PROC_BIND") ?: "not set");
    printf("  OMP_PLACES:           %s\n", getenv("OMP_PLACES") ?: "not set");
    
    // Query OpenBLAS runtime configuration
    actual_threads = openblas_get_num_threads();
    int parallel_mode = openblas_get_parallel();
    
    printf("\nOpenBLAS Runtime Configuration:\n");
    printf("  Actual threads in use: %d\n", actual_threads);
    printf("  Parallel mode: %d\n", parallel_mode);
    
    // Verify thread count matches environment
    char *env_threads = getenv("OPENBLAS_NUM_THREADS");
    if (env_threads != NULL) {
        int expected_threads = atoi(env_threads);
        if (actual_threads != expected_threads) {
            printf("  ⚠️  WARNING: Mismatch! Expected %d threads, using %d\n",
                   expected_threads, actual_threads);
        } else {
            printf("  ✅ Thread count verified: %d\n", actual_threads);
        }
    }
    
    printf("\nMemory before: VmRSS=%.2f MB, VmSize=%.2f MB\n", 
           vmrss_before/1024.0, vmsize_before/1024.0);
    
    // ========================================================================
    // LAPACK Operations with Separate Timing
    // ========================================================================
    
    printf("\n--- Phase 1: LU Factorization (DGETRF) ---\n");
    gettimeofday(&start_total, NULL);
    gettimeofday(&start, NULL);

    //-----------------------------------------------------------------------
    actual_threads = openblas_get_num_threads();
    parallel_mode  = openblas_get_parallel();
    printf("OpenBLAS Runtime Configuration:\n");
    printf("  Actual threads in use: %d\n", actual_threads);
    printf("  Parallel mode: %d\n", parallel_mode);
    //-----------------------------------------------------------------------

    // LU decomposition (Phase 1)
    ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, n, ipiv);
    
    gettimeofday(&finish, NULL);
    dgetrf_time = ((double)(finish.tv_sec - start.tv_sec) * 1000000 + 
                   (double)(finish.tv_usec - start.tv_usec)) / 1000000;
    
    if (ret != 0) {
        printf("ERROR: LAPACKE_dgetrf failed with code %d\n", ret);
        return ret;
    }
    
    // DGETRF performance analysis
    double flops_dgetrf = (2.0/3.0) * (double)n * (double)n * (double)n;
    double gflops_dgetrf = flops_dgetrf / dgetrf_time / 1.0e9;
    
    printf("  Completed in:        %.6f seconds\n", dgetrf_time);
    printf("  FLOPs:               %.2e (%.0f billion)\n", flops_dgetrf, flops_dgetrf/1.0e9);
    printf("  GFLOPS:              %.2f\n", gflops_dgetrf);

    // Matrix inversion (Phase 2)
    printf("\n--- Phase 2: Matrix Inversion (DGETRI) ---\n");
    gettimeofday(&start, NULL);
    
    ret = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, A, n, ipiv);
    
    gettimeofday(&finish, NULL);
    dgetri_time = ((double)(finish.tv_sec - start.tv_sec) * 1000000 + 
                   (double)(finish.tv_usec - start.tv_usec)) / 1000000;
    
    if (ret != 0) {
        printf("ERROR: LAPACKE_dgetri failed with code %d\n", ret);
        return ret;
    }
    
    // DGETRI performance analysis
    double flops_dgetri = (4.0/3.0) * (double)n * (double)n * (double)n;
    double gflops_dgetri = flops_dgetri / dgetri_time / 1.0e9;
    
    printf("  Completed in:        %.6f seconds\n", dgetri_time);
    printf("  FLOPs:               %.2e (%.0f billion)\n", flops_dgetri, flops_dgetri/1.0e9);
    printf("  GFLOPS:              %.2f\n", gflops_dgetri);
    
    // Total duration
    duration = dgetrf_time + dgetri_time;

    // Update global statistics
    g_perf_stats.total_mat_inv_time += duration;
    g_perf_stats.total_mat_inv_calls++;

    get_memory_usage_kb(&vmrss_after, &vmsize_after);
    
    // ========================================================================
    // Performance Summary
    // ========================================================================
    
    // Total FLOPS: (2/3)*n³ + (4/3)*n³ = 2*n³
    double total_flops = flops_dgetrf + flops_dgetri;
    double gflops = total_flops / duration / 1.0e9;
    double gflops_per_thread = gflops / actual_threads;
    
    // Theoretical peak (assuming AVX2, 3.7 GHz)
    // Peak = cores × (ops/cycle) × (FMA factor) × freq
    //      = 1 × 4 (AVX2) × 2 (FMA) × 3.7 GHz = 29.6 GFLOPS per core
    double theoretical_peak_per_core = 29.6;  // GFLOPS
    double efficiency_percent = (gflops_per_thread / theoretical_peak_per_core) * 100.0;
    
    printf("\n=== Performance Summary ===\n");
    printf("Total time:              %.6f seconds\n", duration);
    printf("  - DGETRF (LU):         %.6f s (%.1f%%)\n", 
           dgetrf_time, dgetrf_time/duration*100);
    printf("  - DGETRI (inversion):  %.6f s (%.1f%%)\n", 
           dgetri_time, dgetri_time/duration*100);
    printf("\n");
    printf("FLOPs:\n");
    printf("  - DGETRF:              %.2e (%.0f billion)\n", flops_dgetrf, flops_dgetrf/1.0e9);
    printf("  - DGETRI:              %.2e (%.0f billion)\n", flops_dgetri, flops_dgetri/1.0e9);
    printf("  - Total:               %.2e (%.0f billion)\n", total_flops, total_flops/1.0e9);
    printf("\n");
    printf("Performance:\n");
    printf("  Overall GFLOPS:        %.2f\n", gflops);
    printf("  DGETRF GFLOPS:         %.2f\n", gflops_dgetrf);
    printf("  DGETRI GFLOPS:         %.2f\n", gflops_dgetri);
    printf("  GFLOPS per thread:     %.2f (using %d threads)\n", gflops_per_thread, actual_threads);
    printf("  Efficiency:            %.1f%% of theoretical peak per core\n", efficiency_percent);
    printf("\n");
    printf("Memory:\n");
    printf("  Before: VmRSS=%.2f MB, VmSize=%.2f MB\n",
           vmrss_before/1024.0, vmsize_before/1024.0);
    printf("  After:  VmRSS=%.2f MB, VmSize=%.2f MB\n", 
           vmrss_after/1024.0, vmsize_after/1024.0);
    printf("  Delta:  VmRSS=%.2f MB, VmSize=%.2f MB\n",
           (vmrss_after - vmrss_before)/1024.0, (vmsize_after - vmsize_before)/1024.0);
    printf("=============================================\n\n");
    
    return ret;
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
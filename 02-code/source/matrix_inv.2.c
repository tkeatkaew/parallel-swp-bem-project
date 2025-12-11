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
/* Enhanced mat_inv with timing and memory tracking */
/*----------------------------------------------------------------------------------*/
lapack_int mat_inv(double *A, unsigned n)
{
    int ipiv[n + 1];
    lapack_int ret;
    struct timeval start, finish;
    double duration;
    long vmrss_before, vmsize_before, vmrss_after, vmsize_after;

    get_memory_usage_kb(&vmrss_before, &vmsize_before);
    
    printf("=== LAPACK Matrix Inversion ===\n");
    printf("Matrix size: %u x %u\n", n, n);
    printf("Memory before: VmRSS=%.2f MB, VmSize=%.2f MB\n", 
           vmrss_before/1024.0, vmsize_before/1024.0);

    gettimeofday(&start, NULL);

    // LU decomposition
    ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, n, ipiv);
    if (ret != 0) {
        printf("ERROR: LAPACKE_dgetrf failed with code %d\n", ret);
        return ret;
    }

    // Matrix inversion from LU
    ret = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, A, n, ipiv);
    
    gettimeofday(&finish, NULL);
    
    duration = ((double)(finish.tv_sec - start.tv_sec) * 1000000 + 
                (double)(finish.tv_usec - start.tv_usec)) / 1000000;

    // Update global statistics
    g_perf_stats.total_mat_inv_time += duration;
    g_perf_stats.total_mat_inv_calls++;

    get_memory_usage_kb(&vmrss_after, &vmsize_after);
    
    // Estimate FLOPS for inversion: (2/3)*n^3 for LU + n^3 for inversion
    double flops = (2.0/3.0) * n * n * n + n * n * n;
    double gflops = flops / duration / 1.0e9;
    
    printf("Matrix inversion completed in %.6f seconds (%.2f GFLOPS)\n", duration, gflops);
    printf("Memory after: VmRSS=%.2f MB, VmSize=%.2f MB\n", 
           vmrss_after/1024.0, vmsize_after/1024.0);
    printf("Memory delta: VmRSS=%.2f MB, VmSize=%.2f MB\n\n",
           (vmrss_after - vmrss_before)/1024.0, (vmsize_after - vmsize_before)/1024.0);
    
    return ret;
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

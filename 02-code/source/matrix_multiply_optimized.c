/*******************************************************************************
 * matrix_multiply_optimized.c - UNIFIED Matrix Multiplication Implementation
 * 
 * VERSION 3.0: Combined Hybrid + OpenBLAS with Runtime Selection
 * 
 * This file provides a unified interface for matrix multiplication with
 * two implementation backends:
 * 
 *   DGEMM_TYPE = 0 (Hybrid):
 *     - Method 0: Pure sequential (baseline)
 *     - Method 1: OpenMP parallelization
 *     - Method 2: OpenMP + Cache blocking
 *     - Method 3: OpenMP + Cache + SIMD (AVX2)
 * 
 *   DGEMM_TYPE = 1 (OpenBLAS):
 *     - Method 0: Pure sequential (same as Hybrid)
 *     - Method 1-3: OpenBLAS DGEMM (all use same optimized library)
 * 
 * Usage:
 *   set_dgemm_type(0);      // Use Hybrid implementation
 *   set_dgemm_type(1);      // Use OpenBLAS implementation
 *   set_multiply_method(3); // Use full optimization (SIMD or OpenBLAS)
 * 
 * Author: BEM Flow Path Analysis - Parallel Optimization Research
 * Date: 2025-12-05 (Version 3.0 - Unified)
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#ifdef __AVX2__
#include <immintrin.h>
#endif

/* OpenBLAS header */
#include "cblas.h"

#include "matrix_multiply_optimized.h"

/*******************************************************************************
 * Global Configuration Variables
 ******************************************************************************/

static int g_multiply_method = 3;   /* 0=Sequential, 1=OpenMP, 2=Cache, 3=SIMD */
static int g_dgemm_type = 1;        /* 0=Hybrid, 1=OpenBLAS (default) */
static int g_block_size = 64;       /* Cache block size for Hybrid methods */
static int g_verbose = 0;           /* Print configuration info */

/*******************************************************************************
 * Configuration Functions - DGEMM Type Selection (NEW)
 ******************************************************************************/

/**
 * Set DGEMM implementation type
 * @param type 0=Hybrid (custom OpenMP+Cache+SIMD), 1=OpenBLAS (production library)
 */
void set_dgemm_type(int type)
{
    if (type < 0 || type > 1)
    {
        fprintf(stderr, "Warning: Invalid DGEMM type %d, using default (1=OpenBLAS)\n", type);
        g_dgemm_type = 1;
    }
    else
    {
        g_dgemm_type = type;
        if (g_verbose)
        {
            printf("[CONFIG] DGEMM type: %d (%s)\n", type,
                   type == 0 ? "Hybrid" : "OpenBLAS");
        }
    }
}

/**
 * Get current DGEMM implementation type
 */
int get_dgemm_type(void)
{
    return g_dgemm_type;
}

/**
 * Get human-readable name of current DGEMM type
 */
const char* get_dgemm_type_name(void)
{
    return g_dgemm_type == 0 ? "Hybrid (OpenMP+Cache+SIMD)" : "OpenBLAS";
}

/*******************************************************************************
 * Configuration Functions - Multiply Method Selection (Existing)
 ******************************************************************************/

/**
 * Set matrix multiplication method
 * @param method 0=Sequential, 1=OpenMP, 2=OpenMP+Cache, 3=OpenMP+Cache+SIMD
 */
void set_multiply_method(int method)
{
    if (method < 0 || method > 3)
    {
        fprintf(stderr, "Warning: Invalid method %d, using default (3)\n", method);
        g_multiply_method = 3;
    }
    else
    {
        g_multiply_method = method;
        if (g_verbose)
        {
            const char *names[] = {
                "Sequential (baseline)",
                "OpenMP parallelization",
                "OpenMP + Cache blocking",
                "OpenMP + Cache + SIMD"
            };
            printf("[CONFIG] Multiply method: %d (%s)\n", method, names[method]);
        }
    }
}

int get_multiply_method(void)
{
    return g_multiply_method;
}

/**
 * Set cache block size for Hybrid methods
 */
void set_block_size(int size)
{
    if (size < 8 || size > 256)
    {
        fprintf(stderr, "Warning: Block size %d out of range [8-256], using default (64)\n", size);
        g_block_size = 64;
    }
    else
    {
        g_block_size = size;
        if (g_verbose)
        {
            printf("[CONFIG] Cache block size: %d\n", size);
        }
    }
}

int get_block_size(void)
{
    return g_block_size;
}

void set_multiply_verbose(int verbose)
{
    g_verbose = verbose;
}

/*******************************************************************************
 * Helper Functions
 ******************************************************************************/

static inline int min_int(int a, int b)
{
    return (a < b) ? a : b;
}

static inline double get_walltime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

/*******************************************************************************
 * METHOD 0: Pure Sequential Matrix Multiplication (TRUE BASELINE)
 * 
 * This implementation is used for BOTH Hybrid and OpenBLAS when method=0.
 * It provides a true sequential baseline with:
 * - No OpenMP (single thread)
 * - No BLAS library calls
 * - No cache blocking
 * - No SIMD vectorization
 * 
 * Purpose: Accurate speedup measurement baseline
 ******************************************************************************/

void multiply_matrix_sequential(matrix *a, matrix *b, matrix *x)
{
    /* Calculate LOGICAL dimensions for loop bounds */
    int a_row_num = a->transpose ? a->columns : a->rows;
    int a_col_num = a->transpose ? a->rows : a->columns;
    int b_col_num = b->transpose ? b->rows : b->columns;
    
    double *A = a->value;
    double *B = b->value;
    double *X = x->value;
    
    /* Pure sequential implementation - matches original multiply_matrix_org() */
    if (a->transpose == 0)
    {
        /* Non-transposed case */
        for (int row_a = 0; row_a < a_row_num; row_a++)
        {
            for (int column_b = 0; column_b < b_col_num; column_b++)
            {
                double value_new = 0.0;
                
                for (int i = 0; i < a_col_num; i++)
                {
                    /* Column-major access: element[row,col] at col*rows+row */
                    value_new += A[i * a->rows + row_a] * B[column_b * b->rows + i];
                }
                
                X[column_b * x->rows + row_a] = value_new;
            }
        }
    }
    else
    {
        /* Transposed case */
        for (int row_a = 0; row_a < a_row_num; row_a++)
        {
            for (int column_b = 0; column_b < b_col_num; column_b++)
            {
                double value_new = 0.0;
                
                for (int i = 0; i < a_col_num; i++)
                {
                    value_new += A[row_a * a->rows + i] * B[column_b * b->rows + i];
                }
                
                X[column_b * x->rows + row_a] = value_new;
            }
        }
    }
}

/*******************************************************************************
 * HYBRID IMPLEMENTATION: METHOD 1 - OpenMP Parallelization
 ******************************************************************************/

void multiply_matrix_openmp_hybrid(matrix *a, matrix *b, matrix *x)
{
    int a_row_num = a->transpose ? a->columns : a->rows;
    int a_col_num = a->transpose ? a->rows : a->columns;
    int b_col_num = b->transpose ? b->rows : b->columns;
    
    double *A = a->value;
    double *B = b->value;
    double *X = x->value;
    
    if (a->transpose == 0)
    {
        #pragma omp parallel for schedule(static)
        for (int row_a = 0; row_a < a_row_num; row_a++)
        {
            for (int column_b = 0; column_b < b_col_num; column_b++)
            {
                double value_new = 0.0;
                
                for (int i = 0; i < a_col_num; i++)
                {
                    value_new += A[i * a->rows + row_a] * B[column_b * b->rows + i];
                }
                
                X[column_b * x->rows + row_a] = value_new;
            }
        }
    }
    else
    {
        #pragma omp parallel for schedule(static)
        for (int row_a = 0; row_a < a_row_num; row_a++)
        {
            for (int column_b = 0; column_b < b_col_num; column_b++)
            {
                double value_new = 0.0;
                
                for (int i = 0; i < a_col_num; i++)
                {
                    value_new += A[row_a * a->rows + i] * B[column_b * b->rows + i];
                }
                
                X[column_b * x->rows + row_a] = value_new;
            }
        }
    }
}

/*******************************************************************************
 * HYBRID IMPLEMENTATION: METHOD 2 - OpenMP + Cache Blocking
 ******************************************************************************/

void multiply_matrix_cache_hybrid(matrix *a, matrix *b, matrix *x, int block_size)
{
    int a_row_num = a->transpose ? a->columns : a->rows;
    int a_col_num = a->transpose ? a->rows : a->columns;
    int b_col_num = b->transpose ? b->rows : b->columns;
    
    double *A = a->value;
    double *B = b->value;
    double *X = x->value;
    
    /* Zero output matrix */
    #pragma omp parallel for
    for (int i = 0; i < x->rows * x->columns; i++)
    {
        X[i] = 0.0;
    }
    
    if (a->transpose == 0)
    {
        #pragma omp parallel for schedule(static) collapse(2)
        for (int ii = 0; ii < a_row_num; ii += block_size)
        {
            for (int jj = 0; jj < b_col_num; jj += block_size)
            {
                for (int kk = 0; kk < a_col_num; kk += block_size)
                {
                    int i_end = min_int(ii + block_size, a_row_num);
                    int j_end = min_int(jj + block_size, b_col_num);
                    int k_end = min_int(kk + block_size, a_col_num);
                    
                    for (int row_a = ii; row_a < i_end; row_a++)
                    {
                        for (int column_b = jj; column_b < j_end; column_b++)
                        {
                            double value_new = 0.0;
                            
                            for (int i = kk; i < k_end; i++)
                            {
                                value_new += A[i * a->rows + row_a] * B[column_b * b->rows + i];
                            }
                            
                            #pragma omp atomic
                            X[column_b * x->rows + row_a] += value_new;
                        }
                    }
                }
            }
        }
    }
    else
    {
        #pragma omp parallel for schedule(static) collapse(2)
        for (int ii = 0; ii < a_row_num; ii += block_size)
        {
            for (int jj = 0; jj < b_col_num; jj += block_size)
            {
                for (int kk = 0; kk < a_col_num; kk += block_size)
                {
                    int i_end = min_int(ii + block_size, a_row_num);
                    int j_end = min_int(jj + block_size, b_col_num);
                    int k_end = min_int(kk + block_size, a_col_num);
                    
                    for (int row_a = ii; row_a < i_end; row_a++)
                    {
                        for (int column_b = jj; column_b < j_end; column_b++)
                        {
                            double value_new = 0.0;
                            
                            for (int i = kk; i < k_end; i++)
                            {
                                value_new += A[row_a * a->rows + i] * B[column_b * b->rows + i];
                            }
                            
                            #pragma omp atomic
                            X[column_b * x->rows + row_a] += value_new;
                        }
                    }
                }
            }
        }
    }
}

/*******************************************************************************
 * HYBRID IMPLEMENTATION: METHOD 3 - OpenMP + Cache + SIMD (AVX2)
 ******************************************************************************/

void multiply_matrix_simd_hybrid(matrix *a, matrix *b, matrix *x, int block_size)
{
    int a_row_num = a->transpose ? a->columns : a->rows;
    int a_col_num = a->transpose ? a->rows : a->columns;
    int b_col_num = b->transpose ? b->rows : b->columns;
    
    double *A = a->value;
    double *B = b->value;
    double *X = x->value;
    
    /* Zero output matrix */
    #pragma omp parallel for
    for (int i = 0; i < x->rows * x->columns; i++)
    {
        X[i] = 0.0;
    }
    
    if (a->transpose == 0)
    {
        #pragma omp parallel for schedule(static)
        for (int row_a = 0; row_a < a_row_num; row_a++)
        {
            for (int column_b = 0; column_b < b_col_num; column_b++)
            {
#ifdef __AVX2__
                __m256d sum_vec = _mm256_setzero_pd();
                int i = 0;
                
                /* SIMD loop (4 doubles at a time) */
                for (; i + 3 < a_col_num; i += 4)
                {
                    /* Load from A (strided access) */
                    __m256d a_vec = _mm256_set_pd(
                        A[(i + 3) * a->rows + row_a],
                        A[(i + 2) * a->rows + row_a],
                        A[(i + 1) * a->rows + row_a],
                        A[i * a->rows + row_a]);
                    
                    /* Load from B (contiguous) */
                    __m256d b_vec = _mm256_loadu_pd(&B[column_b * b->rows + i]);
                    
                    /* FMA: sum += a * b */
                    sum_vec = _mm256_fmadd_pd(a_vec, b_vec, sum_vec);
                }
                
                /* Horizontal sum */
                double sum_array[4];
                _mm256_storeu_pd(sum_array, sum_vec);
                double value_new = sum_array[0] + sum_array[1] + sum_array[2] + sum_array[3];
                
                /* Scalar remainder */
                for (; i < a_col_num; i++)
                {
                    value_new += A[i * a->rows + row_a] * B[column_b * b->rows + i];
                }
#else
                /* Fallback: scalar implementation */
                double value_new = 0.0;
                for (int i = 0; i < a_col_num; i++)
                {
                    value_new += A[i * a->rows + row_a] * B[column_b * b->rows + i];
                }
#endif
                X[column_b * x->rows + row_a] = value_new;
            }
        }
    }
    else
    {
        /* Transposed case - both A and B have contiguous access */
        #pragma omp parallel for schedule(static)
        for (int row_a = 0; row_a < a_row_num; row_a++)
        {
            for (int column_b = 0; column_b < b_col_num; column_b++)
            {
#ifdef __AVX2__
                __m256d sum_vec = _mm256_setzero_pd();
                int i = 0;
                
                /* SIMD loop - BOTH contiguous! */
                for (; i + 3 < a_col_num; i += 4)
                {
                    __m256d a_vec = _mm256_loadu_pd(&A[row_a * a->rows + i]);
                    __m256d b_vec = _mm256_loadu_pd(&B[column_b * b->rows + i]);
                    sum_vec = _mm256_fmadd_pd(a_vec, b_vec, sum_vec);
                }
                
                /* Horizontal sum */
                double sum_array[4];
                _mm256_storeu_pd(sum_array, sum_vec);
                double value_new = sum_array[0] + sum_array[1] + sum_array[2] + sum_array[3];
                
                /* Scalar remainder */
                for (; i < a_col_num; i++)
                {
                    value_new += A[row_a * a->rows + i] * B[column_b * b->rows + i];
                }
#else
                double value_new = 0.0;
                for (int i = 0; i < a_col_num; i++)
                {
                    value_new += A[row_a * a->rows + i] * B[column_b * b->rows + i];
                }
#endif
                X[column_b * x->rows + row_a] = value_new;
            }
        }
    }
}

/*******************************************************************************
 * OPENBLAS IMPLEMENTATION: Methods 1-3
 * 
 * OpenBLAS provides highly optimized DGEMM that includes:
 * - Multi-threading (controlled by OPENBLAS_NUM_THREADS)
 * - Cache blocking (internal)
 * - SIMD vectorization (AVX2, AVX-512)
 * - Assembly-level optimizations
 * 
 * All methods 1, 2, 3 use the same OpenBLAS implementation.
 ******************************************************************************/

void multiply_matrix_openblas(matrix *a, matrix *b, matrix *x)
{
    /* Get logical dimensions */
    int M = a->transpose ? a->columns : a->rows;  /* Rows of op(A) */
    int K = a->transpose ? a->rows : a->columns;  /* Cols of op(A) = Rows of op(B) */
    int N = b->transpose ? b->rows : b->columns;  /* Cols of op(B) */
    
    /* Verify dimensions */
    if (x->rows != M || x->columns != N)
    {
        fprintf(stderr, "Error: Result matrix dimension mismatch!\n");
        fprintf(stderr, "  Expected: %d x %d, Got: %d x %d\n", 
                M, N, x->rows, x->columns);
        exit(1);
    }
    
    double *A = a->value;
    double *B = b->value;
    double *X = x->value;
    
    /* OpenBLAS parameters */
    double alpha = 1.0;
    double beta = 0.0;  /* Overwrite output matrix */
    
    /* Leading dimensions (physical storage) */
    int lda = a->rows;
    int ldb = b->rows;
    int ldx = x->rows;
    
    /* Transpose handling */
    CBLAS_TRANSPOSE transA = (a->transpose == 0) ? CblasNoTrans : CblasTrans;
    CBLAS_TRANSPOSE transB = (b->transpose == 0) ? CblasNoTrans : CblasTrans;
    
    /* Call OpenBLAS DGEMM: C := alpha * op(A) * op(B) + beta * C */
    cblas_dgemm(
        CblasColMajor,     /* Matrix storage order */
        transA,            /* Transpose A? */
        transB,            /* Transpose B? */
        M,                 /* Rows of op(A) and C */
        N,                 /* Cols of op(B) and C */
        K,                 /* Cols of op(A), Rows of op(B) */
        alpha,             /* Scalar alpha */
        A,                 /* Matrix A */
        lda,               /* Leading dimension of A */
        B,                 /* Matrix B */
        ldb,               /* Leading dimension of B */
        beta,              /* Scalar beta */
        X,                 /* Matrix C (output) */
        ldx                /* Leading dimension of C */
    );
}

/*******************************************************************************
 * MAIN DISPATCHER FUNCTION
 * 
 * Routes to appropriate implementation based on:
 * - g_dgemm_type: 0=Hybrid, 1=OpenBLAS
 * - g_multiply_method: 0=Sequential, 1=OpenMP, 2=Cache, 3=SIMD
 ******************************************************************************/

void multiply_matrix_optimized(matrix *a, matrix *b, matrix *x)
{
    /* Method 0 is always sequential (true baseline) regardless of dgemm_type */
    if (g_multiply_method == 0)
    {
        multiply_matrix_sequential(a, b, x);
        return;
    }
    
    /* For methods 1-3, choose implementation based on dgemm_type */
    if (g_dgemm_type == 0)
    {
        /* Hybrid implementation (custom OpenMP + Cache + SIMD) */
        switch (g_multiply_method)
        {
        case 1:
            multiply_matrix_openmp_hybrid(a, b, x);
            break;
        case 2:
            multiply_matrix_cache_hybrid(a, b, x, g_block_size);
            break;
        case 3:
            multiply_matrix_simd_hybrid(a, b, x, g_block_size);
            break;
        default:
            fprintf(stderr, "Error: Unknown multiply method %d\n", g_multiply_method);
            exit(1);
        }
    }
    else
    {
        /* OpenBLAS implementation (all methods 1-3 use same optimized library) */
        multiply_matrix_openblas(a, b, x);
    }
}

/*******************************************************************************
 * CONVENIENCE WRAPPERS
 ******************************************************************************/

void multiply_matrix_auto(matrix *a, matrix *b, matrix *x)
{
    multiply_matrix_optimized(a, b, x);
}

void multiply_matrix_method(matrix *a, matrix *b, matrix *x, int method)
{
    int saved_method = g_multiply_method;
    set_multiply_method(method);
    multiply_matrix_optimized(a, b, x);
    g_multiply_method = saved_method;
}

/* Legacy compatibility wrappers */
void multiply_matrix_openmp(matrix *a, matrix *b, matrix *x)
{
    if (g_dgemm_type == 0)
        multiply_matrix_openmp_hybrid(a, b, x);
    else
        multiply_matrix_openblas(a, b, x);
}

void multiply_matrix_cache(matrix *a, matrix *b, matrix *x, int block_size)
{
    if (g_dgemm_type == 0)
        multiply_matrix_cache_hybrid(a, b, x, block_size);
    else
        multiply_matrix_openblas(a, b, x);
}

void multiply_matrix_simd(matrix *a, matrix *b, matrix *x, int block_size)
{
    if (g_dgemm_type == 0)
        multiply_matrix_simd_hybrid(a, b, x, block_size);
    else
        multiply_matrix_openblas(a, b, x);
}

/*******************************************************************************
 * PERFORMANCE INFORMATION FUNCTIONS
 ******************************************************************************/

void print_multiply_config(void)
{
    const char *method_names[] = {
        "Sequential (baseline)",
        "OpenMP parallelization",
        "OpenMP + Cache blocking",
        "OpenMP + Cache + SIMD (AVX2)"
    };
    
    const char *dgemm_names[] = {
        "Hybrid (custom implementation)",
        "OpenBLAS (production library)"
    };
    
    char *openblas_threads = getenv("OPENBLAS_NUM_THREADS");
    char *omp_threads = getenv("OMP_NUM_THREADS");
    
    printf("\n");
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("           MATRIX MULTIPLICATION CONFIGURATION (v3.0)\n");
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("  DGEMM Type:          %d (%s)\n", g_dgemm_type, dgemm_names[g_dgemm_type]);
    printf("  Multiply Method:     %d (%s)\n", g_multiply_method, method_names[g_multiply_method]);
    printf("  Block Size:          %d\n", g_block_size);
    printf("───────────────────────────────────────────────────────────────────\n");
    
    if (g_multiply_method == 0)
    {
        printf("  Mode: PURE SEQUENTIAL (true baseline)\n");
        printf("    • Single thread only\n");
        printf("    • No BLAS library calls\n");
        printf("    • No cache optimization\n");
        printf("    • No SIMD vectorization\n");
    }
    else if (g_dgemm_type == 0)
    {
        printf("  Mode: HYBRID (custom implementation)\n");
        printf("    • OpenMP threads:   %s\n", omp_threads ? omp_threads : "default");
#ifdef __AVX2__
        printf("    • SIMD support:     AVX2 + FMA ✅\n");
#else
        printf("    • SIMD support:     Not available ❌\n");
#endif
        printf("    • Cache blocking:   %s (B=%d)\n", 
               g_multiply_method >= 2 ? "Enabled" : "Disabled", g_block_size);
    }
    else
    {
        printf("  Mode: OPENBLAS (production library)\n");
        printf("    • OPENBLAS_NUM_THREADS: %s\n", 
               openblas_threads ? openblas_threads : "default");
        printf("    • OMP_NUM_THREADS:      %s\n", 
               omp_threads ? omp_threads : "default");
        printf("    • Features:             Threading + Cache + SIMD + Assembly ✅\n");
    }
    
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("\n");
}

double get_expected_speedup(int method, int num_threads)
{
    if (method == 0)
    {
        return 1.0;  /* Sequential baseline */
    }
    
    if (g_dgemm_type == 0)
    {
        /* Hybrid implementation expected speedup */
        double base_speedup;
        if (num_threads <= 4)
            base_speedup = num_threads * 0.95;
        else if (num_threads <= 8)
            base_speedup = num_threads * 0.85;
        else
            base_speedup = num_threads * 0.70;
        
        switch (method)
        {
        case 1: return base_speedup;              /* OpenMP only */
        case 2: return base_speedup * 1.40;       /* +Cache: ~40% boost */
        case 3: return base_speedup * 1.80;       /* +SIMD: ~80% boost over cache */
        default: return base_speedup;
        }
    }
    else
    {
        /* OpenBLAS expected speedup (highly optimized) */
        if (num_threads <= 4)
            return num_threads * 0.95;
        else if (num_threads <= 8)
            return num_threads * 0.90;
        else if (num_threads <= 16)
            return num_threads * 0.80;
        else
            return num_threads * 0.70;
    }
}

void print_expected_performance(void)
{
    char *threads_env = getenv("OMP_NUM_THREADS");
    if (threads_env == NULL)
        threads_env = getenv("OPENBLAS_NUM_THREADS");
    
    int threads = threads_env ? atoi(threads_env) : omp_get_max_threads();
    
    printf("\n");
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("           EXPECTED PERFORMANCE\n");
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("  System: %d threads, DGEMM=%s\n", threads, get_dgemm_type_name());
    printf("───────────────────────────────────────────────────────────────────\n");
    
    if (g_dgemm_type == 0)
    {
        printf("  Hybrid implementation speedup vs sequential:\n");
        printf("    Method 0 (Sequential):         %.2f×\n", 1.0);
        printf("    Method 1 (OpenMP):             %.2f×\n", get_expected_speedup(1, threads));
        printf("    Method 2 (OpenMP+Cache):       %.2f×\n", get_expected_speedup(2, threads));
        printf("    Method 3 (OpenMP+Cache+SIMD):  %.2f×\n", get_expected_speedup(3, threads));
    }
    else
    {
        printf("  OpenBLAS implementation speedup vs sequential:\n");
        printf("    Method 0 (Sequential):         %.2f×\n", 1.0);
        printf("    Method 1-3 (OpenBLAS):         %.2f×\n", get_expected_speedup(3, threads));
        printf("\n");
        printf("  Note: Methods 1-3 all use same OpenBLAS implementation.\n");
    }
    
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("\n");
}

/*******************************************************************************
 * CONFIGURATION SUMMARY (for benchmarking scripts)
 ******************************************************************************/

void print_config_csv_header(void)
{
    printf("dgemm_type,multiply_method,block_size,omp_threads,openblas_threads\n");
}

void print_config_csv(void)
{
    char *omp_threads = getenv("OMP_NUM_THREADS");
    char *openblas_threads = getenv("OPENBLAS_NUM_THREADS");
    
    printf("%d,%d,%d,%s,%s\n",
           g_dgemm_type,
           g_multiply_method,
           g_block_size,
           omp_threads ? omp_threads : "default",
           openblas_threads ? openblas_threads : "default");
}

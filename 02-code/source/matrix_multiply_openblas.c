/*******************************************************************************
 * matrix_multiply_openblas.c - OpenBLAS-based Matrix Multiplication
 * 
 * This implementation replaces custom Level 1-3 optimizations with OpenBLAS
 * DGEMM calls for maximum performance and maintainability.
 * 
 * Key Features:
 * - Method 0: Sequential (baseline, uses naive implementation)
 * - Method 1-3: All use OpenBLAS DGEMM (highly optimized)
 * - Proper handling of transposed matrices
 * - Performance tracking and timing
 * - Multi-threaded BLAS control via environment variables
 * 
 * Author: Optimized for BEM Flow Path Analysis
 * Date: 2025-11-21
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#include "cblas.h"
#include "matrix_multiply_optimized.h"

/*******************************************************************************
 * Global Configuration
 ******************************************************************************/

static int g_multiply_method = 3;  // Default: OpenBLAS (all methods 1-3 use BLAS)
static int g_verbose = 0;          // Print performance info

/*******************************************************************************
 * Configuration Functions
 ******************************************************************************/

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
                "OpenBLAS DGEMM", 
                "OpenBLAS DGEMM (same as Method 1)", 
                "OpenBLAS DGEMM (same as Method 1)"
            };
            printf("[CONFIG] Matrix multiply method: %d (%s)\n", method, names[method]);
        }
    }
}

int get_multiply_method(void)
{
    return g_multiply_method;
}

void set_block_size(int size)
{
    // Block size is not used with OpenBLAS, but keep for API compatibility
    if (g_verbose)
    {
        printf("[CONFIG] Block size parameter ignored (using OpenBLAS)\n");
    }
}

int get_block_size(void)
{
    return 0; // Not applicable for OpenBLAS
}

void set_multiply_verbose(int verbose)
{
    g_verbose = verbose;
}

/*******************************************************************************
 * Helper Functions
 ******************************************************************************/

static inline double get_walltime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

/*******************************************************************************
 * METHOD 0: Sequential Matrix Multiplication (Baseline)
 * 
 * This is a simple triple-loop implementation for baseline comparison.
 * Does NOT use OpenBLAS - pure sequential code for accurate speedup measurement.
 ******************************************************************************/

void multiply_matrix_sequential(matrix *a, matrix *b, matrix *x)
{
    // Calculate LOGICAL dimensions for loop bounds
    int a_row_num = a->transpose ? a->columns : a->rows;
    int a_col_num = a->transpose ? a->rows : a->columns;
    int b_row_num = b->transpose ? b->columns : b->rows;
    int b_col_num = b->transpose ? b->rows : b->columns;
    
    double *A = a->value;
    double *B = b->value;
    double *X = x->value;
    
    // EXACT SAME LOGIC as original multiply_matrix_org()
    if (a->transpose == 0)
    {
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
    { // a->transpose == 1
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
 * METHOD 1-3: OpenBLAS DGEMM (All optimizations use the same implementation)
 * 
 * OpenBLAS provides highly optimized DGEMM that includes:
 * - Multi-threading (controlled by OPENBLAS_NUM_THREADS)
 * - Cache blocking
 * - SIMD vectorization (AVX2, AVX-512)
 * - Assembly-level optimizations
 * 
 * This single implementation is used for methods 1, 2, and 3.
 ******************************************************************************/

void multiply_matrix_openblas(matrix *a, matrix *b, matrix *x)
{
    // Get logical dimensions
    int M = a->transpose ? a->columns : a->rows;  // Rows of A (logical)
    int K = a->transpose ? a->rows : a->columns;  // Cols of A = Rows of B (logical)
    int N = b->transpose ? b->rows : b->columns;  // Cols of B (logical)
    
    // Result matrix dimensions (should match x->rows and x->columns)
    int result_rows = M;
    int result_cols = N;
    
    if (x->rows != result_rows || x->columns != result_cols)
    {
        fprintf(stderr, "Error: Result matrix dimension mismatch!\n");
        fprintf(stderr, "  Expected: %d x %d, Got: %d x %d\n", 
                result_rows, result_cols, x->rows, x->columns);
        exit(1);
    }
    
    double *A = a->value;
    double *B = b->value;
    double *X = x->value;
    
    // OpenBLAS parameters
    double alpha = 1.0;
    double beta = 0.0;  // Overwrite output matrix
    
    // Leading dimensions (physical storage)
    int lda = a->rows;  // Physical row dimension of A
    int ldb = b->rows;  // Physical row dimension of B
    int ldx = x->rows;  // Physical row dimension of X
    
    /*
     * CBLAS_ORDER: We use CblasColMajor because matrices are stored column-major
     * 
     * Matrix storage convention:
     * - Physical storage: column-major (Fortran style)
     * - Element A[i,j] stored at: A[j*rows + i]
     * 
     * Transpose handling:
     * - If a->transpose == 0: Use CblasNoTrans
     * - If a->transpose == 1: Use CblasTrans
     */
    
    CBLAS_TRANSPOSE transA = (a->transpose == 0) ? CblasNoTrans : CblasTrans;
    CBLAS_TRANSPOSE transB = (b->transpose == 0) ? CblasNoTrans : CblasTrans;
    
    // Call OpenBLAS DGEMM
    // C := alpha * op(A) * op(B) + beta * C
    // X := 1.0 * op(A) * op(B) + 0.0 * X
    cblas_dgemm(
        CblasColMajor,     // Matrix storage order
        transA,            // Transpose A?
        transB,            // Transpose B?
        M,                 // Rows of op(A) and C
        N,                 // Cols of op(B) and C
        K,                 // Cols of op(A), Rows of op(B)
        alpha,             // Scalar alpha
        A,                 // Matrix A
        lda,               // Leading dimension of A
        B,                 // Matrix B
        ldb,               // Leading dimension of B
        beta,              // Scalar beta
        X,                 // Matrix C (output)
        ldx                // Leading dimension of C
    );
}

/*******************************************************************************
 * Convenience aliases for methods 1-3
 * (All call the same OpenBLAS implementation)
 ******************************************************************************/

void multiply_matrix_openmp(matrix *a, matrix *b, matrix *x)
{
    multiply_matrix_openblas(a, b, x);
}

void multiply_matrix_cache(matrix *a, matrix *b, matrix *x, int block_size)
{
    // Ignore block_size parameter - OpenBLAS handles blocking internally
    multiply_matrix_openblas(a, b, x);
}

void multiply_matrix_simd(matrix *a, matrix *b, matrix *x, int block_size)
{
    // Ignore block_size parameter - OpenBLAS handles SIMD internally
    multiply_matrix_openblas(a, b, x);
}

/*******************************************************************************
 * Main Dispatcher Function
 ******************************************************************************/

void multiply_matrix_optimized(matrix *a, matrix *b, matrix *x)
{
    // Dispatch to appropriate method
    switch (g_multiply_method)
    {
    case 0:
        multiply_matrix_sequential(a, b, x);
        break;
        
    case 1:
    case 2:
    case 3:
        // All optimized methods use OpenBLAS
        multiply_matrix_openblas(a, b, x);
        break;
        
    default:
        fprintf(stderr, "Error: Unknown multiply method %d\n", g_multiply_method);
        exit(1);
    }
}

/*******************************************************************************
 * Convenience wrappers
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

/*******************************************************************************
 * Performance information
 ******************************************************************************/

void print_multiply_config(void)
{
    const char *method_names[] = {
        "Sequential (baseline)",
        "OpenBLAS DGEMM",
        "OpenBLAS DGEMM (same as Method 1)",
        "OpenBLAS DGEMM (same as Method 1)"
    };
    
    // Get OpenBLAS threading info
    char *openblas_threads = getenv("OPENBLAS_NUM_THREADS");
    char *omp_threads = getenv("OMP_NUM_THREADS");
    
    printf("\n");
    printf("=================================================================\n");
    printf("           MATRIX MULTIPLICATION CONFIGURATION\n");
    printf("=================================================================\n");
    printf("  Method:              %d (%s)\n", g_multiply_method, 
           method_names[g_multiply_method]);
    printf("  Implementation:      %s\n", 
           (g_multiply_method == 0) ? "Sequential loops" : "OpenBLAS cblas_dgemm");
    printf("  OPENBLAS_NUM_THREADS: %s\n", openblas_threads ? openblas_threads : "not set (using default)");
    printf("  OMP_NUM_THREADS:     %s\n", omp_threads ? omp_threads : "not set");
    
    if (g_multiply_method > 0)
    {
        printf("\n");
        printf("  OpenBLAS Features:\n");
        printf("    - Multi-threaded:  ✅\n");
        printf("    - Cache blocking:  ✅\n");
        printf("    - SIMD (AVX2/512): ✅\n");
        printf("    - Assembly opts:   ✅\n");
    }
    
    printf("=================================================================\n");
    printf("\n");
}

double get_expected_speedup(int method, int num_threads)
{
    switch (method)
    {
    case 0:
        return 1.0; // Sequential baseline
        
    case 1:
    case 2:
    case 3:
        // OpenBLAS typically achieves 90-95% of peak theoretical performance
        // with near-linear scaling up to 8-12 threads
        if (num_threads <= 4)
        {
            return num_threads * 0.95; // 95% efficiency
        }
        else if (num_threads <= 8)
        {
            return num_threads * 0.90; // 90% efficiency
        }
        else if (num_threads <= 16)
        {
            return num_threads * 0.80; // 80% efficiency
        }
        else
        {
            return num_threads * 0.70; // 70% efficiency (diminishing returns)
        }
        
    default:
        return 1.0;
    }
}

void print_expected_performance(void)
{
    char *threads_env = getenv("OPENBLAS_NUM_THREADS");
    int threads = threads_env ? atoi(threads_env) : omp_get_max_threads();
    
    printf("\n");
    printf("=================================================================\n");
    printf("           EXPECTED PERFORMANCE (OpenBLAS)\n");
    printf("=================================================================\n");
    printf("  System configuration:\n");
    printf("    CPU cores:         %d\n", threads);
    printf("    BLAS library:      OpenBLAS (multi-threaded)\n");
    printf("\n");
    printf("  Expected speedup vs sequential:\n");
    printf("    Method 0 (Sequential):         %.2f×\n", get_expected_speedup(0, threads));
    printf("    Method 1 (OpenBLAS):           %.2f×\n", get_expected_speedup(1, threads));
    printf("    Method 2 (OpenBLAS):           %.2f×\n", get_expected_speedup(2, threads));
    printf("    Method 3 (OpenBLAS):           %.2f×\n", get_expected_speedup(3, threads));
    printf("\n");
    printf("  Note: Methods 1-3 use identical OpenBLAS implementation.\n");
    printf("        All optimizations (threading, cache, SIMD) are built-in.\n");
    printf("\n");
    printf("  Time reduction (from 2.3 seconds baseline estimate):\n");
    printf("    OpenBLAS: %.2f seconds\n", 2.3 / get_expected_speedup(1, threads));
    printf("\n");
    printf("  Performance (from 1.5 GFLOPS baseline estimate):\n");
    printf("    OpenBLAS: %.1f GFLOPS\n", 1.5 * get_expected_speedup(1, threads));
    printf("=================================================================\n");
    printf("\n");
}

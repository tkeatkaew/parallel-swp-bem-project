/*******************************************************************************
 * matrix_multiply_optimized.h - UNIFIED Matrix Multiplication Header
 * 
 * VERSION 3.0: Combined Hybrid + OpenBLAS with Runtime Selection
 * 
 * Author: BEM Flow Path Analysis - Parallel Optimization Research
 * Date: 2025-12-05
 ******************************************************************************/

#ifndef MATRIX_MULTIPLY_OPTIMIZED_H
#define MATRIX_MULTIPLY_OPTIMIZED_H

#include "matrix_types.h"

/*******************************************************************************
 * DGEMM Type Selection (NEW in v3.0)
 ******************************************************************************/

/**
 * Set DGEMM implementation type
 * @param type 0=Hybrid (custom OpenMP+Cache+SIMD), 1=OpenBLAS (default)
 */
void set_dgemm_type(int type);

/**
 * Get current DGEMM implementation type
 * @return 0=Hybrid, 1=OpenBLAS
 */
int get_dgemm_type(void);

/**
 * Get human-readable name of current DGEMM type
 * @return "Hybrid (OpenMP+Cache+SIMD)" or "OpenBLAS"
 */
const char* get_dgemm_type_name(void);

/*******************************************************************************
 * Multiply Method Selection
 ******************************************************************************/

/**
 * Set matrix multiplication method
 * @param method 0=Sequential, 1=OpenMP, 2=OpenMP+Cache, 3=OpenMP+Cache+SIMD
 * 
 * When DGEMM_TYPE=0 (Hybrid):
 *   Method 0: Pure sequential (baseline)
 *   Method 1: OpenMP parallelization
 *   Method 2: OpenMP + Cache blocking
 *   Method 3: OpenMP + Cache + SIMD (AVX2)
 * 
 * When DGEMM_TYPE=1 (OpenBLAS):
 *   Method 0: Pure sequential (baseline)
 *   Method 1-3: All use OpenBLAS DGEMM (same implementation)
 */
void set_multiply_method(int method);
int get_multiply_method(void);

/**
 * Set cache block size for Hybrid methods
 * @param size Block size in range [8, 256], default=64
 */
void set_block_size(int size);
int get_block_size(void);

/**
 * Enable/disable verbose configuration output
 */
void set_multiply_verbose(int verbose);

/*******************************************************************************
 * Main Multiplication Functions
 ******************************************************************************/

/**
 * Main dispatcher - routes to appropriate implementation
 * based on current dgemm_type and multiply_method settings
 */
void multiply_matrix_optimized(matrix *a, matrix *b, matrix *x);

/**
 * Convenience wrapper - same as multiply_matrix_optimized()
 */
void multiply_matrix_auto(matrix *a, matrix *b, matrix *x);

/**
 * Multiply with temporary method override
 * @param method Method to use for this call only
 */
void multiply_matrix_method(matrix *a, matrix *b, matrix *x, int method);

/*******************************************************************************
 * Direct Implementation Access (for testing/benchmarking)
 ******************************************************************************/

/* Sequential baseline (always the same, regardless of dgemm_type) */
void multiply_matrix_sequential(matrix *a, matrix *b, matrix *x);

/* Hybrid implementations */
void multiply_matrix_openmp_hybrid(matrix *a, matrix *b, matrix *x);
void multiply_matrix_cache_hybrid(matrix *a, matrix *b, matrix *x, int block_size);
void multiply_matrix_simd_hybrid(matrix *a, matrix *b, matrix *x, int block_size);

/* OpenBLAS implementation */
void multiply_matrix_openblas(matrix *a, matrix *b, matrix *x);

/* Legacy wrappers (route based on current dgemm_type) */
void multiply_matrix_openmp(matrix *a, matrix *b, matrix *x);
void multiply_matrix_cache(matrix *a, matrix *b, matrix *x, int block_size);
void multiply_matrix_simd(matrix *a, matrix *b, matrix *x, int block_size);

/*******************************************************************************
 * Performance Information
 ******************************************************************************/

/**
 * Print current configuration (human-readable)
 */
void print_multiply_config(void);

/**
 * Get expected speedup for given method and thread count
 */
double get_expected_speedup(int method, int num_threads);

/**
 * Print expected performance estimates
 */
void print_expected_performance(void);

/**
 * Print configuration in CSV format (for benchmarking)
 */
void print_config_csv_header(void);
void print_config_csv(void);

#endif /* MATRIX_MULTIPLY_OPTIMIZED_H */

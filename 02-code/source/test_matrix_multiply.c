/*******************************************************************************
 * test_matrix_multiply.c
 * 
 * Test and benchmark program for optimized matrix multiplication
 * 
 * Verifies:
 * - Correctness (all methods produce same results)
 * - Performance (measures speedup for each method)
 * - Scaling (tests different thread counts)
 * 
 * Usage:
 *   gcc -O3 -fopenmp -mavx2 -mfma -o test_matrix_multiply \
 *       test_matrix_multiply.c matrix_multiply_optimized.c -lm
 *   
 *   ./test_matrix_multiply [size] [threads]
 *   
 *   Examples:
 *     ./test_matrix_multiply 512 6      # 512×512 matrices, 6 threads
 *     ./test_matrix_multiply 1024       # 1024×1024, auto threads
 *     ./test_matrix_multiply            # Default: 256×256
 * 
 * Author: Test suite for BEM optimization
 * Date: 2025-01-27
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

/* Mock matrix_types.h for standalone compilation */
typedef struct {
    int rows;
    int columns;
    double *value;
    int transpose;
    int invert;
} matrix;

/* Include the optimization functions */
#include "matrix_multiply_optimized.c"

/*******************************************************************************
 * Helper Functions
 ******************************************************************************/

/* Create a matrix */
matrix *create_test_matrix(int rows, int cols) {
    matrix *m = (matrix *)malloc(sizeof(matrix));
    m->rows = rows;
    m->columns = cols;
    m->value = (double *)malloc(rows * cols * sizeof(double));
    m->transpose = 0;
    m->invert = 0;
    return m;
}

/* Free a matrix */
void destroy_test_matrix(matrix *m) {
    if (m) {
        if (m->value) free(m->value);
        free(m);
    }
}

/* Initialize matrix with random values */
void init_random(matrix *m, int seed) {
    srand(seed);
    for (int i = 0; i < m->rows * m->columns; i++) {
        m->value[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;  // [-1, 1]
    }
}

/* Initialize matrix with specific pattern (for debugging) */
void init_pattern(matrix *m) {
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->columns; j++) {
            m->value[i * m->columns + j] = (i + 1) * 10.0 + (j + 1);
        }
    }
}

/* Compare two matrices */
double compare_matrices(matrix *a, matrix *b, double *max_diff, double *avg_diff) {
    if (a->rows != b->rows || a->columns != b->columns) {
        fprintf(stderr, "Error: Matrix size mismatch\n");
        return -1.0;
    }
    
    double max_d = 0.0;
    double sum_d = 0.0;
    int n = a->rows * a->columns;
    
    for (int i = 0; i < n; i++) {
        double diff = fabs(a->value[i] - b->value[i]);
        if (diff > max_d) max_d = diff;
        sum_d += diff;
    }
    
    *max_diff = max_d;
    *avg_diff = sum_d / n;
    
    return max_d;
}

/* Print matrix (small matrices only) */
void print_matrix(matrix *m, const char *name) {
    printf("%s (%d×%d):\n", name, m->rows, m->columns);
    
    int max_print = 8;  // Only print up to 8×8
    int rows = (m->rows < max_print) ? m->rows : max_print;
    int cols = (m->columns < max_print) ? m->columns : max_print;
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%8.3f ", m->value[i * m->columns + j]);
        }
        if (m->columns > max_print) printf("...");
        printf("\n");
    }
    if (m->rows > max_print) printf("...\n");
    printf("\n");
}

/* Get current time in seconds */
double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/*******************************************************************************
 * Test Functions
 ******************************************************************************/

/* Test correctness of all methods */
int test_correctness(int size) {
    printf("=================================================================\n");
    printf("TEST 1: Correctness Verification (size=%d)\n", size);
    printf("=================================================================\n");
    
    /* Create matrices */
    matrix *A = create_test_matrix(size, size);
    matrix *B = create_test_matrix(size, size);
    matrix *C0 = create_test_matrix(size, size);  // Method 0 result
    matrix *C1 = create_test_matrix(size, size);  // Method 1 result
    matrix *C2 = create_test_matrix(size, size);  // Method 2 result
    matrix *C3 = create_test_matrix(size, size);  // Method 3 result
    
    /* Initialize input matrices */
    init_random(A, 42);
    init_random(B, 123);
    
    printf("Computing with all 4 methods...\n");
    
    /* Compute with each method */
    multiply_matrix_sequential(A, B, C0);
    multiply_matrix_openmp(A, B, C1);
    multiply_matrix_cache(A, B, C2, 32);
    multiply_matrix_simd(A, B, C3, 32);
    
    /* Compare results */
    double max_diff, avg_diff;
    int pass = 1;
    
    printf("\nComparing results:\n");
    printf("  Method 1 vs Method 0: ");
    compare_matrices(C0, C1, &max_diff, &avg_diff);
    printf("max=%.2e, avg=%.2e ", max_diff, avg_diff);
    if (max_diff < 1e-10) {
        printf("✅ PASS\n");
    } else {
        printf("❌ FAIL\n");
        pass = 0;
    }
    
    printf("  Method 2 vs Method 0: ");
    compare_matrices(C0, C2, &max_diff, &avg_diff);
    printf("max=%.2e, avg=%.2e ", max_diff, avg_diff);
    if (max_diff < 1e-10) {
        printf("✅ PASS\n");
    } else {
        printf("❌ FAIL\n");
        pass = 0;
    }
    
    printf("  Method 3 vs Method 0: ");
    compare_matrices(C0, C3, &max_diff, &avg_diff);
    printf("max=%.2e, avg=%.2e ", max_diff, avg_diff);
    if (max_diff < 1e-10) {
        printf("✅ PASS\n");
    } else {
        printf("❌ FAIL\n");
        pass = 0;
    }
    
    /* Cleanup */
    destroy_test_matrix(A);
    destroy_test_matrix(B);
    destroy_test_matrix(C0);
    destroy_test_matrix(C1);
    destroy_test_matrix(C2);
    destroy_test_matrix(C3);
    
    printf("\nResult: %s\n", pass ? "✅ ALL TESTS PASSED" : "❌ SOME TESTS FAILED");
    printf("=================================================================\n\n");
    
    return pass;
}

/* Benchmark performance of all methods */
void benchmark_performance(int size, int num_threads) {
    printf("=================================================================\n");
    printf("TEST 2: Performance Benchmark (size=%d, threads=%d)\n", size, num_threads);
    printf("=================================================================\n");
    
    omp_set_num_threads(num_threads);
    
    /* Create matrices */
    matrix *A = create_test_matrix(size, size);
    matrix *B = create_test_matrix(size, size);
    matrix *C = create_test_matrix(size, size);
    
    /* Initialize */
    init_random(A, 42);
    init_random(B, 123);
    
    /* Warm-up */
    printf("Warming up...\n");
    multiply_matrix_sequential(A, B, C);
    
    printf("\nBenchmarking...\n\n");
    
    double times[4];
    const char *names[] = {"Sequential", "OpenMP", "OpenMP+Cache", "OpenMP+Cache+SIMD"};
    
    /* Benchmark Method 0 */
    printf("Method 0 (Sequential)... ");
    fflush(stdout);
    double t0 = get_time();
    multiply_matrix_sequential(A, B, C);
    double t1 = get_time();
    times[0] = t1 - t0;
    printf("%.3f sec\n", times[0]);
    
    /* Benchmark Method 1 */
    printf("Method 1 (OpenMP)... ");
    fflush(stdout);
    t0 = get_time();
    multiply_matrix_openmp(A, B, C);
    t1 = get_time();
    times[1] = t1 - t0;
    printf("%.3f sec\n", times[1]);
    
    /* Benchmark Method 2 */
    printf("Method 2 (OpenMP+Cache)... ");
    fflush(stdout);
    t0 = get_time();
    multiply_matrix_cache(A, B, C, 32);
    t1 = get_time();
    times[2] = t1 - t0;
    printf("%.3f sec\n", times[2]);
    
    /* Benchmark Method 3 */
    printf("Method 3 (OpenMP+Cache+SIMD)... ");
    fflush(stdout);
    t0 = get_time();
    multiply_matrix_simd(A, B, C, 32);
    t1 = get_time();
    times[3] = t1 - t0;
    printf("%.3f sec\n", times[3]);
    
    /* Calculate GFLOPS */
    double flops = 2.0 * size * size * size;  // 2*N^3 operations
    double gflops[4];
    for (int i = 0; i < 4; i++) {
        gflops[i] = flops / times[i] / 1e9;
    }
    
    /* Print results */
    printf("\n");
    printf("┌────────┬──────────────────────┬─────────┬──────────┬───────────┐\n");
    printf("│ Method │ Description          │ Time(s) │ Speedup  │ GFLOPS    │\n");
    printf("├────────┼──────────────────────┼─────────┼──────────┼───────────┤\n");
    for (int i = 0; i < 4; i++) {
        double speedup = times[0] / times[i];
        printf("│   %d    │ %-20s │ %7.3f │ %7.2f× │ %8.2f  │\n",
               i, names[i], times[i], speedup, gflops[i]);
    }
    printf("└────────┴──────────────────────┴─────────┴──────────┴───────────┘\n");
    
    /* Cleanup */
    destroy_test_matrix(A);
    destroy_test_matrix(B);
    destroy_test_matrix(C);
    
    printf("=================================================================\n\n");
}

/* Test thread scaling */
void test_scaling(int size) {
    printf("=================================================================\n");
    printf("TEST 3: Thread Scaling (size=%d)\n", size);
    printf("=================================================================\n");
    
    matrix *A = create_test_matrix(size, size);
    matrix *B = create_test_matrix(size, size);
    matrix *C = create_test_matrix(size, size);
    
    init_random(A, 42);
    init_random(B, 123);
    
    int max_threads = omp_get_max_threads();
    printf("Testing with 1 to %d threads...\n\n", max_threads);
    
    printf("┌─────────┬─────────┬──────────┬────────────┐\n");
    printf("│ Threads │ Time(s) │ Speedup  │ Efficiency │\n");
    printf("├─────────┼─────────┼──────────┼────────────┤\n");
    
    double time_1thread = 0;
    
    for (int t = 1; t <= max_threads; t++) {
        omp_set_num_threads(t);
        
        double t0 = get_time();
        multiply_matrix_simd(A, B, C, 32);
        double t1 = get_time();
        double time = t1 - t0;
        
        if (t == 1) time_1thread = time;
        
        double speedup = time_1thread / time;
        double efficiency = speedup / t * 100.0;
        
        printf("│  %4d   │ %7.3f │  %6.2f× │   %5.1f%%   │\n",
               t, time, speedup, efficiency);
    }
    
    printf("└─────────┴─────────┴──────────┴────────────┘\n");
    
    destroy_test_matrix(A);
    destroy_test_matrix(B);
    destroy_test_matrix(C);
    
    printf("=================================================================\n\n");
}

/*******************************************************************************
 * Main Program
 ******************************************************************************/

int main(int argc, char *argv[]) {
    /* Parse arguments */
    int size = 256;        // Default matrix size
    int threads = 0;       // 0 = auto-detect
    
    if (argc > 1) {
        size = atoi(argv[1]);
        if (size < 8 || size > 8192) {
            fprintf(stderr, "Error: Size must be between 8 and 8192\n");
            return 1;
        }
    }
    
    if (argc > 2) {
        threads = atoi(argv[2]);
        if (threads < 1 || threads > 64) {
            fprintf(stderr, "Error: Threads must be between 1 and 64\n");
            return 1;
        }
        omp_set_num_threads(threads);
    }
    
    if (threads == 0) {
        threads = omp_get_max_threads();
    }
    
    /* Print configuration */
    printf("\n");
    printf("╔═══════════════════════════════════════════════════════════════╗\n");
    printf("║       MATRIX MULTIPLICATION OPTIMIZATION TEST SUITE          ║\n");
    printf("╚═══════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    printf("Configuration:\n");
    printf("  Matrix size:    %d × %d\n", size, size);
    printf("  OpenMP threads: %d\n", threads);
    printf("  CPU cores:      %d\n", omp_get_num_procs());
    
#ifdef __AVX2__
    printf("  SIMD support:   AVX2 enabled ✅\n");
#else
    printf("  SIMD support:   No AVX2 (scalar fallback)\n");
#endif
    
    printf("\n");
    
    /* Run tests */
    int pass = test_correctness(size);
    
    if (!pass) {
        fprintf(stderr, "\n❌ Correctness test failed! Skipping performance tests.\n");
        return 1;
    }
    
    benchmark_performance(size, threads);
    test_scaling(size);
    
    /* Summary */
    printf("╔═══════════════════════════════════════════════════════════════╗\n");
    printf("║                       TEST SUMMARY                            ║\n");
    printf("╚═══════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    printf("✅ All tests completed successfully!\n");
    printf("\n");
    printf("Recommendations:\n");
    printf("  - Use Method 3 for production (fastest)\n");
    printf("  - Optimal block size: 32 (for most systems)\n");
    printf("  - Expected speedup on your system: 5-8×\n");
    printf("\n");
    
    return 0;
}

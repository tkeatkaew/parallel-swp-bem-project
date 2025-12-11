/*******************************************************************************
 * performance_summary.c
 *
 * Comprehensive performance tracking and reporting for BEM solver
 *
 * Key points:
 * - Tracks setup, BEM, and finalization times
 * - Tracks matrix multiply + inversion times and GFLOPS
 * - Provides consistent numbers across:
 *     * TIME BREAKDOWN
 *     * TIMING STATISTICS
 *     * JOURNAL TABLE (LaTeX)
 *     * CSV EXPORT
 *
 * Author: For BEM Flow Path Analysis
 * Date:   2025-10-28
 * Fixed:  2025-11-09
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "performance_summary.h"

/*******************************************************************************
 * Global Performance Summary
 ******************************************************************************/

PerformanceSummary g_perf_summary = {0};

/*******************************************************************************
 * Initialization
 ******************************************************************************/

void init_performance_summary(void) {
    memset(&g_perf_summary, 0, sizeof(PerformanceSummary));
    g_perf_summary.multiply_method   = -1;
    g_perf_summary.inversion_method  = -1;
    g_perf_summary.num_threads       = -1;
    g_perf_summary.block_size        = -1;
}

/*******************************************************************************
 * Update Functions
 ******************************************************************************/

void update_setup_time(double time_sec) {
    g_perf_summary.setup_time += time_sec;
}

void update_bem_time(double time_sec) {
    g_perf_summary.bem_time += time_sec;
}

void update_multiply_time(double time_sec, int rows, int cols, int k) {
    g_perf_summary.matrix_multiply_time += time_sec;
    g_perf_summary.num_multiplications++;

    /* Update total matrix computation time (DGEMM + Inversion) */
    g_perf_summary.matrix_computation_time =
        g_perf_summary.matrix_multiply_time + g_perf_summary.matrix_inversion_time;

    /* Track largest matrix dimensions seen */
    if (rows > g_perf_summary.max_matrix_rows) {
        g_perf_summary.max_matrix_rows = rows;
    }
    if (cols > g_perf_summary.max_matrix_cols) {
        g_perf_summary.max_matrix_cols = cols;
    }

    /* Matrix multiply FLOPs: 2*m*n*k */
    if (time_sec > 0.0) {
        long long flops = 2LL * rows * cols * k;
        double gflops   = (flops / 1e9) / time_sec;

        /* Running average GFLOPS */
        if (g_perf_summary.multiply_gflops == 0.0) {
            g_perf_summary.multiply_gflops = gflops;
        } else {
            g_perf_summary.multiply_gflops =
                (g_perf_summary.multiply_gflops * (g_perf_summary.num_multiplications - 1) + gflops)
                / g_perf_summary.num_multiplications;
        }
    }
}

void update_inversion_time(double time_sec, int n) {
    g_perf_summary.matrix_inversion_time += time_sec;
    g_perf_summary.num_inversions++;

    /* Update total matrix computation time (DGEMM + Inversion) */
    g_perf_summary.matrix_computation_time =
        g_perf_summary.matrix_multiply_time + g_perf_summary.matrix_inversion_time;

    /* Matrix inversion FLOPs: approx (2/3)*n^3 */
    if (time_sec > 0.0) {
        long long flops = (2LL * n * n * n) / 3;
        double gflops   = (flops / 1e9) / time_sec;

        /* Running average GFLOPS */
        if (g_perf_summary.inversion_gflops == 0.0) {
            g_perf_summary.inversion_gflops = gflops;
        } else {
            g_perf_summary.inversion_gflops =
                (g_perf_summary.inversion_gflops * (g_perf_summary.num_inversions - 1) + gflops)
                / g_perf_summary.num_inversions;
        }
    }
}

/*******************************************************************************
 * PATCH: Add this function to performance_summary.c
 * 
 * Location: After the update_inversion_time() function (around line 111)
 * 
 * This is a legacy wrapper - the main stats are already updated by
 * update_inversion_time() which is called earlier in matrix.c
 ******************************************************************************/

/* Legacy compatibility wrapper */
//void update_matrix_inversion_stats(double time_sec) {
    /* This function is called for backward compatibility.
     * The actual statistics are already updated by update_inversion_time()
     * which provides both time and matrix size for GFLOPS calculation.
     * This wrapper does nothing but prevents link errors. */
  //  (void)time_sec;  /* Suppress unused parameter warning */
//}

void update_finalization_time(double time_sec) {
    g_perf_summary.finalization_time += time_sec;
}

void update_memory_usage(long vmrss_kb, long vmsize_kb) {
    if (g_perf_summary.initial_memory_kb == 0) {
        g_perf_summary.initial_memory_kb = vmrss_kb;
    }

    if (vmrss_kb > g_perf_summary.peak_memory_kb) {
        g_perf_summary.peak_memory_kb = vmrss_kb;
    }

    g_perf_summary.final_memory_kb = vmrss_kb;
}

void set_performance_config(int multiply_method, int inversion_method,
                            int num_threads, int block_size) {
    g_perf_summary.multiply_method  = multiply_method;
    g_perf_summary.inversion_method = inversion_method;
    g_perf_summary.num_threads      = num_threads;
    g_perf_summary.block_size       = block_size;
}

void set_problem_parameters(double step, double rm, double dr,
                            int zones, int points) {
    g_perf_summary.step_size  = step;
    g_perf_summary.rm         = rm;
    g_perf_summary.dr         = dr;
    g_perf_summary.num_zones  = zones;
    g_perf_summary.max_points = points;
}

/*******************************************************************************
 * Print Functions
 ******************************************************************************/

void print_performance_summary(void) {
    /* Total wall-clock time for the whole program */
    g_perf_summary.total_time =
        g_perf_summary.setup_time +
        g_perf_summary.bem_time +
        g_perf_summary.finalization_time;

    printf("\n");
    printf("################################################################################\n");
    printf("###                                                                          ###\n");
    printf("###                    PERFORMANCE SUMMARY REPORT                            ###\n");
    printf("###                                                                          ###\n");
    printf("################################################################################\n\n");

    /***** Problem Configuration *****/
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("PROBLEM CONFIGURATION\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("  Step size:              %.3f\n",  g_perf_summary.step_size);
    printf("  Rm:                     %.3f\n",  g_perf_summary.rm);
    printf("  Dr:                     %.6f\n",  g_perf_summary.dr);
    printf("  Number of zones:        %d\n",    g_perf_summary.num_zones);
    printf("  Max points in zone:     %d\n",    g_perf_summary.max_points);
    printf("  Largest matrix size:    %d × %d\n",
           g_perf_summary.max_matrix_rows, g_perf_summary.max_matrix_cols);
    printf("\n");

    /***** System Configuration *****/
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("SYSTEM CONFIGURATION\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n");

    const char *multiply_methods[]  = {"Sequential", "OpenMP", "OpenMP+Cache", "OpenMP+Cache+SIMD"};
    const char *inversion_methods[] = {"Parallel (LAPACK)", "Sequential (Manual)"};

    if (g_perf_summary.multiply_method >= 0 &&
        g_perf_summary.multiply_method <= 3) {
        printf("  Matrix multiply method: %d (%s)\n",
               g_perf_summary.multiply_method,
               multiply_methods[g_perf_summary.multiply_method]);
    }

    if (g_perf_summary.inversion_method >= 0 &&
        g_perf_summary.inversion_method <= 1) {
        printf("  Matrix inversion:       %d (%s)\n",
               g_perf_summary.inversion_method,
               inversion_methods[g_perf_summary.inversion_method]);
    }

    if (g_perf_summary.num_threads > 0) {
        printf("  OpenMP threads:         %d\n", g_perf_summary.num_threads);
    }

    if (g_perf_summary.block_size > 0) {
        printf("  Cache block size:       %d\n", g_perf_summary.block_size);
    }

    printf("\n");

    /***** Time Breakdown *****/
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("TIME BREAKDOWN\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("  Component                      Time (sec)      %% of Total\n");
    printf("  ─────────────────────────────────────────────────────────────────────────────\n");

    double total = g_perf_summary.total_time;
    if (total <= 0.0) total = 1.0; /* Avoid division by zero */

    printf("  Setup & Initialization         %9.4f       %6.2f%%\n",
           g_perf_summary.setup_time,
           (g_perf_summary.setup_time / total) * 100.0);

    printf("  BEM Computation                %9.4f       %6.2f%%\n",
           g_perf_summary.bem_time,
           (g_perf_summary.bem_time / total) * 100.0);

    printf("    ├─ Matrix Multiplication     %9.4f       %6.2f%%\n",
           g_perf_summary.matrix_multiply_time,
           (g_perf_summary.matrix_multiply_time / total) * 100.0);

    printf("    └─ Matrix Inversion          %9.4f       %6.2f%%\n",
           g_perf_summary.matrix_inversion_time,
           (g_perf_summary.matrix_inversion_time / total) * 100.0);

    printf("  Finalization                   %9.4f       %6.2f%%\n",
           g_perf_summary.finalization_time,
           (g_perf_summary.finalization_time / total) * 100.0);

    printf("  ───────────────────────────────────────────────────────────────────────────\n");
    printf("  TOTAL                          %9.4f      100.00%%\n\n", total);

    /***** Timing Statistics (Matrix Ops) *****/
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("TIMING STATISTICS:\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("  Total DGEMM (matrix multiply) time:  %.6f seconds\n",
           g_perf_summary.matrix_multiply_time);
    printf("  Total Matrix Inversion time:         %.6f seconds\n",
           g_perf_summary.matrix_inversion_time);
    printf("  Total computation time:              %.6f seconds\n",
           g_perf_summary.matrix_computation_time);
    printf("\n");

    /***** Operation Counts *****/
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("OPERATION COUNTS:\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("  DGEMM calls:                         %d\n", g_perf_summary.num_multiplications);
    printf("  Matrix Inversion calls:              %d\n", g_perf_summary.num_inversions);
    printf("\n");

    /***** Performance Metrics *****/
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("PERFORMANCE METRICS:\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n");

    if (g_perf_summary.multiply_gflops > 0.0) {
        printf("  DGEMM GFLOPS:                        %.2f\n",
               g_perf_summary.multiply_gflops);
    }

    if (g_perf_summary.num_multiplications > 0) {
        printf("  Average DGEMM time per call:         %.6f seconds\n",
               g_perf_summary.matrix_multiply_time /
               g_perf_summary.num_multiplications);
    }

    if (g_perf_summary.num_inversions > 0) {
        printf("  Average Matrix Inv time per call:    %.6f seconds\n",
               g_perf_summary.matrix_inversion_time /
               g_perf_summary.num_inversions);
    }

    printf("\n");

    /***** Memory Usage *****/
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("MEMORY USAGE:\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("  Peak allocated (tracked):            %.2f MB\n",
           (g_perf_summary.peak_memory_kb -
            g_perf_summary.initial_memory_kb) / 1024.0);
    printf("  VmRSS (resident set):                %.2f MB\n",
           g_perf_summary.peak_memory_kb / 1024.0);
    printf("\n");

    printf("################################################################################\n");
    printf("###                      END OF PERFORMANCE REPORT                           ###\n");
    printf("################################################################################\n\n");
}

void print_journal_table(void) {
    printf("\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n");
    printf("TABLE FOR JOURNAL PAPER (LaTeX Format)\n");
    printf("═══════════════════════════════════════════════════════════════════════════════\n\n");

    const char *multiply_methods[]  = {"Sequential", "OpenMP", "OpenMP+Cache", "OpenMP+Cache+SIMD"};
    const char *inversion_methods[] = {"Parallel (LAPACK)", "Sequential (Manual)"};

    printf("\\begin{table}[htbp]\n");
    printf("\\centering\n");
    printf("\\caption{Performance Results for BEM Surface Water Path Delineation}\n");
    printf("\\label{tab:bem_performance}\n");
    printf("\\begin{tabular}{|l|r|}\n");
    printf("\\hline\n");
    printf("\\textbf{Parameter} & \\textbf{Value} \\\\\n");
    printf("\\hline\n");
    printf("\\hline\n");
    printf("Matrix Size & $%d \\times %d$ \\\\\n",
           g_perf_summary.max_matrix_rows,
           g_perf_summary.max_matrix_cols);
    printf("\\hline\n");

    if (g_perf_summary.multiply_method >= 0) {
        printf("Multiply Method & %s \\\\\n",
               multiply_methods[g_perf_summary.multiply_method]);
    }

    if (g_perf_summary.inversion_method >= 0) {
        printf("Inversion Method & %s \\\\\n",
               inversion_methods[g_perf_summary.inversion_method]);
    }

    printf("\\hline\n");
    printf("Total Time & %.4f s \\\\\n", g_perf_summary.total_time);
    printf("Matrix Operations Time & %.4f s \\\\\n",
           g_perf_summary.matrix_computation_time);
    printf("\\quad DGEMM Time & %.4f s \\\\\n",
           g_perf_summary.matrix_multiply_time);
    printf("\\quad Inversion Time & %.4f s \\\\\n",
           g_perf_summary.matrix_inversion_time);
    printf("\\hline\n");

    if (g_perf_summary.multiply_gflops > 0.0) {
        printf("Multiply Performance & %.2f GFLOPS \\\\\n",
               g_perf_summary.multiply_gflops);
    }

    if (g_perf_summary.inversion_gflops > 0.0) {
        printf("Inversion Performance & %.2f GFLOPS \\\\\n",
               g_perf_summary.inversion_gflops);
    }

    printf("\\hline\n");
    printf("Peak Memory & %.2f MB \\\\\n",
           g_perf_summary.peak_memory_kb / 1024.0);
    printf("\\hline\n");
    printf("\\end{tabular}\n");
    printf("\\end{table}\n\n");

    /* CSV-style inline output for quick copy-paste */
    printf("─────────────────────────────────────────────────────────────────────────────\n");
    printf("CSV FORMAT (for Excel/plotting):\n");
    printf("─────────────────────────────────────────────────────────────────────────────\n\n");

    printf("Parameter,Value\n");
    printf("Matrix_Rows,%d\n",      g_perf_summary.max_matrix_rows);
    printf("Matrix_Cols,%d\n",      g_perf_summary.max_matrix_cols);
    printf("Multiply_Method,%d\n", g_perf_summary.multiply_method);
    printf("Inversion_Method,%d\n",g_perf_summary.inversion_method);
    printf("Num_Threads,%d\n",     g_perf_summary.num_threads);
    printf("Block_Size,%d\n",      g_perf_summary.block_size);
    printf("Total_Time_sec,%.6f\n",        g_perf_summary.total_time);
    printf("Computation_Time_sec,%.6f\n",  g_perf_summary.matrix_computation_time);
    printf("Multiply_Time_sec,%.6f\n",     g_perf_summary.matrix_multiply_time);
    printf("Inversion_Time_sec,%.6f\n",    g_perf_summary.matrix_inversion_time);
    printf("Multiply_GFLOPS,%.2f\n",       g_perf_summary.multiply_gflops);
    printf("Inversion_GFLOPS,%.2f\n",      g_perf_summary.inversion_gflops);
    printf("Peak_Memory_MB,%.2f\n",        g_perf_summary.peak_memory_kb / 1024.0);
    printf("\n");
}

void export_performance_csv(const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Warning: Could not open %s for writing\n", filename);
        return;
    }

    fprintf(fp, "Parameter,Value\n");
    fprintf(fp, "Matrix_Rows,%d\n",      g_perf_summary.max_matrix_rows);
    fprintf(fp, "Matrix_Cols,%d\n",      g_perf_summary.max_matrix_cols);
    fprintf(fp, "Multiply_Method,%d\n", g_perf_summary.multiply_method);
    fprintf(fp, "Inversion_Method,%d\n",g_perf_summary.inversion_method);
    fprintf(fp, "Num_Threads,%d\n",     g_perf_summary.num_threads);
    fprintf(fp, "Block_Size,%d\n",      g_perf_summary.block_size);
    fprintf(fp, "Total_Time_sec,%.6f\n",        g_perf_summary.total_time);
    fprintf(fp, "Setup_Time_sec,%.6f\n",        g_perf_summary.setup_time);
    fprintf(fp, "BEM_Time_sec,%.6f\n",          g_perf_summary.bem_time);
    fprintf(fp, "Computation_Time_sec,%.6f\n",  g_perf_summary.matrix_computation_time);
    fprintf(fp, "Multiply_Time_sec,%.6f\n",     g_perf_summary.matrix_multiply_time);
    fprintf(fp, "Inversion_Time_sec,%.6f\n",    g_perf_summary.matrix_inversion_time);
    fprintf(fp, "Finalization_Time_sec,%.6f\n", g_perf_summary.finalization_time);
    fprintf(fp, "Num_Multiplications,%d\n",     g_perf_summary.num_multiplications);
    fprintf(fp, "Num_Inversions,%d\n",          g_perf_summary.num_inversions);
    fprintf(fp, "Multiply_GFLOPS,%.2f\n",       g_perf_summary.multiply_gflops);
    fprintf(fp, "Inversion_GFLOPS,%.2f\n",      g_perf_summary.inversion_gflops);
    fprintf(fp, "Initial_Memory_MB,%.2f\n",     g_perf_summary.initial_memory_kb / 1024.0);
    fprintf(fp, "Peak_Memory_MB,%.2f\n",        g_perf_summary.peak_memory_kb / 1024.0);
    fprintf(fp, "Final_Memory_MB,%.2f\n",       g_perf_summary.final_memory_kb / 1024.0);

    fclose(fp);

    printf("Performance data exported to: %s\n", filename);
}

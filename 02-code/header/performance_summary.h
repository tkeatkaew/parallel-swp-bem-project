/*******************************************************************************
 * performance_summary.h
 *
 * Header file for comprehensive performance tracking and reporting
 *
 * Author: For BEM Flow Path Analysis
 * Date:   2025-10-28
 * Fixed:  2025-12-05 (Added struct definition)
 ******************************************************************************/

#ifndef PERFORMANCE_SUMMARY_H
#define PERFORMANCE_SUMMARY_H

/*******************************************************************************
 * Performance Summary Structure
 ******************************************************************************/

typedef struct {
    /* Timing data */
    double setup_time;              /* Setup/initialization time (seconds) */
    double bem_time;                /* BEM computation time (seconds) */
    double finalization_time;       /* Finalization/cleanup time (seconds) */
    double total_time;              /* Total execution time (seconds) */
    
    /* Matrix operation timing */
    double matrix_multiply_time;    /* Total matrix multiplication time */
    double matrix_inversion_time;   /* Total matrix inversion time */
    double matrix_computation_time; /* multiply + inversion combined */
    
    /* Operation counts */
    int num_multiplications;        /* Number of DGEMM calls */
    int num_inversions;             /* Number of matrix inversions */
    
    /* Performance metrics */
    double multiply_gflops;         /* Average GFLOPS for multiply */
    double inversion_gflops;        /* Average GFLOPS for inversion */
    
    /* Matrix dimensions */
    int max_matrix_rows;            /* Largest matrix rows */
    int max_matrix_cols;            /* Largest matrix columns */
    
    /* Configuration */
    int multiply_method;            /* 0=Seq, 1=OMP, 2=Cache, 3=SIMD */
    int inversion_method;           /* 0=Parallel, 1=Sequential */
    int num_threads;                /* Number of OpenMP threads */
    int block_size;                 /* Cache block size */
    
    /* Problem parameters */
    double step_size;               /* Step size parameter */
    double rm;                      /* Rm parameter */
    double dr;                      /* Dr parameter */
    int num_zones;                  /* Number of catchment zones */
    int max_points;                 /* Max points in any zone */
    
    /* Memory tracking */
    long initial_memory_kb;         /* Initial memory usage (KB) */
    long peak_memory_kb;            /* Peak memory usage (KB) */
    long final_memory_kb;           /* Final memory usage (KB) */
    
} PerformanceSummary;

/*******************************************************************************
 * Global Performance Summary (defined in performance_summary.c)
 ******************************************************************************/

extern PerformanceSummary g_perf_summary;

/*******************************************************************************
 * Function Prototypes
 ******************************************************************************/

/* Initialization */
void init_performance_summary(void);

/* Update functions */
void update_setup_time(double time_sec);
void update_bem_time(double time_sec);
void update_multiply_time(double time_sec, int rows, int cols, int k);
void update_inversion_time(double time_sec, int n);
void update_finalization_time(double time_sec);
void update_memory_usage(long vmrss_kb, long vmsize_kb);

/* Legacy compatibility - maps to update_inversion_time */
void update_matrix_inversion_stats(double time_sec);

/* Configuration */
void set_performance_config(int multiply_method, int inversion_method, 
                            int num_threads, int block_size);
void set_problem_parameters(double step, double rm, double dr, 
                            int zones, int points);

/* Output functions */
void print_performance_summary(void);
void print_journal_table(void);
void export_performance_csv(const char *filename);

/* Utility */
void get_memory_usage_kb(long *vmrss_kb, long *vmsize_kb);

#endif /* PERFORMANCE_SUMMARY_H */

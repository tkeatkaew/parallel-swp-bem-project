/*
 * zone_processor.h
 * 
 * Thread-safe zone processing module for Level 0 Contour-Layer Parallelism
 * 
 * This module encapsulates the processing of a single elevation zone,
 * including BEM matrix assembly, system solving, and result collection.
 * All functions are designed to be thread-safe for use in OpenMP parallel regions.
 *
 * Author: Generated for BEM SWP Research
 * Date: 2025-11-14
 * Version: 1.0
 */

#ifndef ZONE_PROCESSOR_H
#define ZONE_PROCESSOR_H

#include "boundary_types.h"
#include "matrix_types.h"
#include "co_matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

/*─────────────────────────────────────────────────────────────────────────
 * Data Structures
 *─────────────────────────────────────────────────────────────────────────*/

/**
 * @brief Result structure for a single zone computation
 * 
 * Contains all data produced by processing one elevation zone,
 * including timing information, computed area, and intermediate results.
 * Each thread maintains its own zone_result to avoid data races.
 */
typedef struct {
    /* Zone identification */
    int zone_id;                    /* Zone index (0-based) */
    boundary *b;                    /* Pointer to zone boundary (read-only) */
    
    /* BEM computation data (thread-private) */
    bem_vectors *vectors;           /* BEM vectors (bvv, bcv) */
    matrix *bvv_solved;             /* Solved voltage vector */
    matrix *bcv_solved;             /* Solved current vector */
    
    /* Results */
    double catchment_area;          /* Computed area for this zone */
    int num_streamlines;            /* Number of streamlines computed */
    path **streamlines;             /* Array of streamline paths */
    
    /* Timing statistics */
    double setup_time;              /* Data structure allocation time */
    double assembly_time;           /* BEM matrix assembly time */
    double solve_time;              /* Linear system solve time */
    double streamline_time;         /* Streamline computation time */
    double total_time;              /* Total processing time */
    
    /* Memory usage */
    long peak_memory_kb;            /* Peak memory usage in KB */
    
    /* Validation data */
    double matrix_norm_H;           /* Frobenius norm of H matrix */
    double matrix_norm_G;           /* Frobenius norm of G matrix */
    double residual_norm;           /* Solution residual norm */
    
    /* Thread information */
    int thread_id;                  /* OpenMP thread that processed this zone */
    
} zone_result;

/*─────────────────────────────────────────────────────────────────────────
 * Function Declarations
 *─────────────────────────────────────────────────────────────────────────*/

/**
 * @brief Process a single elevation zone (thread-safe)
 * 
 * This function performs complete BEM analysis for one zone:
 * 1. Allocates thread-private data structures
 * 2. Assembles BEM matrices [H] and [G]
 * 3. Solves the BEM system using LAPACK
 * 4. Computes streamlines (if needed)
 * 5. Calculates zone catchment area
 * 
 * Thread Safety:
 * - All data structures are allocated on the thread's stack or heap
 * - No shared mutable state between threads
 * - Safe to call from OpenMP parallel regions
 * 
 * @param c           Catchment structure (read-only access)
 * @param zone_idx    Zone index to process (0-based)
 * @param step_size   Streamline step size
 * @param rm          BEM parameter rm
 * @param dr          BEM parameter dr
 * @param max_steps   Maximum streamline steps
 * @param max_points  Maximum points in any zone (for buffer allocation)
 * 
 * @return Pointer to zone_result structure (caller must free)
 * @retval NULL on error
 */
zone_result* process_single_zone(
    catchment *c,
    int zone_idx,
    double step_size,
    double rm,
    double dr,
    int max_steps,
    int max_points
);

/**
 * @brief Free memory allocated for zone result
 * 
 * Safely deallocates all memory associated with a zone_result structure,
 * including BEM matrices, streamlines, and the structure itself.
 * 
 * @param result Pointer to zone_result to free (can be NULL)
 */
void free_zone_result(zone_result *result);

/**
 * @brief Aggregate results from multiple zones
 * 
 * Combines results from all processed zones to compute total catchment area
 * and collect overall statistics.
 * 
 * @param results     Array of zone results
 * @param num_zones   Number of zones
 * @param total_area  Output: total catchment area
 * @param total_time  Output: sum of all zone processing times
 */
void aggregate_zone_results(
    zone_result **results,
    int num_zones,
    double *total_area,
    double *total_time
);

/**
 * @brief Print detailed statistics for all zones
 * 
 * Outputs formatted table with per-zone timing, memory usage, and results.
 * Useful for debugging and performance analysis.
 * 
 * @param results    Array of zone results
 * @param num_zones  Number of zones
 */
void print_zone_statistics(zone_result **results, int num_zones);

/**
 * @brief Validate zone results against expected values
 * 
 * Compares computed results with known correct values (from sequential run)
 * to verify correctness of parallel implementation.
 * 
 * @param results           Array of zone results
 * @param num_zones         Number of zones
 * @param expected_areas    Array of expected areas (can be NULL)
 * @param tolerance         Relative error tolerance (e.g., 1e-6)
 * 
 * @return 0 if validation passes, -1 if fails
 */
int validate_zone_results(
    zone_result **results,
    int num_zones,
    const double *expected_areas,
    double tolerance
);

/**
 * @brief Export zone results to CSV file
 * 
 * Writes detailed zone results to CSV for analysis and comparison.
 * 
 * @param results    Array of zone results
 * @param num_zones  Number of zones
 * @param filename   Output CSV filename
 * 
 * @return 0 on success, -1 on error
 */
int export_zone_results_csv(
    zone_result **results,
    int num_zones,
    const char *filename
);

/*─────────────────────────────────────────────────────────────────────────
 * Helper Functions
 *─────────────────────────────────────────────────────────────────────────*/

/**
 * @brief Get number of points in a zone boundary
 * 
 * Thread-safe function to count total boundary points in a zone.
 * 
 * @param b  Boundary structure
 * @return Total number of points
 */
int num_points_in_zone(boundary *b);

/**
 * @brief Check if zone processing should be skipped
 * 
 * Some zones may be empty or invalid - this function checks validity.
 * 
 * @param b  Boundary structure
 * @return 1 if zone is valid, 0 if should skip
 */
int is_zone_valid(boundary *b);

/*─────────────────────────────────────────────────────────────────────────
 * Debug and Diagnostics
 *─────────────────────────────────────────────────────────────────────────*/

/* Debug levels (can be ORed together) */
#define DEBUG_ZONE_NONE        0
#define DEBUG_ZONE_INFO        1
#define DEBUG_ZONE_TIMING      2
#define DEBUG_ZONE_MEMORY      4
#define DEBUG_ZONE_VALIDATION  8
#define DEBUG_ZONE_ALL         15

/**
 * @brief Set debug level for zone processor
 * 
 * Controls verbosity of diagnostic output.
 * 
 * @param level  Debug level (combination of DEBUG_ZONE_* flags)
 */
void set_zone_debug_level(int level);

/**
 * @brief Get current debug level
 * 
 * @return Current debug level
 */
int get_zone_debug_level(void);

#endif /* ZONE_PROCESSOR_H */
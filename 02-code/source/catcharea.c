#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sys/time.h"
#include "time.h"
#include <sys/resource.h>
#include <string.h>
/*--------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

#include "area.h"
#include "catchment.h"
#include "file.h"
#include "memory.h"
#include "path.h"
#include "scan.h"
#include "streamline.h"
#include "trapfloat.h"

#include <omp.h>
#include "catcharea.h"
#include "performance_summary.h"

/* External function declarations */
extern void set_inversion_method(int method); /* 0=Parallel, 1=Sequential */
extern int get_inversion_method(void);
/*--------------------------------------------------------*/
/* External function to print performance summary */
/*--------------------------------------------------------*/
extern void print_matrix_performance_summary();
extern void get_memory_usage_kb(long *vmrss_kb, long *vmsize_kb);

/*--------------------------------------------------------*/
extern void set_multiply_method(int method);
extern void set_block_size(int size);
extern void set_dgemm_type(int type);       /* NEW: 0=Hybrid, 1=OpenBLAS */
extern int get_dgemm_type(void);
extern const char* get_dgemm_type_name(void);
extern void print_expected_performance(void);
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{

  // ═══════════════════════════════════════════════════════════
  // STEP 1: Initialize performance tracking
  // ═══════════════════════════════════════════════════════════
  init_performance_summary();

  struct timeval total_start, total_end;
  gettimeofday(&total_start, NULL);
  // ═══════════════════════════════════════════════════════════

  //---------- OPTIMIZE PATCH ------------------
  int multiply_method = 3; // Default: full optimization
  int block_size = 64;     // Default: 64 for large matrices
  int dgemm_type = 1;      // Default: OpenBLAS (NEW)

  if (argc > 5)
    multiply_method = atoi(argv[5]);
  if (argc > 6)
    block_size = atoi(argv[6]);
  if (argc > 7)
    dgemm_type = atoi(argv[7]);  // NEW: Parse dgemm_type

  // Set methods
  set_multiply_method(multiply_method);
  set_block_size(block_size);
  set_dgemm_type(dgemm_type);  // NEW: Set DGEMM type

  printf("  DGEMM Type:           %d (%s)\n", dgemm_type, get_dgemm_type_name());  // NEW
  printf("  Multiply method:      %d ", multiply_method);
  switch (multiply_method)
  {
  case 0:
    printf("(Sequential)\n");
    break;
  case 1:
    printf("(OpenMP)\n");
    break;
  case 2:
    printf("(OpenMP+Cache, B=%d)\n", block_size);
    break;
  case 3:
    printf("(OpenMP+Cache+SIMD, B=%d)\n", block_size);
    break;
  }

  // Print expected performance
  print_expected_performance();

  //---------------------------------------------

  char title[64];
  char data[] = "catchment.txt"; /* step-01 */

  bem_vectors *vectors;
  catchment *c;
  char *buffer;
  double step_size, SCA, C_area;
  int buf_size, i, max_points, num_zones, max_steps, max_streams;
  int inversion_method = 0; /* 0=Parallel (default), 1=Sequential */
  matrix bvv, bcv;
  path **streamlines;
  section mouth;
  coordinates PA, PB;

  struct timeval all_start, all_finish;
  double all_duration;

  struct timeval phase_start, phase_finish;
  double init_time, bem_time, catchment_time, total_time;

  struct rusage r_usage_start, r_usage_end;
  long vmrss_start, vmsize_start, vmrss_end, vmsize_end;

  printf("\n");
  printf("################################################################################\n");
  printf("###                                                                          ###\n");
  printf("###                   CATCHMENT AREA COMPUTATION PROGRAM                    ###\n");
  printf("###              Boundary Element Method - Optimized Version                ###\n");
  printf("###                                                                          ###\n");
  printf("################################################################################\n");
  printf("\n");

  // Get initial memory state
  getrusage(RUSAGE_SELF, &r_usage_start);
  get_memory_usage_kb(&vmrss_start, &vmsize_start);

  printf("SYSTEM CONFIGURATION:\n");
  printf("---------------------\n");
  printf("  OMP_NUM_THREADS:      %s\n", getenv("OMP_NUM_THREADS") ?: "not set");
  printf("  OPENBLAS_NUM_THREADS: %s\n", getenv("OPENBLAS_NUM_THREADS") ?: "not set");
  printf("  OMP_PROC_BIND:        %s\n", getenv("OMP_PROC_BIND") ?: "not set");
  printf("  OMP_PLACES:           %s\n", getenv("OMP_PLACES") ?: "not set");
  printf("\n");

  printf("INITIAL MEMORY STATE:\n");
  printf("---------------------\n");
  printf("  VmRSS:                %.2f MB\n", vmrss_start / 1024.0);
  printf("  VmSize:               %.2f MB\n", vmsize_start / 1024.0);
  printf("\n");

  printf("================================================================================\n");
  printf("Starting computation at: %s", ctime(&(time_t){time(NULL)}));
  printf("================================================================================\n");
  printf("\n");

  gettimeofday(&all_start, NULL);

  /*--------------------------------------------------------*/
  /* PHASE 1: Initialization */
  /*--------------------------------------------------------*/
  printf("================================================================================\n");
  printf("PHASE 1: Initialization and Setup\n");
  printf("================================================================================\n");
  gettimeofday(&phase_start, NULL);

  trap_floating_errors();
  buf_size = 512 * 12 + 1;
  buffer = (char *)malloc(buf_size * sizeof(char));

  num_zones = catchment_zones(data);
  c = create_catchment(num_zones, 30);
  get_catchment(data, c);
  plot_catchment(c, "catchment.out");
  max_points = max_points_in_any_zone(c);

  printf("  Catchment zones:      %d\n", num_zones);
  printf("  Max points in zone:   %d\n", max_points);

  vectors = create_bem_vectors(&bvv, &bcv, max_points);

  max_steps = 10000;
  step_size = 1.0;

  double rm, dr, seta;
  rm = 99.0;
  dr = 0.001;

  if (argc > 1)
    step_size = atof(argv[1]);
  if (argc > 2)
    rm = atof(argv[2]);
  if (argc > 3)
    dr = atof(argv[3]);
  if (argc > 4)
    inversion_method = atoi(argv[4]); /* NEW: 4th argument */

  // ═══════════════════════════════════════════════════════════
  // STEP 2: Parse arguments and set configuration
  // ═══════════════════════════════════════════════════════════
  set_performance_config( // Set system configuration
      multiply_method,    // 0-3
      inversion_method,   // 0=parallel, 1=sequential
      omp_get_max_threads(),
      block_size);
  // ═══════════════════════════════════════════════════════════

  // ═══════════════════════════════════════════════════════════
  // STEP 3: Start Setup phase
  // ═══════════════════════════════════════════════════════════
  struct timeval setup_start, setup_end;
  gettimeofday(&setup_start, NULL);
  // ═══════════════════════════════════════════════════════════

  /* Set the inversion method */
  set_inversion_method(inversion_method);

  printf("  Step size:            %.3f\n", step_size);
  printf("  Rm:                   %.3f\n", rm);
  printf("  Dr:                   %.6f\n", dr);
  printf("  Max steps:            %d\n", max_steps);
  printf("  Inversion method:     %s\n",
         inversion_method ? "SEQUENTIAL" : "PARALLEL");
  printf("\n");

  // mouth-01
  put_section("P(0) = (581559.0,943674.0)  P(4) = (581743.0,943675.0)", &mouth);
  show_section(&mouth);
  max_streams = mouth.n;

  printf("  Max streams:          %d\n", max_streams);

  streamlines = (path **)malloc(max_streams * sizeof(path *));
  for (i = 0; i < max_streams; i++)
  {
    streamlines[i] = create_path(max_steps, 1, 0);
  }

  gettimeofday(&phase_finish, NULL);
  init_time = ((double)(phase_finish.tv_sec - phase_start.tv_sec) * 1000000 +
               (double)(phase_finish.tv_usec - phase_start.tv_usec)) /
              1000000;

  get_memory_usage_kb(&vmrss_end, &vmsize_end);
  printf("\n  Initialization time:  %.6f seconds\n", init_time);
  printf("  Memory after init:    VmRSS=%.2f MB, VmSize=%.2f MB\n",
         vmrss_end / 1024.0, vmsize_end / 1024.0);
  printf("================================================================================\n\n");

  // ═══════════════════════════════════════════════════════════
  // STEP 3: End Setup phase
  // ═══════════════════════════════════════════════════════════
  gettimeofday(&setup_end, NULL);
  double setup_time = (setup_end.tv_sec - setup_start.tv_sec) +
                      (setup_end.tv_usec - setup_start.tv_usec) / 1000000.0;
  update_setup_time(setup_time);

  // Set problem parameters
  set_problem_parameters(step_size, rm, dr, num_zones, max_points);

  // Track initial memory
  long vmrss_kb, vmsize_kb;
  get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
  update_memory_usage(vmrss_kb, vmsize_kb);
  // ═══════════════════════════════════════════════════════════

  // ═══════════════════════════════════════════════════════════
  // STEP 4: Start BEM computation phase
  // ═══════════════════════════════════════════════════════════
  struct timeval bem_start, bem_end;
  gettimeofday(&bem_start, NULL);
  // ═══════════════════════════════════════════════════════════

  /*--------------------------------------------------------*/
  /* PHASE 2: BEM Computation (Matrix operations happen here) */
  /*--------------------------------------------------------*/
  printf("================================================================================\n");
  printf("PHASE 2: Boundary Element Method Computation\n");
  printf("================================================================================\n");
  printf("This phase includes matrix setup, multiplication, and inversion.\n");
  printf("Detailed timing will be shown below.\n");
  printf("================================================================================\n\n");

  gettimeofday(&phase_start, NULL);

  C_area = catchment_area(c, &mouth, 0, max_steps, step_size, max_streams,
                          streamlines, vectors); // stream down

  gettimeofday(&phase_finish, NULL);
  bem_time = ((double)(phase_finish.tv_sec - phase_start.tv_sec) * 1000000 +
              (double)(phase_finish.tv_usec - phase_start.tv_usec)) /
             1000000;

  get_memory_usage_kb(&vmrss_end, &vmsize_end);

  printf("\n");
  printf("================================================================================\n");
  printf("PHASE 2 COMPLETED\n");
  printf("================================================================================\n");
  printf("  BEM computation time: %.6f seconds\n", bem_time);
  printf("  Memory after BEM:     VmRSS=%.2f MB, VmSize=%.2f MB\n",
         vmrss_end / 1024.0, vmsize_end / 1024.0);
  printf("================================================================================\n\n");

  // ═══════════════════════════════════════════════════════════
  // STEP 4: End BEM computation phase
  // ═══════════════════════════════════════════════════════════
  gettimeofday(&bem_end, NULL);
  bem_time = (bem_end.tv_sec - bem_start.tv_sec) +
             (bem_end.tv_usec - bem_start.tv_usec) / 1000000.0;
  update_bem_time(bem_time);

  // Track peak memory
  get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
  update_memory_usage(vmrss_kb, vmsize_kb);
  // ═══════════════════════════════════════════════════════════

  // ═══════════════════════════════════════════════════════════
  // STEP 5: Start Finalization phase
  // ═══════════════════════════════════════════════════════════
  struct timeval final_start, final_end;
  gettimeofday(&final_start, NULL);
  // ═══════════════════════════════════════════════════════════

  /*--------------------------------------------------------*/
  /* PHASE 3: Post-processing */
  /*--------------------------------------------------------*/
  printf("================================================================================\n");
  printf("PHASE 3: Post-processing and Output\n");
  printf("================================================================================\n");
  gettimeofday(&phase_start, NULL);

  plot_streamlines(c, max_streams, streamlines, "test.out");
  printf("\n  Catchment area:       %.6f\n", C_area);

  for (i = 0; i < max_streams; i++)
  {
    streamlines[i] = destroy_path(streamlines[i]);
  }

  gettimeofday(&phase_finish, NULL);
  catchment_time = ((double)(phase_finish.tv_sec - phase_start.tv_sec) * 1000000 +
                    (double)(phase_finish.tv_usec - phase_start.tv_usec)) /
                   1000000;

  printf("  Post-processing time: %.6f seconds\n", catchment_time);
  printf("================================================================================\n\n");

  /*--------------------------------------------------------*/
  /* Final timing and memory summary */
  /*--------------------------------------------------------*/
  gettimeofday(&all_finish, NULL);
  all_duration = ((double)(all_finish.tv_sec - all_start.tv_sec) * 1000000 +
                  (double)(all_finish.tv_usec - all_start.tv_usec)) /
                 1000000;

  getrusage(RUSAGE_SELF, &r_usage_end);
  get_memory_usage_kb(&vmrss_end, &vmsize_end);

  printf("\n");
  printf("################################################################################\n");
  printf("###                                                                          ###\n");
  printf("###                      OVERALL EXECUTION SUMMARY                          ###\n");
  printf("###                                                                          ###\n");
  printf("################################################################################\n");
  printf("\n");

  printf("TIMING BREAKDOWN:\n");
  printf("-----------------\n");
  printf("  Phase 1 (Initialization):    %10.6f sec (%5.1f%%)\n",
         init_time, init_time / all_duration * 100);
  printf("  Phase 2 (BEM Computation):   %10.6f sec (%5.1f%%)\n",
         bem_time, bem_time / all_duration * 100);
  printf("  Phase 3 (Post-processing):   %10.6f sec (%5.1f%%)\n",
         catchment_time, catchment_time / all_duration * 100);
  printf("  ----------------------------------------------\n");
  printf("  TOTAL EXECUTION TIME:        %10.6f sec\n", all_duration);
  printf("\n");

  printf("MEMORY SUMMARY:\n");
  printf("---------------\n");
  printf("  Initial VmRSS:               %.2f MB\n", vmrss_start / 1024.0);
  printf("  Final VmRSS:                 %.2f MB\n", vmrss_end / 1024.0);
  printf("  Peak VmRSS (delta):          %.2f MB\n", (vmrss_end - vmrss_start) / 1024.0);
  printf("  Initial VmSize:              %.2f MB\n", vmsize_start / 1024.0);
  printf("  Final VmSize:                %.2f MB\n", vmsize_end / 1024.0);
  printf("  Peak VmSize (delta):         %.2f MB\n", (vmsize_end - vmsize_start) / 1024.0);
  printf("  Max RSS (rusage):            %.2f MB\n", r_usage_end.ru_maxrss / 1024.0);
  printf("\n");

  printf("CPU USAGE:\n");
  printf("----------\n");
  printf("  User CPU time:               %.6f sec\n",
         (double)r_usage_end.ru_utime.tv_sec + (double)r_usage_end.ru_utime.tv_usec / 1000000);
  printf("  System CPU time:             %.6f sec\n",
         (double)r_usage_end.ru_stime.tv_sec + (double)r_usage_end.ru_stime.tv_usec / 1000000);
  printf("  Total CPU time:              %.6f sec\n",
         (double)(r_usage_end.ru_utime.tv_sec + r_usage_end.ru_stime.tv_sec) +
             (double)(r_usage_end.ru_utime.tv_usec + r_usage_end.ru_stime.tv_usec) / 1000000);
  printf("\n");

  /* Call the matrix performance summary function */
  print_matrix_performance_summary();

  printf("################################################################################\n");
  printf("Computation completed at: %s", ctime(&(time_t){time(NULL)}));
  printf("################################################################################\n");
  printf("\n");

  /*---------------------------------------------------*/

  //----------- Performance End phase -----------------------
  // gettimeofday(&end, NULL);
  // double setup_time = (end.tv_sec - start.tv_sec) +
  //                    (end.tv_usec - start.tv_usec) / 1000000.0;
  // update_setup_time(setup_time);

  // ═══════════════════════════════════════════════════════════
  // STEP 5: End Finalization phase
  // ═══════════════════════════════════════════════════════════
  gettimeofday(&final_end, NULL);
  double final_time = (final_end.tv_sec - final_start.tv_sec) +
                      (final_end.tv_usec - final_start.tv_usec) / 1000000.0;
  update_finalization_time(final_time);
  // ═══════════════════════════════════════════════════════════

  /*--------------------------------------------------------*/

  destroy_catchment(c);
  destroy_bem_vectors(vectors);
  free((void *)buffer);

  // ═══════════════════════════════════════════════════════════
  // STEP 6: Print summary
  // ═══════════════════════════════════════════════════════════
  printf("\n\n");

  // Print comprehensive summary
  print_performance_summary();

  // Print journal table (LaTeX + CSV)
  print_journal_table();

  // Export to CSV file (optional)
  export_performance_csv("performance_results.csv");
  /* Legacy compatibility for older benchmark scripts:
   echo an easily greppable one-line summary of matrix multiply time. */
  {
    FILE *pf = fopen("performance_results.csv", "r");
    if (pf)
    {
      char line[256];
      double mt = -1.0;
      while (fgets(line, sizeof(line), pf))
      {
        if (sscanf(line, "Multiply_Time_sec,%lf", &mt) == 1)
          break;
      }
      fclose(pf);
      if (mt >= 0.0)
      {
        printf("Total Matrix Multiply time: %.6f s\n", mt);
      }
    }
  }

  return (0);
}
/*--------------------------------------------------------*/

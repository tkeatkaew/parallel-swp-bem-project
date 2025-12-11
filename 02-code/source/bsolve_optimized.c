/*--------------------------------- bsolve.c ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for solving boundary value problem                                      */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sys/time.h"
#include "time.h"
#include <sys/resource.h>
#include <string.h>

/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"

#include "co_matrix.h"
#include "linear_sys.h"
#include "matrix.h"
#include "path.h"
#include "ten_matrix.h"

#include "bsolve.h"

/*----------------------------------------------------------------------------------*/
/* External function to print matrix performance summary */
/*----------------------------------------------------------------------------------*/
extern void print_matrix_performance_summary();

/*----------------------------------------------------------------------------------*/
/* Helper function to get memory usage */
/*----------------------------------------------------------------------------------*/
static void get_memory_usage_kb(long *vmrss_kb, long *vmsize_kb) {
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
/*  make boundary voltage vector */
/*----------------------------------------------------------------------------------*/
void make_boundary_voltage_vector(b,bvv)
     boundary *b;
     matrix *bvv;
{
  int i,paths;
  path *path_i;
  double *result;
  int poi;

  if(b->bvv==(double *)NULL)
    {
      bvv->value=(double *)malloc(get_num_rows(bvv)*sizeof(double));
      b->bvv=bvv->value;
      paths = b->components;
      result = bvv->value;
            for(i=0;i<paths;i++)
	{ 
	  path_i=b->loop[i];      
	  fill_boundary_voltage_vector(path_i, result);
	  poi=path_i->points;
	  result=result+(2*poi);
	}
    }
  else
    {
      bvv->value=b->bvv; /* use array calculated last time */
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill boundary voltage in vector */
/*----------------------------------------------------------------------------------*/
void fill_boundary_voltage_vector(path_i, result)

  path *path_i;
  double *result;
{
  double v1,v2;
  double V, W;
  int j,poi;
  
  poi=path_i->points;
  for(j=0;j<poi;j++)
    {
      v1=get_path_value(path_i,j);
      v2=get_path_value(path_i,j+1);
      V=v2-v1;
      W=(v2+v1)/2.0;
      p2c_2coeff(V,W,&V,&W);
      result[2*j]=V;
      result[2*j+1]=W;
    } 
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*  make boundary current vector (final ?? version) */
/*----------------------------------------------------------------------------------*/
void make_boundary_current_vector(b,V,J)
     matrix *V; /* bvv boundary voltage vector */
     matrix *J; /* bcv boundary current vector */
     boundary *b;
{
  int j,finite;

  finite=0;
  for(j=0;j<b->components;j++)
    {
      if(b->level[j]==0) finite=1; /* the zone is finite; it does not go to infinity */
    }
  if(finite==1)
    {
      make_bcv_use_KCL(b,V,J); /* enforce KCL explicitly */
    }
  else
    {
      make_bcv_no_KCL(b,V,J); /* do not enforce KCL */
    }
}

/*----------------------------------------------------------------------------------*/
/*  make boundary current vector (third version) */
/*----------------------------------------------------------------------------------*/
void make_bcv_no_KCL(b,V,J)
     matrix *V; /* bvv boundary voltage vector */
     matrix *J; /* bcv boundary current vector */
     boundary *b;
{
  matrix A, B, D, DA, DAV;
  matrix BT, BTB, BTDAV;
  double *values;
  int N,j;

  if(b->bcv==(double *)NULL)
    {
      J->value=(double *)malloc(get_num_rows(J)*sizeof(double));
      b->bcv=J->value;
      N=0;
      for(j=0;j<b->components;j++)
	{
	  N=N+b->loop[j]->points;
	}
      values=(double *)malloc(9*N*(4*N+1)*sizeof(double));

      /*-----------------------------*/
      attach_matrix(&A,5*N,2*N,values);
      attach_matrix(&D,5*N,2*N,after_matrix(&A)); 
      attach_matrix(&DA,5*N,2*N,values);            /* write on top of A */
      attach_matrix(&DAV,5*N,1,after_matrix(&D)); 
      
      printf("\nmake voltage matrix A and D\n");
      make_voltage_geometry_matrix(b,&A);
      make_diagonal_matrix(b,&D);  
      add_matrix(&D,&A,&DA); 
      multiply_matrix(&DA,V,&DAV);

      /*-----------------------------*/
      attach_matrix(&B,5*N,4*N,values);            /* write on top of A,D,DA */  
      attach_matrix(&BT,5*N,4*N,values);           /* share data of B */
      attach_matrix(&BTB,4*N,4*N,after_matrix(&DAV));
      attach_matrix(&BTDAV,4*N,1,after_matrix(&BTB));

      printf("make current matrix B\n");
      make_current_geometry_matrix(b,&B);
      transpose_matrix(&BT,&BT);                   /* makes transpose but does not destroy B */
      
      multiply_matrix(&BT,&B,&BTB);
      invert_matrix(&BTB,&BTB);
      multiply_matrix(&BT,&DAV,&BTDAV);
      multiply_matrix(&BTB,&BTDAV,J);
      free((void *)values);
    }
  else
    {
      J->value=b->bcv; /* use array calculated last time */
    }
}

/*----------------------------------------------------------------------------------*/
/*  make boundary current vector (fourth version: imposes KCL) */
/*  ENHANCED VERSION with detailed timing and memory tracking */
/*----------------------------------------------------------------------------------*/
void make_bcv_use_KCL(b,V,J)
     matrix *V; /* bvv boundary voltage vector */
     matrix *J; /* bcv boundary current vector */
     boundary *b;
{
  matrix A, B, D, DA, DAV;
  matrix BT, BTB, BTDAV, KCL;
  double *values;
  int N,j;

  long vmrss, vmsize;
  struct rusage r_usage;

  struct timeval start,finish;
  double duration;

  struct timeval inv_start,inv_finish;
  double inv_duration;

  struct timeval mul_start1,mul_finish1,mul_start2,mul_finish2;
  double mul_duration1, mul_duration2;

  struct timeval all_start,all_finish;
  double all_duration;

  // Phase timing
  struct timeval phase_start, phase_finish;
  double phase1_time, phase2_time, phase3_time, phase4_time;

  if(b->bcv==(double *)NULL)
    {
      printf("\n");
      printf("================================================================================\n");
      printf("           STARTING BOUNDARY CURRENT VECTOR COMPUTATION (KCL)\n");
      printf("================================================================================\n");
      
      gettimeofday(&all_start, NULL);
      
      J->value=(double *)malloc(get_num_rows(J)*sizeof(double));
      b->bcv=J->value;
      N=0;
      for(j=0;j<b->components;j++)
	{
	  N=N+b->loop[j]->points;
	}
      
      printf("\nProblem size: N = %d boundary points\n", N);
      printf("Matrix dimensions:\n");
      printf("  A, D, DA:  %d x %d\n", 5*N+1, 2*N);
      printf("  B, BT:     %d x %d\n", 5*N+1, 4*N);
      printf("  BTB:       %d x %d\n", 4*N, 4*N);
      
      size_t memory_required = (size_t)(9*N+1) * (4*N+1) * sizeof(double);
      printf("\nAllocating %.2f MB for computation matrices\n", memory_required / (1024.0*1024.0));
      
      get_memory_usage_kb(&vmrss, &vmsize);
      printf("Memory before allocation: VmRSS=%.2f MB, VmSize=%.2f MB\n", 
             vmrss/1024.0, vmsize/1024.0);
      
      values=(double *)malloc(memory_required);
      if (values == NULL) {
          printf("ERROR: Failed to allocate memory!\n");
          exit(1);
      }
      
      get_memory_usage_kb(&vmrss, &vmsize);
      printf("Memory after allocation: VmRSS=%.2f MB, VmSize=%.2f MB\n\n", 
             vmrss/1024.0, vmsize/1024.0);

/*---------------------------------------------------*/
printf("================================================================================\n");
printf("PHASE 1: Voltage Geometry Matrix Setup\n");
printf("================================================================================\n");
gettimeofday(&phase_start, NULL);   

      /*-----------------------------*/
      attach_matrix(&A,5*N+1,2*N,values);
      attach_matrix(&D,5*N+1,2*N,after_matrix(&A)); 
      attach_matrix(&DA,5*N+1,2*N,values);            /* write on top of A */
      attach_matrix(&DAV,5*N+1,1,after_matrix(&D)); 
      
      gettimeofday(&start, NULL);
      make_voltage_geometry_matrix(b,&A);
      gettimeofday(&finish, NULL);
      duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
      printf("  make_voltage_geometry_matrix: %.6f sec\n", duration);
      
      zero_last_matrix_row(&A);
      
      gettimeofday(&start, NULL);
      make_diagonal_matrix(b,&D);  
      gettimeofday(&finish, NULL);
      duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
      printf("  make_diagonal_matrix: %.6f sec\n", duration);
      
      zero_last_matrix_row(&D);
      
      gettimeofday(&start, NULL);
      add_matrix(&D,&A,&DA); 
      gettimeofday(&finish, NULL);
      duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
      printf("  add_matrix: %.6f sec\n", duration);
      
      gettimeofday(&start, NULL);
      multiply_matrix(&DA,V,&DAV);
      gettimeofday(&finish, NULL);
      duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
      printf("  multiply_matrix (DA*V): %.6f sec\n", duration);
      
gettimeofday(&phase_finish, NULL);
phase1_time = ((double)(phase_finish.tv_sec-phase_start.tv_sec)*1000000 + (double)(phase_finish.tv_usec-phase_start.tv_usec)) / 1000000;
printf("PHASE 1 Total: %.6f seconds\n\n", phase1_time);

/*---------------------------------------------------*/
printf("================================================================================\n");
printf("PHASE 2: Current Geometry Matrix Setup\n");
printf("================================================================================\n");
gettimeofday(&phase_start, NULL);   

      /*-----------------------------*/
      attach_matrix(&B,5*N+1,4*N,values);            /* write on top of A,D,DA */  
      attach_matrix(&BT,5*N+1,4*N,values);           /* share data of B */
      attach_matrix(&BTB,4*N,4*N,after_matrix(&DAV));
      attach_matrix(&BTDAV,4*N,1,after_matrix(&BTB));
      attach_matrix(&KCL,1,4*N,after_matrix(&BTB));  /* KCL comes after BTB */

gettimeofday(&start, NULL);  
      make_current_geometry_matrix(b,&B);
gettimeofday(&finish, NULL);
duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
printf("  make_current_geometry_matrix: %.6f sec\n", duration);

gettimeofday(&start, NULL);  
      make_kcl_geometry_vector(b,&KCL);
gettimeofday(&finish, NULL);
duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
printf("  make_kcl_geometry_vector: %.6f sec\n", duration);

gettimeofday(&start, NULL);  
      fill_last_matrix_row(&B,&KCL);
gettimeofday(&finish, NULL);
duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
printf("  fill_last_matrix_row: %.6f sec\n", duration);

gettimeofday(&start, NULL);  
      transpose_matrix(&BT,&BT);
gettimeofday(&finish, NULL);
duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
printf("  transpose_matrix: %.6f sec\n", duration);

gettimeofday(&phase_finish, NULL);
phase2_time = ((double)(phase_finish.tv_sec-phase_start.tv_sec)*1000000 + (double)(phase_finish.tv_usec-phase_start.tv_usec)) / 1000000;
printf("PHASE 2 Total: %.6f seconds\n\n", phase2_time);

/*---------------------------------------------------*/
printf("================================================================================\n");
printf("PHASE 3: Matrix Multiplication (BT * B)\n");
printf("================================================================================\n");

  printf("  Matrix B info:\n    ");
  show_matrix_info(&B);
  printf("  Matrix BT info:\n    ");
  show_matrix_info(&BT);

  get_memory_usage_kb(&vmrss, &vmsize);
  printf("  Memory before BTB multiply: VmRSS=%.2f MB, VmSize=%.2f MB\n", 
         vmrss/1024.0, vmsize/1024.0);

gettimeofday(&mul_start1, NULL);        
      multiply_matrix(&BT,&B,&BTB);
gettimeofday(&mul_finish1, NULL);
mul_duration1 = ((double)(mul_finish1.tv_sec-mul_start1.tv_sec)*1000000 + (double)(mul_finish1.tv_usec-mul_start1.tv_usec)) / 1000000;

  get_memory_usage_kb(&vmrss, &vmsize);
  printf("  Memory after BTB multiply: VmRSS=%.2f MB, VmSize=%.2f MB\n", 
         vmrss/1024.0, vmsize/1024.0);

  // Calculate FLOPS for BT*B multiplication
  long long flops_btb = 2LL * (5*N+1) * 4*N * 4*N;
  double gflops_btb = (double)flops_btb / mul_duration1 / 1.0e9;

printf("  multiply_matrix (BT*B): %.6f sec (%.2f GFLOPS)\n", mul_duration1, gflops_btb);
printf("PHASE 3 Total: %.6f seconds\n\n", mul_duration1);

/*---------------------------------------------------*/
printf("================================================================================\n");
printf("PHASE 4: Matrix Inversion\n");
printf("================================================================================\n");

  get_memory_usage_kb(&vmrss, &vmsize);
  printf("  Memory before inversion: VmRSS=%.2f MB, VmSize=%.2f MB\n", 
         vmrss/1024.0, vmsize/1024.0);

gettimeofday(&inv_start, NULL);   
      invert_matrix(&BTB,&BTB);
gettimeofday(&inv_finish, NULL);
inv_duration = ((double)(inv_finish.tv_sec-inv_start.tv_sec)*1000000 + (double)(inv_finish.tv_usec-inv_start.tv_usec)) / 1000000;

  get_memory_usage_kb(&vmrss, &vmsize);
  printf("  Memory after inversion: VmRSS=%.2f MB, VmSize=%.2f MB\n", 
         vmrss/1024.0, vmsize/1024.0);

printf("PHASE 4 Total: %.6f seconds\n\n", inv_duration);

/*---------------------------------------------------*/
printf("================================================================================\n");
printf("PHASE 5: Final Matrix Multiplications\n");
printf("================================================================================\n");

gettimeofday(&mul_start2, NULL);    

      gettimeofday(&start, NULL);
      multiply_matrix(&BT,&DAV,&BTDAV);
      gettimeofday(&finish, NULL);
      duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
      printf("  multiply_matrix (BT*DAV): %.6f sec\n", duration);
      
      gettimeofday(&start, NULL);
      multiply_matrix(&BTB,&BTDAV,J);
      gettimeofday(&finish, NULL);
      duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
      printf("  multiply_matrix (BTB*BTDAV): %.6f sec\n", duration);

gettimeofday(&mul_finish2, NULL);
mul_duration2 = ((double)(mul_finish2.tv_sec-mul_start2.tv_sec)*1000000 + (double)(mul_finish2.tv_usec-mul_start2.tv_usec)) / 1000000;
printf("PHASE 5 Total: %.6f seconds\n\n", mul_duration2);

/*---------------------------------------------------*/
gettimeofday(&all_finish, NULL);
all_duration = ((double)(all_finish.tv_sec-all_start.tv_sec)*1000000 + (double)(all_finish.tv_usec-all_start.tv_usec)) / 1000000;

printf("================================================================================\n");
printf("                    BOUNDARY COMPUTATION SUMMARY\n");
printf("================================================================================\n");
printf("\nTIMING BREAKDOWN:\n");
printf("  Phase 1 (Voltage Setup):        %10.6f sec (%5.1f%%)\n", phase1_time, phase1_time/all_duration*100);
printf("  Phase 2 (Current Setup):        %10.6f sec (%5.1f%%)\n", phase2_time, phase2_time/all_duration*100);
printf("  Phase 3 (BT*B multiply):        %10.6f sec (%5.1f%%)\n", mul_duration1, mul_duration1/all_duration*100);
printf("  Phase 4 (Matrix Inversion):     %10.6f sec (%5.1f%%)\n", inv_duration, inv_duration/all_duration*100);
printf("  Phase 5 (Final multiplies):     %10.6f sec (%5.1f%%)\n", mul_duration2, mul_duration2/all_duration*100);
printf("  -----------------------------------------------------------\n");
printf("  TOTAL:                          %10.6f sec\n\n", all_duration);

printf("OPERATION TOTALS:\n");
printf("  All Matrix Multiplications:     %10.6f sec (%5.1f%%)\n", 
       mul_duration1+mul_duration2, (mul_duration1+mul_duration2)/all_duration*100);
printf("  Matrix Inversion:               %10.6f sec (%5.1f%%)\n", 
       inv_duration, inv_duration/all_duration*100);
printf("  Other Operations:               %10.6f sec (%5.1f%%)\n", 
       all_duration-(mul_duration1+mul_duration2+inv_duration),
       (all_duration-(mul_duration1+mul_duration2+inv_duration))/all_duration*100);

getrusage(RUSAGE_SELF,&r_usage);
get_memory_usage_kb(&vmrss, &vmsize);

printf("\nMEMORY USAGE:\n");
printf("  VmRSS (resident):               %.2f MB\n", vmrss/1024.0);
printf("  VmSize (virtual):               %.2f MB\n", vmsize/1024.0);
printf("  Max RSS:                        %.2f MB\n", r_usage.ru_maxrss/1024.0);
printf("================================================================================\n\n");

      free((void *)values);

    }
  else
    {
      J->value=b->bcv; /* use array calculated last time */
    }
}

/*----------------------------------------------------------------------------------*/
/*  make boundary vector : bvv , bcv */
/*----------------------------------------------------------------------------------*/
void make_boundary_vector(b,bvv,bcv)
     boundary *b;
     matrix *bvv, *bcv;
{
  make_boundary_voltage_vector(b,bvv);
  make_boundary_current_vector(b,bvv,bcv);
}
  
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*  make boundary internal voltage at point P */
/*----------------------------------------------------------------------------------*/
double make_internal_voltage(b,bvv,bcv,P,vgv,cgv)
     boundary *b;
     matrix *bvv, *bcv, *vgv, *cgv;
     coordinates P;
{
  matrix temp1;
  double voltage,v1,v2;

  temp1.transpose=0;
  temp1.invert=0;   
  temp1.rows=1;      
  temp1.columns=1;   
  temp1.value=&voltage;

  make_voltage_geometry_vector(P,b,vgv);
  make_current_geometry_vector(P,b,cgv);

  multiply_matrix(vgv,bvv,&temp1);
  v1=voltage;
  multiply_matrix(cgv,bcv,&temp1);
  v2=voltage;
  voltage= v2-v1;

  return(voltage);
}
/*----------------------------------------------------------------------------------*/
/*  make boundary internal grad_voltage at point P */
/*----------------------------------------------------------------------------------*/
void make_internal_grad_voltage(b,bvv,bcv,P,co_vgv,co_cgv,Gv)
     boundary *b;
     matrix *bvv, *bcv;
     co_matrix *co_vgv, *co_cgv;
     coordinates P,Gv;
{
  co_matrix temp1;
  coordinates voltage, v1 ,v2;

  temp1.transpose=0;
  temp1.invert=0;   
  temp1.rows=1;      
  temp1.columns=1;   
  temp1.value=&voltage;

  /*-------------------------------------*/
  /*------- trying myself ---------------*/   
  make_co_voltage_geometry_vector(P,b,co_vgv);
  make_co_current_geometry_vector(P,b,co_cgv);
  /*-------------------------------------*/

  multiply_co_matrix(co_vgv,bvv,&temp1);
  v1[0]=voltage[0];  v1[1]=voltage[1];

  multiply_co_matrix(co_cgv,bcv,&temp1);
  v2[0]=voltage[0];  v2[1]=voltage[1];

  Gv[0]=v2[0]-v1[0]; Gv[1]=v2[1]-v1[1];  

  Gv[0]=Gv[0]/(2.0*M_PI);
  Gv[1]=Gv[1]/(2.0*M_PI);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*  make boundary internal sec_grad_voltage at point P */
/*----------------------------------------------------------------------------------*/
void make_internal_sec_grad_voltage(b,bvv,bcv,P,ten_vgv,ten_cgv,Gv)
     boundary *b;
     matrix *bvv, *bcv;
     ten_matrix *ten_vgv, *ten_cgv;
     coordinates P;
     tensor Gv;
{
  ten_matrix temp1;
  tensor voltage, v1 ,v2;

  temp1.transpose=0;
  temp1.invert=0;   
  temp1.rows=1;      
  temp1.columns=1;   
  temp1.value=&voltage;

  /*-------------------------------------*/
  /*------- trying myself ---------------*/   
  make_ten_voltage_geometry_vector(P,b,ten_vgv);
  make_ten_current_geometry_vector(P,b,ten_cgv);
  /*-------------------------------------*/

  multiply_ten_matrix(ten_vgv,bvv,&temp1);
  v1[0][0]=voltage[0][0];  v1[0][1]=voltage[0][1];
  v1[1][0]=voltage[1][0];  v1[1][1]=voltage[1][1];

  multiply_ten_matrix(ten_cgv,bcv,&temp1);
  v2[0][0]=voltage[0][0];  v2[0][1]=voltage[0][1];
  v2[1][0]=voltage[1][0];  v2[1][1]=voltage[1][1];

  Gv[0][0]=v2[0][0]-v1[0][0];
  Gv[0][1]=v2[0][1]-v1[0][1];
  Gv[1][0]=v2[1][0]-v1[1][0];
  Gv[1][1]=v2[1][1]-v1[1][1];

  Gv[0][0]=Gv[0][0]/(2.0*M_PI);
  Gv[0][1]=Gv[0][1]/(2.0*M_PI);
  Gv[1][0]=Gv[1][0]/(2.0*M_PI);
  Gv[1][1]=Gv[1][1]/(2.0*M_PI);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

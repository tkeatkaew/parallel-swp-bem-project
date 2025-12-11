/*--------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
/*--------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

#include "scan.h"
#include "streamline.h"

#include "area.h"
/*--------------------------------------------------------*/
/*----------- catchment area loop ------------------------*/
/*--------------------------------------------------------*/
double catchment_area(c,mouth,direction,max_steps,step_size,
	     n_stream,streamline,vectors) 
     catchment *c;
     section *mouth;
     int direction; /* 1 = go to max; 0 = go to min */
     int max_steps; /* +ve = number of steps; -ve = don't check */
     double step_size;
     int n_stream;
     path **streamline;
     bem_vectors *vectors;
{
  bem_results R;
  coordinates P;
  double dx,dy,dw,C_sum;
  double L_old,L_new,s_theta_old,s_theta_new;
  double cosq_theta;  
  int i,n,k;
  
  i=0;
  C_sum=0.0;
  n=mouth->n-1;     /* number of steps = number of points - 1 */
  dx=(mouth->P2[0]-mouth->P1[0])/(double)n;
  dy=(mouth->P2[1]-mouth->P1[1])/(double)n;
  dw=mouth->step;   /* step size across mouth */
  xy_section(mouth,0,P);

  L_old=streamline_loop(P,c,direction,max_steps,step_size,streamline[0],vectors,&R) ;
  cosq_theta=(dx*R.dV[0]+dy*R.dV[1])/dw;
  cosq_theta=cosq_theta*cosq_theta/(R.dV[0]*R.dV[0]+R.dV[1]*R.dV[1]);
  if(cosq_theta>1.0) cosq_theta=1.0;
  s_theta_old=sqrt(1.0-cosq_theta);
  C_sum=0.0;
  k=1;

  for(i=1;i<mouth->n;i++)
    {
      xy_section(mouth,i,P);
      if(i*(n_stream-1)>=k*n){
	L_new=streamline_loop(P,c,direction,max_steps,step_size,streamline[k],vectors,&R);
	k=k+1; }
      else{
	L_new=streamline_loop(P,c,direction,max_steps,step_size,streamline[k],vectors,&R); }
      cosq_theta=(dx*R.dV[0]+dy*R.dV[1])/dw;
      cosq_theta=cosq_theta*cosq_theta/(R.dV[0]*R.dV[0]+R.dV[1]*R.dV[1]);
      if(cosq_theta>1.0) cosq_theta=1.0;
      s_theta_new=sqrt(1.0-cosq_theta);
      C_sum=C_sum+L_old*s_theta_old+L_new*s_theta_new;
      L_old=L_new;
      s_theta_old=s_theta_new;
    }
  C_sum=C_sum*dw/2.0;
  printf("\n dw=%f C_sum=%f",dw,C_sum);
  return(C_sum);
}
/*--------------------------------------------------------*/
/*--------------------- SCA index ------------------------*/
/*--------------------------------------------------------*/
double Cal_SCA(c,mouth,direction,max_steps,step_size,
	     n_stream,streamline,vectors) 
     catchment *c;
     section *mouth;
     int direction; /* 1 = go to max; 0 = go to min */
     int max_steps; /* +ve = number of steps; -ve = don't check */
     double step_size;
     int n_stream;
     path **streamline;
     bem_vectors *vectors;
{
  bem_results R;
  coordinates P;
  double dx,dy,dw,SCA;
  int i,n,k;

  i=0;
  SCA=0.0;

  n=mouth->n-1;
  dx=(mouth->P2[0]-mouth->P1[0])/(double)n;
  dy=(mouth->P2[1]-mouth->P1[1])/(double)n;
  dw=mouth->step;
  xy_section(mouth,0,P);

  //SCA=SCA_loop_GH0_v2(P,c,direction,max_steps,step_size,streamline[0],vectors,&R);
  //SCA=SCA_loop_GH0_v3(P,c,direction,max_steps,step_size,streamline[0],vectors,&R);
  //SCA=SCA_loop_GH0_v4(P,c,direction,max_steps,step_size,streamline[0],vectors,&R);
  SCA=SCA_loop_GH0_v5(P,c,direction,max_steps,step_size,streamline[0],vectors,&R);
  //printf("%f\t%f\t%f---",P[0],P[1],SCA);
  
  return(SCA);
}
/*--------------------------------------------------------*/

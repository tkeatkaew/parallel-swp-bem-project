/*--------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
/*--------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

#include "deep.h"
#include "flow.h"
#include "scan.h"
#include "streamline.h"

#include "mouthflow.h"
/*--------------------------------------------------------*/
/*          flow rate through mouth of catchment          */
/*--------------------------------------------------------*/
double flow_rate(c,mouth,direction,max_steps,step_size,
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
  double dx,dy,dw,dQ_sum;
  double dQ_old,dQ_new,s_theta_old,s_theta_new;
  double cosq_theta,L,d,Q;  
  int i,n,k;
  
  printf("mouth of streamline from (%f, %f) to (%f, %f)\n",
	 mouth->P1[0],mouth->P1[1],mouth->P2[0],mouth->P2[1]);

  dQ_sum=0.0;
  n=mouth->n-1;     /* number of steps = number of points - 1 */
  dx=(mouth->P2[0]-mouth->P1[0])/(double)n;
  dy=(mouth->P2[1]-mouth->P1[1])/(double)n;
  dw=mouth->step;   /* step size across mouth */
  xy_section(mouth,0,P);
  printf("start streamline %2d at (%f, %f)",0,P[0],P[1]);
  L=streamline_loop(P,c,direction,max_steps,step_size,streamline[0],
		    vectors,&R) ;
  d=depth(P,L,R.dV);
  Q=current_density(P,R.dV);
  dQ_old=d*Q;

  cosq_theta=(dx*R.dV[0]+dy*R.dV[1])/dw;
  cosq_theta=cosq_theta*cosq_theta/(R.dV[0]*R.dV[0]+R.dV[1]*R.dV[1]);
  if(cosq_theta>1.0) cosq_theta=1.0;
  s_theta_old=sqrt(1.0-cosq_theta);
  dQ_sum=0.0;
  k=1;

  for(i=1;i<mouth->n;i++)
    {
      xy_section(mouth,i,P);
      printf("\nstart streamline %2d at (%f, %f)",i,P[0],P[1]);
      if(i*(n_stream-1)>=k*n)
	{
	  L=streamline_loop(P,c,direction,max_steps,step_size,streamline[k],
			vectors,&R) ;
	  k=k+1;
	}
      else
	{
	  L=streamline_loop(P,c,direction,max_steps,step_size,(path *)NULL,
			vectors,&R) ;
	}
      d=depth(P,L,R.dV);
      Q=current_density(P,R.dV);
      dQ_new=d*Q;

      cosq_theta=(dx*R.dV[0]+dy*R.dV[1])/dw;
      cosq_theta=cosq_theta*cosq_theta/(R.dV[0]*R.dV[0]+R.dV[1]*R.dV[1]);
      if(cosq_theta>1.0) cosq_theta=1.0;
      s_theta_new=sqrt(1.0-cosq_theta);

      dQ_sum=dQ_sum+dQ_old*s_theta_old+dQ_new*s_theta_new;
      dQ_old=dQ_new;
      s_theta_old=s_theta_new;
    }
  dQ_sum=dQ_sum*dw/2.0;
  return(dQ_sum);
}
/*--------------------------------------------------------*/

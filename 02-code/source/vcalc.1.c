/*---------------------------------- vcalc.c ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for calculating voltage and gradients at point P                        */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

#include "bsolve.h"
#include "catchment.h"
#include "co_matrix.h"
#include "matrix.h"
#include "path.h"
#include "ten_matrix.h"

#include "vcalc.h"
/*----------------------------------------------------------------------------------*/
double voltage_on_path(c,s,segment,this_path)
     catchment *c;
     double s;
     int segment;
     path *this_path;
{
  double v1,v2,vp;

  v1=get_path_value(this_path,segment);
  v2=get_path_value(this_path,segment+1);

  vp=(v2-v1)*s+(v1+v2)*0.5;
  /*
  printf("P on path\n"); 
  printf("vp=%f",vp);
  */
  return(vp);
}

/*----------------------------------------------------------------------------------*/
double voltage_outside_catchment()
{
  /* printf("P not inside this catchment\n"); */
  return(0.0);
}
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
double calculate_in_same_zone(b,P,x,R)
     boundary *b;
     coordinates P;
     bem_vectors *x;
     bem_results *R;
{
  reverse_zone(b);
  R->V=make_internal_voltage(b,x->bvv,x->bcv,P,x->vgv,x->cgv);
  make_internal_grad_voltage(b,x->bvv,x->bcv,P,x->co_vgv,x->co_cgv,R->dV);
  make_internal_sec_grad_voltage(b,x->bvv,x->bcv,P,x->ten_vgv,x->ten_cgv,R->d2V);
  reverse_zone(b);       

  return(R->V);
}
/*----------------------------------------------------------------------------------*/
double calculate_in_new_zone(b,P,x,R)
     boundary *b;
     coordinates P;
     bem_vectors *x;
     bem_results *R;
{
  int N,k;

  N=0;
  for(k=0;k<b->components;k++)  N=N+b->loop[k]->points;

  attach_matrix(x->bvv,2*N,1,(double *)NULL);  
  attach_matrix(x->bcv,4*N,1,(double *)NULL);

  attach_matrix(x->cgv,1,4*N,startof_matrix(x->cgv));           
  attach_matrix(x->vgv,1,2*N,startof_matrix(x->vgv));
  attach_co_matrix(x->co_cgv,1,4*N,startof_co_matrix(x->co_cgv));     
  attach_co_matrix(x->co_vgv,1,2*N,startof_co_matrix(x->co_vgv));
  attach_ten_matrix(x->ten_cgv,1,4*N,startof_ten_matrix(x->ten_cgv));   
  attach_ten_matrix(x->ten_vgv,1,2*N,startof_ten_matrix(x->ten_vgv));

  reverse_zone(b);

  make_boundary_vector(b,x->bvv,x->bcv);

  /*
  show_matrix(x->bvv);
  show_matrix(x->bcv);
  */
  R->V=make_internal_voltage(b,x->bvv,x->bcv,P,x->vgv,x->cgv);
  make_internal_grad_voltage(b,x->bvv,x->bcv,P,x->co_vgv,x->co_cgv,R->dV);
  make_internal_sec_grad_voltage(b,x->bvv,x->bcv,P,x->ten_vgv,x->ten_cgv,R->d2V);

  reverse_zone(b);       

  return(R->V);
}
/*----------------------------------------------------------------------------------*/
/*------calculate on the path-----------------------------*/
/*----------------------------------------------------------------------------------*/
#if 0
void calculate_on_path(state,P,dP,r0)
     int state;
     coordinates P,dP;
     double *r0;
{
  double k=0.0;

  if(state<1)
    {
      printf("\n !warning : the start P is outside catchment\n ");	  
      printf("should to choose the new point P\n");
      exit(0);
    }
  /* return to the last point */
  P[0]=P[0]-dP[0];  
  P[1]=P[1]-dP[1];

  k=*r0;
  *r0=k/2.0;
}
#endif
/*----------------------------------------------------------------------------------*/
/*------calculate inside the catchment--------------------*/
/*----------------------------------------------------------------------------------*/
double calculate_inside_catchment(c,P,vectors,voltage,new_z)
     catchment *c;
     coordinates P;
     int *new_z;
     bem_vectors *vectors;
     bem_results *voltage;
{
  int this_zone,previous_zone;
  double pp;
  boundary *bb;

  this_zone=check_each_zone(c,P);
  previous_zone=c->previous_zone;
  if(this_zone<0)           /* outside catchment */
    {
      (*new_z)=(-1);
      pp=0.0;
      voltage->V=0.0;
      voltage->dV[0]=0.0;
      voltage->dV[1]=0.0;
      voltage->d2V[0][0]=0.0;      voltage->d2V[0][1]=0.0;
      voltage->d2V[1][0]=0.0;      voltage->d2V[1][1]=0.0;
    }
  else                 /* inside catchment */
    {
      bb=c->zones[this_zone];
      if(this_zone==previous_zone) /* same zone */
	{
	  (*new_z)=0;
	  pp=calculate_in_same_zone(bb,P,vectors,voltage);
	}
      else		   /* new zone */
	{
	  c->previous_zone=this_zone;
	  (*new_z)=1;
	  pp=calculate_in_new_zone(bb,P,vectors,voltage);
	}
    }
  return(pp);
}
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------*/
/*-------------------------------- memory.c  ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on memory structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

#include "co_matrix.h"
#include "matrix.h"
#include "ten_matrix.h"

#include "memory.h"
/*----------------------------------------------------------------------------------*/
/* create a bem vector */
/*----------------------------------------------------------------------------------*/
bem_vectors *create_bem_vectors(bvv,bcv,N)
     matrix *bvv, *bcv;
     int N;
{
  bem_vectors *x;

  x=(bem_vectors *)malloc(sizeof(bem_vectors));
  if(x==(bem_vectors *)NULL)
    {
      printf("error allocating memory for bem_vectors\n");
      exit(0);
    }
  x->bvv=bvv;                          /* boundary voltage vector */               
  x->bcv=bcv;			       /* boundary current vector */               
  x->vgv=create_matrix(1,2*N);	       /* voltage geometry vector */               
  x->cgv=create_matrix(1,4*N);	       /* current geometry vector */               
  x->co_vgv=create_co_matrix(1,2*N);   /* voltage geometry vector 1st derivative */
  x->co_cgv=create_co_matrix(1,4*N);   /* current geometry vector 1st derivative */
  x->ten_vgv=create_ten_matrix(1,2*N); /* voltage geometry vector 2nd derivative */
  x->ten_cgv=create_ten_matrix(1,4*N); /* current geometry vector 2nd derivative */
  return(x);
}

/*----------------------------------------------------------------------------------*/
/* destroy bem_vectors */
/*----------------------------------------------------------------------------------*/
bem_vectors *destroy_bem_vectors(x)
     bem_vectors *x;
{
  if(x!=(bem_vectors *)NULL) 
    {
      if(x->vgv!=(matrix *)NULL) destroy_matrix(x->vgv);
      if(x->cgv!=(matrix *)NULL) destroy_matrix(x->cgv);
      if(x->co_vgv!=(co_matrix *)NULL) destroy_co_matrix(x->co_vgv);
      if(x->co_cgv!=(co_matrix *)NULL) destroy_co_matrix(x->co_cgv);
      if(x->ten_vgv!=(ten_matrix *)NULL) destroy_ten_matrix(x->ten_vgv);
      if(x->ten_cgv!=(ten_matrix *)NULL) destroy_ten_matrix(x->ten_cgv);
    }
  free((void *)x);
  return((bem_vectors *)NULL);
}

/*----------------------------------------------------------------------------------*/

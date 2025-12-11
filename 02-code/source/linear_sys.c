/*------------------------------ linear_sys.c ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on path structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"

#include "co_matrix.h"
#include "geometry.h"
#include "matrix.h"
#include "path.h"
#include "ten_matrix.h"
#include "terms.h"

#include "linear_sys.h"
/*----------------------------------------------------------------------------------*/
/*  make voltage geometry matrix */
/*----------------------------------------------------------------------------------*/
void make_voltage_geometry_matrix(b,vgm)
     boundary *b;
     matrix *vgm;
{
  int i,j,paths;
  path *path_i, *path_j;
  int offset_i, offset_j;

  paths = b->components;

  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      offset_i=0;
      for(i=0;i<paths;i++)
	{ 
	  path_i=b->loop[i];
	  fill_voltage_geometry_matrix(offset_i,offset_j,path_i,path_j,vgm);
	  offset_i=offset_i+b->loop[i]->points;
	}
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill voltage geometry matrix */
/*----------------------------------------------------------------------------------*/
void fill_voltage_geometry_matrix(offset_i,offset_j,path_i,path_j,vgm)
     
  matrix *vgm;   
  int offset_i, offset_j;
  path *path_i, *path_j;

{
  int i,j;
  int segment_i, segment_j, points_i, points_j;
  coordinates Qa, Qb, Pa, Pb, Pc, Pd, Pe, Pf;
  double x, y1, y2;
  double V,W;

  offset_i=offset_i*5;
  offset_j=offset_j*2;
  points_i=path_i->points;
  points_j=path_j->points;
  i=0;
  j=0;
 
  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);
      j=segment_j*2;
      for(segment_i=0;segment_i<points_i;segment_i++)
	{  
	  i=segment_i*5;
	  get_path_xy(path_i,segment_i,Pa);
	  get_path_xy(path_i,segment_i+1,Pf);
	  Pb[0]=0.8*Pa[0]+0.2*Pf[0];    Pb[1]=0.8*Pa[1]+0.2*Pf[1];
	  Pc[0]=0.6*Pa[0]+0.4*Pf[0];    Pc[1]=0.6*Pa[1]+0.4*Pf[1];
	  Pd[0]=0.4*Pa[0]+0.6*Pf[0];    Pd[1]=0.4*Pa[1]+0.6*Pf[1];
	  Pe[0]=0.2*Pa[0]+0.8*Pf[0];    Pe[1]=0.2*Pa[1]+0.8*Pf[1];
	  if(path_i!=path_j)  /* Pa, Pb, Pc Pd off segment S */
	    {
	      convert_PQ(Qa, Qb, Pa, &x, &y1, &y2);
	      V=Vterm_PoffS(x,y1,y2);
	      W=Wterm_PoffS(x,y1,y2);
	      p2c_2basis(V,W,&V,&W);
	      put_block_matrix_element(vgm,offset_i,offset_j,i,j,  V);
	      put_block_matrix_element(vgm,offset_i,offset_j,i,j+1,W);

	      convert_PQ(Qa, Qb, Pb, &x, &y1, &y2);
	      V=Vterm_PoffS(x,y1,y2);
	      W=Wterm_PoffS(x,y1,y2);
	      p2c_2basis(V,W,&V,&W);
	      put_block_matrix_element(vgm,offset_i,offset_j,i+1,j,  V);
	      put_block_matrix_element(vgm,offset_i,offset_j,i+1,j+1,W);

	      convert_PQ(Qa, Qb, Pc, &x, &y1, &y2);
	      V=Vterm_PoffS(x,y1,y2);
	      W=Wterm_PoffS(x,y1,y2);
	      p2c_2basis(V,W,&V,&W);
	      put_block_matrix_element(vgm,offset_i,offset_j,i+2,j,  V);
	      put_block_matrix_element(vgm,offset_i,offset_j,i+2,j+1,W);

	      convert_PQ(Qa, Qb, Pd, &x, &y1, &y2);
	      V=Vterm_PoffS(x,y1,y2);
	      W=Wterm_PoffS(x,y1,y2);
	      p2c_2basis(V,W,&V,&W);
	      put_block_matrix_element(vgm,offset_i,offset_j,i+3,j,  V);
	      put_block_matrix_element(vgm,offset_i,offset_j,i+3,j+1,W);

	      convert_PQ(Qa, Qb, Pe, &x, &y1, &y2);
	      V=Vterm_PoffS(x,y1,y2);
	      W=Wterm_PoffS(x,y1,y2);
	      p2c_2basis(V,W,&V,&W);
	      put_block_matrix_element(vgm,offset_i,offset_j,i+4,j,  V);
	      put_block_matrix_element(vgm,offset_i,offset_j,i+4,j+1,W);
	    }
	  else
	    {
	      switch((segment_i-segment_j+points_i)%points_i)
		{
		case 0: /* Pa, Pb, Pc Pd on segment S */
		  convert_PQ(Qa, Qb, Pa, &x, &y1, &y2);
		  V=Vterm_PonS(x,y1,y2);
		  W=Wterm_PonS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i,j+1,W);

		  convert_PQ(Qa, Qb, Pb, &x, &y1, &y2);
		  V=Vterm_PonS(x,y1,y2);
		  W=Wterm_PonS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+1,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+1,j+1,W);
	      
		  convert_PQ(Qa, Qb, Pc, &x, &y1, &y2);
		  V=Vterm_PonS(x,y1,y2);
		  W=Wterm_PonS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+2,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+2,j+1,W);

		  convert_PQ(Qa, Qb, Pd, &x, &y1, &y2);
		  V=Vterm_PonS(x,y1,y2);
		  W=Wterm_PonS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+3,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+3,j+1,W);

		  convert_PQ(Qa, Qb, Pe, &x, &y1, &y2);
		  V=Vterm_PonS(x,y1,y2);
		  W=Wterm_PonS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+4,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+4,j+1,W);
		  break;

		case 1: /* Pa on segment S, Pb, Pc Pd off segment S */
		  convert_PQ(Qa, Qb, Pa, &x, &y1, &y2);
		  V=Vterm_PonS(x,y1,y2);
		  W=Wterm_PonS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i,j+1,W);

		  convert_PQ(Qa, Qb, Pb, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+1,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+1,j+1,W);

		  convert_PQ(Qa, Qb, Pc, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+2,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+2,j+1,W);

		  convert_PQ(Qa, Qb, Pd, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+3,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+3,j+1,W);

		  convert_PQ(Qa, Qb, Pe, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+4,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+4,j+1,W);
		  break;

		default: /* Pa, Pb, Pc Pd off segment S */
		  convert_PQ(Qa, Qb, Pa, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i,j+1,W);

		  convert_PQ(Qa, Qb, Pb, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+1,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+1,j+1,W);

		  convert_PQ(Qa, Qb, Pc, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+2,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+2,j+1,W);

		  convert_PQ(Qa, Qb, Pd, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+3,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+3,j+1,W);

		  convert_PQ(Qa, Qb, Pe, &x, &y1, &y2);
		  V=Vterm_PoffS(x,y1,y2);
		  W=Wterm_PoffS(x,y1,y2);
		  p2c_2basis(V,W,&V,&W);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+4,j,  V);
		  put_block_matrix_element(vgm,offset_i,offset_j,i+4,j+1,W);
		  break;
		}
	    }
	}
    } 
}

/*----------------------------------------------------------------------------------*/
/*  make current geometry matrix */
/*----------------------------------------------------------------------------------*/
void make_current_geometry_matrix(b,cgm)
     boundary *b;
     matrix *cgm;
{
  int i,j,paths;
  path *path_i, *path_j;
  int offset_i, offset_j;

  paths = b->components;

  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      offset_i=0;
      for(i=0;i<paths;i++)
	{ 
	  path_i=b->loop[i];
	  fill_current_geometry_matrix(offset_i,offset_j,path_i,path_j,cgm);
	  offset_i=offset_i+b->loop[i]->points;
	}
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill current geometry matrix */
/*----------------------------------------------------------------------------------*/
void fill_current_geometry_matrix(offset_i,offset_j,path_i,path_j,cgm)
     
  matrix *cgm;   
  int offset_i, offset_j;
  path *path_i, *path_j;

{
  int i,j;
  int segment_i, segment_j, points_i, points_j;
  coordinates Qa, Qb, Pa, Pb, Pc, Pd, Pe, Pf;
  double x, y1, y2;
  double J,K,L,M;

  offset_i=offset_i*5;
  offset_j=offset_j*4;
  points_i=path_i->points;
  points_j=path_j->points;
  i=0;
  j=0;

  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);
      j=segment_j*4;
      for(segment_i=0;segment_i<points_i;segment_i++)
	{  
	  i=segment_i*5;
	  get_path_xy(path_i,segment_i,Pa);
	  get_path_xy(path_i,segment_i+1,Pf);
	  Pb[0]=0.8*Pa[0]+0.2*Pf[0];    Pb[1]=0.8*Pa[1]+0.2*Pf[1];
	  Pc[0]=0.6*Pa[0]+0.4*Pf[0];    Pc[1]=0.6*Pa[1]+0.4*Pf[1];
	  Pd[0]=0.4*Pa[0]+0.6*Pf[0];    Pd[1]=0.4*Pa[1]+0.6*Pf[1];
	  Pe[0]=0.2*Pa[0]+0.8*Pf[0];    Pe[1]=0.2*Pa[1]+0.8*Pf[1];
	 
	  if(path_i!=path_j)  /* Pa, Pb, Pc Pd off segment S */
	    {
	      convert_PQ(Qa, Qb, Pa, &x, &y1, &y2); 
	      J=Jterm_PoffS(x,y1,y2);
	      K=Kterm_PoffS(x,y1,y2);
	      L=Lterm_PoffS(x,y1,y2);
	      M=Mterm_PoffS(x,y1,y2);
	      p2c_4basis(J,K,L,M,&J,&K,&L,&M);
	      put_block_matrix_element(cgm,offset_i,offset_j,i,j,  J);
	      put_block_matrix_element(cgm,offset_i,offset_j,i,j+1,K);
	      put_block_matrix_element(cgm,offset_i,offset_j,i,j+2,L);
	      put_block_matrix_element(cgm,offset_i,offset_j,i,j+3,M);

	      convert_PQ(Qa, Qb, Pb, &x, &y1, &y2); 
	      J=Jterm_PoffS(x,y1,y2); 
	      K=Kterm_PoffS(x,y1,y2);
	      L=Lterm_PoffS(x,y1,y2);
	      M=Mterm_PoffS(x,y1,y2);
	      p2c_4basis(J,K,L,M,&J,&K,&L,&M);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+1,j,  J);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+1,K);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+2,L);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+3,M);

	      convert_PQ(Qa, Qb, Pc, &x, &y1, &y2);  
	      J=Jterm_PoffS(x,y1,y2);
	      K=Kterm_PoffS(x,y1,y2);
	      L=Lterm_PoffS(x,y1,y2);
	      M=Mterm_PoffS(x,y1,y2);
	      p2c_4basis(J,K,L,M,&J,&K,&L,&M);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+2,j,  J);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+1,K);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+2,L);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+3,M);

	      convert_PQ(Qa, Qb, Pd, &x, &y1, &y2);  
	      J=Jterm_PoffS(x,y1,y2);
	      K=Kterm_PoffS(x,y1,y2);
	      L=Lterm_PoffS(x,y1,y2);
	      M=Mterm_PoffS(x,y1,y2);
	      p2c_4basis(J,K,L,M,&J,&K,&L,&M);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+3,j,  J);
 	      put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+1,K);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+2,L);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+3,M);

	      convert_PQ(Qa, Qb, Pe, &x, &y1, &y2);  
	      J=Jterm_PoffS(x,y1,y2);
	      K=Kterm_PoffS(x,y1,y2);
	      L=Lterm_PoffS(x,y1,y2);
	      M=Mterm_PoffS(x,y1,y2);
	      p2c_4basis(J,K,L,M,&J,&K,&L,&M);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+4,j,  J);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+1,K);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+2,L);
	      put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+3,M);
	    }
	  else
	    {
	      switch((segment_i-segment_j+points_i)%points_i)
		{
		case 0: /* Pa, Pb, Pc Pd on segment S */
		  convert_PQ(Qa, Qb, Pa, &x, &y1, &y2); 
		  J=Jterm_PonS(x,y1,y2);
		  K=Kterm_PonS(x,y1,y2);
		  L=Lterm_PonS(x,y1,y2);
		  M=Mterm_PonS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+3,M);

		  convert_PQ(Qa, Qb, Pb, &x, &y1, &y2);  
		  J=Jterm_PonS(x,y1,y2);
		  K=Kterm_PonS(x,y1,y2);
		  L=Lterm_PonS(x,y1,y2);
		  M=Mterm_PonS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+3,M);
	      
		  convert_PQ(Qa, Qb, Pc, &x, &y1, &y2);  
		  J=Jterm_PonS(x,y1,y2);
		  K=Kterm_PonS(x,y1,y2);
		  L=Lterm_PonS(x,y1,y2);
		  M=Mterm_PonS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+3,M);

		  convert_PQ(Qa, Qb, Pd, &x, &y1, &y2);  
		  J=Jterm_PonS(x,y1,y2);
		  K=Kterm_PonS(x,y1,y2);
		  L=Lterm_PonS(x,y1,y2);
		  M=Mterm_PonS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+3,M);

		  convert_PQ(Qa, Qb, Pe, &x, &y1, &y2);  
		  J=Jterm_PonS(x,y1,y2);
		  K=Kterm_PonS(x,y1,y2);
		  L=Lterm_PonS(x,y1,y2);
		  M=Mterm_PonS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+3,M);
		  break;

		case 1: /* Pa on segment S, Pb, Pc Pd off segment S */
		  convert_PQ(Qa, Qb, Pa, &x, &y1, &y2); 
		  J=Jterm_PonS(x,y1,y2);
		  K=Kterm_PonS(x,y1,y2);
		  L=Lterm_PonS(x,y1,y2);
		  M=Mterm_PonS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+3,M);

		  convert_PQ(Qa, Qb, Pb, &x, &y1, &y2);  
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+3,M);

		  convert_PQ(Qa, Qb, Pc, &x, &y1, &y2);  
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+3,M);

		  convert_PQ(Qa, Qb, Pd, &x, &y1, &y2);  
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+3,M);

		  convert_PQ(Qa, Qb, Pe, &x, &y1, &y2);  
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+3,M);
		  break;

		default: /* Pa, Pb, Pc Pd off segment S */
		  convert_PQ(Qa, Qb, Pa, &x, &y1, &y2); 
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i,j+3,M);

		  convert_PQ(Qa, Qb, Pb, &x, &y1, &y2);  
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+1,j+3,M);

		  convert_PQ(Qa, Qb, Pc, &x, &y1, &y2);  
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+2,j+3,M);

		  convert_PQ(Qa, Qb, Pd, &x, &y1, &y2);  
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+3,j+3,M);

		  convert_PQ(Qa, Qb, Pe, &x, &y1, &y2);  
		  J=Jterm_PoffS(x,y1,y2);
		  K=Kterm_PoffS(x,y1,y2);
		  L=Lterm_PoffS(x,y1,y2);
		  M=Mterm_PoffS(x,y1,y2);
		  p2c_4basis(J,K,L,M,&J,&K,&L,&M);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j,  J);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+1,K);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+2,L);
		  put_block_matrix_element(cgm,offset_i,offset_j,i+4,j+3,M);
		  break;
		}
	    }
	}
    } 
}

/*----------------------------------------------------------------------------------*/
/*  make diagonal matrix */
/*----------------------------------------------------------------------------------*/
void make_diagonal_matrix(b,dm)
     boundary *b;
     matrix *dm;
{
  int i,j,paths;
  path *path_i, *path_j;
  int offset_i, offset_j;

  paths = b->components;

  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      offset_i=0;
      for(i=0;i<j;i++)
	{ 
	  path_i=b->loop[i];
	  empty_diagonal_matrix(offset_i,offset_j,path_i,path_j,dm);
	  offset_i=offset_i+b->loop[i]->points;
	}
      i=j;
      path_i=b->loop[i];
      fill_diagonal_matrix(offset_i,offset_j,path_i,path_j,dm);
      offset_i=offset_i+b->loop[i]->points;
      for(i=j+1;i<paths;i++)
	{ 
	  path_i=b->loop[i];
	  empty_diagonal_matrix(offset_i,offset_j,path_i,path_j,dm);
	  offset_i=offset_i+b->loop[i]->points;
	}
      offset_j=offset_j+b->loop[j]->points;
    }
}

void fill_diagonal_matrix(offset_i,offset_j,path_i,path_j,dm)
     path *path_i, *path_j;
     int offset_i,offset_j;
     matrix *dm;
{
  int points_i ,points_j, segment_i ,segment_j ,i ,j;
  double syn_x ,syn_y ,syn;
  double V,W;
  coordinates a, b, c;

  points_i=path_i->points;
  points_j=path_j->points;
  offset_i=offset_i*5;
  offset_j=offset_j*2;

  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      j=segment_j*2;
      for(segment_i=0;segment_i<segment_j;segment_i++)
	{
	  i=segment_i*5;
	  put_block_matrix_element(dm,offset_i,offset_j,i,  j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i,  j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+1,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+1,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+2,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+2,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+3,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+3,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+4,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+4,j+1,0.0);
	}

      segment_i=segment_j;
      i=segment_i*5;
      /*------------ segment_i-1=-1 or last points------------------ */
      
      get_path_xy(path_i,segment_i-1+points_i, a);
      get_path_xy(path_i,segment_i, b);
      get_path_xy(path_i,segment_i+1,c);
      
      /*------------ calculate syn-------------------------------- */
      
      syn_y = ((a[1]-b[1])*(c[0]-b[0]) - (c[1]-b[1])*(a[0]-b[0]));
      syn_x = ((a[0]-b[0])*(c[0]-b[0]) + (c[1]-b[1])*(a[1]-b[1]));
      syn = atan2(syn_y,syn_x);
      syn=syn/(2.0*M_PI);
      if(syn<0.0) { syn = syn + 1.0; }
      
      V=(-0.5)*syn;
      W=syn;
      p2c_2basis(V,W,&V,&W);
      put_block_matrix_element(dm,offset_i,offset_j,i,j,  V);
      put_block_matrix_element(dm,offset_i,offset_j,i,j+1,W);
      V=(-0.3)*0.5;
      W=0.5;
      p2c_2basis(V,W,&V,&W);
      put_block_matrix_element(dm,offset_i,offset_j,i+1,j,  V); 
      put_block_matrix_element(dm,offset_i,offset_j,i+1,j+1,W);
      V=(-0.1)*0.5;
      W=0.5;
      p2c_2basis(V,W,&V,&W);
      put_block_matrix_element(dm,offset_i,offset_j,i+2,j,  V);
      put_block_matrix_element(dm,offset_i,offset_j,i+2,j+1,W);
      V=0.1*0.5;
      W=0.5;
      p2c_2basis(V,W,&V,&W);
      put_block_matrix_element(dm,offset_i,offset_j,i+3,j,  V);
      put_block_matrix_element(dm,offset_i,offset_j,i+3,j+1,W);
      V=0.3*0.5;
      W=0.5;
      p2c_2basis(V,W,&V,&W);
      put_block_matrix_element(dm,offset_i,offset_j,i+4,j,  V);
      put_block_matrix_element(dm,offset_i,offset_j,i+4,j+1,W);

      for(segment_i=segment_j+1;segment_i<points_i;segment_i++)
	{
	  i=segment_i*5;
	  put_block_matrix_element(dm,offset_i,offset_j,i,  j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i,  j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+1,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+1,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+2,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+2,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+3,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+3,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+4,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+4,j+1,0.0);
	}
    }
}

void empty_diagonal_matrix(offset_i,offset_j,path_i,path_j,dm)
     path *path_i, *path_j;
     int offset_i,offset_j;
     matrix *dm;
{
  int points_i ,points_j, segment_i ,segment_j ,i ,j;

  points_i=path_i->points;
  points_j=path_j->points;
  offset_i=offset_i*5;
  offset_j=offset_j*2;

  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      j=segment_j*2;
      for(segment_i=0;segment_i<points_i;segment_i++)
	{
	  i=segment_i*5;
	  put_block_matrix_element(dm,offset_i,offset_j,i,  j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i,  j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+1,j,  0.0); 
	  put_block_matrix_element(dm,offset_i,offset_j,i+1,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+2,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+2,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+3,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+3,j+1,0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+4,j,  0.0);
	  put_block_matrix_element(dm,offset_i,offset_j,i+4,j+1,0.0);
	}
    }
}

/*----------------------------------------------------------------------------------*/
/*  make voltage geometry vector */
/*----------------------------------------------------------------------------------*/
void make_voltage_geometry_vector(P,b,vgv)
     coordinates P;
     boundary *b;
     matrix *vgv;
{
  int j,paths;
  path *path_j;
  int offset_j;

  paths = b->components;

  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      fill_voltage_geometry_vector(P,offset_j,path_j,vgv);
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill voltage geometry vector */
/*----------------------------------------------------------------------------------*/
void fill_voltage_geometry_vector(P,offset_j,path_j,vgv)
     coordinates P;
     matrix *vgv;   
     int offset_j;
     path *path_j;

{
  int j;
  int segment_j, points_j;
  coordinates Qa, Qb;
  double x, y1, y2;
  double V,W;

  offset_j=offset_j*2;
  points_j=path_j->points;
  j=0;
 
  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);
      j=segment_j*2;

      convert_PQ(Qa, Qb, P, &x, &y1, &y2);
      V=Vterm_PoffS(x,y1,y2);
      W=Wterm_PoffS(x,y1,y2);
      p2c_2basis(V,W,&V,&W);
      put_block_matrix_element(vgv,0,offset_j,0,j,  V);
      put_block_matrix_element(vgv,0,offset_j,0,j+1,W);
    } 
}

/*----------------------------------------------------------------------------------*/
/*  make current geometry vector */
/*----------------------------------------------------------------------------------*/
void make_current_geometry_vector(P,b,cgv)

     coordinates P;
     boundary *b;
     matrix *cgv;
{
  int j,paths;
  path *path_j;
  int offset_j;

  paths = b->components;
  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      fill_current_geometry_vector(P, offset_j, path_j, cgv);
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill current geometry vector */
/*----------------------------------------------------------------------------------*/
void fill_current_geometry_vector(P, offset_j, path_j, cgv)
     
     coordinates P;
     matrix *cgv;   
     int  offset_j;
     path *path_j;
{
  int j;
  int segment_j, points_j;
  coordinates Qa, Qb;
  double x, y1, y2;
  double J,K,L,M;

  offset_j=offset_j*4;
  points_j=path_j->points;
  j=0;

  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);

      j=segment_j*4;
	    
      convert_PQ(Qa, Qb, P, &x, &y1, &y2); 
      J=Jterm_PoffS(x,y1,y2);
      K=Kterm_PoffS(x,y1,y2);
      L=Lterm_PoffS(x,y1,y2);
      M=Mterm_PoffS(x,y1,y2);
      p2c_4basis(J,K,L,M,&J,&K,&L,&M);
      put_block_matrix_element(cgv,0,offset_j,0,j,  J);
      put_block_matrix_element(cgv,0,offset_j,0,j+1,K);
      put_block_matrix_element(cgv,0,offset_j,0,j+2,L);
      put_block_matrix_element(cgv,0,offset_j,0,j+3,M);
    }
} 

/*----------------------------------------------------------------------------------*/
/*  make Kirchoff's current law geometry vector */
/*----------------------------------------------------------------------------------*/
void make_kcl_geometry_vector(b,kcl)

     boundary *b;
     matrix *kcl;
{
  int j,paths;
  path *path_j;
  int offset_j;

  paths = b->components;
  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      fill_kcl_geometry_vector(offset_j, path_j, kcl);
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill Kirchoff's current law geometry vector */
/*----------------------------------------------------------------------------------*/
void fill_kcl_geometry_vector(offset_j, path_j, kcl)
     
     matrix *kcl;   
     int  offset_j;
     path *path_j;
{
  int j;
  int segment_j, points_j;
  coordinates P, Qa, Qb;
  double x, y1, y2;
  double J,K,L,M;

  offset_j=offset_j*4;
  points_j=path_j->points;
  j=0;

  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);
      P[0]=(Qa[0]+Qb[0])/2.0;
      P[1]=(Qa[1]+Qb[1])/2.0;

      j=segment_j*4;
	    
      convert_PQ(Qa, Qb, P, &x, &y1, &y2); 
      J=0.0;              /* J term */
      K=(y2-y1)/12.0;     /* K term */
      L=0.0;              /* L term */
      M=(y2-y1);          /* M term */
      p2c_4basis(J,K,L,M,&J,&K,&L,&M);
      put_block_matrix_element(kcl,0,offset_j,0,j,  J);
      put_block_matrix_element(kcl,0,offset_j,0,j+1,K);
      put_block_matrix_element(kcl,0,offset_j,0,j+2,L);
      put_block_matrix_element(kcl,0,offset_j,0,j+3,M);
    }
} 

/*----------------------------------------------------------------------------------*/
/*-------------------- these ones work with coordinates (not scalars) --------------*/
/*----------------------------------------------------------------------------------*/
/*  make coordinates voltage geometry vector */
/*----------------------------------------------------------------------------------*/
void make_co_voltage_geometry_vector(P,b,co_vgv)
     coordinates P;
     boundary *b;
     co_matrix *co_vgv;
{
  int j,paths;
  path *path_j;
  int offset_j;

  paths = b->components;

  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      fill_co_voltage_geometry_vector(P,offset_j,path_j,co_vgv);
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill coordinates voltage geometry vector */
/*----------------------------------------------------------------------------------*/
void fill_co_voltage_geometry_vector(P,offset_j,path_j,co_vgv)
     coordinates P;
     co_matrix *co_vgv;   
     int offset_j;
     path *path_j;
{
  int j;
  int segment_j, points_j;
  coordinates Qa, Qb;
  double x, y1, y2;
  coordinates Vterm,Wterm,A0,A1;

  offset_j=offset_j*2;
  points_j=path_j->points;
  j=0;
 
  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);
      j=segment_j*2;

      convert_PQ(Qa, Qb, P, &x, &y1, &y2);

      V1(x,y1,y2,Vterm);
      W1(x,y1,y2,Wterm);
      rotate_to_PQ(Vterm[0],Vterm[1],Qa,Qb,Vterm);
      rotate_to_PQ(Wterm[0],Wterm[1],Qa,Qb,Wterm);
      p2c_2basis_co(Vterm,Wterm,A0,A1);
      put_block_co_matrix_element(co_vgv,0,offset_j,0,j,  A0);
      put_block_co_matrix_element(co_vgv,0,offset_j,0,j+1,A1);
    } 
}

/*----------------------------------------------------------------------------------*/
/*  make coordinates current geometry vector */
/*----------------------------------------------------------------------------------*/
void make_co_current_geometry_vector(P,b,co_cgv)

     coordinates P;
     boundary *b;
     co_matrix *co_cgv;
{
  int j,paths;
  path *path_j;
  int offset_j;

  paths = b->components;
  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      fill_co_current_geometry_vector(P, offset_j, path_j, co_cgv);
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill coordinates current geometry vector */
/*----------------------------------------------------------------------------------*/
void fill_co_current_geometry_vector(P, offset_j, path_j, co_cgv)
     
     coordinates P;
     co_matrix *co_cgv;   
     int  offset_j;
     path *path_j;
{
  int j;
  int segment_j, points_j;
  coordinates Qa, Qb;
  double x, y1, y2;
  coordinates Jterm,Kterm,Lterm,Mterm,A0,A1,A2,A3;

  offset_j=offset_j*4;
  points_j=path_j->points;
  j=0;

  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);

      j=segment_j*4;
	    
      convert_PQ(Qa, Qb, P, &x, &y1, &y2); 
      J1(x,y1,y2,Jterm);
      K1(x,y1,y2,Kterm);
      L1(x,y1,y2,Lterm);
      M1(x,y1,y2,Mterm);
      rotate_to_PQ(Jterm[0],Jterm[1],Qa,Qb,Jterm);
      rotate_to_PQ(Kterm[0],Kterm[1],Qa,Qb,Kterm);
      rotate_to_PQ(Lterm[0],Lterm[1],Qa,Qb,Lterm);
      rotate_to_PQ(Mterm[0],Mterm[1],Qa,Qb,Mterm);
      p2c_4basis_co(Jterm,Kterm,Lterm,Mterm,A0,A1,A2,A3);
      put_block_co_matrix_element(co_cgv,0,offset_j,0,j,  A0);
      put_block_co_matrix_element(co_cgv,0,offset_j,0,j+1,A1);
      put_block_co_matrix_element(co_cgv,0,offset_j,0,j+2,A2);
      put_block_co_matrix_element(co_cgv,0,offset_j,0,j+3,A3);
    }
} 

/*----------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------*/
/*-------------------- these ones work with tensor----------------------------------*/
/*----------------------------------------------------------------------------------*/
/*  make tensor voltage geometry vector */
/*----------------------------------------------------------------------------------*/
void make_ten_voltage_geometry_vector(P,b,ten_vgv)
     coordinates P;
     boundary *b;
     ten_matrix *ten_vgv;
{
  int j,paths;
  path *path_j;
  int offset_j;

  paths = b->components;

  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      fill_ten_voltage_geometry_vector(P,offset_j,path_j,ten_vgv);
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill tensor voltage geometry vector */
/*----------------------------------------------------------------------------------*/
void fill_ten_voltage_geometry_vector(P,offset_j,path_j,ten_vgv)
     coordinates P;
     ten_matrix *ten_vgv;   
     int offset_j;
     path *path_j;

{
  int j;
  int segment_j, points_j;
  coordinates Qa, Qb;
  double x, y1, y2;
  tensor Vterm,Wterm,A0,A1;

  offset_j=offset_j*2;
  points_j=path_j->points;
  j=0;
 
  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);
      j=segment_j*2;

      convert_PQ(Qa, Qb, P, &x, &y1, &y2);
      V2(x,y1,y2,Vterm);
      W2(x,y1,y2,Wterm);
      double_rotate_to_PQ(Vterm[0][0],Vterm[0][1],Vterm[1][0],Vterm[1][1],Qa,Qb,Vterm);
      double_rotate_to_PQ(Wterm[0][0],Wterm[0][1],Wterm[1][0],Wterm[1][1],Qa,Qb,Wterm);
      p2c_2basis_ten(Vterm,Wterm,A0,A1);
      put_block_ten_matrix_element(ten_vgv,0,offset_j,0,j,  A0);
      put_block_ten_matrix_element(ten_vgv,0,offset_j,0,j+1,A1);
    } 
}

/*----------------------------------------------------------------------------------*/
/*  make tensor current geometry vector */
/*----------------------------------------------------------------------------------*/
void make_ten_current_geometry_vector(P,b,ten_cgv)

     coordinates P;
     boundary *b;
     ten_matrix *ten_cgv;
{
  int j,paths;
  path *path_j;
  int offset_j;

  paths = b->components;
  offset_j=0;  
  for(j=0;j<paths;j++)
    { 
      path_j=b->loop[j];
      fill_ten_current_geometry_vector(P, offset_j, path_j, ten_cgv);
      offset_j=offset_j+b->loop[j]->points;
    }
}

/*----------------------------------------------------------------------------------*/
/*  fill tensor current geometry vector */
/*----------------------------------------------------------------------------------*/
void fill_ten_current_geometry_vector(P, offset_j, path_j, ten_cgv)
     
     coordinates P;
     ten_matrix *ten_cgv;   
     int  offset_j;
     path *path_j;
{
  int j;
  int segment_j, points_j;
  coordinates Qa, Qb;
  double x, y1, y2;
  tensor Jterm,Kterm,Lterm,Mterm,A0,A1,A2,A3;

  offset_j=offset_j*4;
  points_j=path_j->points;
  j=0;

  for(segment_j=0;segment_j<points_j;segment_j++)
    {
      get_path_xy(path_j,segment_j,Qa);
      get_path_xy(path_j,segment_j+1,Qb);

      j=segment_j*4;
	    
      convert_PQ(Qa, Qb, P, &x, &y1, &y2); 

      J2(x,y1,y2,Jterm);
      K2(x,y1,y2,Kterm);
      L2(x,y1,y2,Lterm);
      M2(x,y1,y2,Mterm);
      double_rotate_to_PQ(Jterm[0][0],Jterm[0][1],Jterm[1][0],Jterm[1][1],Qa,Qb,Jterm);
      double_rotate_to_PQ(Kterm[0][0],Kterm[0][1],Kterm[1][0],Kterm[1][1],Qa,Qb,Kterm);
      double_rotate_to_PQ(Lterm[0][0],Lterm[0][1],Lterm[1][0],Lterm[1][1],Qa,Qb,Lterm);
      double_rotate_to_PQ(Mterm[0][0],Mterm[0][1],Mterm[1][0],Mterm[1][1],Qa,Qb,Mterm);
      p2c_4basis_ten(Jterm,Kterm,Lterm,Mterm,A0,A1,A2,A3);
      put_block_ten_matrix_element(ten_cgv,0,offset_j,0,j,  A0);
      put_block_ten_matrix_element(ten_cgv,0,offset_j,0,j+1,A1);
      put_block_ten_matrix_element(ten_cgv,0,offset_j,0,j+2,A2);
      put_block_ten_matrix_element(ten_cgv,0,offset_j,0,j+3,A3);
    }
} 
/*----------------------------------------------------------------------------------*/
/*---------- routines to convert from polynomial basis to Chebyshev basis ----------*/
/*----------------------------------------------------------------------------------*/
void p2c_2coeff(v,w,a0,a1)
     double v,w;
     double *a0,*a1;
{
  (*a0)=w;
  (*a1)=v/4.0;
}

/*----------------------------------------------------------------------------------*/
void p2c_4coeff(j,k,l,m,a0,a1,a2,a3)
     double j,k,l,m;
     double *a0,*a1,*a2,*a3;
{
  (*a0)=m+k/16.0;
  (*a1)=l/4.0+j/32.0;
  (*a2)=k/16.0;
  (*a3)=j/64.0;
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void p2c_2basis(Vterm,Wterm,A0,A1)
     double Vterm,Wterm;
     double *A0,*A1;
{
  (*A0)=Wterm;
  (*A1)=4.0*Vterm;
}

/*----------------------------------------------------------------------------------*/
void p2c_4basis(Jterm,Kterm,Lterm,Mterm,A0,A1,A2,A3)
     double Jterm,Kterm,Lterm,Mterm;
     double *A0,*A1,*A2,*A3;
{
  (*A0)=Mterm;
  (*A1)=4.0*Lterm;
  (*A2)=16.0*Kterm-Mterm;
  (*A3)=64.0*Jterm-8.0*Lterm;
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void p2c_2basis_co(Vterm,Wterm,A0,A1)
     coordinates Vterm,Wterm;
     coordinates A0,A1;
{
  A0[0]=Wterm[0];      A0[1]=Wterm[1];
  A1[0]=4.0*Vterm[0];  A1[1]=4.0*Vterm[1];
}

/*----------------------------------------------------------------------------------*/
void p2c_4basis_co(Jterm,Kterm,Lterm,Mterm,A0,A1,A2,A3)
     coordinates Jterm,Kterm,Lterm,Mterm;
     coordinates A0,A1,A2,A3;
{
  A0[0]=Mterm[0];                    A0[1]=Mterm[1];
  A1[0]=4.0*Lterm[0];                A1[1]=4.0*Lterm[1];
  A2[0]=16.0*Kterm[0]-Mterm[0];      A2[1]=16.0*Kterm[1]-Mterm[1];
  A3[0]=64.0*Jterm[0]-8.0*Lterm[0];  A3[1]=64.0*Jterm[1]-8.0*Lterm[1];
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void p2c_2basis_ten(Vterm,Wterm,A0,A1)
     tensor Vterm,Wterm;
     tensor A0,A1;
{
  A0[0][0]=Wterm[0][0];      A0[0][1]=Wterm[0][1];
  A0[1][0]=Wterm[1][0];      A0[1][1]=Wterm[1][1];
  A1[0][0]=4.0*Vterm[0][0];  A1[0][1]=4.0*Vterm[0][1];
  A1[1][0]=4.0*Vterm[1][0];  A1[1][1]=4.0*Vterm[1][1];
}

/*----------------------------------------------------------------------------------*/
void p2c_4basis_ten(Jterm,Kterm,Lterm,Mterm,A0,A1,A2,A3)
     tensor Jterm,Kterm,Lterm,Mterm;
     tensor A0,A1,A2,A3;
{
  A0[0][0]=Mterm[0][0];                       A0[0][1]=Mterm[0][1];
  A0[1][0]=Mterm[1][0];                       A0[1][1]=Mterm[1][1];
  A1[0][0]=4.0*Lterm[0][0];                   A1[0][1]=4.0*Lterm[0][1];
  A1[1][0]=4.0*Lterm[1][0];                   A1[1][1]=4.0*Lterm[1][1];
  A2[0][0]=16.0*Kterm[0][0]-Mterm[0][0];      A2[0][1]=16.0*Kterm[0][1]-Mterm[0][1];
  A2[1][0]=16.0*Kterm[1][0]-Mterm[1][0];      A2[1][1]=16.0*Kterm[1][1]-Mterm[1][1];
  A3[0][0]=64.0*Jterm[0][0]-8.0*Lterm[0][0];  A3[0][1]=64.0*Jterm[0][1]-8.0*Lterm[0][1];
  A3[1][0]=64.0*Jterm[1][0]-8.0*Lterm[1][0];  A3[1][1]=64.0*Jterm[1][1]-8.0*Lterm[1][1];
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

/*-------------------------------------------------*/
/*      routines to set up scanning parameters     */
/*-------------------------------------------------*/
#include <stdio.h>
#include <math.h>
/*-------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

#include "scan.h"
/*-------------------------------------------------*/
void put_raster(data,ras)
     char *data;
     raster *ras;
{
  int i1,j1,i2,j2;
  coordinates first,last;

  sscanf(data,"P(%d, %d) = (%lf, %lf) P(%d, %d) = (%lf, %lf)",
	 &i1,&j1,&first[0],&first[1],&i2,&j2,&last[0],&last[1]);
  ras->nx=i2-i1+1;
  ras->ny=j2-j1+1;
  ras->P1[0]=first[0];
  ras->P1[1]=first[1];
  ras->P2[0]=last[0];
  ras->P2[1]=last[1];
}

/*-------------------------------------------------*/
void show_raster(ras)
     raster *ras;
{
  printf("P(0, 0) = (%f, %f)    P(%d, %d) = (%f, %f)\n",
	 ras->P1[0],ras->P1[1],ras->nx-1,ras->ny-1,
	 ras->P2[0],ras->P2[1]);
}

/*-------------------------------------------------*/
double x_raster(ras,i)
     raster *ras;
     int i;
{
  double x;

  if(i==0) { 
    x=ras->P1[0];}
  else { 
    if(i==ras->nx-1) {
      x=ras->P2[0]; }
    else {
      x=(ras->P1[0]*(ras->nx-1-i)+ras->P2[0]*i)/(ras->nx-1); }}
  return(x);
}

/*-------------------------------------------------*/
double y_raster(ras,j)
     raster *ras;
     int j;
{
  double y;

  if(j==0) { 
    y=ras->P1[1];}
  else { 
    if(j==ras->ny-1) {
      y=ras->P2[1]; }
    else {
      y=(ras->P1[1]*(ras->ny-1-j)+ras->P2[1]*j)/(ras->ny-1); }}
  return(y);
}

/*-------------------------------------------------*/
void put_section(data,sec)
     char *data;
     section *sec;
{
  int i1,i2;
  coordinates first,last;

  sscanf(data,"P(%d) = (%lf, %lf) P(%d) = (%lf, %lf)",
	 &i1,&first[0],&first[1],&i2,&last[0],&last[1]);
  sec->n=i2-i1+1;
  sec->P1[0]=first[0];
  sec->P1[1]=first[1];
  sec->P2[0]=last[0];
  sec->P2[1]=last[1];
  first[0]=(last[0]-first[0]);
  first[1]=(last[1]-first[1]);
  sec->step=sqrt(first[0]*first[0]+first[1]*first[1])/(i2-i1);
}
/*-------------------------------------------------*/
void put_sectionV2(Nseg,PA,PB,sec)
     int Nseg;
     coordinates PA, PB;
     section *sec;
{
  coordinates dP;
  
  sec->n=Nseg+1;
  sec->P1[0]=PA[0];
  sec->P1[1]=PA[1];
  sec->P2[0]=PB[0];
  sec->P2[1]=PB[1];
  dP[0]=(PB[0]-PA[0]);
  dP[1]=(PB[1]-PA[1]);
  sec->step=sqrt(dP[0]*dP[0]+dP[1]*dP[1])/(Nseg);
}

/*-------------------------------------------------*/
void show_section(sec)
     section *sec;
{
  printf("P(0) = (%f, %f) P(%d) = (%f, %f) ",
	 sec->P1[0],sec->P1[1],sec->n-1,sec->P2[0],sec->P2[1]);
  printf("dW = %f N_SWP = %d ",sec->step,sec->n);
}

/*-------------------------------------------------*/
void xy_section(sec,i,xy)
     section *sec;
     int i;
     coordinates xy;
{
  if(i==0) { 
    xy[0]=sec->P1[0];
    xy[1]=sec->P1[1];}
  else { 
    if(i==sec->n-1) {
      xy[0]=sec->P2[0]; 
      xy[1]=sec->P2[1]; }
    else {
      xy[0]=(sec->P1[0]*(sec->n-1-i)+sec->P2[0]*i)/(sec->n-1); 
      xy[1]=(sec->P1[1]*(sec->n-1-i)+sec->P2[1]*i)/(sec->n-1); }}
}

/*-------------------------------------------------*/
/*-------------------------------------------------*/
void put_interval(data,inter)
     char *data;
     interval *inter;
{
  int i1,i2;
  double first,last;

  sscanf(data,"t(%d) = (%lf) t(%d) = (%lf)",&i1,&first,&i2,&last);
  inter->nt=i2-i1+1;
  inter->t1=first;
  inter->t2=last;
}

/*-------------------------------------------------*/
void show_interval(inter)
     interval *inter;
{
  printf("t(0) = (%f)    t(%d) = (%f)\n",
	 inter->t1,inter->nt-1,inter->t2);
}

/*-------------------------------------------------*/
double t_interval(inter,i)
     interval *inter;
     int i;
{
  double t;

  if(i==0) { 
    t=inter->t1;}
  else { 
    if(i==inter->nt-1) {
      t=inter->t2; }
    else {
      t=(inter->t1*(inter->nt-1-i)+inter->t2*i)/(inter->nt-1); }}
  return(t);
}

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

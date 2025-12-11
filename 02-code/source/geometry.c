/*------------------------------ geometry.c -----------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all geometric operations */
/*----------------------------------------------------------------------------------*/
#include <math.h>
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"

#include "geometry.h"
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*  convert coordinate system from global to local */
/*----------------------------------------------------------------------------------*/
void convert_PQ(Qa, Qb, P, x, y1, y2)
     coordinates Qa, Qb, P;
     double *x, *y1, *y2;
{
  double yu, yv, xu, xv, d;

  yu = Qb[0] - Qa[0];    yv = Qb[1] - Qa[1];
  d = sqrt(yu*yu + yv*yv);
  yu = yu/d;
  yv = yv/d;
  xu = yv;
  xv = -yu;
  *y1 = (Qa[0] - P[0])*yu + (Qa[1] - P[1])*yv;
  *y2 = (Qb[0] - P[0])*yu + (Qb[1] - P[1])*yv;
  *x  = (Qa[0] - P[0])*xu + (Qa[1] - P[1])*xv;
}

/*----------------------------------------------------------------------------------*/
/*  rotate coordinate system from local to global */
/*----------------------------------------------------------------------------------*/
void rotate_to_PQ(x, y, Qa, Qb, R)
     coordinates Qa, Qb, R;
     double x, y;
{
  double yu, yv, xu, xv, d;

  yu = Qb[0] - Qa[0];    yv = Qb[1] - Qa[1];
  d = sqrt(yu*yu + yv*yv);
  yu = yu/d;
  yv = yv/d;
  xu = yv;
  xv = -yu;
  R[0]=x*xu+y*yu;
  R[1]=x*xv+y*yv;
}

/*----------------------------------------------------------------------------------*/
/*  rotate tensor system from local to global */
/*----------------------------------------------------------------------------------*/
void double_rotate_to_PQ(a, b, c, d, Qa, Qb, R)
     coordinates Qa, Qb;
     tensor R;
     double a,b,c,d;
{
  double yu, yv, dd , alphasq, alphabeta, betasq;

  yu = Qb[0] - Qa[0];    yv = Qb[1] - Qa[1];
  dd = sqrt(yu*yu + yv*yv);
  yu = yu/dd;
  yv = yv/dd;
  /* xu = yv , xv = -yu is rotated about 90 degree */
  alphasq   =  yv*yv;
  alphabeta = -yu*yv;
  betasq    =  yu*yu;
  /*
  printf("\n  V[0][0]=%f V[0][1]=%f",a,b);
  printf("\n  V[1][0]=%f V[1][1]=%f",c,d);
  printf("\n");
  */
  R[0][0] = a*alphasq - (b+c)*alphabeta + d*betasq;
  R[0][1] = b*alphasq + (a-d)*alphabeta - c*betasq;
  R[1][0] = c*alphasq + (a-d)*alphabeta - b*betasq;
  R[1][1] = d*alphasq + (b+c)*alphabeta + a*betasq;
}

/*----------------------------------------------------------------------------------*/
/*  calculate atan(y2/x)-atan(y1/x) */
/*----------------------------------------------------------------------------------*/
double atan3(y2, y1, x)
     double y2, y1, x;
{
  double a;
      a = atan2( x*(y2-y1) , (x*x)+(y1*y2) );
      return (a);
}

/*----------------------------------------------------------------------------------*/
/*  convert coordinate system  */
/*----------------------------------------------------------------------------------*/
double atanv(Q1, Q2, P)
     coordinates Q1, Q2, P;
{
  double a, x1, x2, y1, y2;
      
  x1=Q1[0]-P[0];
  y1=Q1[1]-P[1];
  x2=Q2[0]-P[0];
  y2=Q2[1]-P[1];

  a = atan2( (x1*y2-y1*x2) , (x1*x2+y1*y2) );
  return (a);
}
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/


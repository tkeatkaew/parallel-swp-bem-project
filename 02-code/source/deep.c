/*---------------------------------------------------*/
/* depth model */
/*---------------------------------------------------*/
#include <math.h>
/*---------------------------------------------------*/
#include "boundary_types.h"

#include "flow.h"
#include "rain.h"

#include "deep.h"
/*---------------------------------------------------*/
double depth(P,L,grad_h)
     coordinates P,grad_h;
     double L;
{
  double d,v;

  v=velocity(P,grad_h);
  if(v>0.0)
    {
      d=rainfall(P)*L/v;
    }
  else
    {
      d=0.0;
    }
  return(d);
}

/*---------------------------------------------------*/

/*---------------------------------------------------*/
/* flow model */
/*---------------------------------------------------*/
#include <math.h>
/*---------------------------------------------------*/
#include "boundary_types.h"

#include "impedance.h"

#include "flow.h"
/*---------------------------------------------------*/
/* density of water: rho*/
/*---------------------------------------------------*/
double density(void)
{
  double rho;

  rho=997.0; /* kg/m^3 at 25 degrees Celcius */
  return(rho);
}

/*---------------------------------------------------*/
/* velocity v */
/*---------------------------------------------------*/
double velocity(P,grad_h)
     coordinates P,grad_h;
{
  double gradh_sq,v;

  gradh_sq=(grad_h[0]*grad_h[0]+grad_h[1]*grad_h[1]);
  v=conductivity(P)*sqrt(gradh_sq/(1.0+gradh_sq));
  return(v);
}

/*---------------------------------------------------*/
/* current density Q */
/*---------------------------------------------------*/
double current_density(P,grad_h)
     coordinates P,grad_h;
{
  double Q;

  Q=density()*velocity(P,grad_h);
  return(Q);
}

/*---------------------------------------------------*/

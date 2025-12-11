#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*--------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

#include "catchment.h"
#include "file.h"
#include "memory.h"
#include "scan.h"
#include "trapfloat.h"
#include "vcalc.h"
#include "flow.h"

#include "ssection.h"
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main()
{
  char title[64];
  char data[]="c_valley.txt"; 

  FILE *output1, *output2;
  bem_results voltage;
  bem_vectors *vectors;
  catchment *c;
  char *buffer[2];
  coordinates P;
  double pp,v,Q;
  int buf_size,i,column,max_points,num_zones,new_z;
  matrix bvv, bcv;
  section mouth;

  trap_floating_errors();
  buf_size=512*13+1;
  buffer[0]=(char *)malloc(2*buf_size*sizeof(char));
  buffer[1]=buffer[0]+buf_size;

  num_zones=catchment_zones(data);
  c=create_catchment(num_zones,16);
  get_catchment(data,c);

  max_points=max_points_in_any_zone(c);
  printf("maximum points in any zone is %d\n",max_points);
  vectors=create_bem_vectors(&bvv,&bcv,max_points);

  put_section("P(0) = (4.0,5.0) P(100) = (5.0,4.0)",&mouth);
  show_section(&mouth);
  printf("step size across mouth is %f\n",mouth.step);

  output1=open_file(0,"vsection.out","w");
  output2=open_file(0,"qsection.out","w");
  column=0;
/*------------------------------------------------------*/
/* this loop scans across the mouth */
/*------------------------------------------------------*/
  for(i=0;i<mouth.n;i++)
    {
      xy_section(&mouth,i,P);    	  
      pp=calculate_inside_catchment(c,P,vectors,&voltage,&new_z);
      if(new_z==(-1)) printf("point is outside catchment\n");

      v=velocity(P,voltage.dV);
      Q=current_density(P,voltage.dV);

      put_buffer(buf_size,buffer[0],column,"%14.5e ",i*mouth.step);
      column=put_buffer(buf_size,buffer[1],column,"%14.5e ",i*mouth.step);
      put_buffer(buf_size,buffer[0],column,"%14.5e\n",v);
      column=put_buffer(buf_size,buffer[1],column,"%14.5e\n",Q);
    }
  put_next_line(output1,buffer[0]);
  put_next_line(output2,buffer[1]);

  sprintf(title,"cross-section from P1=(%f,%f) to P2=(%f,%f)",
	  mouth.P1[0],mouth.P1[1],mouth.P2[0],mouth.P2[1]);
  make_gpl2_file("vsection.out",title,"distance [m]","velocity [m/s]");

  sprintf(title,"cross-section from P1=(%f,%f) to P2=(%f,%f)",
	  mouth.P1[0],mouth.P1[1],mouth.P2[0],mouth.P2[1]);
  make_gpl2_file("qsection.out",title,"distance [m]","current density [(kg/s)/m^2]");
/*------------------------------------------------------*/

  printf("\n");
  fclose(output1);
  fclose(output2);
  destroy_catchment(c);
  destroy_bem_vectors(vectors);
  free((void *)buffer[0]);
  return(0);
}
/*--------------------------------------------------------*/

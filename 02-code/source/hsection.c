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

#include "hsection.h"
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main()
{
  char title[64];
  char data[]="c_valley.txt"; 

  FILE *output;
  bem_results voltage;
  bem_vectors *vectors;
  catchment *c;
  char *buffer;
  coordinates P;
  double pp;
  int buf_size,i,column,max_points,num_zones,new_z;
  matrix bvv, bcv;
  section mouth;

  trap_floating_errors();
  buf_size=512*13+1;
  buffer=(char *)malloc(buf_size*sizeof(char));

  num_zones=catchment_zones(data);
  c=create_catchment(num_zones,16);
  get_catchment(data,c);

  max_points=max_points_in_any_zone(c);
  printf("maximum points in any zone is %d\n",max_points);
  vectors=create_bem_vectors(&bvv,&bcv,max_points);

  put_section("P(0) = (4.0,5.0) P(100) = (5.0,4.0)",&mouth);
  show_section(&mouth);
  printf("step size across mouth is %f\n",mouth.step);

  output=open_file(0,"hsection.out","w");
  column=0;
/*------------------------------------------------------*/
/* this loop scans across the mouth */
/*------------------------------------------------------*/
  for(i=0;i<mouth.n;i++)
    {
      xy_section(&mouth,i,P);    	  
      pp=calculate_inside_catchment(c,P,vectors,&voltage,&new_z);
      if(new_z==(-1)) printf("point is outside catchment\n");
      
      column=put_buffer(buf_size,buffer,column,"%14.5e ",i*mouth.step);
      column=put_buffer(buf_size,buffer,column,"%14.5e\n",pp);
    }
  put_next_line(output,buffer);

  sprintf(title,"cross-section from P1=(%f,%f) to P2=(%f,%f)",
	  mouth.P1[0],mouth.P1[1],mouth.P2[0],mouth.P2[1]);
  make_gpl2_file("hsection.out",title,"distance [m]","height [m]");
/*------------------------------------------------------*/

  printf("\n");
  fclose(output);
  destroy_catchment(c);
  destroy_bem_vectors(vectors);
  free((void *)buffer);
  return(0);
}
/*--------------------------------------------------------*/

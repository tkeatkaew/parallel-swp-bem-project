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
#include "path.h"
#include "scan.h"
#include "streamline.h"
#include "trapfloat.h"
#include "mouthflow.h"

#include "flowrate.h"
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main()
{
  char title[64];
  char data[]="c_valley2.txt"; 

  bem_vectors *vectors;
  catchment *c;
  char *buffer;
  double step_size,I;
  int buf_size,i,max_points,num_zones,max_steps,max_streams;
  matrix bvv, bcv;
  path **streamlines;
  section mouth;

  trap_floating_errors();
  buf_size=512*12+1;
  buffer=(char *)malloc(buf_size*sizeof(char));

  num_zones=catchment_zones(data);
  c=create_catchment(num_zones,16);
  get_catchment(data,c);

  plot_catchment(c,"catchment.out");

  max_points=max_points_in_any_zone(c);
  printf("maximum points in any zone is %d\n",max_points);
  vectors=create_bem_vectors(&bvv,&bcv,max_points);

  max_steps=500;
  step_size=0.1;

  put_section("P(0) = (2.0,3.0) P(40) = (3.0,2.0)",&mouth);
  show_section(&mouth);
  printf("step size across mouth is %f\n",mouth.step);

  max_streams=mouth.n;
  streamlines=(path **)malloc(max_streams*sizeof(path *));
  for(i=0;i<max_streams;i++){
    streamlines[i]=create_path(max_steps,1,0);}

/*------------------------------------------------------*/
/* this loop scans a line of x,y */
/*------------------------------------------------------*/
  I=flow_rate(c,&mouth,1,max_steps,step_size,max_streams,
			streamlines,vectors);

  plot_streamlines(c,max_streams,streamlines,"flowrate.out");

  sprintf(title,"flow rate [%f kg/s]",I);
  make_gpl_file("flowrate.out",title,"[0:20]","[0:20]");
/*------------------------------------------------------*/

  printf("\nflow rate is %f\n",I);
  printf("destroy %d streamlines:\n",max_streams);
 
  for(i=0;i<max_streams;i++){
    streamlines[i]=destroy_path(streamlines[i]);}

  destroy_catchment(c);
  destroy_bem_vectors(vectors);
  free((void *)buffer);
  return(0);
}
/*--------------------------------------------------------*/

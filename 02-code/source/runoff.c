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

#include "runoff.h"
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main()
{
  /*  char data[]="c_valley.txt"; */

  char data[]="catchment1.txt";

  FILE *output;
  bem_results voltage;
  bem_vectors *vectors;
  catchment *c;
  coordinates P, Pstart;
  double L,step_size;
  int i,j,k,max_points,num_zones,max_steps;
  matrix bvv, bcv;
  path *streamup, *streamdown;
  raster ras;

  trap_floating_errors();

  num_zones=catchment_zones(data);
  c=create_catchment(num_zones,16);
  get_catchment(data,c);

  plot_catchment(c,"catchment.out");

  max_points=max_points_in_any_zone(c);
  printf("maximum points in any zone is %d\n",max_points);
  vectors=create_bem_vectors(&bvv,&bcv,max_points);

  max_steps=500;
  step_size=0.01;
  streamup=create_path(max_steps,1,0);
  streamdown=create_path(max_steps,1,0);
 
  output=open_file(0,"runoff.out","w");
  
  put_raster("P(0,0)=(0.2,0.0) P(5,5)=(0.4,0.0)",&ras);
  show_raster(&ras);
  k=0;
/*------------------------------------------------------*/
/* this loop scans a grid of x and y */
/*------------------------------------------------------*/
  for(j=0;j<ras.ny;j++)
   {
      P[1]=y_raster(&ras,j);
      for(i=0;i<ras.nx;i++)
	{
	  P[0]=x_raster(&ras,i);    	  
 	  Pstart[0]=P[0];
	  Pstart[1]=P[1];
 	  streamup->points=max_steps;
 	  L=streamline_loop(Pstart,c,1,max_steps,step_size,
 			    streamup,vectors,&voltage) ;
  	  Pstart[0]=P[0];
	  Pstart[1]=P[1];
 	  streamdown->points=max_steps;
 	  L=streamline_loop(Pstart,c,0,max_steps,step_size,
 			    streamdown,vectors,&voltage) ;
 	  plot_1_streamline(c,streamup,output);
	  /*
 	  plot_1_streamline(c,streamdown,output);
	  */
	}
      if(k%50==0) printf("\n");
      printf("#"); fflush(stdout);
      k=k+1;
   }
/*------------------------------------------------------*/

  printf("\n");
  fclose(output);
  destroy_path(streamup);
  destroy_path(streamdown);
  destroy_catchment(c);
  destroy_bem_vectors(vectors);
  return(0);
}
/*--------------------------------------------------------*/





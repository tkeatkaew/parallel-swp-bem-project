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
#include "streamline.h"
#include "flow.h"

#include "risk.h"
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main()
{
  char data[]="c_valley2.txt"; 

  FILE *output1, *output2, *output3, *output4;
  bem_results voltage;
  bem_vectors *vectors;
  catchment *c;
  char *buffer[3];
  coordinates P, Pstart;
  double L,step_size,v,risk;
  int buf_size,i,j,k,column,max_points,num_zones,max_steps;
  matrix bvv, bcv;
  raster ras;

  trap_floating_errors();
  buf_size=512*13+1;
  buffer[0]=(char *)malloc((2*buf_size+64)*sizeof(char));
  buffer[1]=buffer[0]+buf_size;
  buffer[2]=buffer[1]+buf_size;

  num_zones=catchment_zones(data);
  c=create_catchment(num_zones,16);
  get_catchment(data,c);

  plot_catchment(c,"catchment.out");

  max_points=max_points_in_any_zone(c);
  printf("maximum points in any zone is %d\n",max_points);
  vectors=create_bem_vectors(&bvv,&bcv,max_points);

  max_steps=300;
  step_size=0.1;
 
  put_raster("P(0,0)=(0.0,0.0) P(50,50)=(20.0,20.0)",&ras); /* slow */
  put_raster("P(0,0)=(0.0,0.0) P(20,20)=(20.0,20.0)",&ras);
  show_raster(&ras);

  output1=open_file(0,"dc_area.out","w");
  output2=open_file(0,"dc_area2.out","w");
  output3=open_file(0,"risk.out","w");
  output4=open_file(0,"risk2.out","w");
  column=0;
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
 	  L=streamline_loop(Pstart,c,1,max_steps,step_size,
 			    (path *)NULL,vectors,&voltage) ;
 	  v=velocity(P,voltage.dV);
	  if(v>0.0) { risk=L/v; } else { risk=0.0; }

 	  put_buffer(buf_size,buffer[0],column,"%14.5e",L);
  	  column=put_buffer(buf_size,buffer[1],column,"%14.5e",risk);
	  put_buffer(64,buffer[2],
	      put_buffer(64,buffer[2],
		  put_buffer(64,buffer[2],0,"%14.5e",P[0]),
				      " %14.5e",P[1]),
		                  " %14.5e",L);
	  put_next_line(output2,buffer[2]);
	  put_buffer(64,buffer[2],
	      put_buffer(64,buffer[2],
		  put_buffer(64,buffer[2],0,"%14.5e",P[0]),
				      " %14.5e",P[1]),
		                  " %14.5e",risk);
	  put_next_line(output4,buffer[2]);
	}
      column=0;
      put_next_line(output1,buffer[0]);
      put_next_line(output2,"");
      put_next_line(output3,buffer[1]);
      put_next_line(output4,"");
      if(k%50==0) printf("\n");
      printf("#"); fflush(stdout);
      k=k+1;
   }
/*------------------------------------------------------*/

  printf("\n");
  fclose(output1);
  fclose(output2);
  fclose(output3);
  fclose(output4);
  destroy_catchment(c);
  destroy_bem_vectors(vectors);
  free((void *)buffer[0]);
  return(0);
}
/*--------------------------------------------------------*/

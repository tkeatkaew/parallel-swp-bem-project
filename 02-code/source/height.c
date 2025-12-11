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

#include "height.h"
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main()
{
  /*  char data[]="c_valley.txt"; */
  /*  char data[]="catchment27.txt"; */
  /*  char data[]="c_valley2.txt"; */
  /*  char data[]="catchment1.txt"; */

  char data[]="c_valley.txt"; 

  FILE *output1, *output2;
  bem_results voltage;
  bem_vectors *vectors;
  catchment *c;
  char *buffer[2];
  coordinates P;
  double pp;
  int buf_size,i,j,k,column,max_points,num_zones,new_z;
  matrix bvv, bcv;
  raster ras;

  trap_floating_errors();
  buf_size=512*13+1;
  buffer[0]=(char *)malloc((buf_size+64)*sizeof(char));
  buffer[1]=buffer[0]+buf_size;

  num_zones=catchment_zones(data);
  c=create_catchment(num_zones,16);
  get_catchment(data,c);

  plot_catchment(c,"catchment.out");

  max_points=max_points_in_any_zone(c);
  printf("maximum points in any zone is %d\n",max_points);
  vectors=create_bem_vectors(&bvv,&bcv,max_points);
  /*
  put_raster("P(0,0)=(0.0,0.0) P(100,100)=(20.0,20.0)",&ras);
  put_raster("P(0,0)=(0.0,0.0) P(50,50)=(20.0,20.0)",&ras);
  put_raster("P(0,0)=(-1.2,-0.1) P(50,50)=(1.2,1.1)",&ras); 
  put_raster("P(0,0)=(0.1,0.1) P(9,9)=(0.9,0.9)",&ras);
  */

  put_raster("P(0,0)=(0.0,0.0) P(20,20)=(20.0,20.0)",&ras);
  show_raster(&ras);

  printf("\n------after scan----\n");

  output1=open_file(0,"height.out","w");
  output2=open_file(0,"height2.out","w");
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
	  pp=calculate_inside_catchment(c,P,vectors,&voltage,&new_z);
	  if(new_z==(-1)) printf("point is outside catchment\n");

	  column=put_buffer(buf_size,buffer[0],column,"%14.5e",pp);

	  put_buffer(64,buffer[1],
	      put_buffer(64,buffer[1],
		  put_buffer(64,buffer[1],0,"%14.5e",P[0]),
				      " %14.5e",P[1]),
		                  " %14.5e",pp);
	  put_next_line(output2,buffer[1]);
	}
      column=0;
      put_next_line(output1,buffer[0]);
      put_next_line(output2,"");
      if(k%50==0) printf("\n");
      printf("#"); fflush(stdout);
      k=k+1;
   }
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

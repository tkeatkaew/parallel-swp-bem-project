#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*--------------------------------------------------------*/
#include "boundary_types.h"
#include "co_matrix_types.h"
#include "matrix_types.h"
#include "ten_matrix_types.h"
#include "memory_types.h"

#include "area.h"
#include "catchment.h"
#include "file.h"
#include "memory.h"
#include "path.h"
#include "scan.h"
#include "streamline.h"
#include "trapfloat.h"

#include "catcharea.h"
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
  char title[64];
  char data[]="catchment.txt"; /* step-01 */

  bem_vectors *vectors;
  catchment *c;
  char *buffer;
  double step_size,SCA;
  int buf_size,i,max_points,num_zones,max_steps,max_streams;
  matrix bvv, bcv;
  path **streamlines;
  section mouth;
  coordinates PA,PB;

  trap_floating_errors();
  buf_size=512*12+1;
  buffer=(char *)malloc(buf_size*sizeof(char));

  num_zones=catchment_zones(data);
  c=create_catchment(num_zones,16);
  get_catchment(data,c);
  plot_catchment(c,"catchment.out");
  max_points=max_points_in_any_zone(c);
  vectors=create_bem_vectors(&bvv,&bcv,max_points);

  max_steps=10000;
  step_size=1.0;

  double rm, dr, seta;
  rm=99.0;
  dr=0.001;

  step_size=atof(argv[1]);
  rm=atof(argv[2]);
  dr=atof(argv[3]);

  PA[0]=   1.0;
  PA[1]=1000.0;
  PB[0]=   1.0;
  PB[1]=1001.0;  
  put_sectionV2(1,PA,PB,&mouth);
  max_streams=mouth.n;
  streamlines=(path **)malloc(max_streams*sizeof(path *));
  for(i=0;i<max_streams;i++){
    streamlines[i]=create_path(max_steps,1,0);} 
  SCA=Cal_SCA(c,&mouth,1,max_steps,step_size,max_streams,
	      streamlines,vectors); //stream up  
  plot_streamlines(c,max_streams,streamlines,"catcharea-01.out");
  printf("%f\t%f\t",rm,SCA);

  PA[0]=250.0;
  PA[1]=  1.0;
  PB[0]=250.0;
  PB[1]=  2.0;  
  put_sectionV2(1,PA,PB,&mouth);
  max_streams=mouth.n;
  streamlines=(path **)malloc(max_streams*sizeof(path *));
  for(i=0;i<max_streams;i++){
    streamlines[i]=create_path(max_steps,1,0);} 
  SCA=Cal_SCA(c,&mouth,1,max_steps,step_size,max_streams,
	      streamlines,vectors); //stream up  
  plot_streamlines(c,max_streams,streamlines,"catcharea-02.out");
  printf("%f\t",SCA);

  PA[0]=1250.0;
  PA[1]=   1.0;
  PB[0]=1250.0;
  PB[1]=   2.0;  
  put_sectionV2(1,PA,PB,&mouth);
  max_streams=mouth.n;
  streamlines=(path **)malloc(max_streams*sizeof(path *));
  for(i=0;i<max_streams;i++){
    streamlines[i]=create_path(max_steps,1,0);} 
  SCA=Cal_SCA(c,&mouth,1,max_steps,step_size,max_streams,
	      streamlines,vectors); //stream up  
  plot_streamlines(c,max_streams,streamlines,"catcharea-03.out");
  printf("%f\n",SCA);
  
  for(i=0;i<max_streams;i++){
    streamlines[i]=destroy_path(streamlines[i]);}
  
  destroy_catchment(c);
  destroy_bem_vectors(vectors);
  free((void *)buffer);
  return(0);
}
/*--------------------------------------------------------*/

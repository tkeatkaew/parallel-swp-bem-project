/*----------------------------------------------------------*/
/*    calculate difference between two raster files a.out b.out  */
/*----------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "trapfloat.h"
#include "image.h"

#include "diffras.h"
/*----------------------------------------------------------*/
int main(argv,argc)
     int argv;
     char **argc;
{
  char file[64],file2[64];
  double rmserror;
  float *in_image, *in_image2;
  int nx,ny,nx2,ny2,len;

  trap_floating_errors();
  if(argv<3)
    {
      printf("usage: difference file_a.out file_b.out\n");
      exit(0);
    }
  if(strlen(argc[1])>59)
    {
      printf("file name too long: %s\n",argc[1]);
      exit(0);
    }
  strcpy(file,argc[1]);
  if(strlen(argc[2])>59)
    {
      printf("file name too long: %s\n",argc[2]);
      exit(0);
    }
  strcpy(file2,argc[2]);

  measure_image(file,&nx,&ny);
  measure_image(file2,&nx2,&ny2);
  if(nx2!=nx||ny2!=ny)
    {
      printf("images are not the same size\n");
      exit(0);
    }

  in_image=(float *)malloc(nx*ny*sizeof(float));
  in_image2=(float *)malloc(nx2*ny2*sizeof(float));
  
  load_image(file,nx,ny,in_image);
  load_image(file2,nx2,ny2,in_image2);
  subtractfrom_image(nx,ny,in_image,in_image2);
  rmserror=rmsof_image(nx,ny,in_image);
  printf("rms difference is %e\n",rmserror);

  len=strcspn(file,".");
  strcpy(file+len,".err");  /* for output file name.err */
  write_image_ras(file,nx,ny,in_image);

  return(0);
}

/*----------------------------------------------------------*/
/*----------------------------------------------------------*/

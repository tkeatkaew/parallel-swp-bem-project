/*----------------------------------------------------------*/
/*    convert an ascii rater scan file to pgm image file    */
/*----------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "trapfloat.h"
#include "image.h"

#include "makepgm.h"
/*----------------------------------------------------------*/
int main(argv,argc)
     int argv;
     char **argc;
{
  char file[64];
  double maxval;
  float *in_image, *out_image;
  int nx,ny,len,invert,quantize;

  trap_floating_errors();
  if(argv<2)
    {
      printf("usage: makepgm file_name.out\n");
      exit(0);
    }
  if(strlen(argc[1])>59)
    {
      printf("file name too long: %s\n",argc[1]);
      exit(0);
    }
  strcpy(file,argc[1]);
  invert=0;
  if(argv>2)
    {
      sscanf(argc[2],"%d",&invert);
    }
  quantize=0;
  if(argv>3)
    {
      sscanf(argc[3],"%d",&quantize);
    }

  measure_image(file,&nx,&ny);

  in_image=(float *)malloc(nx*ny*sizeof(float));
  out_image=(float *)malloc(512*512*sizeof(float));
  
  load_image(file,nx,ny,in_image);
  maxval=maxof_image(nx,ny,in_image);
  rescale_image(nx,ny,in_image,maxval,256.0);
  if(invert==(-1))
    {
      addto_image(nx,ny,in_image,-128.0);
      invert_image(nx,ny,in_image);
      addto_image(nx,ny,in_image,128.0);
    }
  enlarge_image(nx,ny,in_image,512,512,out_image);
  if(quantize>0)
    {
      quantize_image(512,512,out_image,256.0,quantize);
    }

  len=strcspn(file,".");
  strcpy(file+len,".pgm");  /* for output file name.pgm */
  write_image_pgm(file,512,512,out_image);
  return(0);
}

/*----------------------------------------------------------*/
/*----------------------------------------------------------*/

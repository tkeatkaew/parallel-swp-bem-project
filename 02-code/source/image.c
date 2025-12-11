/*----------------------------------------------------------------------------------*/
/*-------------------------------- image.c -----------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing simple image operations */
/*----------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*----------------------------------------------------------------------------------*/

#include "image.h"
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void measure_image(file,nx,ny)
     char *file;
     int *nx, *ny;
{
  FILE *input;
  int status;
  int this_line,last_line,start_line,end_line;
  int this_item,last_item,start_item,end_item;
  unsigned char byte;

  input=fopen(file,"r");
  if(input==(FILE *)NULL) 
    {
      printf("cannot find image file: %s\n",file);
      exit(0);
    }
  this_line=0;  start_line=0;  end_line=0;
  this_item=0;  start_item=0;  end_item=0;
  status=1;
  while(status==1)
    {
      last_line=this_line;
      last_item=this_item;
      status=fread(&byte,1,1,input);
      if(status==1)
	{
	  switch (byte)
	    {
	    case '\n':
	      this_line=0;
	    case ' ':
	    case '':
	      this_item=0;
	      break;
	    default:
	      this_line=1;
	      this_item=1;
	    }
	  if(last_item==0&&this_item==1) start_item=start_item+1;
	  if(last_item==1&&this_item==0) end_item=end_item+1;
	  if(last_line==0&&this_line==1) start_line=start_line+1;
	  if(last_line==1&&this_line==0) end_line=end_line+1;
	}
    }
  (*ny)=start_line;
  (*nx)=start_item/start_line;
  if((*nx)*(*ny)!=start_item)
    {
      printf("lines not the same length in image file: %s\n",file);
      fclose(input);
      exit(0);
    }
  printf("image file: %s has (%dx%d) pixels\n",file,*nx,*ny);
  fclose(input);
}

/*----------------------------------------------------------------------------------*/
void measure_image_xyz(file,nx,ny)
     char *file;
     int *nx, *ny;
{
  FILE *input;
  int status;
  int this_block,last_block,start_block,end_block;
  int this_line,last_line,start_line,end_line;
  int this_item,last_item,start_item,end_item;
  unsigned char byte;

  input=fopen(file,"r");
  if(input==(FILE *)NULL) 
    {
      printf("cannot find image file: %s\n",file);
      exit(0);
    }
  this_block=0;  start_block=0;  end_block=0;
  this_line=0;  start_line=0;  end_line=0;
  this_item=0;  start_item=0;  end_item=0;
  status=1;
  while(status==1)
    {
      last_block=this_block;
      last_line=this_line;
      last_item=this_item;
      status=fread(&byte,1,1,input);
      if(status==1)
	{
	  switch (byte)
	    {
	    case '\n':
	      if(this_line==0) this_block=0;
	      this_line=0;
	    case ' ':
	    case '':
	      this_item=0;
	      break;
	    default:
	      this_block=1;
	      this_line=1;
	      this_item=1;
	    }
	  if(last_item==0&&this_item==1) start_item=start_item+1;
	  if(last_item==1&&this_item==0) end_item=end_item+1;
	  if(last_line==0&&this_line==1) start_line=start_line+1;
	  if(last_line==1&&this_line==0) end_line=end_line+1;
	  if(last_block==0&&this_block==1) start_block=start_block+1;
	  if(last_block==1&&this_block==0) end_block=end_block+1;
	}
    }
  (*ny)=start_block;
  (*nx)=start_line/start_block;
  if((*nx)*(*ny)*3!=start_item)
    {
      printf("lines not the same length in image file: %s\n",file);
      fclose(input);
      exit(0);
    }
  printf("image file: %s has (%dx%d) pixels\n",file,*nx,*ny);
  fclose(input);
}

/*----------------------------------------------------------------------------------*/
void load_image(file,nx,ny,image)
     char *file;
     int nx,ny;
     float *image;
{
  FILE *input;
  int i,n,status,count;

  input=fopen(file,"r");
  if(input==(FILE *)NULL) 
    {
      printf("cannot find image file: %s\n",file);
      exit(0);
    }
  count=0;
  n=nx*ny;
  for(i=0;i<n;i++)
    {
      status=fscanf(input," %f",&image[i]);
      if(status!=1) 
	{
	  printf("error reading file: %s after %d/%d values\n",file,count,n);
	  fclose(input);
	  exit(0);
	}
      count=count+1;
    }
  printf("loaded image file: %s\n",file);
  fclose(input);
}

/*----------------------------------------------------------------------------------*/
void load_image_xyz(file,nx,ny,x,y,image)
     char *file;
     int nx,ny;
     float *x, *y, *image;
{
  FILE *input;
  int i,n,status,count;

  input=fopen(file,"r");
  if(input==(FILE *)NULL) 
    {
      printf("cannot find image file: %s\n",file);
      exit(0);
    }
  count=0;
  n=nx*ny;
  for(i=0;i<n;i++)
    {
      status=fscanf(input," %f %f %f",&x[i],&y[i],&image[i]);
      if(status!=3) 
	{
	  printf("error reading file: %s after %d/%d values\n",file,count,n);
	  fclose(input);
	  exit(0);
	}
      count=count+1;
    }
  printf("loaded image file: %s\n",file);
  fclose(input);
}

/*----------------------------------------------------------------------------------*/
double maxof_image(nx,ny,image)
     int nx,ny;
     float *image;
{
  int i,n;
  double maxval,thisval;

  n=nx*ny;
  maxval=image[0];
  for(i=1;i<n;i++)
    {
      thisval=image[i];
      if(thisval>maxval) maxval=thisval;
    }
  return(maxval);
}

/*----------------------------------------------------------------------------------*/
double rmsof_image(nx,ny,image)
     int nx,ny;
     float *image;
{
  int i,n;
  double rms,thisval;

  n=nx*ny;
  thisval=image[0];
  rms=thisval*thisval;
  for(i=1;i<n;i++)
    {
      thisval=image[i];
      rms=rms+thisval*thisval;
    }
  rms=sqrt(rms/((double)n));
  return(rms);
}

/*----------------------------------------------------------------------------------*/
void addto_image(nx,ny,image,value)
     int nx,ny;
     float *image;
     double value;
{
  int i,n;

  n=nx*ny;
  for(i=0;i<n;i++)
    {
      image[i]=image[i]+value;
    }
}

/*----------------------------------------------------------------------------------*/
void invert_image(nx,ny,image)
     int nx,ny;
     float *image;
{
  int i,n;

  n=nx*ny;
  for(i=0;i<n;i++)
    {
      image[i]=(-image[i]);
    }
}

/*----------------------------------------------------------------------------------*/
void rescale_image(nx,ny,image,oldmax,newmax)
     int nx,ny;
     float *image;
     double oldmax,newmax;
{
  int i,n;
  double scale;

  n=nx*ny;
  scale=newmax/oldmax;
  for(i=0;i<n;i++)
    {
      image[i]=image[i]*scale;
    }
}

/*----------------------------------------------------------------------------------*/
void subtractfrom_image(nx,ny,image,image2)
     int nx,ny;
     float *image, *image2;
{
  int i,n;

  n=nx*ny;
  for(i=0;i<n;i++)
    {
      image[i]=image[i]-image2[i];
    }
}

/*----------------------------------------------------------------------------------*/
void quantize_image(nx,ny,image,maxval,n_levels)
     int nx,ny,n_levels;
     float *image;
     double maxval;
{
  int i,n;
  double scale,scale2,val;

  n=nx*ny;
  scale=n_levels/maxval;
  scale2=(n_levels-1)/maxval;
  for(i=0;i<n;i++)
    {
      val=floor(image[i]*scale)/scale2;
      if(val>maxval) val=maxval;
      image[i]=val;
    }
}

/*----------------------------------------------------------------------------------*/
void enlarge_image(nx,ny,image,kx,ky,new_image)
     int nx,ny;
     float *image;
     int kx,ky;
     float *new_image;
{
  int i,j,y0,y1,x0,x1;
  double x,y,xscale,yscale,a0,a1;

  xscale=(double)(nx-1)/(double)(kx-1);
  yscale=(double)(ny-1)/(double)(ky-1);

  for(j=0;j<ky;j++)
    {
      y=j*yscale;
      y0=floor(y);      if(y0<0) y0=0;      if(y0>ny-2) y0=ny-2;
      y1=y0+1;
      for(i=0;i<kx;i++)
	{
	  x=i*xscale;
	  x0=floor(x);      if(x0<0) x0=0;      if(x0>nx-2) x0=nx-2;
	  x1=x0+1;
	  a0=image[x0+nx*y0]*(x1-x)+image[x1+nx*y0]*(x-x0);
	  a1=image[x0+nx*y1]*(x1-x)+image[x1+nx*y1]*(x-x0);
	  x =a0*(y1-y)+a1*(y-y0);
	  new_image[i+kx*j]=x;
	}
    }
}

/*----------------------------------------------------------------------------------*/
void write_image_pgm(file,nx,ny,image)
     char *file;
     int nx,ny;
     float *image;
{
  FILE *output;
  int i,j,val,count,status;
  unsigned char byte;

  output=fopen(file,"w");
  if(output==(FILE *)NULL) 
    {
      printf("cannot create image file: %s\n",file);
      exit(0);
    }
  fprintf(output,"P5\n%d %d\n255\n",nx,ny);

  count=0;
  for(j=ny-1;j>=0;j--)
    {
      for(i=0;i<nx;i++)
	{
	  val=image[i+nx*j]; /* truncate to integer */
	  if (val<0) val=0;
	  if (val>255) val=255;
	  byte=val;
	  status=fwrite(&byte,1,1,output);
	  if(status!=1) 
	    {
	      printf("error writing file: %s after %d/%d values\n",file,count,nx*ny);
	      fclose(output);
	      exit(0);
	    }
	  count=count+1;
	}
    }
  printf("wrote image file: %s\n",file);
  fclose(output);
}

/*----------------------------------------------------------------------------------*/
void write_image_ras(file,nx,ny,image)
     char *file;
     int nx,ny;
     float *image;
{
  FILE *output;
  int i,j,count,status;
  double val;
  unsigned char byte;

  output=fopen(file,"w");
  if(output==(FILE *)NULL) 
    {
      printf("cannot create image file: %s\n",file);
      exit(0);
    }

  count=0;
  for(j=0;j<ny;j++)
    {
      for(i=0;i<nx;i++)
	{
	  val=image[count];
	  fprintf(output,"%14.5e ",val);
	  count=count+1;
	}
      fprintf(output,"\n");
    }
  printf("wrote image file: %s\n",file);
  fclose(output);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void write_image_xyz(file,nx,ny,x,y,image)
     char *file;
     int nx,ny;
     float *x,*y,*image;
{
  FILE *output;
  int i,j,count,status;
  double val;
  unsigned char byte;

  output=fopen(file,"w");
  if(output==(FILE *)NULL) 
    {
      printf("cannot create image file: %s\n",file);
      exit(0);
    }

  count=0;
  for(j=0;j<ny;j++)
    {
      for(i=0;i<nx;i++)
	{
	  val=image[count];
	  fprintf(output,"%14.5e %14.5e %14.5e\n",x[count],y[count],val);
	  count=count+1;
	}
      fprintf(output,"\n");
    }

  printf("wrote image file: %s\n",file);
  fclose(output);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/


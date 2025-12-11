/*----------------------------------------------------------------------------------*/
/*-------------------------------- - path.c  ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on path structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"
#include "file.h"
#include "path.h"
/*----------------------------------------------------------------------------------*/
/* create a path */
/*----------------------------------------------------------------------------------*/
path *create_path(points,make_xy,make_value)
     int points;
     int make_xy, make_value;
{
  path *p;
  double *values;
  coordinates *coordinate;

  p=(path *)malloc(sizeof(path));
  if(p==(path *)NULL)
    {
      printf("error allocating memory for path\n");
      exit(0);
    }
  
  values=(double *)NULL;  /* give initial pointer */
  if(make_value==1)
    {
      
      values=(double *)malloc(points*sizeof(double));
      if(values==(double *)NULL)
	{
	  printf("error allocating memory for path\n");
	  exit(0);
	}
    }  
  
  coordinate=(coordinates  *)NULL;  /* give initial pointer */
  if(make_xy==1)
    {
      coordinate=(coordinates *)malloc(points*sizeof(coordinates));
      if(coordinate==(coordinates *)NULL)
	{
	  printf("error allocating memory for path\n");
	  exit(0);
	}
    }  

  p->links=0;
  p->close=0;
  p->reverse=0;
  p->value=values;
  p->xy=coordinate;
  p->points=points;
  return(p);
}

/*----------------------------------------------------------------------------------*/
/* destroy a path */
/*----------------------------------------------------------------------------------*/

path *destroy_path(x)
     path *x;
{
  //printf("address: %x, links: %d\n",(unsigned int)x,x->links); //old version
  //printf("address: %p, links: %d\n",x,x->links); //new version
  if(x->links>0)
    {
      x->links=x->links-1;
    }
else
  {
    if(x->value!=(double *)NULL)   free((void *)x->value);
    if(x->xy!=(coordinates *)NULL) free((void *)x->xy);
    if(x!=(path *)NULL) free((void *)x);
  }

  return((void *)NULL);
}

/*----------------------------------------------------------------------------------*/
/* check memory to see if allocated already */
/*----------------------------------------------------------------------------------*/
void check_value_memory(p)
     path *p;
{
  if(p->value==(double *)NULL)
    {
      printf("you have forgotten to allocate memory for the values\n");
      exit(0);
    }
}

void check_xy_memory(p)
     path *p;
{
  if(p->xy==(coordinates *)NULL)
    {
      printf("you have forgotten to allocate memory for the xy coordinates\n");
      exit(0);
    }
}

/*----------------------------------------------------------------------------------*/
/* check index */
/*----------------------------------------------------------------------------------*/
void check_path_index(p,i)
     path *p;
     int i;
{
  if(p->close==0)
    {
      if(0>i||i>=p->points)
	{
	  printf("path index (%d) out of bounds (0,%d)\n",i,p->points-1);
	  exit(0);
	}
    }
  else
    {
      if(i<0)
	{
	  printf("path index (%d) out of bounds (0,*)\n",i);
	  exit(0);      
	}
    }
}

/*----------------------------------------------------------------------------------*/
/* reverse a path */
/*----------------------------------------------------------------------------------*/
void reverse_path(p)
     path *p;
{
  p->reverse=1-p->reverse;
}

/*----------------------------------------------------------------------------------*/
/* close a path */
/*----------------------------------------------------------------------------------*/
void close_path(p)
     path *p;
{
  p->close=1;
}

/*----------------------------------------------------------------------------------*/
/* open a path */
/*----------------------------------------------------------------------------------*/
void open_path(p)
     path *p;
{
  p->close=0;
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* get value from path */
/*----------------------------------------------------------------------------------*/
double get_path_value(p,i)
     path *p;
     int i;
{
  double value;
  
  check_value_memory(p);
  check_path_index(p,i);
  i=i%p->points;
  if(p->reverse==0)
    {
      value=p->value[i];
    }
  /*
  else
    {
      if(i==0)  value=p->value[i];      
      else     value=p->value[p->points-i];      
    }
  */
  
  else
    {
      value=p->value[p->points-1-i];      
    }
  
  return(value);
}

/*----------------------------------------------------------------------------------*/
/* get xy from path */
/*----------------------------------------------------------------------------------*/

void get_path_xy(p,i,xy)
     path *p;
     int i;
     coordinates xy;
{
  check_xy_memory(p);
  check_path_index(p,i);
  i=i%p->points;
  if(p->reverse==0)
    {
      xy[0]=p->xy[i][0]; 
      xy[1]=p->xy[i][1]; 
    }
  /*
  else
    {
      if(i==0)  
	{
	        xy[0]=p->xy[i][0];
		xy[1]=p->xy[i][1];
	}
      else 
	{
	        xy[0]=p->xy[p->points-i][0];
		xy[1]=p->xy[p->points-i][1];
	}
    }
  */
  
  else
    {
      xy[0]=p->xy[p->points-1-i][0];
      xy[1]=p->xy[p->points-1-i][1];
    }
  

}

/*----------------------------------------------------------------------------------*/
/* put value into path */
/*----------------------------------------------------------------------------------*/

void put_path_value(p,i,val)
     path *p;
     int i;
     double val;
{
  check_value_memory(p);
  check_path_index(p,i);
  i=i%p->points;
  if(p->reverse==0)
    {
      p->value[i]=val;
    }
  else
    {
      p->value[p->points-1-i]=val;
    }
}

/*----------------------------------------------------------------------------------*/
/* put coordinates into path */
/*----------------------------------------------------------------------------------*/
void put_path_xy(p,i,xy)
     path *p;
     int i;
     coordinates xy;
{
  check_xy_memory(p);
  check_path_index(p,i);
  i=i%p->points;
  if(p->reverse==0)
    {
      p->xy[i][0]=xy[0]; 
      p->xy[i][1]=xy[1];       
    }
  else
    {
      p->xy[p->points-1-i][0]=xy[0]; 
      p->xy[p->points-1-i][1]=xy[1];       
    }
}

/*----------------------------------------------------------------------------------*/
int num_points_in_path(p)
     path *p;
{
  return(p->points);
}
/*----------------------------------------------------------------------------------*/
/* show path information */
/*----------------------------------------------------------------------------------*/

void show_path_info(p)
     path *p;
{
  //  printf("links=%d, close=%d, reverse=%d, points=%d, xy=0x%0X, value=0x%0X\n",
  //	 p->links,p->close,p->reverse,p->points,(unsigned int)p->xy,
  //	 (unsigned int)p->value);

  //printf("links=%d, close=%d, reverse=%d, points=%d, xy=0x%0X, value=0x%0X\n",
	// p->links,p->close,p->reverse,p->points,p->xy,p->value);

  printf("links=%d, close=%d, reverse=%d, points=%d",
	 p->links,p->close,p->reverse,p->points);
}

/*----------------------------------------------------------------------------------*/
/* show a path */
/*----------------------------------------------------------------------------------*/

void show_path(p)
     path *p;
{
  int i;
  int point;
  coordinates a;

  point=p->points;
  
  if(p->value!=(double *)NULL)
  { 
    printf("show value \n");
    for(i=0;i<point;i++)
	{
    	  printf("%7.6f \n",get_path_value(p,i));
	}
    printf("\n");
  }
  
  if(p->xy!=(coordinates *)NULL)
  { 
    printf("show coordinates(x,y)\n");
     
    for(i=0;i<point;i++)
	{
	  get_path_xy(p,i,a);
	  printf("(%7.6f %7.6f)\n",a[0],a[1]);
	}
    printf("\n");
  } 
   
}

/*----------------------------------------------------------------------------------*/
/* path file operations */
/*----------------------------------------------------------------------------------*/
int path_length(file)
     unsigned char *file;
{
  int length;

  length=count_lines(1,file);
  return(length);
}

/*----------------------------------------------------------------------------------*/
#define NBYTES 96
/*----------------------------------------------------------------------------------*/
void get_path(file,p)
     unsigned char *file;
     path *p;
{
  unsigned char buffer[NBYTES];
  FILE *input;
  int i,n,status,data_missing;
  double value;
  coordinates xy;

  data_missing=0;
  input=open_file(1,file,"r");
  n=p->points;
  for(i=0;i<n;i++)
    {
      status=get_next_line_verbose(input,1,NBYTES,buffer);
      if(sscanf(buffer," %lf %lf %lf",&xy[0],&xy[1],&value)!=3)
	{
	  data_missing=1;
	}
      put_path_value(p,i,value);
      put_path_xy(p,i,xy);
    }
  if(status==0)
    {
      printf("Fewer than %d lines in file '%s'\n",n,file);
      exit(0);
    }
  if(data_missing==1)
    {
      printf("Fewer than 3 data values/line in file '%s'\n",file);
      exit(0);
    }
}

/*----------------------------------------------------------------------------------*/
void plot_path(p,file)
     path *p;
     unsigned char *file;
{
  unsigned char buffer[NBYTES];
  FILE *output;
  int i,n,offset;
  coordinates xy;

  output=open_file(0,file,"w");
  n=p->points;
  for(i=0;i<n;i++)
    {
      get_path_xy(p,i,xy);
      offset=put_buffer(NBYTES,buffer,0,"%f",xy[0]);
      offset=put_buffer(NBYTES,buffer,offset," %f",xy[1]);
      put_next_line(output,buffer);
    }
  if(offset>NBYTES-1)
    {
      printf("Output buffer not big enough (size=%d)\n",NBYTES);
      exit(0);
    }
  fclose(output);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

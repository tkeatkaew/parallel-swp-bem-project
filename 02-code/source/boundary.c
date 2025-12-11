/*----------------------------------------------------------------------------------*/
/*------------------------------ boundary.c ----------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on boundary structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"

#include "file.h"
#include "path.h"

#include "boundary.h"
/*----------------------------------------------------------------------------------*/
int boundary_loops(file)
     unsigned char *file;
{
  int length;

  length=count_lines(1,file);
  return(length);
}

/*----------------------------------------------------------------------------------*/
boundary *create_boundary(n)
     int n;
{
  int i;
  boundary *b;
  int *level;
  path **loops;

  b=(boundary *)malloc(sizeof(boundary));
  if(b==(boundary *)NULL)
    {
      printf("Failed to allocate memory for boundary\n");
      exit(0);
    }

  level=(int *)malloc(n*sizeof(int));
  if(level==(int *)NULL)
    {
      printf("Failed to allocate memory for array of levels for paths\n");
      exit(0);
    }

  loops=(path **)malloc(n*sizeof(path *));
  if(loops==(path **)NULL)
    {
      printf("Failed to allocate memory for list of pointers to paths\n");
      exit(0);
    }

  b->curve=0;
  b->components=n;
  b->level=level;
  b->loop=loops;
  for(i=0;i<n;i++)
    {
      level[i]=0;
      loops[i]=(path *)NULL;
    }
  b->bvv=(double *)NULL;
  b->bcv=(double *)NULL;
  return(b);
}

/*----------------------------------------------------------------------------------*/
boundary *destroy_boundary(b)
     boundary *b;
{
  int i,n;
  int *level;
  path **loops;

  n=b->components;
  level=b->level;
  loops=b->loop;
  if(level!=(int *)NULL)
    {
      free((void *)level);
    }
  if(loops!=(path **)NULL)
    {
      for(i=0;i<n;i++)
	{
	  if(loops[i]!=(path *)NULL)
	    {
	      destroy_path(loops[i]);
	    }
	}
      free((void *)loops);
    }
  if(b->bvv!=(double *)NULL) free((void *)b->bvv);
  if(b->bcv!=(double *)NULL) free((void *)b->bcv);
  free((void *)b);
  return((boundary *)NULL);
}

/*----------------------------------------------------------------------------------*/
boundary *destroy_boundary_ignore_paths(b)
     boundary *b;
{
  int n;
  int *level;
  path **loops;

  n=b->components;
  level=b->level;
  loops=b->loop;
  if(level!=(int *)NULL)
    {
      free((void *)level);
    }
  if(loops!=(path **)NULL)
    {
      free((void *)loops);
    }
  free((void *)b);
  return((boundary *)NULL);
}

/*----------------------------------------------------------------------------------*/
int num_points_in_zone(b)
     boundary *b;
{
  int i,n,n_points;
  path **loops;

  n=b->components;
  loops=b->loop;
  n_points=0;
  for(i=0;i<n;i++)
    {
      if(loops[i]!=(path *)NULL)
	{
	  n_points=n_points+num_points_in_path(loops[i]);
	}
    }
  return(n_points);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void get_boundary(file,b)
     unsigned char *file;
     boundary *b;
{
  unsigned char path_file[64];
  FILE *input;
  int i,n,status,points;
  path **loops;

  input=open_file(1,file,"r");
  n=b->components;
  loops=b->loop;

  for(i=0;i<n;i++)
    {
      status=get_next_line(input,64,path_file); /* gets the name of the path file */
      points=path_length(path_file);
      loops[i]=create_path(points,1,1);
      get_path(path_file,loops[i]);
    }
  if(status==0)
    {
      printf("Fewer than %d lines in file '%s'\n",n,file);
      exit(0);
    }
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void plot_boundary(b,file)
     boundary *b;
     unsigned char *file;
{
  unsigned char buffer[64];
  FILE *output;
  path **loops;
  path *p;
  int i,j,n,offset,n_paths;
  coordinates xy;

  n_paths=b->components;
  loops=b->loop;

  output=open_file(0,file,"w");
  for(j=0;j<n_paths;j++)
    {
      p=loops[j];
      n=p->points;
      for(i=0;i<n;i++)
	{
	  get_path_xy(p,i,xy);
	  offset=put_buffer(64,buffer,0,"%f",xy[0]);
	  offset=put_buffer(64,buffer,offset," %f",xy[1]);
	  put_next_line(output,buffer);
	}
      if(offset>63)
	{
	  printf("Output buffer not big enough (size=%d)\n",64);
	  exit(0);
	}
    }
  fclose(output);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/






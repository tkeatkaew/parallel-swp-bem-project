/*----------------------------------------------------------------------------------*/
/*------------------------------ catchment.c ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on catchment structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"

#include "boundary.h"
#include "file.h"
#include "geometry.h"
#include "path.h"
#include "path_list.h"

#include "catchment.h"
/*----------------------------------------------------------------------------------*/
int catchment_zones(file)
    unsigned char *file;
{
  int length;

  length=count_lines(1,file);
  return(length);
}

/*----------------------------------------------------------------------------------*/
catchment *create_catchment(zones,paths)
     int zones,paths;
{
  int i;
  catchment *c;
  boundary **zone_list;
  path_link *p_link;

  c=(catchment *)malloc(sizeof(catchment));
  if(c==(catchment *)NULL) {
      printf("Failed to allocate memory for catchment\n");
      exit(0); }

  zone_list=(boundary **)malloc(zones*sizeof(boundary *));
  if(zone_list==(boundary **)NULL) {
      printf("Failed to allocate memory for list of pointers to boundaries\n");
      exit(0); }

  p_link=(path_link *)malloc(paths*sizeof(path_link));
  if(p_link==(path_link *)NULL) {
      printf("Failed to allocate memory for list of links to paths\n");
      exit(0); }

  for(i=0;i<zones;i++)
    {
      zone_list[i]=(boundary *)NULL;
    }
  for(i=0;i<paths;i++)
    {
      p_link[i].path_p=(path *)NULL;
      p_link[i].name[0]='\0';
    }

  c->num_zones=0;
  c->max_zones=zones;
  c->previous_zone=(-1);
  c->zones=zone_list;
  c->num_paths=0;
  c->max_paths=paths;
  c->path_list=p_link;
  return(c);
}

/*----------------------------------------------------------------------------------*/
catchment *destroy_catchment(c)
     catchment *c;
{
  int i,n;
  boundary **zones;
  path_link *links;
  
  n=c->num_paths;
  links=c->path_list;
  destroy_path_list(n,links);

  n=c->num_zones;
  zones=c->zones;
  if(zones!=(boundary **)NULL)
    {
      for(i=0;i<n;i++)
	{
	  destroy_boundary_ignore_paths(zones[i]);  /* paths destroyed already */
	}
      free((void *)zones);
    }

  free((void *)c);
  return((catchment *)NULL);
}

/*----------------------------------------------------------------------------------*/
int max_points_in_any_zone(c)
     catchment *c;
{
  int i,n,temp,max;
  boundary **zones;
  
  n=c->num_zones;
  zones=c->zones;
  max=0;
  for(i=0;i<n;i++)
    {
      temp=num_points_in_zone(zones[i]);
      if(temp>max) max=temp;
    }
  return(max);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
#define NBYTES 96
/*----------------------------------------------------------------------------------*/
void get_catchment(file,c)
     unsigned char *file;
     catchment *c;
{
  unsigned char boundary_file[NBYTES];
  unsigned char path_file[NBYTES];
  FILE *input, *path_input;
  int i,j,index,n,status,n_paths;
  boundary **zones;
  path_link *p_list;
  path *new_path;
  
  input=open_file(1,file,"r");

  n=c->max_zones;
  zones=c->zones;
  p_list=c->path_list;
  for(i=0;i<n;i++)
    {                                    /* gets the name of the boundary file */
      status=get_next_line_verbose(input,1,NBYTES,boundary_file);
      n_paths=boundary_loops(boundary_file);
      zones[i]=create_boundary(n_paths);

      path_input=open_file(1,boundary_file,"r");
      for(j=0;j<n_paths;j++)
	{
	  status=get_next_line_verbose(path_input,1,NBYTES,path_file);
	  index=search_path_list(path_file,c->num_paths,p_list);
	  if(index<0)
	    {
	      if(c->num_paths<c->max_paths)
		{
		  index=load_path_list(path_file,c->num_paths,p_list);
		  c->num_paths=c->num_paths+1;
		}
	      else
		{
		  printf("error :- only %d paths reserved for catchment\n",c->max_paths);
		  printf("         but trying to load more than %d\n",c->max_paths);
		  exit(0);
		}
	    }
	  new_path=get_path_list(index,p_list);
	  zones[i]->loop[j]=new_path;
	}
      fclose(path_input);

      mark_curve(zones[i]);
      /* show_curve(zones,i); */
      mark_paths(zones[i]);
      /* show_paths(zones,i); */
      
      c->num_zones=c->num_zones+1;
      if(status==0)
	{
	  printf("Fewer than %d lines in file '%s'\n",n,file);
	  exit(0);
	}
    }
  fclose(input);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void plot_catchment(c,file)
     catchment *c;
     unsigned char *file;
{
  unsigned char buffer[NBYTES];
  FILE *output;
  path **loops;
  path *p;
  int i,j,k,n,offset,n_paths,nz;
  coordinates xy;
  boundary *b;

  output=open_file(0,file,"w");

  nz=c->num_zones;
  for(k=0;k<nz;k++)
    {
      b=c->zones[k];
      n_paths=b->components;
      loops=b->loop;
      for(j=0;j<n_paths;j++)
	{
	  p=loops[j];
	  n=p->points;
	  for(i=0;i<n+1;i++) /* plot first point at beginning and end (2 times!) */
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
	  put_next_line(output,"\0");
	}
    }
  fclose(output);
}
/*----------------------------------------------------------------------------*/
/*---------------tell point p locate at th-zone (second version) -------------*/ 
/*----------------------------------------------------------------------------*/
int check_each_zone(c,P)
     catchment *c;
     coordinates P;
{
  int k;
  boundary *b;
  int nz, not_this_zone;

  nz=c->num_zones;
  not_this_zone=1;
  for(k=0; k<nz && not_this_zone; k++) /* each zone */
    {
      b=c->zones[k];
      if(check_zone(b,P)==1) not_this_zone=0;
    }
  if(not_this_zone==1) k=0;  /* outside this catchment */
  return(k-1); //return index of the zone, return -1 when point p is outside catchment
}   
/*----------------------------------------------------------------------------------*/
/*----------------------- check_zone -----------------------------------------------*/
/*----------------------------------------------------------------------------------*/
int check_zone(b,P)
     boundary *b;
     coordinates P;
{
  int n_path,j,count;
  path **loops;
  double s,d;
  int segment;
  
  n_path=b->components;
  loops=b->loop;
  count=0;

  reverse_zone(b);
  for(j=0;j<n_path;j++){      /* each path */
    if(distance_to_path(P,loops[j],&d,&s,&segment)==1){
      count=count+1;}}
  reverse_zone(b);

  j=0;                        /* P is not in this zone */
  if(count==n_path) j=1;      /* P is in this zone (inside all paths) */
  return(j);
}

/*----------------------------------------------------------------------------------*/
/*----------------------- reverse_zone----------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void reverse_zone(b)
     boundary *b;
{
  int j,zone_type,path_type;
  path *this_path;

  zone_type=b->curve;

  for(j=0;j<b->components;j++)
    {
      this_path=b->loop[j];
      path_type=b->level[j];
      if((zone_type==0 && path_type==1) || (zone_type==1 && path_type==0)) 
	{
	  reverse_path(this_path);
	}
      /*show_path(b->loop[j]);*/
    }
}	
/*----------------------------------------------------------------------------------*/
/*----------------------- reverse all paths in zone --------------------------------*/
/*----------------------------------------------------------------------------------*/
void reverse_all_paths(b)
     boundary *b;
{
  int j;
  
  for(j=0;j<b->components;j++)
    {
      reverse_path(b->loop[j]);
    }
}	
/*----------------------------------------------------------------------------------*/
/*----------------------- distance to path------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* returns 0 if point outside; 1 if point inside */
/* points at infinity are outside anticlockwise paths, and inside clockwise paths */
int distance_to_path(P,this_path,d,s,segment)
     coordinates P;
     double *s,*d;
     path *this_path;
     int *segment;
{
  int i,n_segment,imin,new_value;
  double dmin,dsq,x,y,y1,y2;
  double PminusQdotN;
  coordinates Qa,Qb; 

  n_segment=this_path->points;

  get_path_xy(this_path,0,Qa);
  x=Qa[0]-P[0];     y=Qa[1]-P[1];
  dmin=x*x+y*y;
  imin=0;
  for(i=1;i<n_segment;i++)
    {
      get_path_xy(this_path,i,Qa);
      x=Qa[0]-P[0];     y=Qa[1]-P[1];
      dsq=x*x+y*y;
      if(dsq<dmin)
	{
	  dmin=dsq;
	  imin=i;
	}
    }
  dmin=sqrt(dmin);
  (*s)=-0.5;
  (*d)=dmin;
  (*segment)=imin;

  new_value=0;
  for(i=0;i<n_segment;i++)
    {
      get_path_xy(this_path,i,Qa);
      get_path_xy(this_path,i+1,Qb);
      convert_PQ(Qa,Qb,P,&x,&y1,&y2);
      if(y1<=0.0 && y2>=0.0)
	{
	  x=fabs(x);
	  if(x<dmin)
	    {
	      new_value=1;
	      dmin=x;
	      imin=i;
	    }
	}
    }

  if(new_value==0)
    {
      get_path_xy(this_path,imin+n_segment-1,Qa);
      get_path_xy(this_path,imin,Qb);
      convert_PQ(Qa,Qb,P,&x,&y1,&y2);
      PminusQdotN=(-x);      
      get_path_xy(this_path,imin,Qa);
      get_path_xy(this_path,imin+1,Qb);
      convert_PQ(Qa,Qb,P,&x,&y1,&y2);
      PminusQdotN=PminusQdotN-x;
    }
  else
    {
      get_path_xy(this_path,imin,Qa);
      get_path_xy(this_path,imin+1,Qb);
      convert_PQ(Qa,Qb,P,&x,&y1,&y2);
      PminusQdotN=(-x);
      (*s)=-(y1+y2)/2.0/(y2-y1);
      (*d)=dmin;
      (*segment)=imin;
    }
  i=0;                      /* outside boundary */
  if(PminusQdotN<0.0) i=1;  /* inside boundary */
  return(i);
}
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* mark paths */
/* paths can be two types: outside the zone or inside the zone */
/* an outside path separates the zone from infinity      (value of level = 0) */
/* an inside path is separated from infinity by the zone (value of level = 1) */

void mark_paths(b)
     boundary *b;
{
  int n_paths, i;
  int bpath, outside, inside;

  /*
  printf("mark_paths\n");
  */
  n_paths=b->components;

  for(i=0;i<n_paths;i++) {        /* assume all paths are inside zone */
    b->level[i]=1; }
            
  bpath=count_paths(b,&outside,&inside);
  if(outside==1) {
    b->level[bpath]=0; }         /* 1 path outside zone */
  else {
    if(outside>1) {
      printf("error :- more than 1 paths outside zone\n");
      exit(0);}}

  if(outside+inside!=n_paths) {
    printf("error :- outside + inside paths for zone not same as total\n");
    exit(0); }
}

/*----------------------------------------------------------------------------------*/
/* mark the style of zone, according to all of the paths that make up its boundary */
/* all paths anti-clockwise -> 0  */
/* all paths      clockwise -> 1  */
void mark_curve(b)
     boundary *b;
{
  int n_paths, i;
  int clockwise;
  /*
  printf("mark_curve\n");
  */
  clockwise=0;
  n_paths=b->components;

  for(i=0;i<n_paths;i++) {
    clockwise=clockwise+find_orientation(b->loop[i]); }

  if(clockwise==n_paths) {
    b->curve=1; }
  else {
    if(clockwise==0) {
      b->curve=0; }
    else {
      printf("error : zone has mixed clockwise and anti-clockwise paths\n");
      exit(0); }}
}

/*----------------------------------------------------------------------------------*/
/* find orientation of path */
/* 0 = anti-clockwise;      1 = clockwise */
int find_orientation(this_path)
     path *this_path;
{
  int location;
  double d,s;
  int orient;
  coordinates Pmin,Pmax;

  find_limits(this_path,Pmin,Pmax);
  Pmin[0]=(3.0*Pmin[0]-Pmax[0])/2.0;
  Pmin[1]=(3.0*Pmin[1]-Pmax[1])/2.0; /* this point is in the same region as infinity */

  location=distance_to_path(Pmin,this_path,&d,&s,&orient);
  if(location==0) /* point is excluded, path must be normal (anti-clockwise) */
    {
      orient=0;
    }
  else /* point is included, path must be reversed (clockwise) */
    {
      orient=1; 
    }
  return(orient);
}

/*----------------------------------------------------------------------------------*/
void find_limits(this_path,min,max)
     path *this_path;
     coordinates min,max;
{
  coordinates P;
  int i,n_segment;

  get_path_xy(this_path,0,P);
  min[0]=P[0];      max[0]=P[0];
  min[1]=P[1];      max[1]=P[1];

  n_segment=this_path->points;
  for(i=1;i<n_segment;i++)
    {
      get_path_xy(this_path,i,P);
      if(P[0]<min[0])                       { min[0]=P[0]; }
      else	      { if(P[0]>max[0])	    { max[0]=P[0]; }}
      if(P[1]<min[1])	                    { min[1]=P[1]; }
      else	      { if(P[1]>max[1])	    { max[1]=P[1]; }}
    }
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* count the number of paths outside the zone, and inside the zone            */
/* return also: the index of the outside path (if there is one)               */
/*            : -1 (if there is no path outside the zone)                     */
/* if the zone is badly formed the sum of inside and outside paths            */
/*    will be less than the total (paths not on the boundary are not counted) */
int count_paths(b,outside_zone,inside_zone)
     boundary *b;
     int *outside_zone, *inside_zone;
{
  int n_paths,zone_type, index_path;
  int paths_enclosed,i,j;
  int inside,outside,segment;
  path *this_path;
  coordinates P;
  double d,s;

  outside=0;
  inside=0;
  index_path=(-1);
  zone_type=b->curve;
  n_paths=b->components;

  if(n_paths==1)
    {
      if(zone_type==0) { outside=1; index_path=0; } 
      else             { inside=1;                } 
    }

  else
    {
      if(zone_type==1) reverse_all_paths(b);
      for(j=0;j<n_paths;j++)
	{
	  this_path=b->loop[j]; /* test this path (j) */
	  paths_enclosed=0;
	  for(i=0;i<n_paths;i++)
	    {
	      if(i!=j)
		{
		  get_path_xy(b->loop[i],0,P); /* P is on path i */
		  if(distance_to_path(P,this_path,&d,&s,&segment)==1)
		    {
		      paths_enclosed=paths_enclosed+1;
		    }
		}
	    }
	  if(paths_enclosed==0)
	    {
	      inside=inside+1; /* this path (j) is inside zone */
	    }
	  if(paths_enclosed==n_paths-1)
	    {
	      outside=outside+1; /* this path (j) is outside zone */
	      index_path=j;
	    }
	}
      if(zone_type==1) reverse_all_paths(b);
    }
      
  if(outside>1)
    {
      printf("error :- found more than 1 path outside zone\n");
      exit(0);
    }

  (*outside_zone)=outside;
  (*inside_zone)=inside;
  return(index_path);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void show_curve(zones,i)
     boundary **zones;
     int i;
{
  printf("zone = %d, curve = %d, ",i,zones[i]->curve);
  if(zones[i]->curve==0)
    {
      printf("anti-clockwise, left , standard\n");
    }
  else
    {
      if(zones[i]->curve==1)
	{
	  printf("clockwise,      right, reversed\n");
	}
      else
	{
	  printf("bad value\n");
	}
    }
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void show_paths(zones,i)
     boundary **zones;
     int i;
{
  boundary *b;
  int n_paths,j;

  printf("zone = %d\n",i);

  b=zones[i];
  n_paths=b->components;
  for(j=0;j<n_paths;j++)
    {
      printf("  path = %3d, ",j);
      if(b->level[j]==0)
	{
	  printf(" outside zone, bounding path, separates zone from infinity\n");
	}
      else
	{
	  if(b->level[j]==1)
	    {
	  printf(" inside zone,  hole path,     separated from infinity by zone\n");
	    }
	  else
	    {
	      printf(" bad value\n");
	    }
	}
    }
}
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*--------------------------------------------------------*/
/* return minimum distance to closest path */
void check_each_path(c,P,d,s,segment,this_path)
     catchment *c;
     coordinates P;
     double *d,*s;
     int *segment;
     path **this_path;
{
  int num_path, index, test_segment;
  path_link *p_list;
  path *test_path;
  double test_d, d_min, test_s;
   
  num_path=c->num_paths;   
  /*printf(" num_path=%d",num_path);*/
  p_list=c->path_list;

  test_path=get_path_list(0,p_list);
  distance_to_path(P,test_path,&d_min,&test_s,&test_segment);
 
  (*s)=test_s;
  (*segment)=test_segment;
  (*this_path)=test_path;
  
  for(index=1;index<num_path;index++)
    {
      test_path=get_path_list(index,p_list);
      distance_to_path(P,test_path,&test_d,&test_s,&test_segment);
 
      if(test_d<d_min)
	{
	  d_min=test_d;
	  (*s)=test_s;
	  (*segment)=test_segment;
	  (*this_path)=test_path;
	}
    }
  (*d)=d_min;
}

/*--------------------------------------------------------*/

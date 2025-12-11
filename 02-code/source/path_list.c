/*----------------------------------------------------------------------------------*/
/*------------------------------ path_list.c  --------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on path list structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h"

#include "path.h"

#include "path_list.h"
/*----------------------------------------------------------------------------------*/
path_link *create_path_list(n)
     int n;
{
  int i;
  path_link *path_list;

  path_list=(path_link *)malloc(n*sizeof(path_link));
  if(path_list==(path_link *)NULL)
    {
      printf("Failed to allocate memory for path list\n");
      exit(0);
    }
  for(i=0;i<n;i++)
    {
      path_list[i].path_p=(path *)NULL;
      path_list[i].name[0]='\0';
    }
  return(path_list);
}

/*----------------------------------------------------------------------------------*/
path_link *destroy_path_list(n,path_list)
     int n;
     path_link *path_list;
{
  int i;

  if(path_list!=(path_link *)NULL)
    {
      for(i=0;i<n;i++)
	{
	  if(path_list[i].path_p!=(path *)NULL)
	    {
	      //printf("destroy path %d: ",i+1);
	      destroy_path(path_list[i].path_p);
	    }
	}
      free((void *)path_list);
    }
  return((path_link *)NULL);
}

/*----------------------------------------------------------------------------------*/
/* returns -1 if cannot find in list, or index of position in list */

int search_path_list(file_name,n,path_list)
     int n;
     unsigned char *file_name;
     path_link *path_list;
{
  int i,index,compare;

  index=-1;

  i=0;
  compare=1;
  while(compare!=0&&i<n)
    {
      compare=strncmp(file_name,path_list[i].name,31);
      if(compare==0) { index=i; }
      i=i+1;
    }
  return(index);
}
/*----------------------------------------------------------------------------------*/
/* loads path and returns index of position in list */
int load_path_list(file_name,index,path_list)
     int index;
     unsigned char *file_name;
     path_link *path_list;
{
  int n;
  path *this_path;

  n=path_length(file_name);
  this_path=create_path(n,1,1);
  close_path(this_path);
  get_path(file_name,this_path);
  path_list[index].path_p=this_path;
  strncpy(path_list[index].name,file_name,32);
  path_list[index].name[31]='\0';
  return(index);
}

/*----------------------------------------------------------------------------------*/
/* gets pointer to path from place in path list */
path *get_path_list(index,path_list)
     int index;
     path_link *path_list;
{
  path *this_path;

  this_path=path_list[index].path_p;
  return(this_path);
}


/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

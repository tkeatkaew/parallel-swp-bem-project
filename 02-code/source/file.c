/*----------------------------------------------------------------------------------*/
/*-------------------------------- - file.c  ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing simple file input/output */
/*----------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*----------------------------------------------------------------------------------*/

#include "file.h"
/*----------------------------------------------------------------------------------*/
static unsigned char lbuffer[128];
/*----------------------------------------------------------------------------------*/
int put_buffer(n,buffer,offset,format,value) 
    int n,offset;
     unsigned char *buffer, *format;
     double value;
{
  int status;

  status=snprintf(NULL,0,format,value); /* length of output string */
  if(status+1>n-offset)
    {
    //  printf("size of line + null longer than output buffer (%d)\n",n);
    //  printf("offset into buffer = %d, free buffer = %d, space required = %d\n",
	  //   offset,n-offset,status+1);
      status=0; /* nothing put into buffer */
      //exit(0);
    }
  else
    {
      status=sprintf(buffer+offset,format,value);
    }
  return(offset+status);
}

/*----------------------------------------------------------------------------------*/
/* opens the file 
   if use_path = 0, no path is attached
   if use_path = 1, the path is attached using the environment variable CATCHMENT */
 
FILE *open_file(use_path,file,mode)
     int use_path;
     char *file, *mode;
{
  FILE *input;
  char file_name[128];
  int length;
  int nn;
  nn=128;
  
  if(use_path==0)
    {
      strncpy(file_name,file,128);
      file_name[127]='\0';
    }
  else
    {
      //catchment_path(128,file_name);
      catchment_path(nn,file_name);
      length=strlen(file_name);
      strncat(file_name,file,128-length);
      file_name[127]='\0';
    }

  input=fopen(file_name,mode);
  if(input==(FILE *)NULL)
    {
      printf("Cannot open file: '%s' for access mode '%s'\n",file_name,mode);
      exit(0);
    }
  return(input);
}

/*----------------------------------------------------------------------------------*/
int count_lines(use_path,file)
     int use_path;
     char *file;
{
  FILE *input;
  int status,n;

  input=open_file(use_path,file,"r");
  
  n=0;
  status=get_next_line(input,128,lbuffer);
  while(status>0)
    {
      if(status==1) n=n+1;             /* count lines which are not comments */
      status=get_next_line(input,128,lbuffer);
    }
  return(n);
  fclose(input);
}

/*----------------------------------------------------------------------------------*/
int get_next_line(input,n,buffer)
     FILE *input;
     int n;
     unsigned char *buffer;
{
  int i,count,status;
  unsigned char byte;

  n=n-1;
  count=1;
  i=0;
  status=fread(&byte,1,1,input);
  while(status==1&&byte!='\n')
    {
      count=count+1;
      if(i<n)
	{
	  buffer[i]=byte;
	  i=i+1;
	}
      status=fread(&byte,1,1,input);
    }
  buffer[i]='\0';
  if(count>n)
    {
      printf("size of line + null (%d) longer than input buffer (%d)\n",count,n);
      exit(0);
    }
  if(status==0 && i==0)     {status=0;}   /* no more lines in file */
  else { if(buffer[0]!='#') {status=1;}   /* line does not start with comment character */
         else               {status=2;}}  /* line starts with comment character */ 
  return(status);
}

/*----------------------------------------------------------------------------------*/
int get_next_line_verbose(input,k,n,buffer)
     FILE *input;
     int n,k;
     unsigned char *buffer;
{
  int status,count;

  count=0;
  status=get_next_line(input,n,buffer);
  while(status==2) 
    {
      if(count<k)
	{
	  printf("%s\n",buffer);
	}
      count=count+1;
      status=get_next_line(input,n,buffer);
    } 
  return(status);
}

/*----------------------------------------------------------------------------------*/
void put_next_line(output,buffer)
     FILE *output;
     unsigned char *buffer;
{
  int n,status;

  n=strlen(buffer);
  status=fprintf(output,"%s\n",buffer);
  if(status<n)
    {
      printf("Failed to write data to output file\n");
      exit(0);
    }
}

/*----------------------------------------------------------------------------------*/
void catchment_path(n,c_path)
     int n;
     unsigned char *c_path;
{
  char *env;

  env=getenv("CATCHMENT");
  if(env==(char *)NULL)
    {
      printf("Cannot find environment variable: CATCHMENT\n");
      exit(0);
    }
  strncpy(c_path,env,n);
  c_path[n-1]='\0';
}

/*----------------------------------------------------------------------------------*/
void make_gpl_file(datafile,title,xrange,yrange)
     char *datafile, *title, *xrange, *yrange;
{
  FILE *output;
  int len;

  len=strcspn(datafile,".");
  strcpy(lbuffer,datafile);
  strcpy(lbuffer+len,".gpl");  /* for output file name.gpl */
  output=fopen(lbuffer,"w");
  printf("-----------------------------------------------\n");
  printf("to plot results type: gnuplot %s\n",lbuffer);
  printf("-----------------------------------------------\n");

  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"set nokey\n");
  fprintf(output,"set data style lines\n");
  fprintf(output,"set size ratio -1\n");
  fprintf(output,"set size square\n");
  fprintf(output,"set xrange %s\n",xrange);
  fprintf(output,"set yrange %s\n",yrange);
  fprintf(output,"set title '%s'\n",title);
  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"# plot to screen\n");
  fprintf(output,"set multiplot\n");
  fprintf(output,"plot '%s'\n",datafile);
  fprintf(output,"plot 'catchment.out'\n");
  fprintf(output,"set nomultiplot\n");
  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"# plot to postscript file\n");
  strcpy(lbuffer+len,".ps");           /* for output file name.ps */
  fprintf(output,"set terminal postscript\n");
  fprintf(output,"set output '%s'\n",lbuffer);
  fprintf(output,"set multiplot\n");
  fprintf(output,"plot '%s'\n",datafile);
  fprintf(output,"plot 'catchment.out'\n");
  fprintf(output,"set nomultiplot\n");
  fprintf(output,"set output\n");
  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"# plot to pbm file\n");
  strcpy(lbuffer+len,".pbm");           /* for output file name.pbm */
  fprintf(output,"set bmargin 0\n");
  fprintf(output,"set lmargin 0\n");
  fprintf(output,"set rmargin 0\n");
  fprintf(output,"set tmargin 0\n");
  fprintf(output,"set size 0.8,1.06666666\n");
  fprintf(output,"set terminal pbm\n");
  fprintf(output,"set output '%s'\n",lbuffer);
  fprintf(output,"set multiplot\n");
  fprintf(output,"plot '%s'\n",datafile);
  fprintf(output,"plot 'catchment.out'\n");
  fprintf(output,"set nomultiplot\n");
  fprintf(output,"set output\n");
  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"# wait for interactive user\n");
  fprintf(output,"if(pi>3) pause -1\n");
  fprintf(output,"if(pi<3) pause 5\n");
  fprintf(output,"#-----------------------------------------------\n");

  fclose(output);
}
/*----------------------------------------------------------------------------------*/
void make_gpl2_file(datafile,title,xlabel,ylabel)
     char *datafile, *title, *xlabel, *ylabel;
{
  FILE *output;
  int len;

  len=strcspn(datafile,".");
  strcpy(lbuffer,datafile);
  strcpy(lbuffer+len,".gpl");  /* for output file name.gpl */
  output=fopen(lbuffer,"w");
  printf("-----------------------------------------------\n");
  printf("to plot results type: gnuplot %s\n",lbuffer);
  printf("-----------------------------------------------\n");

  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"set nokey\n");
  fprintf(output,"set data style lines\n");
  fprintf(output,"set xlabel '%s'\n",xlabel);
  fprintf(output,"set ylabel '%s'\n",ylabel);
  fprintf(output,"set title '%s'\n",title);
  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"# plot to screen\n");
  fprintf(output,"plot '%s'\n",datafile);
  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"# plot to postscript file\n");
  strcpy(lbuffer+len,".ps");           /* for output file name.ps */
  fprintf(output,"set terminal postscript\n");
  fprintf(output,"set output '%s'\n",lbuffer);
  fprintf(output,"plot '%s'\n",datafile);
  fprintf(output,"set output\n");
  fprintf(output,"#-----------------------------------------------\n");
  fprintf(output,"# wait for interactive user\n");
  fprintf(output,"if(pi>3) pause -1\n");
  fprintf(output,"if(pi<3) pause 5\n");
  fprintf(output,"#-----------------------------------------------\n");

  fclose(output);
}
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

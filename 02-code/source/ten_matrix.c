/*----------------------------------------------------------------------------------*/
/*--------------------------------- ten_matrix.c ------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on matrix structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h> 
#include <stdio.h> 
/*----------------------------------------------------------------------------------*/
#include "boundary_types.h" 
#include "matrix_types.h" 
#include "ten_matrix_types.h" 

#include "matrix.h" 

#include "ten_matrix.h" 
/*----------------------------------------------------------------------------------*/
/*------------------------ tensor matrices -----------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* create a matrix */
/*----------------------------------------------------------------------------------*/
ten_matrix *create_ten_matrix(rows,columns) 
     int rows, columns; 
{ 
 ten_matrix *x; 
 tensor *data; 
 
 x=(ten_matrix *)malloc(sizeof(ten_matrix)); 
 if(x==(ten_matrix *)NULL) 
   { 
 printf("error allocating memory for matrix\n"); 
 exit(0); 
   } 
 data=(tensor *)malloc(rows*columns*sizeof(tensor)); 
 if(data==(tensor *)NULL) 
   { 
 printf("error allocating memory for matrix\n"); 
 exit(0); 
   } 
 x->transpose=0; 
 x->invert=0; 
 x->rows=rows; 
 x->columns=columns; 
 x->value=data; 
 
 return(x); 
} 

/*----------------------------------------------------------------------------------*/
/* attach a matrix to pre-defined data */
/*----------------------------------------------------------------------------------*/
void attach_ten_matrix(x,rows,columns,data) 
     ten_matrix *x; 
     int rows, columns; 
     tensor *data; 
{ 
 x->transpose=0; 
 x->invert=0; 
 x->rows=rows; 
 x->columns=columns; 
 x->value=data; 
} 

/*----------------------------------------------------------------------------------*/
/* destroy a matrix */
/*----------------------------------------------------------------------------------*/
ten_matrix *destroy_ten_matrix(x) 
     ten_matrix *x; 
{ 
  if(x!=(ten_matrix *)NULL)
    {
      if(x->value!=(tensor *)NULL) free((void *)x->value); 
      free((void *)x); 
    }
  return((void *)NULL); 
} 

/*----------------------------------------------------------------------------------*/
/* check index to see if it is out of bounds */
/*----------------------------------------------------------------------------------*/
void check_ten_matrix_index(x,i,j) 
     ten_matrix *x; 
     int i,j; 
{ 
  int rows,cols; 
  
  rows=get_ten_num_rows(x); 
  cols=get_ten_num_columns(x); 
  
  if(0>i||i>=rows) { 
    printf("first array index (%d) out of bounds (0,%d)\n",i,rows-1); 
    exit(0); } 
  if(0>j||j>=cols) { 
    printf("second array index (%d) out of bounds (0,%d)\n",i,cols-1); 
    exit(0); } 
} 

/*----------------------------------------------------------------------------------*/
/* check matrix to see if need to invert */
/*----------------------------------------------------------------------------------*/
void check_ten_invert(x) 
     ten_matrix *x; 
{ 
  if(x->invert!=0) { 
    printf("you need to invert the matrix first\n"); 
    exit(0); } 
} 
/*----------------------------------------------------------------------------------*/
/* check to see if matrix has memory allocated yet */
/*----------------------------------------------------------------------------------*/
void check_ten_memory(x) 
     ten_matrix *x; 
{ 
  if(x->value==(tensor *)NULL) 
    { 
      printf("you have forgotten to provide memory for the result\n"); 
      exit(0); 
    } 
} 
/*----------------------------------------------------------------------------------*/
/*   check to see if matrices have shapes which can be multiplied */
/*----------------------------------------------------------------------------------*/
void check_ten_multiply_shape(a,b) 
     ten_matrix *a; 
     matrix *b; 
{ 
  int row_a,row_b,col_a,col_b; 
  
  row_a=get_ten_num_rows(a); col_a=get_ten_num_columns(a); 
  row_b=get_num_rows(b); col_b=get_num_columns(b); 
  if(col_a!=row_b) 
    { 
      printf("cannot multiply matrix shape (%dx%d) by matrix shape (%dx%d)\n", 
	     row_a,col_a,row_b,col_b); 
      exit(0); 
    } 
} 
/*----------------------------------------------------------------------------------*/
/* check to see if output matrix has size O.K. for product */
/*----------------------------------------------------------------------------------*/
void check_ten_multiply_size(a,b,x) 
     matrix *b; 
     ten_matrix *a, *x; 
{ 
  int row_a,col_b,row_x,col_x; 
  
  row_a=get_ten_num_rows(a); col_b=get_num_columns(b); 
  row_x=get_ten_num_rows(x); col_x=get_ten_num_columns(x); 
  if(row_a*col_b!=row_x*col_x) 
    { 
      printf("cannot put matrix product (%dx%d) into matrix shape (%dx%d)\n", 
	     row_a,col_b,row_x,col_x); 
      exit(0); 
    } 
  if(a==x) 
    { 
      printf("matrix for result must be different from input\n"); 
      exit(0); 
    } 
} 
/*----------------------------------------------------------------------------------*/
/*   get number of columns */
/*----------------------------------------------------------------------------------*/
int get_ten_num_columns(x) 
     ten_matrix *x; 
{ 
  int num_columns; 
  if(x->transpose==0) { num_columns=x->columns; } 
  else                { num_columns=x->rows;    } 
  return(num_columns); 
} 
/*----------------------------------------------------------------------------------*/
/* get number of rows */
/*----------------------------------------------------------------------------------*/
int get_ten_num_rows(x) 
     ten_matrix *x; 
{ 
  int num_rows; 
  if(x->transpose==0) { num_rows=x->rows;    } 
  else                { num_rows=x->columns; } 
  return(num_rows); 
} 
/*----------------------------------------------------------------------------------*/
/* get data area */
/*----------------------------------------------------------------------------------*/
tensor *startof_ten_matrix(x) 
     ten_matrix *x; 
{ 
  tensor *data; 
  
  data=x->value; 
  return(data); 
} 

/*----------------------------------------------------------------------------------*/
/*   put element in matrix */
/*----------------------------------------------------------------------------------*/
void put_ten_matrix_element(x,i,j,value) 
     ten_matrix *x; 
     int i,j; 
     tensor value; 
{ 
  check_ten_invert(x); 
  check_ten_matrix_index(x,i,j); 
  if(x->transpose==0) 
    { 
      x->value[j*x->rows+i][0][0]=value[0][0]; 
      x->value[j*x->rows+i][0][1]=value[0][1]; 
      x->value[j*x->rows+i][1][0]=value[1][0]; 
      x->value[j*x->rows+i][1][1]=value[1][1]; 
    } 
  else 
    { 
      x->value[i*x->rows+j][0][0]=value[0][0]; 
      x->value[i*x->rows+j][0][1]=value[0][1]; 
      x->value[i*x->rows+j][1][0]=value[1][0]; 
      x->value[i*x->rows+j][1][1]=value[1][1]; 
    } 
} 

/*----------------------------------------------------------------------------------*/
/* get element from matrix */
/*----------------------------------------------------------------------------------*/
void get_ten_matrix_element(x,i,j,value) 
     ten_matrix *x; 
     int i,j; 
     tensor value; 
{ 
  check_ten_invert(x); 
  check_ten_matrix_index(x,i,j); 
  if(x->transpose==0) 
    { 
      value[0][0]=x->value[j*x->rows+i][0][0]; 
      value[0][1]=x->value[j*x->rows+i][0][1]; 
      value[1][0]=x->value[j*x->rows+i][1][0]; 
      value[1][1]=x->value[j*x->rows+i][1][1]; 
    } 
  else 
    { 
      value[0][0]=x->value[i*x->rows+j][0][0]; 
      value[0][1]=x->value[i*x->rows+j][0][1]; 
      value[1][0]=x->value[i*x->rows+j][1][0]; 
      value[1][1]=x->value[i*x->rows+j][1][1]; 
    } 
} 

/*----------------------------------------------------------------------------------*/
/* put element in matrix of blocks */
/*----------------------------------------------------------------------------------*/
void put_block_ten_matrix_element(x,offset_i,offset_j,i,j,value) 
     ten_matrix *x; 
     int offset_i,offset_j,i,j; 
     tensor value; 
{ 
  i=offset_i+i; 
  j=offset_j+j; 
  check_ten_invert(x); 
  check_ten_matrix_index(x,i,j); 
  if(x->transpose==0) 
    { 
      x->value[j*x->rows+i][0][0]=value[0][0]; 
      x->value[j*x->rows+i][0][1]=value[0][1]; 
      x->value[j*x->rows+i][1][0]=value[1][0]; 
      x->value[j*x->rows+i][1][1]=value[1][1]; 
    } 
  else 
    { 
      x->value[i*x->rows+j][0][0]=value[0][0]; 
      x->value[i*x->rows+j][0][1]=value[0][1]; 
      x->value[i*x->rows+j][1][0]=value[1][0]; 
      x->value[i*x->rows+j][1][1]=value[1][1]; 
    } 
} 

/*----------------------------------------------------------------------------------*/
/* get element from matrix of blocks */
/*----------------------------------------------------------------------------------*/
void get_block_ten_matrix_element(x,offset_i,offset_j,i,j,value) 
     ten_matrix *x; 
     int offset_i,offset_j,i,j; 
     tensor value; 
{ 
  i=offset_i+i; 
  j=offset_j+j; 
  check_ten_invert(x); 
  check_ten_matrix_index(x,i,j); 
  if(x->transpose==0) 
    { 
      value[0][0]=x->value[j*x->rows+i][0][0]; 
      value[0][1]=x->value[j*x->rows+i][0][1]; 
      value[1][0]=x->value[j*x->rows+i][1][0]; 
      value[1][1]=x->value[j*x->rows+i][1][1]; 
    } 
  else 
    { 
      value[0][0]=x->value[i*x->rows+j][0][0]; 
      value[0][1]=x->value[i*x->rows+j][0][1]; 
      value[1][0]=x->value[i*x->rows+j][1][0]; 
      value[1][1]=x->value[i*x->rows+j][1][1]; 
    } 
} 
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void multiply_ten_matrix(a,b,x)
     ten_matrix *a,*x;
     matrix *b;
{
  int row_a;

  check_ten_multiply_shape(a,b);
  check_ten_multiply_size(a,b,x);
  check_ten_memory(x);
  
  if(a->invert==1) {printf("no routine to invert ten_matrix\n"); exit(0);}
  if(b->invert==1) invert_this_matrix(b);

  x->transpose=0;
  x->invert=0;
  x->rows=get_ten_num_rows(a);
  x->columns=get_num_columns(b);

  for(row_a=0;row_a<get_ten_num_rows(a);row_a++)
    {
     multiply_ten_element_row(a,b,row_a,x);
    }
}

/*----------------------------------------------------------------------------------*/
void multiply_ten_element_row(a,b,row_a,x)
     ten_matrix *a,*x;
     matrix *b;
     int row_a;
{
  int column_b;
  for(column_b=0;column_b<get_num_columns(b);column_b++)
    { 
      multiply_ten_element_column(a,b,row_a,column_b,x);
    }
}

/*----------------------------------------------------------------------------------*/
void multiply_ten_element_column(a,b,row_a,column_b,x)
     ten_matrix *a,*x;
     matrix *b;
     int row_a,column_b;
{   
  int i;
  tensor value_new,xy;
  double scalar;

  value_new[0][0]=0.0;  value_new[0][1]=0.0;
  value_new[1][0]=0.0;  value_new[1][1]=0.0;
  for(i=0;i<get_ten_num_columns(a);i++)
    {
      get_ten_matrix_element(a,row_a,i,xy);
      /*
      printf("\n----i=%d---",i);
      printf("\nxy[0][0]=%f,xy[0][1]=%f", xy[0][0] , xy[0][1]);
      printf("\nxy[1][0]=%f,xy[1][1]=%f", xy[1][0] , xy[1][1]);
      */
      scalar=get_matrix_element(b,i,column_b);

      value_new[0][0]=value_new[0][0]+xy[0][0]*scalar;
      value_new[0][1]=value_new[0][1]+xy[0][1]*scalar;
      value_new[1][0]=value_new[1][0]+xy[1][0]*scalar;
      value_new[1][1]=value_new[1][1]+xy[1][1]*scalar;
    }  

  /*
  printf("\n value_new[0][0]=%f,value_new[0][1]=%f", value_new[0][0] , value_new[0][1]);
  printf("\n value_new[1][0]=%f,value_new[1][1]=%f", value_new[1][0] , value_new[1][1]);
  */
  put_ten_matrix_element(x,row_a,column_b,value_new);
}
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

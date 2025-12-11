/*----------------------------------------------------------------------------------*/
/*--------------------------------- matrix.c ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on matrix structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include "sys/time.h"
#include "time.h"
#include <string.h>
#include <omp.h>

// #include "cblas.h"
// #include "lapacke.h"
#include <x86intrin.h>
/*----------------------------------------------------------------------------------*/
#include "matrix_types.h"

#include "file.h"

#include "matrix.h"

#include "matrix_inv.h"

#include "performance_summary.h"

/* External function declarations */
extern void dgemm_(char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
extern int mat_inv(double *A, unsigned n);
extern void update_multiply_matrix_stats(double duration, long long flops);
extern void get_memory_usage_kb(long *vmrss_kb, long *vmsize_kb);

/* External declarations for optimized matrix multiply */
extern void multiply_matrix_optimized(matrix *a, matrix *b, matrix *x);
extern void set_multiply_method(int method);
extern int get_multiply_method(void);

/*----------------------------------------------------------------------------------*/

/* Global variable to control inversion method */
static int use_sequential_inversion = 0; /* 0=Parallel, 1=Sequential */

/* Function to set inversion method */
void set_inversion_method(int method)
{
  use_sequential_inversion = method;
  if (method == 0)
  {
    printf("\n[CONFIG] Matrix inversion method: PARALLEL (LAPACK)\n");
  }
  else
  {
    printf("\n[CONFIG] Matrix inversion method: SEQUENTIAL (manual)\n");
  }
}

/* Function to get current inversion method */
int get_inversion_method(void)
{
  return use_sequential_inversion;
}

/* Helper function for memory tracking */
/*
void get_memory_usage_kb(long *vmrss_kb, long *vmsize_kb) {
    FILE* file = fopen("/proc/self/status", "r");
    char line[128];

    *vmrss_kb = 0;
    *vmsize_kb = 0;

    if (file) {
        while (fgets(line, 128, file) != NULL) {
            if (strncmp(line, "VmRSS:", 6) == 0) {
                sscanf(line + 6, "%ld", vmrss_kb);
            }
            if (strncmp(line, "VmSize:", 7) == 0) {
                sscanf(line + 7, "%ld", vmsize_kb);
            }
        }
        fclose(file);
    }
}
*/

/*----------------------------------------------------------------------------------*/
#define CHECK_MATRIX 0
/*----------------------------------------------------------------------------------*/
/* create a matrix */
/*----------------------------------------------------------------------------------*/
matrix *create_matrix(rows, columns)
int rows, columns;
{
  matrix *x;
  double *data;

  x = (matrix *)malloc(sizeof(matrix));
  if (x == (matrix *)NULL)
  {
    printf("error allocating memory for matrix\n");
    exit(0);
  }
  data = (double *)malloc(rows * columns * sizeof(double));
  if (data == (double *)NULL)
  {
    printf("error allocating memory for matrix\n");
    exit(0);
  }
  x->transpose = 0;
  x->invert = 0;
  x->rows = rows;
  x->columns = columns;
  x->value = data;

  return (x);
}

/*----------------------------------------------------------------------------------*/
/* attach a matrix to predefined memory */
/*----------------------------------------------------------------------------------*/
void attach_matrix(x, rows, columns, data)
    matrix *x;
int rows, columns;
double *data;
{
  x->transpose = 0;
  x->invert = 0;
  x->rows = rows;
  x->columns = columns;
  x->value = data;
}

/*----------------------------------------------------------------------------------*/
/* destroy a matrix */
/*----------------------------------------------------------------------------------*/
matrix *destroy_matrix(x)
matrix *x;
{
  if (x != (matrix *)NULL)
  {
    if (x->value != (double *)NULL)
      free((void *)x->value);
    free((void *)x);
  }
  return ((void *)NULL);
}

/*----------------------------------------------------------------------------------*/
/* check index to see if it is out of bounds */
/*----------------------------------------------------------------------------------*/
void check_matrix_index(x, i, j)
    matrix *x;
int i, j;
{
  int rows, cols;

  rows = get_num_rows(x);
  cols = get_num_columns(x);

  if (0 > i || i >= rows)
  {
    printf("first array index (%d) out of bounds (0,%d)\n", i, rows - 1);
    exit(0);
  }
  if (0 > j || j >= cols)
  {
    printf("second array index (%d) out of bounds (0,%d)\n", j, cols - 1);
    exit(0);
  }
}

/*----------------------------------------------------------------------------------*/
/* check matrix to see if need to invert */
/*----------------------------------------------------------------------------------*/
void check_invert(x)
    matrix *x;
{
  if (x->invert != 0)
  {
    printf("you need to invert the matrix first\n");
    exit(0);
  }
}

/*----------------------------------------------------------------------------------*/
/* check to see if two matrices are same size */
/*----------------------------------------------------------------------------------*/
void check_size(a, b)
    matrix *a,
    *b;
{
  int size_a, size_b;

  size_a = a->rows * a->columns;
  size_b = b->rows * b->columns;
  if (size_a != size_b)
  {
    printf("cannot put matrix size (%d=%dx%d) into matrix size (%d=%dx%d)\n",
           size_a, get_num_rows(a), get_num_columns(a), size_b, get_num_rows(b),
           get_num_columns(b));
    exit(0);
  }
}

/*----------------------------------------------------------------------------------*/
/* check to see if matrix and row vector are same size */
/*----------------------------------------------------------------------------------*/
void check_size_row(a, b)
    matrix *a,
    *b;
{
  int size_a, size_b;

  size_a = get_num_columns(a);
  size_b = get_num_columns(b);
  if (size_a != size_b || get_num_rows(b) != 1)
  {
    printf("cannot put matrix size (%dx%d) into last row of matrix size (%dx%d)\n",
           get_num_rows(b), size_b, get_num_rows(a), size_a);
    exit(0);
  }
}

/*----------------------------------------------------------------------------------*/
/* check to see if second matrix has size of first minus one row */
/*----------------------------------------------------------------------------------*/
void check_collapse_size(a, b)
    matrix *a,
    *b;
{
  int size_a, size_b;

  size_a = a->rows * a->columns;
  size_b = b->rows * b->columns;
  if (size_a - get_num_columns(a) != size_b)
  {
    printf("cannot collapse matrix size (%d=%dx%d) into matrix size (%d=%dx%d)\n",
           size_a, get_num_rows(a), get_num_columns(a), size_b, get_num_rows(b),
           get_num_columns(b));
    exit(0);
  }
}

/*----------------------------------------------------------------------------------*/
/* check to see if matrix has memory allocated yet */
/*----------------------------------------------------------------------------------*/
void check_memory(x)
    matrix *x;
{
  if (x->value == (double *)NULL)
  {
    printf("you have forgotten to provide memory for the result\n");
    exit(0);
  }
}

/*----------------------------------------------------------------------------------*/
/* check to see if matrices have shapes which can be added */
/*----------------------------------------------------------------------------------*/
void check_add_shape(a, b)
    matrix *a,
    *b;
{
  int row_a, row_b, col_a, col_b;

  row_a = get_num_rows(a);
  col_a = get_num_columns(a);
  row_b = get_num_rows(b);
  col_b = get_num_columns(b);
  if (row_a != row_b || col_a != col_b)
  {
    printf("cannot add matrix shape (%dx%d) to matrix shape (%dx%d)\n",
           row_a, col_a, row_b, col_b);
    exit(0);
  }
}
/*----------------------------------------------------------------------------------*/
/*   check to see if matrices have shapes which can be multiplied */
/*----------------------------------------------------------------------------------*/
void check_multiply_shape(a, b)
    matrix *a,
    *b;
{
  int row_a, row_b, col_a, col_b;

  row_a = get_num_rows(a);
  col_a = get_num_columns(a);
  row_b = get_num_rows(b);
  col_b = get_num_columns(b);
  if (col_a != row_b)
  {
    printf("cannot multiply matrix shape (%dx%d) by matrix shape (%dx%d)\n",
           row_a, col_a, row_b, col_b);
    exit(0);
  }
}
/*----------------------------------------------------------------------------------*/
/* check to see if output matrix has size O.K. for product */
/*----------------------------------------------------------------------------------*/
void check_multiply_size(a, b, x)
    matrix *a,
    *b, *x;
{
  int row_a, col_b, row_x, col_x;

  row_a = get_num_rows(a);
  col_b = get_num_columns(b);
  row_x = get_num_rows(x);
  col_x = get_num_columns(x);
  if (row_a * col_b != row_x * col_x)
  {
    printf("cannot put matrix product (%dx%d) into matrix shape (%dx%d)\n",
           row_a, col_b, row_x, col_x);
    exit(0);
  }
  if (a == x || b == x)
  {
    printf("matrix for result must be different from input\n");
    exit(0);
  }
}
/*----------------------------------------------------------------------------------*/
/* check to see if matrix has square shape */
/*----------------------------------------------------------------------------------*/
void check_invert_shape(a)
    matrix *a;
{
  int row_a, col_a;

  row_a = get_num_rows(a);
  col_a = get_num_columns(a);
  if (row_a != col_a)
  {
    printf("cannot invert matrix shape (%dx%d)\n", row_a, col_a);
    exit(0);
  }
}
/*----------------------------------------------------------------------------------*/
/*   get number of columns */
/*----------------------------------------------------------------------------------*/
int get_num_columns(x)
matrix *x;
{
  int num_columns;
  if (x->transpose == 0)
  {
    num_columns = x->columns;
  }
  else
  {
    num_columns = x->rows;
  }
  return (num_columns);
}

/*----------------------------------------------------------------------------------*/
/* get number of rows */
/*----------------------------------------------------------------------------------*/
int get_num_rows(x)
matrix *x;
{
  int num_rows;
  if (x->transpose == 0)
  {
    num_rows = x->rows;
  }
  else
  {
    num_rows = x->columns;
  }
  return (num_rows);
}

/*----------------------------------------------------------------------------------*/
/* get data area */
/*----------------------------------------------------------------------------------*/
double *startof_matrix(x)
matrix *x;
{
  double *data;

  data = x->value;
  return (data);
}

/*----------------------------------------------------------------------------------*/
/* get address beyond data area */
/*----------------------------------------------------------------------------------*/
double *after_matrix(x)
matrix *x;
{
  double *data;

  data = x->value + (x->rows * x->columns);
  return (data);
}

/*----------------------------------------------------------------------------------*/
/*   put element in matrix */
/*----------------------------------------------------------------------------------*/
void put_matrix_element(x, i, j, value)
    matrix *x;
int i, j;
double value;
{
#if CHECK_MATRIX
  check_invert(x);
  check_matrix_index(x, i, j);
#endif
  if (x->transpose == 0)
  {
    x->value[j * x->rows + i] = value;
  }
  else
  {
    x->value[i * x->rows + j] = value;
  }
}

/*----------------------------------------------------------------------------------*/
/* get element from matrix */
/*----------------------------------------------------------------------------------*/
double get_matrix_element(x, i, j)
matrix *x;
int i, j;
{
  double value;

#if CHECK_MATRIX
  check_invert(x);
  check_matrix_index(x, i, j);
#endif
  if (x->transpose == 0)
  {
    value = x->value[j * x->rows + i];
  }
  else
  {
    value = x->value[i * x->rows + j];
  }
  return (value);
}

/*----------------------------------------------------------------------------------*/
/* put element in matrix of blocks */
/*----------------------------------------------------------------------------------*/
void put_block_matrix_element(x, offset_i, offset_j, i, j, value)
    matrix *x;
int offset_i, offset_j, i, j;
double value;
{
  i = offset_i + i;
  j = offset_j + j;
#if CHECK_MATRIX
  check_invert(x);
  check_matrix_index(x, i, j);
#endif
  if (x->transpose == 0)
  {
    x->value[j * x->rows + i] = value;
  }
  else
  {
    x->value[i * x->rows + j] = value;
  }
}

/*----------------------------------------------------------------------------------*/
/* get element from matrix of blocks */
/*----------------------------------------------------------------------------------*/
double get_block_matrix_element(x, offset_i, offset_j, i, j)
matrix *x;
int offset_i, offset_j, i, j;
{
  double value;

  i = offset_i + i;
  j = offset_j + j;
#if CHECK_MATRIX
  check_invert(x);
  check_matrix_index(x, i, j);
#endif
  if (x->transpose == 0)
  {
    value = x->value[j * x->rows + i];
  }
  else
  {
    value = x->value[i * x->rows + j];
  }
  return (value);
}

/*----------------------------------------------------------------------------------*/
/* transpose a matrix */
/*----------------------------------------------------------------------------------*/
void transpose_matrix(a, x)
    matrix *a,
    *x;
{
  check_size(a, x);
  check_memory(x);
  a->transpose = 1 - a->transpose;
  if (a != x)
  {
    copy_matrix(a, x);
    a->transpose = 1 - a->transpose;
  }
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* add two matrices */
/*----------------------------------------------------------------------------------*/
void add_matrix(a, b, x)
    matrix *a,
    *b, *x;
{
  double value_a, value_b;
  int rows, columns;
  int i, j;

  check_add_shape(a, b);
  check_size(a, x);
  check_memory(x);

  if (a->invert == 1)
    invert_this_matrix(a);
  if (b->invert == 1)
    invert_this_matrix(b);

  x->transpose = 0;
  if (a == x)
    x->transpose = a->transpose;
  if (b == x)
    x->transpose = b->transpose;
  x->invert = 0;

  rows = get_num_rows(a);
  columns = get_num_columns(a);
  if (x->transpose == 0)
  {
    x->rows = rows;
    x->columns = columns;
  }
  else
  {
    x->rows = columns;
    x->columns = rows;
  }

  for (j = 0; j < columns; j++) /* different columns of matrix A */
  {
    for (i = 0; i < rows; i++) /*different rows of matrix A */
    {
      value_a = get_matrix_element(a, i, j);
      value_b = get_matrix_element(b, i, j);
      put_matrix_element(x, i, j, value_a + value_b);
    }
  }
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* subtract second matrix from first matrix */
/*----------------------------------------------------------------------------------*/
void subtract_matrix(a, b, x)
    matrix *a,
    *b, *x;
{
  double value_a, value_b;
  int rows, columns;
  int i, j;

  check_add_shape(a, b);
  check_size(a, x);
  check_memory(x);

  if (a->invert == 1)
    invert_this_matrix(a);
  if (b->invert == 1)
    invert_this_matrix(b);

  x->transpose = 0;
  if (a == x)
    x->transpose = a->transpose;
  if (b == x)
    x->transpose = b->transpose;
  x->invert = 0;

  rows = get_num_rows(a);
  columns = get_num_columns(a);
  if (x->transpose == 0)
  {
    x->rows = rows;
    x->columns = columns;
  }
  else
  {
    x->rows = columns;
    x->columns = rows;
  }

  for (j = 0; j < columns; j++) /* different columns of matrix A */
  {
    for (i = 0; i < rows; i++) /*different rows of matrix A */
    {
      value_a = get_matrix_element(a, i, j);
      value_b = get_matrix_element(b, i, j);
      put_matrix_element(x, i, j, value_a - value_b);
    }
  }
}

/*----------------------------------------------------------------------------------*/
/* scale a matrix by a constant value */
/*----------------------------------------------------------------------------------*/
void scale_matrix(s, a, x) double s;
matrix *a, *x;
{
  double *data;
  int i;

  check_size(a, x);
  check_memory(x);

  if (a == x)
  {
    if (a->invert == 1)
      s = 1.0 / s;
    data = startof_matrix(a);
  }
  else
  {
    copy_matrix(a, x);
    data = startof_matrix(x);
  }

  for (i = 0; i < a->rows * a->columns; i++) /* different rows and columns of matrix A */
  {
    data[i] = data[i] * s;
  }
}

/*----------------------------------------------------------------------------------*/
/* put zero in last row of matrix */
/*----------------------------------------------------------------------------------*/
void zero_last_matrix_row(a)
    matrix *a;
{
  int i, j, rows, columns;

  check_memory(a);
  rows = get_num_rows(a);
  columns = get_num_columns(a);

  i = rows - 1;
  for (j = 0; j < columns; j++) /* different columns of matrix A */
  {
    put_matrix_element(a, i, j, 0.0);
  }
}

/*----------------------------------------------------------------------------------*/
/* put values in last row of matrix */
/*----------------------------------------------------------------------------------*/
void fill_last_matrix_row(a, x)
    matrix *a,
    *x;
{
  int i, j, rows, columns;
  double value;

  check_size_row(a, x);
  check_memory(a);
  check_memory(x);
  rows = get_num_rows(a);
  columns = get_num_columns(a);

  i = rows - 1;
  for (j = 0; j < columns; j++) /* different columns of matrix A */
  {
    value = get_matrix_element(x, 0, j);
    put_matrix_element(a, i, j, value);
  }
}

/*----------------------------------------------------------------------------------*/
/* copy a matrix into another matrix */
/*----------------------------------------------------------------------------------*/
void copy_matrix(a, x)
    matrix *a,
    *x;
{
  double value_a;
  int rows, columns;
  int i, j;

  check_size(a, x);
  check_memory(x);

  if (a->invert == 1)
    invert_this_matrix(a);

  if (a != x)
  {
    x->transpose = 0;
    x->invert = 0;
    rows = get_num_rows(a);
    x->rows = rows;
    columns = get_num_columns(a);
    x->columns = columns;

    for (j = 0; j < columns; j++) /* different columns of matrix A */
    {
      for (i = 0; i < rows; i++) /* different rows of matrix A */
      {
        value_a = get_matrix_element(a, i, j);
        put_matrix_element(x, i, j, value_a);
      }
    }
  }
}

/*----------------------------------------------------------------------------------*/

void transpose(double *A, double *B, int n)
{
  int i, j;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      B[j * n + i] = A[i * n + j];
    }
  }
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

// Function to multiply two matrices in squential.
void matmul(a, b, x) //(pa, pb, pmul) --> (a, b, x)
    matrix *a,
    *b, *x;
{
  int row_a;
  int column_b;
  int i;
  double value, value1, value2, value0, value_new, value_new2;
  double cofactor;

  int a_col_num = get_num_columns(a);
  int b_col_num = get_num_columns(b);

  int a_row_num = get_num_rows(a);
  int b_row_num = get_num_rows(b);

  check_multiply_shape(a, b);
  check_multiply_size(a, b, x);
  check_memory(x);

  if (a->invert == 1)
  {
    printf("\nmultiply_matrix() with (a->invert==1)\n");
    invert_this_matrix(a);
  }
  if (b->invert == 1)
  {
    printf("\nmultiply_matrix() with (b->invert==1)\n");
    invert_this_matrix(b);
  }

  x->transpose = 0;
  x->invert = 0;
  x->rows = a_row_num;
  x->columns = b_col_num;

  printf("\nPlease wait... Processing two matrices multiplication in squential.....\n");

  if (x->rows > 1 && x->columns > 1)
  {
    printf("\n------- Multiply1  with (x->rows > 1 && x->columns > 1) ------- \n");
    printf("1--- multiply matrix A ; result=(%dx%d)\n", a->rows, a->columns);
    printf("1--- multiply matrix B ; result=(%dx%d)\n", b->rows, b->columns);
    printf("1--- multiply matrix X ; result=(%dx%d)\n", x->rows, x->columns);

    // if ((a->rows == b->rows) && (a->columns == b->columns))
    if ((a->rows == a->columns) && (b->rows == b->columns))
    {
      printf("\n@@@@@@@@@@@@@ multiply square matrix A*B ; result=(%dx%d)\n\n", a->rows, a->columns);
    }

    if (a->transpose == 0)
    {
      printf("\n------- (a->transpose == 0) ------- \n");

      /*
          for (int i = 0; i < n * n; i++)
          {
              cofactor = 0.0;
              for (int j = 0; j < n; j++)
              {
                  cofactor = cofactor + *(pa + (i / n) * n + j) * (*(pb + (i % n) + n * j));

              }
              *(pmul + i) = cofactor;
          }
      */

      for (int row_a = 0; row_a < a_row_num; row_a++)
      {
        for (int column_b = 0; column_b < b_col_num; column_b++)
        {
          value_new = 0.0;

          for (int i = 0; i < a_col_num; i++)
          {
            value_new = value_new + (a->value[i * a->rows + row_a] * b->value[column_b * b_row_num + i]);
          }

          x->value[column_b * x->rows + row_a] = value_new;
        }
      }
    }
    else
    {
      for (int row_a = 0; row_a < a_row_num; row_a++)
      {
        for (int column_b = 0; column_b < b_col_num; column_b++)
        {
          value_new = 0.0;
          value_new2 = 0.0;

          for (int i = 0; i < a_col_num; i++)
          {
            value_new = value_new + (a->value[row_a * a->rows + i] * b->value[column_b * b_row_num + i]);
            // value_new2 = value_new2 + (a->value[row_a * a->rows + (i+1)] * b->value[column_b * b_row_num + (i+1)]);
          }

          x->value[column_b * x->rows + row_a] = value_new;
          // x->value[column_b * x->rows + (row_a+1)] = value_new2;
        }
      }
    }
  }
  else
  {
    /*
        printf("\n------- Multiply2  ------- \n");
        printf("2--- multiply matrix A ; result=(%dx%d)\n", a->rows, a->columns);
        printf("2--- multiply matrix B ; result=(%dx%d)\n", b->rows, b->columns);
        printf("2--- multiply matrix X ; result=(%dx%d)\n", x->rows, x->columns);
   */

    if (a->transpose == 0)
    {

      for (int row_a = 0; row_a < a_row_num; row_a++)
      {

        for (int column_b = 0; column_b < b_col_num; column_b++)
        {

          value_new = 0.0;

          for (int i = 0; i < a_col_num; i++)
          {
            value_new = value_new + (a->value[i * a->rows + row_a] * b->value[column_b * b_row_num + i]);
          }

          x->value[column_b * x->rows + row_a] = value_new;
        }
      }
    }
    else
    {

      for (int row_a = 0; row_a < a_row_num; row_a++)
      {

        for (int column_b = 0; column_b < b_col_num; column_b++)
        {

          value_new = 0.0;

          for (int i = 0; i < a_col_num; i++)
          {
            value_new = value_new + (a->value[row_a * a->rows + i] * b->value[column_b * b_row_num + i]);
          }

          x->value[column_b * x->rows + row_a] = value_new;
        }
      }
    }
  }

  printf("\n");
}

//----------------------------------------------------------------------------------*/
void multiply_matrix(a, b, x)
    matrix *a,
    *b, *x;
{
  int row_a;
  int column_b;
  int i;
  double value, value1, value2, value0, value_new, value_new2;

  int a_col_num = get_num_columns(a);
  int b_col_num = get_num_columns(b);
  int a_row_num = get_num_rows(a);
  int b_row_num = get_num_rows(b);

  check_multiply_shape(a, b);
  check_multiply_size(a, b, x);
  check_memory(x);

//double inv_start_time = omp_get_wtime();

  // Handle invert flags
  if (a->invert == 1)
  {
    printf("\nmultiply_matrix() with (a->invert==1)\n");
    invert_this_matrix(a);
  }
  if (b->invert == 1)
  {
    printf("\nmultiply_matrix() with (b->invert==1)\n");
    invert_this_matrix(b);
  }

//double inv_time = omp_get_wtime() - inv_start_time;
//update_inversion_time(inv_time, a->rows);

  // Set output matrix properties
  x->transpose = 0;
  x->invert = 0;
  x->rows = a_row_num;
  x->columns = b_col_num;

  // Use optimized multiplication for large matrices
  if (x->rows > 1 && x->columns > 1)
  {
    printf("\n------- Using Optimized Multiply (Method %d) -------\n",
           get_multiply_method());
    printf("Matrix A: (%dx%d), transpose=%d\n", a->rows, a->columns, a->transpose);
    printf("Matrix B: (%dx%d), transpose=%d\n", b->rows, b->columns, b->transpose);
    printf("Matrix X: (%dx%d)\n", x->rows, x->columns);

    if ((a->rows == a->columns) && (b->rows == b->columns))
    {
      printf("Square matrix multiplication\n");
    }

    // Before multiplication
    double start_time = omp_get_wtime();

    // Call optimized multiplication
    multiply_matrix_optimized(a, b, x);

    // After multiplication
    double mult_time = omp_get_wtime() - start_time;
    update_multiply_time(mult_time, a->rows, b->columns, a->columns);

  }
  else
  {
    // Small matrix - use original sequential code
    printf("\n------- Using Sequential (small matrix) -------\n");

    if (a->transpose == 0)
    {
      for (row_a = 0; row_a < a_row_num; row_a++)
      {
        for (column_b = 0; column_b < b_col_num; column_b++)
        {
          value_new = 0.0;

          for (i = 0; i < a_col_num; i++)
          {
            value_new += a->value[i * a->rows + row_a] *
                         b->value[column_b * b_row_num + i];
          }

          x->value[column_b * x->rows + row_a] = value_new;
        }
      }
    }
    else
    {
      for (row_a = 0; row_a < a_row_num; row_a++)
      {
        for (column_b = 0; column_b < b_col_num; column_b++)
        {
          value_new = 0.0;

          for (i = 0; i < a_col_num; i++)
          {
            value_new += a->value[row_a * a->rows + i] *
                         b->value[column_b * b_row_num + i];
          }

          x->value[column_b * x->rows + row_a] = value_new;
        }
      }
    }
  }
}

// void multiply_matrix_before_optimize(a, b, x)
void multiply_matrix_org(a, b, x)
    matrix *a,
    *b, *x;
{
  int row_a;
  int column_b;
  int i;
  double value, value1, value2, value0, value_new, value_new2;

  int a_col_num = get_num_columns(a);
  int b_col_num = get_num_columns(b);

  int a_row_num = get_num_rows(a);
  int b_row_num = get_num_rows(b);

  check_multiply_shape(a, b);
  check_multiply_size(a, b, x);
  check_memory(x);

  // printf("\nmultiply matrix size=(%dx%dx%d)\n", a->rows, a->columns, b->columns);

  if (a->invert == 1)
  {
    printf("\nmultiply_matrix() with (a->invert==1)\n");
    invert_this_matrix(a);
  }
  if (b->invert == 1)
  {
    printf("\nmultiply_matrix() with (b->invert==1)\n");
    invert_this_matrix(b);
  }

  x->transpose = 0;
  x->invert = 0;
  x->rows = a_row_num;
  x->columns = b_col_num;

  // printf("\nmultiplying matrix; result=(%dx%d)\n",x->rows,x->columns);

  if (x->rows > 1 && x->columns > 1)
  {
    printf("\n------- Multiply1  with (x->rows > 1 && x->columns > 1) ------- \n");
    printf("1--- multiply matrix A ; result=(%dx%d)\n", a->rows, a->columns);
    printf("1--- multiply matrix B ; result=(%dx%d)\n", b->rows, b->columns);
    printf("1--- multiply matrix X ; result=(%dx%d)\n", x->rows, x->columns);

    // if ((a->rows == b->rows) && (a->columns == b->columns))
    if ((a->rows == a->columns) && (b->rows == b->columns))
    {
      printf("\n@@@@@@@@@@@@@ multiply square matrix A*B ; result=(%dx%d)\n\n", a->rows, a->columns);
    }

    // int numberofthread = 4;

    // omp_set_num_threads(2);

    if (a->transpose == 0)
    {

      printf("\n------- (a->transpose == 0) ------- \n");

#pragma omp parallel for schedule(dynamic, 50) collapse(2) private(row_a, column_b, i, value_new, value_new2)

      for (int row_a = 0; row_a < a_row_num; row_a++)
      {
        for (int column_b = 0; column_b < b_col_num; column_b++)
        {
          value_new = 0.0;
          value_new2 = 0.0;

          for (int i = 0; i < a_col_num; i += 2)
          {

            value_new = value_new + (a->value[i * a->rows + row_a] * b->value[column_b * b_row_num + i]);
            value_new2 = value_new2 + (a->value[(i + 1) * a->rows + row_a] * b->value[column_b * b_row_num + (i + 1)]);
          }

          x->value[column_b * x->rows + row_a] = value_new;
          x->value[column_b * x->rows + (row_a + 1)] = value_new2;
        }
      }
    }
    else
    {
      printf("\n------- (a->transpose == 0) ------- \n");
#pragma omp parallel for schedule(dynamic, 50) collapse(2) private(row_a, column_b, i, value_new)

      for (int row_a = 0; row_a < a_row_num; row_a++)
      {
        for (int column_b = 0; column_b < b_col_num; column_b++)
        {
          value_new = 0.0;
          value_new2 = 0.0;

          for (int i = 0; i < a_col_num; i++)
          {
            value_new = value_new + (a->value[row_a * a->rows + i] * b->value[column_b * b_row_num + i]);
            // value_new2 = value_new2 + (a->value[row_a * a->rows + (i+1)] * b->value[column_b * b_row_num + (i+1)]);
          }

          x->value[column_b * x->rows + row_a] = value_new;
          // x->value[column_b * x->rows + (row_a+1)] = value_new2;
        }
      }
    }
  }
  else
  {
    /*
        printf("\n------- Multiply2  ------- \n");
        printf("2--- multiply matrix A ; result=(%dx%d)\n", a->rows, a->columns);
        printf("2--- multiply matrix B ; result=(%dx%d)\n", b->rows, b->columns);
        printf("2--- multiply matrix X ; result=(%dx%d)\n", x->rows, x->columns);
    */

    if (a->transpose == 0)
    {

      for (int row_a = 0; row_a < a_row_num; row_a++)
      {

        for (int column_b = 0; column_b < b_col_num; column_b++)
        {

          value_new = 0.0;

          for (int i = 0; i < a_col_num; i++)
          {
            value_new = value_new + (a->value[i * a->rows + row_a] * b->value[column_b * b_row_num + i]);
          }

          x->value[column_b * x->rows + row_a] = value_new;
        }
      }
    }
    else
    {

      for (int row_a = 0; row_a < a_row_num; row_a++)
      {

        for (int column_b = 0; column_b < b_col_num; column_b++)
        {

          value_new = 0.0;

          for (int i = 0; i < a_col_num; i++)
          {
            value_new = value_new + (a->value[row_a * a->rows + i] * b->value[column_b * b_row_num + i]);
          }

          x->value[column_b * x->rows + row_a] = value_new;
        }
      }
    }
  }
}

void multiply_matrix_(a, b, x) // Original matrix multiplication function
    matrix *a,
    *b, *x;
{
  int row_a;

  check_multiply_shape(a, b);
  check_multiply_size(a, b, x);
  check_memory(x);

  if (a->invert == 1)
  {
    printf("\nmultiply_matrix() with (a->invert==1)\n");
    invert_this_matrix(a);
  }
  if (b->invert == 1)
  {
    printf("\nmultiply_matrix() with (b->invert==1)\n");
    invert_this_matrix(b);
  }

  x->transpose = 0;
  x->invert = 0;
  x->rows = get_num_rows(a);
  x->columns = get_num_columns(b);

  if (x->rows > 1 && x->columns > 1)
  {
    printf("\n--- multiply matrix A ; result=(%dx%d)\n", a->rows, a->columns);
    printf("--- multiply matrix B ; result=(%dx%d)\n", b->rows, b->columns);
    printf("--- multiply matrix X ; result=(%dx%d)\n", x->rows, x->columns);

    // printf("\nmultiply matrix; result=(%dx%d)",x->rows,x->columns);
    for (row_a = 0; row_a < get_num_rows(a); row_a++)
    {
      // if(row_a%50==0) printf("\n");
      // printf("*"); fflush(stdout);
      multiply_element_row(a, b, row_a, x);
    }
    // printf("\n");
  }
  else
  {
    for (row_a = 0; row_a < get_num_rows(a); row_a++)
    {
      multiply_element_row(a, b, row_a, x);
    }
  }
}

/*----------------------------------------------------------------------------------*/
void multiply_element_row(a, b, row_a, x)
    matrix *a,
    *b, *x;
int row_a;
{
  int column_b;
  for (column_b = 0; column_b < get_num_columns(b); column_b++)
  {
    multiply_element_column(a, b, row_a, column_b, x);
  }
}

/*----------------------------------------------------------------------------------*/
void multiply_element_column(a, b, row_a, column_b, x)
    matrix *a,
    *b, *x;
int row_a, column_b;
{
  int i;
  double value, value_new;

  value_new = 0.0;
  for (i = 0; i < get_num_columns(a); i++)
  {
    value = get_matrix_element(a, row_a, i) * get_matrix_element(b, i, column_b);
    value_new = value_new + value;
  }
  put_matrix_element(x, row_a, column_b, value_new);
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void invert_matrix(a, x)
    matrix *a,
    *x;
{
  check_size(a, x);
  check_invert_shape(a);
  check_memory(x);

  a->invert = 1 - a->invert;
  if (a != x)
  {
    copy_matrix(a, x);
  }
}

/*----------------------------------------------------------------------------------*/
/*========inverse_matrix===============*/
/* COMPLETE FIXED VERSION - Replace lines 1152-1214 in matrix.c */
/*----------------------------------------------------------------------------------*/

void invert_this_matrix(a)
    matrix *a;
{
  int n, j;
  struct timeval inv_start, inv_finish;
  double inv_duration;
  long vmrss_before, vmsize_before, vmrss_after, vmsize_after;

  a->invert = 0;
  n = get_num_columns(a);

/*----------------------------------------------------------------------------------*/
  if (a->value == NULL) {
      printf("ERROR: invert_this_matrix: a->value is NULL (n=%d)\n", n);
      exit(1);
  }
/*----------------------------------------------------------------------------------*/

  /* Print method being used */
  if (use_sequential_inversion == 0)
  {
    printf("\n[PARALLEL INVERSION] Inverting (%dx%d) using LAPACK\n", n, n);
  }
  else
  {
    printf("\n[SEQUENTIAL INVERSION] Inverting (%dx%d) using Gauss-Jordan\n", n, n);
  }

  /* Get memory before */
  get_memory_usage_kb(&vmrss_before, &vmsize_before);
  printf("Memory before: VmRSS=%.2f MB\n", vmrss_before / 1024.0);

double inv_start_time = omp_get_wtime();

  printf("Inverting the matrix, Please wait..\n");
  gettimeofday(&inv_start, NULL);

  if (use_sequential_inversion == 0)
  {
    /*-------------- Parallel (LAPACK) -----------*/
    printf("Using LAPACK mat_inv() - parallel LU decomposition\n");
    mat_inv(a->value, n);
  }
  else
  {
    /*-------------- Sequential (Manual) -----------*/
    printf("Using manual Gauss-Jordan - sequential\n");
    for (j = 0; j < n; j++)
    {
      if (j % 50 == 0 && j > 0)
      {
        printf("Progress: %d/%d rows (%.1f%%)\n", j, n, (100.0 * j) / n);
        fflush(stdout);
      }
      scale_row(a, j);
      reduce_column(a, j);
    }
    printf("Progress: %d/%d rows (100.0%%)\n", n, n);
  }

  gettimeofday(&inv_finish, NULL);
  inv_duration = ((double)(inv_finish.tv_sec - inv_start.tv_sec) * 1000000 +
                  (double)(inv_finish.tv_usec - inv_start.tv_usec)) /
                 1000000;

double inv_time = omp_get_wtime() - inv_start_time;
update_inversion_time(inv_time, a->rows);

  /*----------------------------------------------------------------------------------*/
  /* ⭐⭐⭐ CRITICAL FIX - ADD THIS LINE ⭐⭐⭐ */
  /*----------------------------------------------------------------------------------*/
  /* Update global statistics so bsolve.c and summary can see the timing */
  update_matrix_inversion_stats(inv_duration);
  /*----------------------------------------------------------------------------------*/

  /* Get memory after */
  get_memory_usage_kb(&vmrss_after, &vmsize_after);

  /* Calculate GFLOPS */
  double flops = (2.0 / 3.0) * n * n * n + n * n * n;
  double gflops = flops / inv_duration / 1.0e9;

  printf("\n========== MATRIX INVERSION COMPLETE ==========\n");
  printf("Method:     %s\n", use_sequential_inversion ? "SEQUENTIAL" : "PARALLEL");
  printf("Matrix:     %d x %d\n", n, n);
  printf("Time:       %.6f seconds\n", inv_duration);
  printf("GFLOPS:     %.2f\n", gflops);
  printf("Memory:     VmRSS=%.2f MB (delta: %.2f MB)\n",
         vmrss_after / 1024.0, (vmrss_after - vmrss_before) / 1024.0);
  printf("==============================================\n\n");
}

/*----------------------------------------------------------------------------------*/
/* END OF FIXED FUNCTION */
/*----------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------*/
/*                           WHAT CHANGED                                          */
/*----------------------------------------------------------------------------------*/
/*
 * ADDED ONE LINE (marked with ⭐ above):
 *
 *   update_matrix_inversion_stats(inv_duration);
 *
 * This line is placed RIGHT AFTER calculating inv_duration (line after gettimeofday)
 * and BEFORE getting memory after.
 *
 * This single line passes the timing data to matrix_inv.c where global statistics
 * live, so that:
 *   - bsolve.c can show correct Phase 4 timing
 *   - Summary can show correct "Total Matrix Inversion time"
 *
 * WITHOUT this line:
 *   - Inversion happens ✅
 *   - Timing is calculated ✅
 *   - Timing is printed ✅
 *   - Global stats NOT updated ❌ (Phase 4 shows 0.000000)
 *
 * WITH this line:
 *   - Inversion happens ✅
 *   - Timing is calculated ✅
 *   - Timing is printed ✅
 *   - Global stats updated ✅ (Phase 4 shows actual time!)
 */
/*----------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------*/
/*
void invert_this_matrix(a)
    matrix *a;
{
  int n, j;

  struct timeval inv_start, inv_finish;
  double inv_duration;

  a->invert = 0;
  n = get_num_columns(a);
  printf("Inverting matrix (%dx%d) \n", n, n);

  //---------------------------------------------------
  printf("\n-----------[matrix.c]-invert_this_matrix -----------\n");
  printf("Inverting the matrix, Please wait..\n");
  gettimeofday(&inv_start, NULL);
  //---------------------------------------------------

  //-------------- Parallel -----------
  mat_inv(a->value, n);

  //-------------- Sequential -----------


  //  for(j=0;j<n;j++)
  //    {
  //      //if(j%50==0) printf("\n");
  //      //printf("*"); fflush(stdout);
  //      scale_row(a,j);
  //      reduce_column(a,j);
  //   }


  gettimeofday(&inv_finish, NULL);
  inv_duration = ((double)(inv_finish.tv_sec - inv_start.tv_sec) * 1000000 + (double)(inv_finish.tv_usec - inv_start.tv_usec)) / 1000000;
  printf("############# Matrix inversion took: %lf seconds\n", inv_duration);
  printf("------ (%dx%d) Matrix inversion is done. -------- \n", n, n);


  // printf("\n\n----------- End -----------\n");
}
*/

/*------------------------------------------------*/
/*------------reduce_column----------------*/
void reduce_column(a, col)
    matrix *a;
int col;
{
  int n, i;
  n = get_num_rows(a);
  for (i = col + 1; i < col + n; i++)
  {
    reduce_row(a, i % n, col);
  }
}
/*-------------------------------------------*/
/*-------------reduce_row------------------*/
void reduce_row(a, row, pivot)
    matrix *a;
int row, pivot;
{
  double scale, value;
  int n, j;

  n = get_num_columns(a);
  scale = get_matrix_element(a, row, pivot);
  value = -get_matrix_element(a, pivot, pivot) * scale;
  put_matrix_element(a, row, pivot, value);
  for (j = pivot + 1; j < pivot + n; j++)
  {
    value = get_matrix_element(a, row, j % n) - get_matrix_element(a, pivot, j % n) * scale;
    put_matrix_element(a, row, j % n, value);
  }
}

void scale_row(a, row)
    matrix *a;
int row;
{
  double scale, value;
  int n, j;

  n = get_num_columns(a);
  scale = get_matrix_element(a, row, row);
  /*
  printf("scale = %e\n",scale);
  */
  scale = 1.0 / scale;
  put_matrix_element(a, row, row, scale);
  for (j = row + 1; j < row + n; j++)
  {
    value = get_matrix_element(a, row, j % n);
    put_matrix_element(a, row, j % n, value * scale);
  }
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* print the shape and other details for the matrix */
/*----------------------------------------------------------------------------------*/
void show_matrix_info(x)
    matrix *x;
{
  // printf("transpose=%d, invert=%d, rows=%d, columns=%d, value=0x%X\n",
  //  x->transpose,x->invert,x->rows,x->columns,x->value);
  printf("transpose=%d, invert=%d, rows=%d, columns=%d\n",
         x->transpose, x->invert, x->rows, x->columns);
}

/*----------------------------------------------------------------------------------*/
/* print the values in the matrix */
/*----------------------------------------------------------------------------------*/
/*
void show_matrix(x)
     matrix *x;
{
  int i,j;
  int row,col;
  int column;
  double k;
  FILE *output;
  char *buffer;
  int buf_size;

  check_invert(x);
  row=get_num_rows(x);
  col=get_num_columns(x);

  buf_size=16*4*12;
  buffer=(char *)malloc(buf_size*sizeof(char));

  output=open_file(0,"test1.out","w");
  column=0;

  for(i=0;i<row;i++)
    {
      printf("\nrow_i=%2d ",i);
      for(j=0;j<col;j++)
  {
    printf("%9.2e ",get_matrix_element(x,i,j));

    k=get_matrix_element(x,i,j);
    column=put_buffer(buf_size,buffer,column,"%9.2e ",k);
  }
      put_next_line(output,buffer);
      column=0;
    }
  printf("\n");
  fclose(output);
}
*/
void show_matrix(x)
    matrix *x;
{
  int i, j;
  int row, col;
  int column;
  double k;
  FILE *output;
  char *buffer;
  int buf_size;

  FILE *output2;

  check_invert(x);
  row = get_num_rows(x);
  col = get_num_columns(x);

  buf_size = 16 * 4 * 12;
  buffer = (char *)malloc(buf_size * sizeof(char));

  output = open_file(0, "test2.out", "w");
  output2 = fopen("before_inv_mat.out", "a+");
  column = 0;

  printf("\n");
  printf("\n");

  for (i = 0; i < row; i++)
  {
    // printf("\nrow_i=%2d ",i);
    for (j = 0; j < col; j++)
    {
      // printf("%9.2e ",get_matrix_element(x,i,j));
      // printf("%9.17f ",get_matrix_element(x,i,j));

      fprintf(output2, "%.1f ", get_matrix_element(x, i, j));

      // printf("###");
      k = get_matrix_element(x, i, j);
      // column=put_buffer(buf_size,buffer,column,"%9.2e ",k);
      column = put_buffer(buf_size, buffer, column, "%9.0f ", k);
    }
    fprintf(output2, " \n\n");

    put_next_line(output, buffer);
    column = 0;
  }

  printf("\n");

  fclose(output);
  fclose(output2);
}

/*----------------------------------------------------------------------------------*/
/* print the special values in the matrix */
/*----------------------------------------------------------------------------------*/
void show_matrix1(x)
    matrix *x;
{
  int i, j;
  int row, col;

  check_invert(x);
  row = get_num_rows(x);
  col = get_num_columns(x);

  /*for(i=0;i<row;i++)*/
  for (i = 0; i < 4; i++)
  {
    /*for(j=0;j<col;j++)*/
    for (j = 0; j < 4; j++)
    {
      /* printf("%7.6f ",get_matrix_element(x,i,j)); */
      printf("%5.4f ", get_matrix_element(x, i, j));
    }
    printf("\n row i= %d \n", i);
  }
  printf("\n");
}

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
void show_matrix2_(x)
    matrix *x;
{
  int i, j;
  int row, col;
  int column;
  double k;
  FILE *output;
  char *buffer;
  int buf_size;

  FILE *output2;

  check_invert(x);
  row = get_num_rows(x);
  col = get_num_columns(x);

  buf_size = 16 * 4 * 12;
  buffer = (char *)malloc(buf_size * sizeof(char));

  output = open_file(0, "test_outout2.out", "w");
  output2 = fopen("output2.out", "w");
  column = 0;

  printf("\n");
  printf("\n");

  for (i = 0; i < row; i++)
  {
    // printf("\nrow_i=%2d ",i);
    for (j = 0; j < col; j++)
    {
      // printf("%9.2e ",get_matrix_element(x,i,j));
      // printf("%9.17f ",get_matrix_element(x,i,j));

      fprintf(output2, "%.1f ", get_matrix_element(x, i, j));

      // printf("###");
      k = get_matrix_element(x, i, j);
      // column=put_buffer(buf_size,buffer,column,"%9.2e ",k);
      column = put_buffer(buf_size, buffer, column, "%9.0f ", k);
    }
    fprintf(output2, " \n\n");

    put_next_line(output, buffer);
    column = 0;
  }

  printf("\n");

  fclose(output);
  fclose(output2);
}

void show_matrix2(x)
    matrix *x;
{
  int i, j;
  double k;
  int row, col;
  FILE *fptr;

  fptr = fopen("matrixA.txt", "w");

  check_invert(x);
  row = get_num_rows(x);
  col = get_num_columns(x);

  /*for(i=0;i<row;i++)*/
  for (i = 0; i < row; i++)
  {
    /*for(j=0;j<col;j++)*/
    for (j = 0; j < col; j++)
    {
      /* printf("%7.6f ",get_matrix_element(x,i,j)); */
      fprintf(fptr, "%5.4f ", get_matrix_element(x, i, j));
    }
    fprintf(fptr, "\n");
  }
  printf("\n");
  fclose(fptr);
}

/*
void show_matrix2(x)
    matrix *x;
{
  int i, j;
  double k;
  int row, col;
  FILE *fptr;
  row = get_num_rows(x);
  col = get_num_columns(x);

  fptr = fopen("matrixA.txt", "w");
  for (int i = 0; i < row; i++)
  {
    for (int j = 0; j < col; j++)
    {
      k = get_matrix_element(x, i, j);
       fprintf(fptr, "%5.4f ", k);
       //fprintf(fptr, "%d ", rand() % 10);
    }
    fprintf(fptr, "\n");
  }
  fclose(fptr);
}
*/
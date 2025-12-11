/*----------------------------------------------------------------------------------*/
/* INSTRUMENTED multiply_matrix() function - ADD TO matrix.c */
/*----------------------------------------------------------------------------------*/

/* Add these includes at the top of matrix.c if not already present */
#include "sys/time.h"
#include <sys/resource.h>

/* Add this external variable declaration near the top of matrix.c */
extern void update_multiply_matrix_stats(double duration, long long flops);

/* REPLACE the existing multiply_matrix() function with this instrumented version */

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

  struct timeval mul_start, mul_finish;
  double mul_duration;
  long long flops;
  long vmrss_before, vmsize_before, vmrss_after, vmsize_after;

  check_multiply_shape(a, b);
  check_multiply_size(a, b, x);
  check_memory(x);

  /* START TIMING */
  gettimeofday(&mul_start, NULL);

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

  if (x->rows > 1 && x->columns > 1)
  {
    printf("\n------- Multiply1  with (x->rows > 1 && x->columns > 1) ------- \n");
    printf("1--- multiply matrix A ; result=(%dx%d)\n", a->rows, a->columns);
    printf("1--- multiply matrix B ; result=(%dx%d)\n", b->rows, b->columns);
    printf("1--- multiply matrix X ; result=(%dx%d)\n", x->rows, x->columns);

    if ((a->rows == a->columns) && (b->rows == b->columns))
    {
      printf("\n@@@@@@@@@@@@@ multiply square matrix A*B ; result=(%dx%d)\n\n", a->rows, a->columns);
    }

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

          x->value[column_b * x->rows + row_a] = value_new + value_new2;
        }
      }
    }
    else
    {
      printf("\n------- (a->transpose == 1) ------- \n");

#pragma omp parallel for schedule(dynamic, 50) collapse(2) private(row_a, column_b, i, value_new, value_new2)
      for (int row_a = 0; row_a < a_row_num; row_a++)
      {
        for (int column_b = 0; column_b < b_col_num; column_b++)
        {
          value_new = 0.0;
          value_new2 = 0.0;

          for (int i = 0; i < a_col_num; i += 2)
          {
            value_new = value_new + (a->value[row_a * a->columns + i] * b->value[column_b * b_row_num + i]);
            value_new2 = value_new2 + (a->value[row_a * a->columns + (i + 1)] * b->value[column_b * b_row_num + (i + 1)]);
          }

          x->value[column_b * x->rows + row_a] = value_new + value_new2;
        }
      }
    }
  }
  else /* Vector operations - keep original code */
  {
    for (row_a = 0; row_a < a_row_num; row_a++)
    {
      for (column_b = 0; column_b < b_col_num; column_b++)
      {
        value_new = 0.0;
        for (i = 0; i < a_col_num; i++)
        {
          value = get_matrix_element(a, row_a, i) * get_matrix_element(b, i, column_b);
          value_new = value_new + value;
        }
        put_matrix_element(x, row_a, column_b, value_new);
      }
    }
  }

  /* END TIMING */
  gettimeofday(&mul_finish, NULL);
  mul_duration = ((double)(mul_finish.tv_sec - mul_start.tv_sec) * 1000000 + 
                  (double)(mul_finish.tv_usec - mul_start.tv_usec)) / 1000000;

  /* Calculate FLOPS: 2 * m * n * k for matrix multiply */
  flops = 2LL * (long long)a_row_num * (long long)b_col_num * (long long)a_col_num;
  
  /* Update global statistics through external function */
  update_multiply_matrix_stats(mul_duration, flops);

  /* Print timing for this multiplication */
  double gflops = (double)flops / mul_duration / 1.0e9;
  printf("==> multiply_matrix(%dx%d * %dx%d) took %.6f sec (%.2f GFLOPS)\n",
         a_row_num, a_col_num, b_row_num, b_col_num, mul_duration, gflops);
}

/*----------------------------------------------------------------------------------*/

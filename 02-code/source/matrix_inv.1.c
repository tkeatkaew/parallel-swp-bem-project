/*----------------------------------------------------------------------------------*/
/*--------------------------------- matrix_inv.c ---------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* routines for performing all operations on matrix structures */
/*----------------------------------------------------------------------------------*/
#include <stdlib.h> 
#include <stdio.h> 
/*----------------------------------------------------------------------------------*/
#include "matrix_types.h" 

#include <time.h>
#include <omp.h>
#include <x86intrin.h>

#include "cblas.h"
#include "lapacke.h"
#include "sys/time.h"

extern void dgemm_(char*, char*, int*, int*,int*, double*, double*, int*, double*, int*, double*, double*, int*);

/*
for (row_a = 0; row_a < a_row_num; row_a++)
  for (column_b = 0; column_b < b_col_num; column_b++)
    for (i = 0; i < a_col_num; i++)
          value = get_matrix_element(a, row_a, i) * get_matrix_element(b, i, column_b);
          value_new = value_new + value;
*/

/*----------------------------------------------------------------------------------*/
/* get element from matrix */
/*----------------------------------------------------------------------------------*/
double get_matrix_element1(x, i, j,rows)
double *x;
int i, j,rows;
{
  double value;

    value = x[j * rows + i];


  return (value);
}

void dgemm(double * a, double * b, double * x, int a_row_num, int a_col_num, int b_row_num, int b_col_num)
{
    double value,value_new;
    int i;
  //x->rows=get_num_rows(a);
  //x->columns=get_num_columns(b);


    for (int row_a = 0;row_a < a_row_num; row_a++)
    {
        for (int column_b = 0; column_b < b_col_num; column_b++)
        {
            value_new = 0.0;
            for (i = 0; i < a_col_num; i++)
            {

               // value = get_matrix_element(a, row_a, i) * get_matrix_element(b, i, column_b);
               // value_new=value_new+value;


                //double get_matrix_element(x, i, j)
                //value = a[i * n + k] * b[k * n + j];
                
 /*       
 double get_matrix_element(x, i, j)     
    if (x->transpose == 0)
    {
        value = x->value[j * x->rows + i]; // only for matrix_multiply
    }
    else
    {
        value = x->value[i * x->rows + j];
    }
*/
                //get_matrix_element(a, row_a, i)*get_matrix_element(b, i, column_b);

                value = a[i * a_row_num + row_a] * b[column_b * b_row_num + i];

                value_new=value_new+value;

               //   c[i * n + j] += a[i * n + k] * b[k * n + j];
            }

/*
void put_matrix_element(x, i, j, value)
    x->value[j * x->rows + i] = value;

*/

            //put_matrix_element(x, row_a, column_b, value_new);
            // x->value[column_b * x->rows + row_a] = value_new;

            x[column_b * a_row_num + row_a] = value_new;

        }
    }
}

/*
void dgemm(double * a, double * b, double * c, int a_row_num, int a_col_num, int b_row_num, int b_col_num)
{
    double value,value_new;

    for (int row_a = 0;row_a < a_row_num; row_a++)
    {
        for (int column_b = 0; column_b < b_col_num; column_b++)
        {
            value_new = 0.0;
            for (int i = 0; i < a_col_num; i++)
            {

               // value=get_matrix_element1(a,row_a,i,a_row_num)*get_matrix_element1(b,i,column_b,b_row_num);
               // value_new=value_new+value;

                //double get_matrix_element(x, i, j)

                //value = a[i * n + k] * b[k * n + j];
                
                //value = a[i * a_row_num + row_a] *  b[column_b * b_row_num + i];
                //value_new=value_new+value;

               // c[i * n + j] += a[i * n + k] * b[k * n + j];
            }

            #if CHECK_MATRIX
                    check_invert(x);
                    check_matrix_index(x, i, j);
            #endif

            c[column_b * a_row_num + row_a] = value_new;

        }
    }
}
*/

/*
void dgemm(double *restrict a, double *restrict b, double *restrict c, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < n; j++)
            {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
}
*/


void dgemm2(double *a, double *b, double *c, int n) 
{ 
  for(int i=0; i<n; i++) { 
    for(int j=0; j<n; j++) {  // row(k,C) 
      for(int k=0; k<n; k++) { 
        c[i*n + j] = a[i*n + k]*b[k*n + j] + c[i*n + j];
      } 
    } 
  }
}

void dgemm_tpl(double *restrict a, double *restrict b, double *restrict c, int n)
{
 // #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < n; j++)
            {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
}


void gemm_tlp_simd(double *a, double *b, double *c, int n)
{
//#pragma omp parallel for schedule(dynamic)
//#pragma omp parallel for schedule(static)
// pragma omp parallel for private(m0,m1,m2,m3)
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    { // D1
        for (int j = 0; j < n; j += 8)
        { // D2
            __m256 m0 = _mm256_setzero_ps();
            for (int k = 0; k < n; k++)
            { // D3
                __m256 m1 = _mm256_broadcast_ss(a + i * n + k);
                __m256 m2 = _mm256_load_ps((b + k * n + j));
                __m256 m3 = _mm256_mul_ps(m1, m2);
                m0 = _mm256_add_ps(m0, m3);
            }
            _mm256_store_ps(c + i * n + j, m0);
        }
    }
}


/*
void gemm_tlp_simd2(double *restrict a, double *restrict b, double *restrict c, int n)
{
    int mb = 256;
    int kb = 32;
    int nb = 32;

    int ii, jj, kk;
    //#pragma omp parallel for schedule(dynamic)
    //#pragma omp parallel for schedule(static)
    // pragma omp parallel for private(m0,m1,m2,m3)
    //#pragma omp parallel for

    for (ii = 0; ii < n; ii += nb)
    {
        for (jj = 0; jj < n; jj += kb)
        {

            //#pragma omp parallel for schedule(static)
#pragma omp parallel for
            for (int i = ii; i < ii + nb; i++)
            { // D1
                for (int j = jj; j < jj + kb; j += 8)
                { // D2
                    __m256 m0 = _mm256_setzero_ps();
                    for (int k = 0; k < n; k++)
                    // for (int k = i; k < i+nb; k++)
                    // for (int k = j; k < j+kb; k++)
                    { // D3
                        __m256 m1 = _mm256_broadcast_ss(a + i * n + k);
                        __m256 m2 = _mm256_load_ps((b + k * n + j));
                        __m256 m3 = _mm256_mul_ps(m1, m2);
                        m0 = _mm256_add_ps(m0, m3);
                    }
                    _mm256_store_ps(c + i * n + j, m0);
                }
            }
        }
    }
}

*/


//----------------------------------------------

//void mat_mul(double *A,double *B,double *X,int rows,int columns)
void mat_mul(double *A,double *B,double *X,int m,int k,int n)
{

  int i;
  //printf("test!\n");
  /*
  if(argc<4){
    printf("Input Error\n");
    return 1;
  }
*/
//int N=2048;
/*
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int k = atoi(argv[3]);
  */

  //int m = rows;
  //int n = columns;
  //int k = columns;
  
  //int sizeofa = m * k;
  //int sizeofb = k * n;
  int sizeofx = m * n;
  char ta = 'N';
  char tb = 'N';
  double alpha = 1.0;
  double beta = 0.0;

  struct timeval start,finish;
  double duration;


  //#if 0
  printf("m=%d,n=%d,k=%d,alpha=%lf,beta=%lf,sizeofc=%d\n",m,n,k,alpha,beta,sizeofx);
  
  gettimeofday(&start, NULL);

 //dgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, X, &m);


// cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, k,n,m, 1.,A,k, B,n,0., X,n);

cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, 
           m, n, k, alpha, A, k, B, n, beta, X, n);

//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
//                m, n, k, alpha, A, k, B, n, beta, X, n);

  gettimeofday(&finish, NULL);

  duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  double gflops = 2.0 * m *n*k;
  gflops = gflops/duration*1.0e-6;

 // FILE *fp;
  //fp = fopen("timeDGEMM.txt", "a");
 // fprintf(fp, "%dx%dx%d\t%lf s\t%lf MFLOPS\n", m, n, k, duration, gflops);
 // fclose(fp);


  //return 0;
}


/*
extern void dgetrf_(const int* M, const int* N, double* A, const int* LDA,
    int* IPIV, int* INFO);
*/

/*
using namespace std;

extern "C"
{
    // LU decomoposition of a general matrix
    //  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    //  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    void openblas_set_num_threads(int num_threads);
    void goto_set_num_threads(int num_threads);

    int openblas_get_num_procs(void);
    int openblas_get_num_threads(void);
    int openblas_get_parallel(void);
    // int openblas_get_config(void);
}
*/

lapack_int mat_inv(double *A, unsigned n)
{
    int ipiv[n + 1];
    lapack_int ret;

    ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR,
                         n,
                         n,
                         A,
                         n,
                         ipiv);

    if (ret != 0)
        return ret;

    ret = LAPACKE_dgetri(LAPACK_COL_MAJOR,
                         n,
                         A,
                         n,
                         ipiv);
    return ret;
}

/*
// void inverse(vector<vector<double>> A, int N)
void inverse(double *A, int N)
{
    int *IPIV = new int[N];
    int LWORK = N * N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N, &N, A, &N, IPIV, &INFO);
    dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

    delete[] IPIV;
    delete[] WORK;
}
*/
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------*/
/*------------------------------- matrix_types.h -----------------------------------*/
/*----------------------------------------------------------------------------------*/
/* structure for holding matrices and row- and column-vectors of scalars */

typedef struct{
  int transpose; /* =1 means the transpose must be taken to get correct values */
  int invert;    /* =1 means the inverse must be taken to get correct values */
  int rows;      /* number of rows; =1 means a row-vector */
  int columns;   /* number of columns; =1 means a column-vector */
  double *value; /* pointer to one-dimensional array of matrix values */
} matrix;

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

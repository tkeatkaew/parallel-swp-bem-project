/*----------------------------------------------------------------------------------*/
/*-------------------------------- memory_types.h ----------------------------------*/
/*----------------------------------------------------------------------------------*/
/* structure for holding the vectors we need in memory */

typedef struct {
  matrix *bvv;         /* boundary voltage vector */
  matrix *bcv;         /* boundary current vector */
  matrix *vgv;         /* voltage geometry vector */
  matrix *cgv;         /* current geometry vector */
  co_matrix *co_vgv;   /* voltage geometry vector 1st derivative */
  co_matrix *co_cgv;   /* current geometry vector 1st derivative */
  ten_matrix *ten_vgv; /* voltage geometry vector 2nd derivative */
  ten_matrix *ten_cgv; /* current geometry vector 2nd derivative */
} bem_vectors;

/*----------------------------------------------------------------------------------*/
/* structure for holding the final results from the bem calculations */

typedef struct {
  double V;           /* voltage at point P */
  coordinates dV;     /* gradient of voltage */
  tensor  d2V;        /* second derivatives of voltage */
} bem_results;

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/* structure for holding definition of raster over which to do calculations */

typedef struct {
  int nx;           /* number of points in x-direction */
  int ny;           /* number of points in y-direction */
  coordinates P1;   /* x,y coordinates of first point */
  coordinates P2;   /* x,y coordinates of last point */
} raster;

/*----------------------------------------------------------------------------------*/
/* structure for holding definition of section over which to do calculations */

typedef struct {
  int n;            /* number of points */
  coordinates P1;   /* x,y coordinates of first point */
  coordinates P2;   /* x,y coordinates of last point */
  double step;      /* distance between points */
} section;

/*----------------------------------------------------------------------------------*/
/* structure for holding definition of interval over which to do calculations */

typedef struct {
  int nt;           /* number of instants */
  double t1;        /* time of first instant */
  double t2;        /* time of last instant */
} interval;

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/

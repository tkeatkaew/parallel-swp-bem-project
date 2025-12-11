/*----------------------------------------------------------------------------------*/
/*------------------------------ boundary_types.h ----------------------------------*/
/*----------------------------------------------------------------------------------*/
/* structure for holding Cartesian coordinates of points in 2 dimensions */
typedef double coordinates[2]; /* array: 1st and 2nd elements = x and y coordinates */

/* structure for holding Cartesian tensor of points in 2 dimensions */
typedef double tensor[2][2]; /* array: elements = xx,xy,yx,yy components of tensor */

/*----------------------------------------------------------------------------------*/
/* structure for holding open or closed plane curves (i.e. in 2 dimesions) */

typedef struct{
  int links;       /* number of extra links to this path */
  int close;       /* =1 means a closed curve; last point joins to first */
  int reverse;     /* =1 means the curve must be reversed to get correct orientation */
  int points;      /* number of points on the curve */
  coordinates *xy; /* pointer to one-dimensional array of coordinates for points */
  double *value;   /* pointer to one-dimensional array of values at each point */
} path;

/*----------------------------------------------------------------------------------*/
/* structure for holding boundary of region, with holes, in 2 dimesions */
/* the first closed path is the outermost path */

typedef struct{
  int curve;      /* 0 = anticlockwise (left); 1 = clockwise (right) */ 
  int components; /* number of closed paths for the boundary */
  int *level;     /* array of flags; flag 0 = path outside zone; 1 = path inside zone */ 
  path **loop;    /* pointer to one-dimensional array of pointers to closed paths */
  double *bvv;    /* pointer to boundary voltage vector calculated from data */
  double *bcv;    /* pointer to boundary current vector found by matrix solution */
} boundary;

/*----------------------------------------------------------------------------------*/
/* structure for holding a name (of a file) and a link to a path */

typedef struct{
  path *path_p;
  char name[32];
} path_link;

/*----------------------------------------------------------------------------------*/
/* structure for holding catchment region */
/* there is a list of zones; 1 zone = 1 region between contours = 1 boundary */
/* there is a list of paths; this is for all paths in the catchment region */

typedef struct{
  int num_zones;        /* number of zones (boundaries) in the catchment region */
  int max_zones;        /* maximum number of zones */
  int previous_zone;    /* calculation done in this zone last time */
  boundary **zones;     /* pointer to list of pointers to boundaries */
  int num_paths;        /* number of paths in the catchement region */
  int max_paths;        /* maximum number of paths */
  path_link *path_list; /* pointer to list of links to paths */
} catchment;

/*----------------------------------------------------------------------------------*/


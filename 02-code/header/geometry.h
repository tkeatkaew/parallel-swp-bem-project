/* ../source/geometry.c */
void convert_PQ(coordinates Qa, coordinates Qb, coordinates P, double *x, double *y1, double *y2);
void rotate_to_PQ(double x, double y, coordinates Qa, coordinates Qb, coordinates R);
void double_rotate_to_PQ(double a, double b, double c, double d, coordinates Qa, coordinates Qb, tensor R);
double atan3(double y2, double y1, double x);
double atanv(coordinates Q1, coordinates Q2, coordinates P);

/* ../source/streamline.c */
int quadratic(double a[3], double z[2][2]);
int cubic(double a[4], double z[3][2]);
int quartic(double a[5], double z[4][2]);
int check_gv(int n, double v[4], double gv[4]);
double choose_max_or_min(int direction, double first[2], double second[2], double angle[4], int n, double r);
int myquartic(double *c, double *x);
int myquadratic(double *c, double *x);
void edit1_my_follow_stream(int direction, coordinates P, coordinates dV, tensor d2V, coordinates dP, double r);
void plot_streamlines(catchment *c, int n_streams, path **streamlines, unsigned char *file);
void plot_streamlines_v2(catchment *c, int n_streams, path **streamlines, unsigned char *file);
void plot_1_streamline(catchment *c, path *streamline, FILE *output);
double streamline_loop(coordinates P, catchment *c, int direction, int max_steps, double step_size, path *streamline, bem_vectors *vectors, bem_results *v1);
double SCA_loop_GH0_v2(coordinates P, catchment *c, int direction, int max_steps, double step_size, path *streamline, bem_vectors *vectors, bem_results *v1);
double SCA_loop_GH0_v3(coordinates P, catchment *c, int direction, int max_steps, double step_size, path *streamline, bem_vectors *vectors, bem_results *v1);
double SCA_loop_GH0_v4(coordinates P, catchment *c, int direction, int max_steps, double step_size, path *streamline, bem_vectors *vectors, bem_results *v1);
double SCA_loop_GH0_v5(coordinates P, catchment *c, int direction, int max_steps, double step_size, path *streamline, bem_vectors *vectors, bem_results *v1);

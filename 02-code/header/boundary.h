/* ../source/boundary.c */
int boundary_loops(unsigned char *file);
boundary *create_boundary(int n);
boundary *destroy_boundary(boundary *b);
boundary *destroy_boundary_ignore_paths(boundary *b);
int num_points_in_zone(boundary *b);
void get_boundary(unsigned char *file, boundary *b);
void plot_boundary(boundary *b, unsigned char *file);

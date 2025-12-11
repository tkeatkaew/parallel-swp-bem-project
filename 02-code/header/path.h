/* ../source/path.c */
path *create_path(int points, int make_xy, int make_value);
path *destroy_path(path *x);
void check_value_memory(path *p);
void check_xy_memory(path *p);
void check_path_index(path *p, int i);
void reverse_path(path *p);
void close_path(path *p);
void open_path(path *p);
double get_path_value(path *p, int i);
void get_path_xy(path *p, int i, coordinates xy);
void put_path_value(path *p, int i, double val);
void put_path_xy(path *p, int i, coordinates xy);
int num_points_in_path(path *p);
void show_path_info(path *p);
void show_path(path *p);
int path_length(unsigned char *file);
void get_path(unsigned char *file, path *p);
void plot_path(path *p, unsigned char *file);

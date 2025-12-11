/* ../source/vcalc.c */
double voltage_on_path(catchment *c, double s, int segment, path *this_path);
double voltage_outside_catchment(void);
double calculate_in_same_zone(boundary *b, coordinates P, bem_vectors *x, bem_results *R);
double calculate_in_new_zone(boundary *b, coordinates P, bem_vectors *x, bem_results *R);
double calculate_inside_catchment(catchment *c, coordinates P, bem_vectors *vectors, bem_results *voltage, int *new_z);

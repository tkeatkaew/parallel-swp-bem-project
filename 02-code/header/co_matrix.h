/* ../source/co_matrix.c */
co_matrix *create_co_matrix(int rows, int columns);
void attach_co_matrix(co_matrix *x, int rows, int columns, coordinates *data);
co_matrix *destroy_co_matrix(co_matrix *x);
void check_co_matrix_index(co_matrix *x, int i, int j);
void check_co_invert(co_matrix *x);
void check_co_memory(co_matrix *x);
void check_co_multiply_shape(co_matrix *a, matrix *b);
void check_co_multiply_size(co_matrix *a, matrix *b, co_matrix *x);
int get_co_num_columns(co_matrix *x);
int get_co_num_rows(co_matrix *x);
coordinates *startof_co_matrix(co_matrix *x);
void put_co_matrix_element(co_matrix *x, int i, int j, coordinates value);
void get_co_matrix_element(co_matrix *x, int i, int j, coordinates value);
void put_block_co_matrix_element(co_matrix *x, int offset_i, int offset_j, int i, int j, coordinates value);
void get_block_co_matrix_element(co_matrix *x, int offset_i, int offset_j, int i, int j, coordinates value);
void multiply_co_matrix(co_matrix *a, matrix *b, co_matrix *x);
void multiply_co_element_row(co_matrix *a, matrix *b, int row_a, co_matrix *x);
void multiply_co_element_column(co_matrix *a, matrix *b, int row_a, int column_b, co_matrix *x);

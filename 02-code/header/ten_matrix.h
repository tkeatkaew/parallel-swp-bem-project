/* ../source/ten_matrix.c */
ten_matrix *create_ten_matrix(int rows, int columns);
void attach_ten_matrix(ten_matrix *x, int rows, int columns, tensor *data);
ten_matrix *destroy_ten_matrix(ten_matrix *x);
void check_ten_matrix_index(ten_matrix *x, int i, int j);
void check_ten_invert(ten_matrix *x);
void check_ten_memory(ten_matrix *x);
void check_ten_multiply_shape(ten_matrix *a, matrix *b);
void check_ten_multiply_size(ten_matrix *a, matrix *b, ten_matrix *x);
int get_ten_num_columns(ten_matrix *x);
int get_ten_num_rows(ten_matrix *x);
tensor *startof_ten_matrix(ten_matrix *x);
void put_ten_matrix_element(ten_matrix *x, int i, int j, tensor value);
void get_ten_matrix_element(ten_matrix *x, int i, int j, tensor value);
void put_block_ten_matrix_element(ten_matrix *x, int offset_i, int offset_j, int i, int j, tensor value);
void get_block_ten_matrix_element(ten_matrix *x, int offset_i, int offset_j, int i, int j, tensor value);
void multiply_ten_matrix(ten_matrix *a, matrix *b, ten_matrix *x);
void multiply_ten_element_row(ten_matrix *a, matrix *b, int row_a, ten_matrix *x);
void multiply_ten_element_column(ten_matrix *a, matrix *b, int row_a, int column_b, ten_matrix *x);

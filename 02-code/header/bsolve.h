/* ../source/bsolve.c */
void make_boundary_voltage_vector(boundary *b, matrix *bvv);
void fill_boundary_voltage_vector(path *path_i, double *result);
void make_boundary_current_vector(boundary *b, matrix *V, matrix *J);
void make_bcv_no_KCL(boundary *b, matrix *V, matrix *J);
void make_bcv_use_KCL(boundary *b, matrix *V, matrix *J);
void make_boundary_vector(boundary *b, matrix *bvv, matrix *bcv);
double make_internal_voltage(boundary *b, matrix *bvv, matrix *bcv, coordinates P, matrix *vgv, matrix *cgv);
void make_internal_grad_voltage(boundary *b, matrix *bvv, matrix *bcv, coordinates P, co_matrix *co_vgv, co_matrix *co_cgv, coordinates Gv);
void make_internal_sec_grad_voltage(boundary *b, matrix *bvv, matrix *bcv, coordinates P, ten_matrix *ten_vgv, ten_matrix *ten_cgv, tensor Gv);

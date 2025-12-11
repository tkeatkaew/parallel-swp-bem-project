/* ../source/file.c */
int put_buffer(int n, unsigned char *buffer, int offset, unsigned char *format, double value);
FILE *open_file(int use_path, char *file, char *mode);
int count_lines(int use_path, char *file);
int get_next_line(FILE *input, int n, unsigned char *buffer);
int get_next_line_verbose(FILE *input, int k, int n, unsigned char *buffer);
void put_next_line(FILE *output, unsigned char *buffer);
void catchment_path(int n, unsigned char *c_path);
void make_gpl_file(char *datafile, char *title, char *xrange, char *yrange);
void make_gpl2_file(char *datafile, char *title, char *xlabel, char *ylabel);

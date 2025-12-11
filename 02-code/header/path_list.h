/* ../source/path_list.c */
path_link *create_path_list(int n);
path_link *destroy_path_list(int n, path_link *path_list);
int search_path_list(unsigned char *file_name, int n, path_link *path_list);
int load_path_list(unsigned char *file_name, int index, path_link *path_list);
path *get_path_list(int index, path_link *path_list);

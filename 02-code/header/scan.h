/* ../source/scan.c */
void put_raster(char *data, raster *ras);
void show_raster(raster *ras);
double x_raster(raster *ras, int i);
double y_raster(raster *ras, int j);
void put_section(char *data, section *sec);
void put_sectionV2(int Nseg, coordinates PA, coordinates PB, section *sec);
void show_section(section *sec);
void xy_section(section *sec, int i, coordinates xy);
void put_interval(char *data, interval *inter);
void show_interval(interval *inter);
double t_interval(interval *inter, int i);

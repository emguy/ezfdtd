/* NOTE: This program is free software; you can redistribute it 
 * and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either 
 * version 3, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * Bugs can be reported to Yu Zhang <zhang235@mcmaster.ca>.
 *
 *     File Name :   
 * Last Modified :
 */



#ifndef H5IO_H
#define H5IO_H
void ***h5_load3(const char *file_name, const char *dset_name, size_t xmax, size_t ymax, size_t zmax);
void **h5_load2(const char *file_name, const char *dset_name, size_t xmax, size_t ymax);
void **h5_load1(const char *file_name, const char *dset_name, size_t xmax);
int h5_get_attr(const char *file_name, const char *grp_name, const char *attr_name, void *value);
int h5_set_file(const char *file_name);
int h5_write2(const char *file_name, const char *dset_name, int type, void **data, size_t xmax, size_t ymax);
int h5_empty3(const char *file_name, const char *dset_name, int type, size_t xmax, size_t ymax, size_t zmax);
int h5_slab(const char *file_name, const char *dset_name, int type, void **buffer, size_t xmax, size_t ymax, int index);
#endif

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
 *     File Name : domain.h
 * Last Modified : Fri 12 Oct 2012 02:56:16 PM EDT
 */

#ifndef DOMAIN_H
#define DOMAIN_H

extern double ***epsilon;

extern double d_tx;
extern double d_ty;
extern double d_tz;

int setup_domain (char *file_name);
int setup_fields (char *file_name);
int get_partition(int x, int y, int z);
int in_partition_x0(int x, int y, int z);
int in_partition_x1(int x, int y, int z);
int in_partition_y0(int x, int y, int z);
int in_partition_y1(int x, int y, int z);
int in_partition_z0(int x, int y, int z);
int in_partition_z1(int x, int y, int z);
int in_partition_main(int x, int y, int z);


#endif

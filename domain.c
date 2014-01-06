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
 *     File Name : domain.c
 * Last Modified : Fri 12 Oct 2012 02:56:13 PM EDT
 */

#include <string.h>
#include "tools.h"
#include "h5io.h"
#include "mem.h"
#include "domain.h"

double ***epsilon;
double d_tx;
double d_ty;
double d_tz;

int setup_domain (char *file_name) 
{
    int partition_index;
    char attr_name[30];
    int boundary_type;
    int status;

    status = h5_get_attr(file_name, "settings", "domain_size_x", &main_length_x);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "domain_size_y", &main_length_y);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "domain_size_z", &main_length_z);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "boundaries", "abc_thickness", &abc_size);
    inspect(status, "fail to get h5 attributes");

    /* fill region data */
    for (partition_index = 1; partition_index < 7; partition_index++)
    {
        sprintf(attr_name, "boundary_%01d", partition_index);
        status = h5_get_attr(file_name, "boundaries", attr_name, &boundary_type);
        inspect(status, "fail to get h5 attributes");
        partition_data[partition_index].boundary_type = boundary_type;

        if (boundary_type == boundary_pml)
            partition_data[partition_index].thickness = abc_size;
        else
            partition_data[partition_index].thickness = 0;
    }

    /* domain size */
    total_length_x = main_length_x + partition_data[partition_x0].thickness + partition_data[partition_x1].thickness;
    total_length_y = main_length_y + partition_data[partition_y0].thickness + partition_data[partition_y1].thickness;
    total_length_z = main_length_z + partition_data[partition_z0].thickness + partition_data[partition_z1].thickness;

    /* main grid */
    partition_data[partition_main].x_start = partition_data[partition_x0].thickness;
    partition_data[partition_main].x_stop  = total_length_x - partition_data[partition_x1].thickness;
    partition_data[partition_main].y_start = partition_data[partition_y0].thickness;
    partition_data[partition_main].y_stop  = total_length_y - partition_data[partition_y1].thickness;
    partition_data[partition_main].z_start = partition_data[partition_z0].thickness;
    partition_data[partition_main].z_stop  = total_length_z - partition_data[partition_z1].thickness;
    partition_data[partition_main].size_x = partition_data[partition_main].x_stop - partition_data[partition_main].x_start;
    partition_data[partition_main].size_y = partition_data[partition_main].y_stop - partition_data[partition_main].y_start;
    partition_data[partition_main].size_z = partition_data[partition_main].z_stop - partition_data[partition_main].z_start;

    /* region 1 z0 */
    partition_data[partition_z0].x_start = 0;
    partition_data[partition_z0].x_stop  = total_length_x;
    partition_data[partition_z0].y_start = 0;
    partition_data[partition_z0].y_stop  = total_length_y;
    partition_data[partition_z0].z_start = 0;
    partition_data[partition_z0].z_stop  = partition_data[partition_z0].thickness;
    partition_data[partition_z0].size_x = partition_data[partition_z0].x_stop - partition_data[partition_z0].x_start;
    partition_data[partition_z0].size_y = partition_data[partition_z0].y_stop - partition_data[partition_z0].y_start;
    partition_data[partition_z0].size_z = partition_data[partition_z0].z_stop - partition_data[partition_z0].z_start;

    /* region 2 z1 */
    partition_data[partition_z1].x_start = 0;
    partition_data[partition_z1].x_stop  = total_length_x;
    partition_data[partition_z1].y_start = 0;
    partition_data[partition_z1].y_stop  = total_length_y;
    partition_data[partition_z1].z_start = total_length_z - partition_data[partition_z1].thickness;
    partition_data[partition_z1].z_stop  = total_length_z;
    partition_data[partition_z1].size_x = partition_data[partition_z1].x_stop - partition_data[partition_z1].x_start;
    partition_data[partition_z1].size_y = partition_data[partition_z1].y_stop - partition_data[partition_z1].y_start;
    partition_data[partition_z1].size_z = partition_data[partition_z1].z_stop - partition_data[partition_z1].z_start;

    /* region 3 x0 */
    partition_data[partition_x0].x_start = 0;
    partition_data[partition_x0].x_stop  = partition_data[partition_x0].thickness;
    partition_data[partition_x0].y_start = 0;
    partition_data[partition_x0].y_stop  = total_length_y;
    partition_data[partition_x0].z_start = 0;
    partition_data[partition_x0].z_stop  = total_length_z;
    partition_data[partition_x0].size_x = partition_data[partition_x0].x_stop - partition_data[partition_x0].x_start;
    partition_data[partition_x0].size_y = partition_data[partition_x0].y_stop - partition_data[partition_x0].y_start;
    partition_data[partition_x0].size_z = partition_data[partition_x0].z_stop - partition_data[partition_x0].z_start;

    /* region 4 x1 */
    partition_data[partition_x1].x_start = total_length_x - partition_data[partition_x1].thickness;
    partition_data[partition_x1].x_stop  = total_length_x;
    partition_data[partition_x1].y_start = 0;
    partition_data[partition_x1].y_stop  = total_length_y;
    partition_data[partition_x1].z_start = 0;
    partition_data[partition_x1].z_stop  = total_length_z;
    partition_data[partition_x1].size_x = partition_data[partition_x1].x_stop - partition_data[partition_x1].x_start;
    partition_data[partition_x1].size_y = partition_data[partition_x1].y_stop - partition_data[partition_x1].y_start;
    partition_data[partition_x1].size_z = partition_data[partition_x1].z_stop - partition_data[partition_x1].z_start;

    /* region 5 y0 */
    partition_data[partition_y0].x_start = 0;
    partition_data[partition_y0].x_stop  = total_length_x;
    partition_data[partition_y0].y_start = 0;
    partition_data[partition_y0].y_stop  = partition_data[partition_y0].thickness;
    partition_data[partition_y0].z_start = 0;
    partition_data[partition_y0].z_stop  = total_length_z;
    partition_data[partition_y0].size_x = partition_data[partition_y0].x_stop - partition_data[partition_y0].x_start;
    partition_data[partition_y0].size_y = partition_data[partition_y0].y_stop - partition_data[partition_y0].y_start;
    partition_data[partition_y0].size_z = partition_data[partition_y0].z_stop - partition_data[partition_y0].z_start;

    /* region 6 y1 */
    partition_data[partition_y1].x_start = 0;
    partition_data[partition_y1].x_stop  = total_length_x;
    partition_data[partition_y1].y_start = total_length_y - partition_data[partition_y1].thickness;
    partition_data[partition_y1].y_stop  = total_length_y;
    partition_data[partition_y1].z_start = 0;
    partition_data[partition_y1].z_stop  = total_length_z;
    partition_data[partition_y1].size_x = partition_data[partition_y1].x_stop - partition_data[partition_y1].x_start;
    partition_data[partition_y1].size_y = partition_data[partition_y1].y_stop - partition_data[partition_y1].y_start;
    partition_data[partition_y1].size_z = partition_data[partition_y1].z_stop - partition_data[partition_y1].z_start;

    /*  
    printf("x1 start: %d\n", partition_data[partition_x1].x_start);
    printf("x1 stop: %d\n",  partition_data[partition_x1].x_stop);
    printf("x0 start: %d\n", partition_data[partition_x0].x_start);
    printf("x0 stop: %d\n",  partition_data[partition_x0].x_stop);
    */
    return 1;
}

int setup_fields (char* file_name) 
{
    int x_main, y_main, z_main;
    int x, y, z;
    int status;

    status = h5_get_attr(file_name, "settings", "mode", &mode);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "EPSILON0", &EPSILON0);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "MU0", &MU0);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "C0", &C0);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "d_x", &d_x);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "d_y", &d_y);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "d_z", &d_z);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "settings", "d_t", &d_t);
    inspect(status, "fail to get h5 attributes");

    d_tx = d_t / d_x;
    d_ty = d_t / d_y;
    d_tz = d_t / d_z;

    switch (mode)
    {
        case (mode_full):
            ex = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ex, "fail to allocate memory for field ex");
            ey = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ey, "fail to allocate memory for field ey");
            ez = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ez, "fail to allocate memory for field ez");
            hx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hx, "fail to allocate memory for field hx");
            hy = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hy, "fail to allocate memory for field hy");
            hz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hz, "fail to allocate memory for field hz");
            dx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dx, "fail to allocate memory for field dx");
            dy = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dy, "fail to allocate memory for field dy");
            dz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dz, "fail to allocate memory for field dz");
            bx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(bx, "fail to allocate memory for field bx");
            by = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(by, "fail to allocate memory for field by");
            bz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(bz, "fail to allocate memory for field bz");
            dipole_ex = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ex, "fail to allocate memory for field dipole_ex");
            dipole_ey = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ey, "fail to allocate memory for field dipole_ey");
            dipole_ez = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ez, "fail to allocate memory for field dipole_ez");
            dipole_hx = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hx, "fail to allocate memory for field dipole_hx");
            dipole_hy = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hy, "fail to allocate memory for field dipole_hy");
            dipole_hz = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hz, "fail to allocate memory for field dipole_hz");
            break;
        case (mode_tmx):
            ex = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ex, "fail to allocate memory for field ex");
            hy = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hy, "fail to allocate memory for field hy");
            hz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hz, "fail to allocate memory for field hz");
            dx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dx, "fail to allocate memory for field dx");
            by = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(by, "fail to allocate memory for field by");
            bz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(bz, "fail to allocate memory for field bz");
            dipole_ex = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ex, "fail to allocate memory for field dipole_ex");
            dipole_hy = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hy, "fail to allocate memory for field dipole_hy");
            dipole_hz = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hz, "fail to allocate memory for field dipole_hz");
            break;
        case (mode_tmy):
            ey = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ey, "fail to allocate memory for field ey");
            hx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hx, "fail to allocate memory for field hx");
            hz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hz, "fail to allocate memory for field hz");
            dy = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dy, "fail to allocate memory for field dy");
            bx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(bx, "fail to allocate memory for field bx");
            bz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(bz, "fail to allocate memory for field bz");
            dipole_ey = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ey, "fail to allocate memory for field dipole_ey");
            dipole_hx = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hx, "fail to allocate memory for field dipole_hx");
            dipole_hz = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hz, "fail to allocate memory for field dipole_hz");
            break;
        case (mode_tmz):
            ez = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ez, "fail to allocate memory for field ez");
            hx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hx, "fail to allocate memory for field hx");
            hy = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hy, "fail to allocate memory for field hy");
            dz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dz, "fail to allocate memory for field dz");
            bx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(bx, "fail to allocate memory for field bx");
            by = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(by, "fail to allocate memory for field by");
            dipole_ez = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ez, "fail to allocate memory for field dipole_ez");
            dipole_hx = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hx, "fail to allocate memory for field dipole_hx");
            dipole_hy = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hy, "fail to allocate memory for field dipole_hy");
            break;
        case (mode_tex):
            hx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hx, "fail to allocate memory for field hx");
            ey = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ey, "fail to allocate memory for field ey");
            ez = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ez, "fail to allocate memory for field ez");
            bx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(bx, "fail to allocate memory for field bx");
            dy = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dy, "fail to allocate memory for field dy");
            dz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dz, "fail to allocate memory for field dz");
            dipole_hx = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hx, "fail to allocate memory for field dipole_hx");
            dipole_ey = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ey, "fail to allocate memory for field dipole_ey");
            dipole_ez = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ez, "fail to allocate memory for field dipole_ez");
            break;
        case (mode_tey):
            hy = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hy, "fail to allocate memory for field hy");
            ex = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ex, "fail to allocate memory for field ex");
            ez = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ez, "fail to allocate memory for field ez");
            by = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(by, "fail to allocate memory for field by");
            dx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dx, "fail to allocate memory for field dx");
            dz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dz, "fail to allocate memory for field dz");
            dipole_hy = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hy, "fail to allocate memory for field dipole_hy");
            dipole_ex = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ex, "fail to allocate memory for field dipole_ex");
            dipole_ez = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ez, "fail to allocate memory for field dipole_ez");
            break;
        case (mode_tez):
            hz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(hz, "fail to allocate memory for field hz");
            ex = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ex, "fail to allocate memory for field ex");
            ey = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(ey, "fail to allocate memory for field ey");
            bz = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(bz, "fail to allocate memory for field bz");
            dx = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dx, "fail to allocate memory for field dx");
            dy = (double ***)mem3(type_double, total_length_x + 1, total_length_y + 1, total_length_z + 1);
            inspect(dy, "fail to allocate memory for field dy");
            dipole_hz = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_hz, "fail to allocate memory for field dipole_hz");
            dipole_ex = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ex, "fail to allocate memory for field dipole_ex");
            dipole_ey = (double ***)mem3(type_double, main_length_x, main_length_y, main_length_z);
            inspect(dipole_ey, "fail to allocate memory for field dipole_ey");
            break;
    }

    epsilon = (double ***)h5_load3(file_name, "/materials/epsilon",  total_length_x, total_length_y, total_length_z);
    inspect(epsilon, "fail to load hdf5 dataset");
    for (x = 0; x < total_length_x; x++)
        for (y = 0; y < total_length_y; y++)
            for (z = 0; z < total_length_z; z++)
            {
                if (in_partition_main(x, y, z))
                {
                    x_main = x - partition_data[partition_main].x_start;
                    y_main = y - partition_data[partition_main].y_start;
                    z_main = z - partition_data[partition_main].z_start;
                    epsilon[x][y][z] = EPSILON0 * epsilon[x_main][y_main][z_main];
                }
                else epsilon[x][y][z] = EPSILON0;
            }

    return 1;
}

int get_partition(int x, int y, int z)
{
    if (in_partition_y0(x, y, z)) return partition_y0;
    if (in_partition_y1(x, y, z)) return partition_y1;
    if (in_partition_x0(x, y, z)) return partition_x0;
    if (in_partition_x1(x, y, z)) return partition_x1;
    if (in_partition_z0(x, y, z)) return partition_z0;
    if (in_partition_z1(x, y, z)) return partition_z1;

    return partition_main;
}

int in_partition_x0(int x, int y, int z)
{
    if (   x >= partition_data[partition_x0].x_start && x < partition_data[partition_x0].x_stop
        && y >= partition_data[partition_x0].y_start && y < partition_data[partition_x0].y_stop
        && z >= partition_data[partition_x0].z_start && z < partition_data[partition_x0].z_stop)
        return 1;
    return 0;
}


int in_partition_x1(int x, int y, int z)
{
    if (   x >= partition_data[partition_x1].x_start && x < partition_data[partition_x1].x_stop
        && y >= partition_data[partition_x1].y_start && y < partition_data[partition_x1].y_stop
        && z >= partition_data[partition_x1].z_start && z < partition_data[partition_x1].z_stop)
        return 1;
    return 0;
}


int in_partition_y0(int x, int y, int z)
{
    if (   x >= partition_data[partition_y0].x_start && x < partition_data[partition_y0].x_stop
        && y >= partition_data[partition_y0].y_start && y < partition_data[partition_y0].y_stop
        && z >= partition_data[partition_y0].z_start && z < partition_data[partition_y0].z_stop)
        return 1;
    return 0;
}

int in_partition_y1(int x, int y, int z)
{
    if (   x >= partition_data[partition_y1].x_start && x < partition_data[partition_y1].x_stop
        && y >= partition_data[partition_y1].y_start && y < partition_data[partition_y1].y_stop
        && z >= partition_data[partition_y1].z_start && z < partition_data[partition_y1].z_stop)
        return 1;
    return 0;
}


int in_partition_z0(int x, int y, int z)
{
    if (   x >= partition_data[partition_z0].x_start && x < partition_data[partition_z0].x_stop
        && y >= partition_data[partition_z0].y_start && y < partition_data[partition_z0].y_stop
        && z >= partition_data[partition_z0].z_start && z < partition_data[partition_z0].z_stop)
        return 1;
    return 0;
}


int in_partition_z1(int x, int y, int z)
{
    if (   x >= partition_data[partition_z1].x_start && x < partition_data[partition_z1].x_stop
        && y >= partition_data[partition_z1].y_start && y < partition_data[partition_z1].y_stop
        && z >= partition_data[partition_z1].z_start && z < partition_data[partition_z1].z_stop)
        return 1;
    return 0;
}

int in_partition_main(int x, int y, int z)
{
    if (   x >= partition_data[partition_main].x_start && x < partition_data[partition_main].x_stop
        && y >= partition_data[partition_main].y_start && y < partition_data[partition_main].y_stop
        && z >= partition_data[partition_main].z_start && z < partition_data[partition_main].z_stop)
        return 1;
    return 0;
}


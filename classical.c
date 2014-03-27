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
 *     File Name : step.c
 * Last Modified : Fri 07 Dec 2012 10:30:15 PM EST
 */

#include "tools.h"
#include "h5io.h"
#include "domain.h"
#include "excitation.h"
#include "mem.h"

int x, y, z;
int classical;
double*** ca;
double*** cb_x;
double*** cb_y;
double*** cb_z;

int setup_classical (char *file_name)
{
    double tmp;
    double*** sigma;
    int status;

    status = h5_get_attr(file_name, "materials", "classical_mode", &classical);
    inspect(status, "fail to get h5 attributes");
    if (classical != 1) return 1;

    sigma = (double ***)h5_load3(file_name, "/materials/sigma",  total_x, total_y, total_z);
    inspect(sigma, "fail to load hdf5 dataset");

    ca = (double ***)mem3(type_double, total_x, total_y, total_z);
    if (mode != mode_tex && mode != mode_tmx)
    {
        cb_x = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(cb_x, "fail to allocate memory for cb_x");
    }
    if (mode != mode_tey && mode != mode_tmy)
    {
        cb_y = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(cb_y, "fail to allocate memory for cb_y");
    }
    if (mode != mode_tez && mode != mode_tmz)
    {
        cb_z = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(cb_z, "fail to allocate memory for cb_z");
    }

    for (z = 0; z < total_z; z++) 
        for (y = 0; y < total_y; y++) 
            for (x = 0; x < total_x; x++) 
            {
                tmp = sigma[x][y][z] * d_t / (2 * epsilon[x][y][z]);
                ca[x][y][z] = (1 - tmp) / (1 + tmp);
                if (cb_x)
                    cb_x[x][y][z] = d_tx / epsilon[x][y][z] / (1 + tmp);
                if (cb_y)
                    cb_y[x][y][z] = d_ty / epsilon[x][y][z] / (1 + tmp);
                if (cb_z)
                    cb_z[x][y][z] = d_tz / epsilon[x][y][z] / (1 + tmp);
            }

    free_mem3((void ***)sigma, total_x, total_y);
    return 1;
}

void classical_e()
{
    for (z = 0; z < total_z; z++) 
    {
        for (y = 0; y < total_y; y++) 
        {
            for (x = 0; x < total_x; x++) 
            {
                if (!in_partition_main(x, y, z))
                    continue;

                switch (mode)
                {
                    case 0:
                        ex[x][y][z] = ca[x][y][z] * ex[x][y][z] - d_t * dipole_ex[x][y][z] / epsilon[x][y][z];
                        ey[x][y][z] = ca[x][y][z] * ey[x][y][z] - d_t * dipole_ey[x][y][z] / epsilon[x][y][z];
                        ez[x][y][z] = ca[x][y][z] * ez[x][y][z] - d_t * dipole_ez[x][y][z] / epsilon[x][y][z];
                        if (y != 0)
                            ex[x][y][z] = ex[x][y][z] + cb_y[x][y][z] * (hz[x][y][z] - hz[x][y-1][z]);
                        if (z != 0)
                            ex[x][y][z] = ex[x][y][z] + cb_z[x][y][z] * (hy[x][y][z-1] - hy[x][y][z]);
                        if (z != 0)
                            ey[x][y][z] = ey[x][y][z] + cb_z[x][y][z] * (hx[x][y][z] - hx[x][y][z-1]);
                        if (x != 0)
                            ey[x][y][z] = ey[x][y][z] - cb_x[x][y][z] * (hz[x][y][z] - hz[x-1][y][z]);
                        if (x != 0)
                            ez[x][y][z] = ez[x][y][z] + cb_x[x][y][z] * (hy[x][y][z] - hy[x-1][y][z]);
                        if (y != 0)
                            ez[x][y][z] = ez[x][y][z] - cb_y[x][y][z] * (hx[x][y][z] - hx[x][y-1][z]);
                        break;
                    case 1: ex[x][y][z] = ca[x][y][z] * ex[x][y][z] - d_t * dipole_ex[x][y][z] / epsilon[x][y][z]; 
                        if (y != 0) 
                            ex[x][y][z] = ex[x][y][z] + cb_y[x][y][z] * (hz[x][y][z] - hz[x][y-1][z]); 
                        if (z != 0) 
                            ex[x][y][z] = ex[x][y][z] + cb_z[x][y][z] * (hy[x][y][z] - hy[x][y][z-1]); 
                            break; 
                    case 2: ey[x][y][z] = ca[x][y][z] * ey[x][y][z] - d_t * dipole_ey[x][y][z] / epsilon[x][y][z]; 
                        if (z != 0)
                            ey[x][y][z] = ey[x][y][z] + cb_z[x][y][z] * (hx[x][y][z] - hx[x][y][z-1]);
                        if (x != 0)
                            ey[x][y][z] = ey[x][y][z] - cb_x[x][y][z] * (hz[x][y][z] - hz[x-1][y][z]);
                        break;
                    case 3:
                        ez[x][y][z] = ca[x][y][z] * ez[x][y][z] - d_t * dipole_ez[x][y][z] / epsilon[x][y][z];
                        if (x != 0)
                            ez[x][y][z] = ez[x][y][z] + cb_x[x][y][z] * (hy[x][y][z] - hy[x-1][y][z]);
                        if (y != 0)
                            ez[x][y][z] = ez[x][y][z] - cb_y[x][y][z] * (hx[x][y][z] - hx[x][y-1][z]);
                        break;
                    case 4:
                        ey[x][y][z] = ca[x][y][z] * ey[x][y][z] - d_t * dipole_ey[x][y][z] / epsilon[x][y][z];
                        ez[x][y][z] = ca[x][y][z] * ez[x][y][z] - d_t * dipole_ez[x][y][z] / epsilon[x][y][z];
                        if (z != 0)
                            ey[x][y][z] = ey[x][y][z] + cb_z[x][y][z] * (hx[x][y][z] - hx[x][y][z-1]);
                        if (y != 0)
                            ez[x][y][z] = ez[x][y][z] - cb_y[x][y][z] * (hx[x][y][z] - hx[x][y-1][z]);
                        break;
                    case 5:
                        ex[x][y][z] = ca[x][y][z] * ex[x][y][z] - d_t * dipole_ex[x][y][z] / epsilon[x][y][z];
                        ez[x][y][z] = ca[x][y][z] * ez[x][y][z] - d_t * dipole_ez[x][y][z] / epsilon[x][y][z];
                        if (z != 0)
                            ex[x][y][z] = ex[x][y][z] - cb_z[x][y][z] * (hy[x][y][z] - hy[x][y][z-1]);
                        if (x != 0)
                            ez[x][y][z] = ez[x][y][z] + cb_x[x][y][z] * (hy[x][y][z] - hy[x-1][y][z]);
                        break;
                    case 6:
                        ex[x][y][z] = ca[x][y][z] * ex[x][y][z] - d_t * dipole_ex[x][y][z] / epsilon[x][y][z];
                        ey[x][y][z] = ca[x][y][z] * ey[x][y][z] - d_t * dipole_ey[x][y][z] / epsilon[x][y][z];
                        if (y != 0)
                            ex[x][y][z] = ex[x][y][z] + cb_y[x][y][z] * (hz[x][y][z] - hz[x][y-1][z]);
                        if (x != 0)
                            ey[x][y][z] = ey[x][y][z] - cb_x[x][y][z] * (hz[x][y][z] - hz[x-1][y][z]);
                        break;
                }
            }
        }
    }
}


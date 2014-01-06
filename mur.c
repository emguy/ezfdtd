/* NOTE: This program is free software; you can redistribute it 
 * and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either 
 * version 3, or (at your option) any later version.  *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * Bugs can be reported to Yu Zhang <zhang235@mcmaster.ca>.
 *
 *     File Name : mur.c
 * Last Modified : Sat 04 Jan 2014 02:03:02 PM EST
 */

#include "tools.h"
#include <math.h>
#include "mem.h"
#include "domain.h"
#include "h5io.h"
#include "pml.h"

double **ey_x0;
double **ez_x0;
double **ey_x1;
double **ez_x1;

double **ex_y0;
double **ez_y0;
double **ex_y1;
double **ez_y1;

double **ex_z0;
double **ey_z0;
double **ex_z1;
double **ey_z1;

double **ey_x0_d;
double **ez_x0_d;
double **ey_x1_d;
double **ez_x1_d;
             
double **ex_y0_d;
double **ez_y0_d;
double **ex_y1_d;
double **ez_y1_d;
            
double **ex_z0_d;
double **ey_z0_d;
double **ex_z1_d;
double **ey_z1_d;

double **k1_x0;
double **k2_x0;
double **ky_x0;
double **kz_x0;
double **k1_x1;
double **k2_x1;
double **ky_x1;
double **kz_x1;

double **k1_y0;
double **k2_y0;
double **kx_y0;
double **kz_y0;
double **k1_y1;
double **k2_y1;
double **kx_y1;
double **kz_y1;

double **k1_z0;
double **k2_z0;
double **kx_z0;
double **ky_z0;
double **k1_z1;
double **k2_z1;
double **kx_z1;
double **ky_z1;

int setup_mur (char* file_name)
{
    int x, y, z;
    int c_dt;

    if (partition_data[partition_x0].boundary_type == boundary_mur)
    {
        if (ey)
        {
            ey_x0   = (double **)mem2(type_double, total_length_y, total_length_z);
            ey_x0_d = (double **)mem2(type_double, total_length_y, total_length_z);
            inspect(ey_x0,   "fail to allocate memory for field");
            inspect(ey_x0_d, "fail to allocate memory for field");
        }
        if (ez)
        {
            ez_x0   = (double **)mem2(type_double, total_length_y, total_length_z);
            ez_x0_d = (double **)mem2(type_double, total_length_y, total_length_z);
            inspect(ez_x0,   "fail to allocate memory for field");
            inspect(ez_x0_d, "fail to allocate memory for field");
        }

        k1_x0   = (double **)mem2(type_double, total_length_y, total_length_z);
        k2_x0   = (double **)mem2(type_double, total_length_y, total_length_z);
        ky_x0   = (double **)mem2(type_double, total_length_y, total_length_z);
        kz_x0   = (double **)mem2(type_double, total_length_y, total_length_z);
        inspect(k1_x0, "fail to allocate memory for field");
        inspect(k2_x0, "fail to allocate memory for field");
        inspect(ky_x0, "fail to allocate memory for field");
        inspect(kz_x0, "fail to allocate memory for field");

        for (y = 0; y < total_length_y; y++)
            for (z = 0; z < total_length_z; z++)
            {
                c_dt = d_t / sqrt(epsilon[0][y][z] * MU0);
                k1_x0[y][z] = (c_dt - d_x) / (c_dt + d_x);
                k2_x0[y][z] = 2 * d_x / (c_dt + d_x);
                ky_x0[y][z] = c_dt * c_dt * d_x / (2 * d_y * d_y * (c_dt + d_x));
                kz_x0[y][z] = c_dt * c_dt * d_x / (2 * d_z * d_z * (c_dt + d_x));
            }
    }

    if (partition_data[partition_x1].boundary_type == boundary_mur)
    {
        k1_x1   = (double **)mem2(type_double, total_length_y, total_length_z);
        k2_x1   = (double **)mem2(type_double, total_length_y, total_length_z);
        ky_x1   = (double **)mem2(type_double, total_length_y, total_length_z);
        kz_x1   = (double **)mem2(type_double, total_length_y, total_length_z);
        inspect(k1_x1, "fail to allocate memory for field");
        inspect(k2_x1, "fail to allocate memory for field");
        inspect(ky_x1, "fail to allocate memory for field");
        inspect(kz_x1, "fail to allocate memory for field");

        if (ey)
        {
            ey_x1   = (double **)mem2(type_double, total_length_y, total_length_z);
            ey_x1_d = (double **)mem2(type_double, total_length_y, total_length_z);
            inspect(ey_x1, "fail to allocate memory for field");
            inspect(ey_x1_d, "fail to allocate memory for field");
        }
        if (ez)
        {
            ez_x1   = (double **)mem2(type_double, total_length_y, total_length_z);
            ez_x1_d = (double **)mem2(type_double, total_length_y, total_length_z);
            inspect(ez_x1, "fail to allocate memory for field");
            inspect(ez_x1_d, "fail to allocate memory for field");
        }

        for (y = 0; y < total_length_y; y++)
            for (z = 0; z < total_length_z; z++)
            {
                c_dt = d_t / sqrt(epsilon[total_length_x - 1][y][z] * MU0);
                k1_x1[y][z] = (c_dt - d_x) / (c_dt + d_x);
                k2_x1[y][z] = 2 * d_x / (c_dt + d_x);
                ky_x1[y][z] = c_dt * c_dt * d_x / (2 * d_y * d_y * (c_dt + d_x));
                kz_x1[y][z] = c_dt * c_dt * d_x / (2 * d_z * d_z * (c_dt + d_x));
            }
    }

    if (partition_data[partition_y0].boundary_type == boundary_mur)
    {
        k1_y0   = (double **)mem2(type_double, total_length_x, total_length_z);
        k2_y0   = (double **)mem2(type_double, total_length_x, total_length_z);
        kx_y0   = (double **)mem2(type_double, total_length_x, total_length_z);
        kz_y0   = (double **)mem2(type_double, total_length_x, total_length_z);
        inspect(k1_y0, "fail to allocate memory for field");
        inspect(k2_y0, "fail to allocate memory for field");
        inspect(kx_y0, "fail to allocate memory for field");
        inspect(kz_y0, "fail to allocate memory for field");
        if (ex)
        {
            ex_y0   = (double **)mem2(type_double, total_length_x, total_length_z);
            ex_y0_d = (double **)mem2(type_double, total_length_x, total_length_z);
            inspect(ex_y0, "fail to allocate memory for field");
            inspect(ex_y0_d, "fail to allocate memory for field");
        }
        if (ez)
        {
            ez_y0   = (double **)mem2(type_double, total_length_x, total_length_z);
            ez_y0_d = (double **)mem2(type_double, total_length_x, total_length_z);
            inspect(ez_y0, "fail to allocate memory for field");
            inspect(ez_y0_d, "fail to allocate memory for field");
        }

        for (x = 0; x < total_length_x; x++)
            for (z = 0; z < total_length_z; z++)
            {
                c_dt = d_t / sqrt(epsilon[x][1][z] * MU0);
                k1_y0[x][z] = (c_dt - d_y) / (c_dt + d_y);
                k2_y0[x][z] = 2 * d_y / (c_dt + d_y);
                kx_y0[x][z] = c_dt * c_dt * d_y / (2 * d_x * d_x * (c_dt + d_y));
                kz_y0[x][z] = c_dt * c_dt * d_y / (2 * d_z * d_z * (c_dt + d_y));
            }
    }

    if (partition_data[partition_y1].boundary_type == boundary_mur)
    {
        k1_y1   = (double **)mem2(type_double, total_length_x, total_length_z);
        k2_y1   = (double **)mem2(type_double, total_length_x, total_length_z);
        kx_y1   = (double **)mem2(type_double, total_length_x, total_length_z);
        kz_y1   = (double **)mem2(type_double, total_length_x, total_length_z);
        inspect(k1_y1, "fail to allocate memory for field");
        inspect(k2_y1, "fail to allocate memory for field");
        inspect(kx_y1, "fail to allocate memory for field");
        inspect(kz_y1, "fail to allocate memory for field");
        if (ex)
        {
            ex_y1   = (double **)mem2(type_double, total_length_x, total_length_z);
            ex_y1_d = (double **)mem2(type_double, total_length_x, total_length_z);
            inspect(ex_y1, "fail to allocate memory for field");
            inspect(ex_y1_d, "fail to allocate memory for field");
        }
        if (ez)
        {
            ez_y1   = (double **)mem2(type_double, total_length_x, total_length_z);
            ez_y1_d = (double **)mem2(type_double, total_length_x, total_length_z);
            inspect(ez_y1, "fail to allocate memory for field");
            inspect(ez_y1_d, "fail to allocate memory for field");
        }

        for (x = 0; x < total_length_x; x++)
            for (z = 0; z < total_length_z; z++)
            {
                c_dt = d_t / sqrt(epsilon[x][total_length_y - 1][z] * MU0);
                k1_y1[x][z] = (c_dt - d_y) / (c_dt + d_y);
                k2_y1[x][z] = 2 * d_y / (c_dt + d_y);
                kx_y1[x][z] = c_dt * c_dt * d_y / (2 * d_x * d_x * (c_dt + d_y));
                kz_y1[x][z] = c_dt * c_dt * d_y / (2 * d_z * d_z * (c_dt + d_y));
            }
    }

    if (partition_data[partition_z0].boundary_type == boundary_mur)
    {
        k1_z0   = (double **)mem2(type_double, total_length_x, total_length_y);
        k2_z0   = (double **)mem2(type_double, total_length_x, total_length_y);
        kx_z0   = (double **)mem2(type_double, total_length_x, total_length_y);
        ky_z0   = (double **)mem2(type_double, total_length_x, total_length_y);
        inspect(k1_z0, "fail to allocate memory for field");
        inspect(k2_z0, "fail to allocate memory for field");
        inspect(kx_z0, "fail to allocate memory for field");
        inspect(ky_z0, "fail to allocate memory for field");
        if (ex)
        {
            ex_z0   = (double **)mem2(type_double, total_length_x, total_length_y);
            ex_z0_d = (double **)mem2(type_double, total_length_x, total_length_y);
            inspect(ex_z0, "fail to allocate memory for field");
            inspect(ex_z0_d, "fail to allocate memory for field");
        }
        if (ey)
        {
            ey_z0   = (double **)mem2(type_double, total_length_x, total_length_y);
            ey_z0_d = (double **)mem2(type_double, total_length_x, total_length_y);
            inspect(ey_z0, "fail to allocate memory for field");
            inspect(ey_z0_d, "fail to allocate memory for field");
        }

        for (y = 0; y < total_length_y; y++)
            for (z = 0; z < total_length_z; z++)
            {
                c_dt = d_t / sqrt(epsilon[1][y][z] * MU0);
                k1_z0[y][z] = (c_dt - d_z) / (c_dt + d_z);
                k2_z0[y][z] = 2 * d_z / (c_dt + d_z);
                kx_z0[y][z] = c_dt * c_dt * d_z / (2 * d_x * d_x * (c_dt + d_z));
                ky_z0[y][z] = c_dt * c_dt * d_z / (2 * d_y * d_y * (c_dt + d_z));
            }
    }

    if (partition_data[partition_z1].boundary_type == boundary_mur)
    {
        k1_z1   = (double **)mem2(type_double, total_length_x, total_length_y);
        k2_z1   = (double **)mem2(type_double, total_length_x, total_length_y);
        kx_z1   = (double **)mem2(type_double, total_length_x, total_length_y);
        ky_z1   = (double **)mem2(type_double, total_length_x, total_length_y);
        inspect(k1_z1, "fail to allocate memory for field");
        inspect(k2_z1, "fail to allocate memory for field");
        inspect(kx_z1, "fail to allocate memory for field");
        inspect(ky_z1, "fail to allocate memory for field");
        if (ex)
        {
            ex_z1   = (double **)mem2(type_double, total_length_x, total_length_y);
            ex_z1_d = (double **)mem2(type_double, total_length_x, total_length_y);
            inspect(ex_z1, "fail to allocate memory for field");
            inspect(ex_z1_d, "fail to allocate memory for field");
        }
        if (ey)
        {
            ey_z1   = (double **)mem2(type_double, total_length_x, total_length_y);
            ey_z1_d = (double **)mem2(type_double, total_length_x, total_length_y);
            inspect(ey_z1, "fail to allocate memory for field");
            inspect(ey_z1_d, "fail to allocate memory for field");
        }

        for (y = 0; y < total_length_y; y++)
            for (z = 0; z < total_length_z; z++)
            {
                c_dt = d_t / sqrt(epsilon[total_length_x - 1][y][z] * MU0);
                k1_z1[y][z] = (c_dt - d_z) / (c_dt + d_z);
                k2_z1[y][z] = 2 * d_z / (c_dt + d_z);
                kx_z1[y][z] = c_dt * c_dt * d_z / (2 * d_x * d_x * (c_dt + d_z));
                ky_z1[y][z] = c_dt * c_dt * d_z / (2 * d_y * d_y * (c_dt + d_z));
            }
    }

    return 1;
}
    





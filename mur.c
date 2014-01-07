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

static void update_ey_x0 ();
static void update_ez_x0 ();
static void update_ey_x1 ();
static void update_ez_x1 ();

static void update_ex_y0 ();
static void update_ez_y0 ();
static void update_ex_y1 ();
static void update_ez_y1 ();

static void update_ex_z0 ();
static void update_ey_z0 ();
static void update_ex_z1 ();
static void update_ey_z1 ();

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
    
void update_mur ()
{
    if (partition_data[partition_x0].boundary_type == boundary_mur)
    {
        if (mode == mode_full && mode == mode_tmy && mode == mode_tez) update_ey_x0();
        if (mode == mode_full && mode == mode_tmz && mode == mode_tey) update_ez_x0();
    }
    if (partition_data[partition_x1].boundary_type == boundary_mur)
    {
        if (mode == mode_full && mode == mode_tmy && mode == mode_tez) update_ey_x1();
        if (mode == mode_full && mode == mode_tmz && mode == mode_tey) update_ez_x1();
    }
    if (partition_data[partition_y0].boundary_type == boundary_mur)
    {
        if (mode == mode_full && mode == mode_tmx && mode == mode_tez) update_ex_y0();
        if (mode == mode_full && mode == mode_tmz && mode == mode_tex) update_ez_y0();
    }
    if (partition_data[partition_y1].boundary_type == boundary_mur)
    {
        if (mode == mode_full && mode == mode_tmx && mode == mode_tez) update_ex_y1();
        if (mode == mode_full && mode == mode_tmz && mode == mode_tex) update_ez_y1();
    }
    if (partition_data[partition_z0].boundary_type == boundary_mur)
    {
        if (mode == mode_full && mode == mode_tmx && mode == mode_tey) update_ex_z0();
        if (mode == mode_full && mode == mode_tmy && mode == mode_tex) update_ey_z0();
    }
    if (partition_data[partition_z1].boundary_type == boundary_mur)
    {
        if (mode == mode_full && mode == mode_tmx && mode == mode_tey) update_ex_z1();
        if (mode == mode_full && mode == mode_tmy && mode == mode_tex) update_ey_z1();
    }
}

static void update_ey_x0 ()
{
    double temp;
    int y, z;

    for (y = 0; y < total_length_y; y++)
        for (z = 0; z < total_length_z; z++)
        {
            temp = ey[0][y][z];

            ey[0][y][z] = ey_x0[y][z] - ey_x0_d[y][z] + k1_x0[y][z] * (ey[1][y][z] + ey_x0_d[y][z]);

            ey_x0_d[y][z] = temp;

            ey_x0[y][z] = k2_x0[y][z] * (ey[0][y][z] + ey[1][y][z]);
            if (y > 0)
            {
                ey_x0[y][z] += ky_x0[y][z] * (ey[0][y-1][z] - ey[0][y][z]);
                ey_x0[y][z] += ky_x0[y][z] * (ey[1][y-1][z] - ey[1][y][z]);
            }
            if (y < total_length_y)
            {
                ey_x0[y][z] += ky_x0[y][z] * (ey[0][y+1][z] - ey[0][y][z]);
                ey_x0[y][z] += ky_x0[y][z] * (ey[1][y+1][z] - ey[1][y][z]);
            }
            if (z > 0)
            {
                ey_x0[y][z] += kz_x0[y][z] * (ey[0][y][z-1] - ey[0][y][z]);
                ey_x0[y][z] += kz_x0[y][z] * (ey[1][y][z-1] - ey[1][y][z]);
            }
            if (z < total_length_z)
            {
                ey_x0[y][z] += kz_x0[y][z] * (ey[0][y][z+1] - ey[0][y][z]);
                ey_x0[y][z] += kz_x0[y][z] * (ey[1][y][z+1] - ey[1][y][z]);
            }
        }
}

static void update_ez_x0 ()
{
    double temp;
    int y, z;

    for (y = 0; y < total_length_y; y++)
        for (z = 0; z < total_length_z; z++)
        {
            temp = ez[0][y][z];

            ez[0][y][z] = ez_x0[y][z] - ez_x0_d[y][z] + k1_x0[y][z] * (ez[1][y][z] + ez_x0_d[y][z]);

            ez_x0_d[y][z] = temp;

            ez_x0[y][z] = k2_x0[y][z] * (ez[0][y][z] + ez[1][y][z]);
            if (y > 0)
            {
                ez_x0[y][z] += ky_x0[y][z] * (ez[0][y-1][z] - ez[0][y][z]);
                ez_x0[y][z] += ky_x0[y][z] * (ez[1][y-1][z] - ez[1][y][z]);
            }
            if (y < total_length_y)
            {
                ez_x0[y][z] += ky_x0[y][z] * (ez[0][y+1][z] - ez[0][y][z]);
                ez_x0[y][z] += ky_x0[y][z] * (ez[1][y+1][z] - ez[1][y][z]);
            }
            if (z > 0)
            {
                ez_x0[y][z] += kz_x0[y][z] * (ez[0][y][z-1] - ez[0][y][z]);
                ez_x0[y][z] += kz_x0[y][z] * (ez[1][y][z-1] - ez[1][y][z]);
            }
            if (z < total_length_z)
            {
                ez_x0[y][z] += kz_x0[y][z] * (ez[0][y][z+1] - ez[0][y][z]);
                ez_x0[y][z] += kz_x0[y][z] * (ez[1][y][z+1] - ez[1][y][z]);
            }
        }
}

static void update_ey_x1 ()
{
    double temp;
    int y, z;

    for (y = 0; y < total_length_y; y++)
        for (z = 0; z < total_length_z; z++)
        {
            temp = ey[total_length_x][y][z];

            ey[total_length_x][y][z] = ey_x1[y][z] - ey_x1_d[y][z] + k1_x1[y][z] * (ey[total_length_x - 1][y][z] + ey_x1_d[y][z]);

            ey_x1_d[y][z] = temp;

            ey_x1[y][z] = k2_x1[y][z] * (ey[total_length_x][y][z] + ey[total_length_x - 1][y][z]);
            if (y > 0)
            {
                ey_x1[y][z] += ky_x1[y][z] * (ey[total_length_x][y-1][z] - ey[total_length_x][y][z]);
                ey_x1[y][z] += ky_x1[y][z] * (ey[total_length_x - 1][y-1][z] - ey[total_length_x - 1][y][z]);
            }
            if (y < total_length_y)
            {
                ey_x1[y][z] += ky_x1[y][z] * (ey[total_length_x][y+1][z] - ey[total_length_x][y][z]);
                ey_x1[y][z] += ky_x1[y][z] * (ey[total_length_x - 1][y+1][z] - ey[total_length_x - 1][y][z]);
            }
            if (z > 0)
            {
                ey_x1[y][z] += kz_x1[y][z] * (ey[total_length_x][y][z-1] - ey[total_length_x][y][z]);
                ey_x1[y][z] += kz_x1[y][z] * (ey[total_length_x - 1][y][z-1] - ey[total_length_x - 1][y][z]);
            }
            if (z < total_length_z)
            {
                ey_x1[y][z] += kz_x1[y][z] * (ey[total_length_x][y][z+1] - ey[total_length_x][y][z]);
                ey_x1[y][z] += kz_x1[y][z] * (ey[total_length_x - 1][y][z+1] - ey[total_length_x - 1][y][z]);
            }
        }
}

static void update_ez_x1 ()
{
    double temp;
    int y, z;

    for (y = 0; y < total_length_y; y++)
        for (z = 0; z < total_length_z; z++)
        {
            temp = ey[total_length_x][y][z];

            ez[total_length_x][y][z] = ez_x1[y][z] - ez_x1_d[y][z] + k1_x1[y][z] * (ez[total_length_x - 1][y][z] + ez_x1_d[y][z]);

            ez_x1_d[y][z] = temp;

            ez_x1[y][z] = k2_x1[y][z] * (ez[total_length_x][y][z] + ez[total_length_x - 1][y][z]);
            if (y > 0)
            {
                ez_x1[y][z] += ky_x1[y][z] * (ez[total_length_x][y-1][z] - ez[total_length_x][y][z]);
                ez_x1[y][z] += ky_x1[y][z] * (ez[total_length_x - 1][y-1][z] - ez[total_length_x - 1][y][z]);
            }
            if (y < total_length_y)
            {
                ez_x1[y][z] += ky_x1[y][z] * (ez[total_length_x][y+1][z] - ez[total_length_x][y][z]);
                ez_x1[y][z] += ky_x1[y][z] * (ez[total_length_x - 1][y+1][z] - ez[total_length_x - 1][y][z]);
            }
            if (z > 0)
            {
                ez_x1[y][z] += kz_x1[y][z] * (ez[total_length_x][y][z-1] - ez[total_length_x][y][z]);
                ez_x1[y][z] += kz_x1[y][z] * (ez[total_length_x - 1][y][z-1] - ez[total_length_x - 1][y][z]);
            }
            if (z < total_length_z)
            {
                ez_x1[y][z] += kz_x1[y][z] * (ez[total_length_x][y][z+1] - ez[total_length_x][y][z]);
                ez_x1[y][z] += kz_x1[y][z] * (ez[total_length_x - 1][y][z+1] - ez[total_length_x - 1][y][z]);
            }
        }
}

static void update_ex_y0 ()
{
    double temp;
    int x, z;

    for (x = 0; x < total_length_x; x++)
        for (z = 0; z < total_length_z; z++)
        {
            temp = ex[x][0][z];

            ex[x][0][z] = ex_y0[x][z] - ex_y0_d[x][z] + k1_y0[x][z] * (ex[x][1][z] + ex_y0_d[x][z]);

            ex_y0_d[x][z] = temp;

            ex_y0[x][z] = k2_y0[x][z] * (ex[x][0][z] + ex[x][1][z]);
            if (x > 0)
            {
                ex_y0[x][z] += kx_y0[x][z] * (ex[x-1][0][z] - ex[x][0][z]);
                ex_y0[x][z] += kx_y0[x][z] * (ex[x-1][1][z] - ex[x][1][z]);
            }
            if (x < total_length_x)
            {
                ex_y0[x][z] += kx_y0[x][z] * (ex[x+1][0][z] - ex[x][0][z]);
                ex_y0[x][z] += kx_y0[x][z] * (ex[x+1][1][z] - ex[x][1][z]);
            }
            if (z > 0)
            {
                ex_y0[x][z] += kz_y0[x][z] * (ex[x][0][z-1] - ex[x][0][z]);
                ex_y0[x][z] += kz_y0[x][z] * (ex[x][1][z-1] - ex[x][1][z]);
            }
            if (z < total_length_z)
            {
                ex_y0[x][z] += kz_y0[x][z] * (ex[x][0][z+1] - ex[x][0][z]);
                ex_y0[x][z] += kz_y0[x][z] * (ex[x][1][z+1] - ex[x][1][z]);
            }
        }
}

static void update_ez_y0 ()
{
    double temp;
    int x, z;

    for (x = 0; x < total_length_x; x++)
        for (z = 0; z < total_length_z; z++)
        {
            temp = ez[x][0][z];

            ez[x][0][z] = ez_y0[x][z] - ez_y0_d[x][z] + k1_y0[x][z] * (ez[x][1][z] + ez_y0_d[x][z]);

            ez_y0_d[x][z] = temp;

            ez_y0[x][z] = k2_y0[x][z] * (ez[x][0][z] + ez[x][1][z]);
            if (x > 0)
            {
                ez_y0[x][z] += kx_y0[x][z] * (ez[x-1][0][z] - ez[x][0][z]);
                ez_y0[x][z] += kx_y0[x][z] * (ez[x-1][1][z] - ez[x][1][z]);
            }
            if (x < total_length_x)
            {
                ez_y0[x][z] += kx_y0[x][z] * (ez[x+1][0][z] - ez[x][0][z]);
                ez_y0[x][z] += kx_y0[x][z] * (ez[x+1][1][z] - ez[x][1][z]);
            }
            if (z > 0)
            {
                ez_y0[x][z] += kz_y0[x][z] * (ez[x][0][z-1] - ez[x][0][z]);
                ez_y0[x][z] += kz_y0[x][z] * (ez[x][1][z-1] - ez[x][1][z]);
            }
            if (z < total_length_z)
            {
                ez_y0[x][z] += kz_y0[x][z] * (ez[x][0][z+1] - ez[x][0][z]);
                ez_y0[x][z] += kz_y0[x][z] * (ez[x][1][z+1] - ez[x][1][z]);
            }
        }
}

static void update_ex_y1 ()
{
    double temp;
    int x, z;

    for (x = 0; x < total_length_x; x++)
        for (z = 0; z < total_length_z; z++)
        {
            temp = ex[x][total_length_y][z];

            ex[x][total_length_y][z] = ex_y1[x][z] - ex_y1_d[x][z] + k1_y1[x][z] * (ex[x][total_length_y - 1][z] + ex_y1_d[x][z]);

            ex_y1_d[x][z] = temp;

            ex_y1[x][z] = k2_y1[x][z] * (ex[x][total_length_y][z] + ex[x][total_length_y - 1][z]);
            if (x > 0)
            {
                ex_y1[x][z] += kx_y1[x][z] * (ex[x-1][total_length_y][z] - ex[x][total_length_y][z]);
                ex_y1[x][z] += kx_y1[x][z] * (ex[x-1][total_length_y - 1][z] - ex[x][total_length_y - 1][z]);
            }
            if (x < total_length_x)
            {
                ex_y1[x][z] += kx_y1[x][z] * (ex[x+1][total_length_y][z] - ex[x][total_length_y][z]);
                ex_y1[x][z] += kx_y1[x][z] * (ex[x+1][total_length_y - 1][z] - ex[x][total_length_y - 1][z]);
            }
            if (z > 0)
            {
                ex_y1[x][z] += kz_y1[x][z] * (ex[x][total_length_y][z-1] - ex[x][total_length_y][z]);
                ex_y1[x][z] += kz_y1[x][z] * (ex[x][total_length_y - 1][z-1] - ex[x][total_length_y - 1][z]);
            }
            if (z < total_length_z)
            {
                ex_y1[x][z] += kz_y1[x][z] * (ex[x][total_length_y][z+1] - ex[x][total_length_y][z]);
                ex_y1[x][z] += kz_y1[x][z] * (ex[x][total_length_y - 1][z+1] - ex[x][total_length_y - 1][z]);
            }
        }
}

static void update_ez_y1 ()
{
    double temp;
    int x, z;

    for (x = 0; x < total_length_x; x++)
        for (z = 0; z < total_length_z; z++)
        {
            temp = ez[x][total_length_y][z];

            ez[x][total_length_y][z] = ez_y1[x][z] - ez_y1_d[x][z] + k1_y1[x][z] * (ez[x][total_length_y - 1][z] + ez_y1_d[x][z]);

            ez_y1_d[x][z] = temp;

            ez_y1[x][z] = k2_y1[x][z] * (ez[x][total_length_y][z] + ez[x][total_length_y - 1][z]);
            if (x > 0)
            {
                ez_y1[x][z] += kx_y1[x][z] * (ez[x-1][total_length_y][z] - ez[x][total_length_y][z]);
                ez_y1[x][z] += kx_y1[x][z] * (ez[x-1][total_length_y - 1][z] - ez[x][total_length_y - 1][z]);
            }
            if (x < total_length_x)
            {
                ez_y1[x][z] += kx_y1[x][z] * (ez[x+1][total_length_y][z] - ez[x][total_length_y][z]);
                ez_y1[x][z] += kx_y1[x][z] * (ez[x+1][total_length_y - 1][z] - ez[x][total_length_y - 1][z]);
            }
            if (z > 0)
            {
                ez_y1[x][z] += kz_y1[x][z] * (ez[x][total_length_y][z-1] - ez[x][total_length_y][z]);
                ez_y1[x][z] += kz_y1[x][z] * (ez[x][total_length_y - 1][z-1] - ez[x][total_length_y - 1][z]);
            }
            if (z < total_length_z)
            {
                ez_y1[x][z] += kz_y1[x][z] * (ez[x][total_length_y][z+1] - ez[x][total_length_y][z]);
                ez_y1[x][z] += kz_y1[x][z] * (ez[x][total_length_y - 1][z+1] - ez[x][total_length_y - 1][z]);
            }
        }
}

static void update_ex_z0 ()
{
    double temp;
    int x, y;

    for (x = 0; x < total_length_x; x++)
        for (y = 0; y < total_length_y; y++)
        {
            temp = ex[x][y][0];

            ex[x][y][0] = ex_z0[x][y] - ex_z0_d[x][y] + k1_z0[x][y] * (ex[x][y][1] + ex_z0_d[x][y]);

            ex_z0_d[x][y] = temp;

            ex_z0[x][y] = k2_z0[x][y] * (ex[x][y][0] + ex[x][y][1]);
            if (x > 0)
            {
                ex_z0[x][y] += kx_z0[x][y] * (ex[x-1][y][0] - ex[x][y][0]);
                ex_z0[x][y] += kx_z0[x][y] * (ex[x-1][y][1] - ex[x][y][1]);
            }
            if (x < total_length_x)
            {
                ex_z0[x][y] += kx_z0[x][y] * (ex[x+1][y][0] - ex[x][y][0]);
                ex_z0[x][y] += kx_z0[x][y] * (ex[x+1][y][1] - ex[x][y][1]);
            }
            if (y > 0)
            {
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y-1][0] - ex[x][y][0]);
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y-1][1] - ex[x][y][1]);
            }
            if (y < total_length_y)
            {
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y+1][0] - ex[x][y][0]);
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y+1][1] - ex[x][y][1]);
            }
        }
}

static void update_ey_z0 ()
{
    double temp;
    int x, y;

    for (x = 0; x < total_length_x; x++)
        for (y = 0; y < total_length_y; y++)
        {
            temp = ey[x][y][0];

            ey[x][y][0] = ey_z0[x][y] - ey_z0_d[x][y] + k1_z0[x][y] * (ey[x][y][1] + ey_z0_d[x][y]);

            ey_z0_d[x][y] = temp;

            ey_z0[x][y] = k2_z0[x][y] * (ey[x][y][0] + ey[x][y][1]);
            if (x > 0)
            {
                ey_z0[x][y] += kx_z0[x][y] * (ey[x-1][y][0] - ey[x][y][0]);
                ey_z0[x][y] += kx_z0[x][y] * (ey[x-1][y][1] - ey[x][y][1]);
            }
            if (x < total_length_x)
            {
                ey_z0[x][y] += kx_z0[x][y] * (ey[x+1][y][0] - ey[x][y][0]);
                ey_z0[x][y] += kx_z0[x][y] * (ey[x+1][y][1] - ey[x][y][1]);
            }
            if (y > 0)
            {
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y-1][0] - ey[x][y][0]);
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y-1][1] - ey[x][y][1]);
            }
            if (y < total_length_y)
            {
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y+1][0] - ey[x][y][0]);
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y+1][1] - ey[x][y][1]);
            }
        }
}

static void update_ex_z1 ()
{
    double temp;
    int x, y;

    for (x = 0; x < total_length_x; x++)
        for (y = 0; y < total_length_y; y++)
        {
            temp = ex[x][y][total_length_z];

            ex[x][y][total_length_z] = ex_z0[x][y] - ex_z0_d[x][y] + k1_z0[x][y] * (ex[x][y][total_length_z - 1] + ex_z0_d[x][y]);

            ex_z0_d[x][y] = temp;

            ex_z0[x][y] = k2_z0[x][y] * (ex[x][y][total_length_z] + ex[x][y][total_length_z - 1]);
            if (x > 0)
            {
                ex_z0[x][y] += kx_z0[x][y] * (ex[x-1][y][total_length_z] - ex[x][y][total_length_z]);
                ex_z0[x][y] += kx_z0[x][y] * (ex[x-1][y][total_length_z - 1] - ex[x][y][total_length_z - 1]);
            }
            if (x < total_length_x)
            {
                ex_z0[x][y] += kx_z0[x][y] * (ex[x+1][y][total_length_z] - ex[x][y][total_length_z]);
                ex_z0[x][y] += kx_z0[x][y] * (ex[x+1][y][total_length_z - 1] - ex[x][y][total_length_z - 1]);
            }
            if (y > 0)
            {
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y-1][total_length_z] - ex[x][y][total_length_z]);
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y-1][total_length_z - 1] - ex[x][y][total_length_z - 1]);
            }
            if (y < total_length_y)
            {
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y+1][total_length_z] - ex[x][y][total_length_z]);
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y+1][total_length_z - 1] - ex[x][y][total_length_z - 1]);
            }
        }
}

static void update_ey_z1 ()
{
    double temp;
    int x, y;

    for (x = 0; x < total_length_x; x++)
        for (y = 0; y < total_length_y; y++)
        {
            temp = ey[x][y][total_length_z];

            ey[x][y][total_length_z] = ey_z0[x][y] - ey_z0_d[x][y] + k1_z0[x][y] * (ey[x][y][total_length_z - 1] + ey_z0_d[x][y]);

            ey_z0_d[x][y] = temp;

            ey_z0[x][y] = k2_z0[x][y] * (ey[x][y][total_length_z] + ey[x][y][total_length_z - 1]);
            if (x > 0)
            {
                ey_z0[x][y] += kx_z0[x][y] * (ey[x-1][y][total_length_z] - ey[x][y][total_length_z]);
                ey_z0[x][y] += kx_z0[x][y] * (ey[x-1][y][total_length_z - 1] - ey[x][y][total_length_z - 1]);
            }
            if (x < total_length_x)
            {
                ey_z0[x][y] += kx_z0[x][y] * (ey[x+1][y][total_length_z] - ey[x][y][total_length_z]);
                ey_z0[x][y] += kx_z0[x][y] * (ey[x+1][y][total_length_z - 1] - ey[x][y][total_length_z - 1]);
            }
            if (y > 0)
            {
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y-1][total_length_z] - ey[x][y][total_length_z]);
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y-1][total_length_z - 1] - ey[x][y][total_length_z - 1]);
            }
            if (y < total_length_y)
            {
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y+1][total_length_z] - ey[x][y][total_length_z]);
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y+1][total_length_z - 1] - ey[x][y][total_length_z - 1]);
            }
        }
}

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

double **ey_x0_dd;
double **ez_x0_dd;
double **ey_x1_dd;
double **ez_x1_dd;

double **ex_y0_dd;
double **ez_y0_dd;
double **ex_y1_dd;
double **ez_y1_dd;

double **ex_z0_dd;
double **ey_z0_dd;
double **ex_z1_dd;
double **ey_z1_dd;

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

int x, y, z;

int setup_mur ()
{
    double c_dt;

    if (partition_data[partition_x0].boundary_type == boundary_mur)
    {
        if (ey)
        {
            ey_x0    = (double **)mem2(type_double, total_y + 1, total_z + 1);
            ey_x0_d  = (double **)mem2(type_double, total_y + 1, total_z + 1);
            ey_x0_dd = (double **)mem2(type_double, total_y + 1, total_z + 1);
            inspect(ey_x0,   "fail to allocate memory for field");
            inspect(ey_x0_d, "fail to allocate memory for field");
            inspect(ey_x0_dd, "fail to allocate memory for field");
        }
        if (ez)
        {
            ez_x0    = (double **)mem2(type_double, total_y + 1, total_z + 1);
            ez_x0_d  = (double **)mem2(type_double, total_y + 1, total_z + 1);
            ez_x0_dd = (double **)mem2(type_double, total_y + 1, total_z + 1);
            inspect(ez_x0,   "fail to allocate memory for field");
            inspect(ez_x0_d, "fail to allocate memory for field");
            inspect(ez_x0_dd, "fail to allocate memory for field");
        }

        k1_x0   = (double **)mem2(type_double, total_y, total_z);
        k2_x0   = (double **)mem2(type_double, total_y, total_z);
        ky_x0   = (double **)mem2(type_double, total_y, total_z);
        kz_x0   = (double **)mem2(type_double, total_y, total_z);
        inspect(k1_x0, "fail to allocate memory for field");
        inspect(k2_x0, "fail to allocate memory for field");
        inspect(ky_x0, "fail to allocate memory for field");
        inspect(kz_x0, "fail to allocate memory for field");

        for (z = 0; z < total_z; z++)
            for (y = 0; y < total_y; y++)
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
        k1_x1   = (double **)mem2(type_double, total_y, total_z);
        k2_x1   = (double **)mem2(type_double, total_y, total_z);
        ky_x1   = (double **)mem2(type_double, total_y, total_z);
        kz_x1   = (double **)mem2(type_double, total_y, total_z);
        inspect(k1_x1, "fail to allocate memory for field");
        inspect(k2_x1, "fail to allocate memory for field");
        inspect(ky_x1, "fail to allocate memory for field");
        inspect(kz_x1, "fail to allocate memory for field");

        if (ey)
        {
            ey_x1    = (double **)mem2(type_double, total_y + 1, total_z + 1);
            ey_x1_d  = (double **)mem2(type_double, total_y + 1, total_z + 1);
            ey_x1_dd = (double **)mem2(type_double, total_y + 1, total_z + 1);
            inspect(ey_x1, "fail to allocate memory for field");
            inspect(ey_x1_d, "fail to allocate memory for field");
            inspect(ey_x1_dd, "fail to allocate memory for field");
        }
        if (ez)
        {
            ez_x1    = (double **)mem2(type_double, total_y + 1, total_z + 1);
            ez_x1_d  = (double **)mem2(type_double, total_y + 1, total_z + 1);
            ez_x1_dd = (double **)mem2(type_double, total_y + 1, total_z + 1);
            inspect(ez_x1, "fail to allocate memory for field");
            inspect(ez_x1_d, "fail to allocate memory for field");
            inspect(ez_x1_dd, "fail to allocate memory for field");
        }

        for (y = 0; y < total_y; y++)
            for (z = 0; z < total_z; z++)
            {
                c_dt = d_t / sqrt(epsilon[total_x - 1][y][z] * MU0);
                k1_x1[y][z] = (c_dt - d_x) / (c_dt + d_x);
                k2_x1[y][z] = 2 * d_x / (c_dt + d_x);
                ky_x1[y][z] = c_dt * c_dt * d_x / (2 * d_y * d_y * (c_dt + d_x));
                kz_x1[y][z] = c_dt * c_dt * d_x / (2 * d_z * d_z * (c_dt + d_x));
            }
    }

    if (partition_data[partition_y0].boundary_type == boundary_mur)
    {
        k1_y0   = (double **)mem2(type_double, total_x, total_z);
        k2_y0   = (double **)mem2(type_double, total_x, total_z);
        kx_y0   = (double **)mem2(type_double, total_x, total_z);
        kz_y0   = (double **)mem2(type_double, total_x, total_z);
        inspect(k1_y0, "fail to allocate memory for field");
        inspect(k2_y0, "fail to allocate memory for field");
        inspect(kx_y0, "fail to allocate memory for field");
        inspect(kz_y0, "fail to allocate memory for field");
        if (ex)
        {
            ex_y0    = (double **)mem2(type_double, total_x, total_z + 1);
            ex_y0_d  = (double **)mem2(type_double, total_x, total_z + 1);
            ex_y0_dd = (double **)mem2(type_double, total_x, total_z + 1);
            inspect(ex_y0, "fail to allocate memory for field");
            inspect(ex_y0_d, "fail to allocate memory for field");
            inspect(ex_y0_dd, "fail to allocate memory for field");
        }
        if (ez)
        {
            ez_y0    = (double **)mem2(type_double, total_x, total_z + 1);
            ez_y0_d  = (double **)mem2(type_double, total_x, total_z + 1);
            ez_y0_dd = (double **)mem2(type_double, total_x, total_z + 1);
            inspect(ez_y0, "fail to allocate memory for field");
            inspect(ez_y0_d, "fail to allocate memory for field");
            inspect(ez_y0_dd, "fail to allocate memory for field");
        }

        for (x = 0; x < total_x; x++)
            for (z = 0; z < total_z; z++)
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
        k1_y1   = (double **)mem2(type_double, total_x, total_z);
        k2_y1   = (double **)mem2(type_double, total_x, total_z);
        kx_y1   = (double **)mem2(type_double, total_x, total_z);
        kz_y1   = (double **)mem2(type_double, total_x, total_z);
        inspect(k1_y1, "fail to allocate memory for field");
        inspect(k2_y1, "fail to allocate memory for field");
        inspect(kx_y1, "fail to allocate memory for field");
        inspect(kz_y1, "fail to allocate memory for field");
        if (ex)
        {
            ex_y1    = (double **)mem2(type_double, total_x + 1, total_z + 1);
            ex_y1_d  = (double **)mem2(type_double, total_x + 1, total_z + 1);
            ex_y1_dd = (double **)mem2(type_double, total_x + 1, total_z + 1);
            inspect(ex_y1, "fail to allocate memory for field");
            inspect(ex_y1_d, "fail to allocate memory for field");
            inspect(ex_y1_dd, "fail to allocate memory for field");
        }
        if (ez)
        {
            ez_y1    = (double **)mem2(type_double, total_x + 1, total_z + 1);
            ez_y1_d  = (double **)mem2(type_double, total_x + 1, total_z + 1);
            ez_y1_dd = (double **)mem2(type_double, total_x + 1, total_z + 1);
            inspect(ez_y1, "fail to allocate memory for field");
            inspect(ez_y1_d, "fail to allocate memory for field");
            inspect(ez_y1_dd, "fail to allocate memory for field");
        }

        for (x = 0; x < total_x; x++)
            for (z = 0; z < total_z; z++)
            {
                c_dt = d_t / sqrt(epsilon[x][total_y - 1][z] * MU0);
                k1_y1[x][z] = (c_dt - d_y) / (c_dt + d_y);
                k2_y1[x][z] = 2 * d_y / (c_dt + d_y);
                kx_y1[x][z] = c_dt * c_dt * d_y / (2 * d_x * d_x * (c_dt + d_y));
                kz_y1[x][z] = c_dt * c_dt * d_y / (2 * d_z * d_z * (c_dt + d_y));
            }
    }

    if (partition_data[partition_z0].boundary_type == boundary_mur)
    {
        k1_z0   = (double **)mem2(type_double, total_x, total_y);
        k2_z0   = (double **)mem2(type_double, total_x, total_y);
        kx_z0   = (double **)mem2(type_double, total_x, total_y);
        ky_z0   = (double **)mem2(type_double, total_x, total_y);
        inspect(k1_z0, "fail to allocate memory for field");
        inspect(k2_z0, "fail to allocate memory for field");
        inspect(kx_z0, "fail to allocate memory for field");
        inspect(ky_z0, "fail to allocate memory for field");
        if (ex)
        {
            ex_z0    = (double **)mem2(type_double, total_x + 1, total_y + 1);
            ex_z0_d  = (double **)mem2(type_double, total_x + 1, total_y + 1);
            ex_z0_dd = (double **)mem2(type_double, total_x + 1, total_y + 1);
            inspect(ex_z0, "fail to allocate memory for field");
            inspect(ex_z0_d, "fail to allocate memory for field");
            inspect(ex_z0_dd, "fail to allocate memory for field");
        }
        if (ey)
        {
            ey_z0    = (double **)mem2(type_double, total_x + 1, total_y + 1);
            ey_z0_d  = (double **)mem2(type_double, total_x + 1, total_y + 1);
            ey_z0_dd = (double **)mem2(type_double, total_x + 1, total_y + 1);
            inspect(ey_z0, "fail to allocate memory for field");
            inspect(ey_z0_d, "fail to allocate memory for field");
            inspect(ey_z0_dd, "fail to allocate memory for field");
        }
        for (x = 0; x < total_x; x++)
            for (y = 0; y < total_y; y++)
            {
                c_dt = d_t / sqrt(epsilon[x][y][0] * MU0);
                k1_z0[x][y] = (c_dt - d_z) / (c_dt + d_z);
                k2_z0[x][y] = 2 * d_z / (c_dt + d_z);
                kx_z0[x][y] = c_dt * c_dt * d_z / (2 * d_x * d_x * (c_dt + d_z));
                ky_z0[x][y] = c_dt * c_dt * d_z / (2 * d_y * d_y * (c_dt + d_z));
                //printf("(%d %d) k1: %e k2:%e  \n",x, y, k1_z0[x][y], k2_z0[x][y]);
            }
    }

    if (partition_data[partition_z1].boundary_type == boundary_mur)
    {
        k1_z1   = (double **)mem2(type_double, total_x, total_y);
        k2_z1   = (double **)mem2(type_double, total_x, total_y);
        kx_z1   = (double **)mem2(type_double, total_x, total_y);
        ky_z1   = (double **)mem2(type_double, total_x, total_y);
        inspect(k1_z1, "fail to allocate memory for field");
        inspect(k2_z1, "fail to allocate memory for field");
        inspect(kx_z1, "fail to allocate memory for field");
        inspect(ky_z1, "fail to allocate memory for field");
        if (ex)
        {
            ex_z1    = (double **)mem2(type_double, total_x + 1, total_y + 1);
            ex_z1_d  = (double **)mem2(type_double, total_x + 1, total_y + 1);
            ex_z1_dd = (double **)mem2(type_double, total_x + 1, total_y + 1);
            inspect(ex_z1, "fail to allocate memory for field");
            inspect(ex_z1_d, "fail to allocate memory for field");
            inspect(ex_z1_dd, "fail to allocate memory for field");
        }
        if (ey)
        {
            ey_z1    = (double **)mem2(type_double, total_x + 1, total_y + 1);
            ey_z1_d  = (double **)mem2(type_double, total_x + 1, total_y + 1);
            ey_z1_dd = (double **)mem2(type_double, total_x + 1, total_y + 1);
            inspect(ey_z1, "fail to allocate memory for field");
            inspect(ey_z1_d, "fail to allocate memory for field");
            inspect(ey_z1_dd, "fail to allocate memory for field");
        }

        for (x = 0; x < total_x; x++)
            for (y = 0; y < total_y; y++)
            {
                c_dt = d_t / sqrt(epsilon[x][y][total_z - 1] * MU0);
                k1_z1[x][y] = (c_dt - d_z) / (c_dt + d_z);
                k2_z1[x][y] = 2 * d_z / (c_dt + d_z);
                kx_z1[x][y] = c_dt * c_dt * d_z / (2 * d_x * d_x * (c_dt + d_z));
                ky_z1[x][y] = c_dt * c_dt * d_z / (2 * d_y * d_y * (c_dt + d_z));
            }
    }

    return 1;
}
    
void update_mur ()
{
    if (partition_data[partition_x0].boundary_type == boundary_mur)
    {
        if (mode == mode_full || mode == mode_tmy || mode == mode_tez) update_ey_x0();
        if (mode == mode_full || mode == mode_tmz || mode == mode_tey) update_ez_x0();
    }
    if (partition_data[partition_x1].boundary_type == boundary_mur)
    {
        if (mode == mode_full || mode == mode_tmy || mode == mode_tez) update_ey_x1();
        if (mode == mode_full || mode == mode_tmz || mode == mode_tey) update_ez_x1();
    }
    if (partition_data[partition_y0].boundary_type == boundary_mur)
    {
        if (mode == mode_full || mode == mode_tmx || mode == mode_tez) update_ex_y0();
        if (mode == mode_full || mode == mode_tmz || mode == mode_tex) update_ez_y0();
    }
    if (partition_data[partition_y1].boundary_type == boundary_mur)
    {
        if (mode == mode_full || mode == mode_tmx || mode == mode_tez) update_ex_y1();
        if (mode == mode_full || mode == mode_tmz || mode == mode_tex) update_ez_y1();
    }
    if (partition_data[partition_z0].boundary_type == boundary_mur)
    {
        if (mode == mode_full || mode == mode_tmx || mode == mode_tey) update_ex_z0();
        if (mode == mode_full || mode == mode_tmy || mode == mode_tex) update_ey_z0();
    }
    if (partition_data[partition_z1].boundary_type == boundary_mur)
    {
        if (mode == mode_full || mode == mode_tmx || mode == mode_tey) update_ex_z1();
        if (mode == mode_full || mode == mode_tmy || mode == mode_tex) update_ey_z1();
    }
}

static void update_ey_x0 ()
{
    for (y = 0; y < total_y; y++)
        for (z = 0; z < total_z; z++)
        {
            //printf("(%d %d) k1: %e k2:%e ey:%e \n",y, z, k1_x0[y][z], k2_x0[y][z], ey[0][y][z]);
            ey[0][y][z] = ey_x0[y][z] - ey_x0_d[y][z] + k1_x0[y][z] * (ey[1][y][z] + ey_x0_d[y][z]);

            ey_x0_dd[y][z] = ey_x0_d[y][z];
            ey_x0_d[y][z] = ey[0][y][z];

            ey_x0[y][z] = k2_x0[y][z] * (ey[0][y][z] + ey[1][y][z]);

            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_x0[y][z] += ky_x0[y][z] * ey[0][y-1][z] - ky_x0[y][z] * ey[0][y][z];
                ey_x0[y][z] += ky_x0[y][z] * ey[1][y-1][z] - ky_x0[y][z] * ey[1][y][z];
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_x0[y][z] += ky_x0[y][z] * ey[0][y+1][z] - ky_x0[y][z] * ey[0][y][z];
                ey_x0[y][z] += ky_x0[y][z] * ey[1][y+1][z] - ky_x0[y][z] * ey[1][y][z];
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ey_x0[y][z] += kz_x0[y][z] * ey[0][y][z-1] - kz_x0[y][z] * ey[0][y][z];
                ey_x0[y][z] += kz_x0[y][z] * ey[1][y][z-1] - kz_x0[y][z] * ey[1][y][z];
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ey_x0[y][z] += kz_x0[y][z] * ey[0][y][z+1] - kz_x0[y][z] * ey[0][y][z];
                ey_x0[y][z] += kz_x0[y][z] * ey[1][y][z+1] - kz_x0[y][z] * ey[1][y][z];
            }

            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_x0[y][z] += ky_x0[y][z] * (ey[0][y-1][z] + ey[0][y+1][z]) - 2 * ky_x0[y][z] * ey[0][y][z];
                ey_x0[y][z] += ky_x0[y][z] * (ey[1][y-1][z] + ey[1][y+1][z]) - 2 * ky_x0[y][z] * ey[1][y][z];
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ey_x0[y][z] += kz_x0[y][z] * (ey[0][y][z-1] + ey[0][y][z+1]) - 2 * kz_x0[y][z] * ey[0][y][z];
                ey_x0[y][z] += kz_x0[y][z] * (ey[1][y][z-1] + ey[1][y][z+1]) - 2 * kz_x0[y][z] * ey[1][y][z];
            }
        }
}

static void update_ez_x0 ()
{
    for (y = 0; y < total_y; y++)
        for (z = 0; z < total_z; z++)
        {
            ez[0][y][z] = ez_x0[y][z] - ez_x0_d[y][z] + k1_x0[y][z] * (ez[1][y][z] + ez_x0_d[y][z]);

            ez_x0_dd[y][z] = ez_x0_d[y][z];
            ez_x0_d[y][z] = ez[0][y][z];

            ez_x0[y][z] = k2_x0[y][z] * (ez[0][y][z] + ez[1][y][z]);
            //printf("%e %e %e\n", ez[0][y][z],  k1_x0[y][z], k2_x0[y][z]);

            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ez_x0[y][z] += ky_x0[y][z] * (ez[0][y-1][z] - ez[0][y][z]);
                ez_x0[y][z] += ky_x0[y][z] * (ez[1][y-1][z] - ez[1][y][z]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ez_x0[y][z] += ky_x0[y][z] * (ez[0][y+1][z] - ez[0][y][z]);
                ez_x0[y][z] += ky_x0[y][z] * (ez[1][y+1][z] - ez[1][y][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ez_x0[y][z] += kz_x0[y][z] * (ez[0][y][z-1] - ez[0][y][z]);
                ez_x0[y][z] += kz_x0[y][z] * (ez[1][y][z-1] - ez[1][y][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ez_x0[y][z] += kz_x0[y][z] * (ez[0][y][z+1] - ez[0][y][z]);
                ez_x0[y][z] += kz_x0[y][z] * (ez[1][y][z+1] - ez[1][y][z]);
            }
        }
}

static void update_ey_x1 ()
{
    for (y = 0; y < total_y; y++)
        for (z = 0; z < total_z; z++)
        {
            //printf("(%d %d) k1: %e k2:%e ey:%e \n",y, z, k1_x1[y][z], k2_x1[y][z], ey[0][y][z]);
            ey[total_x][y][z] = ey_x1[y][z] - ey_x1_d[y][z] + k1_x1[y][z] * (ey[total_x - 1][y][z] + ey_x1_d[y][z]);

            ey_x1_dd[y][z] = ey_x1_d[y][z];
            ey_x1_d[y][z] = ey[total_x][y][z];

            ey_x1[y][z] = k2_x1[y][z] * (ey[total_x][y][z] + ey[total_x - 1][y][z]);

            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_x1[y][z] += ky_x1[y][z] * (ey[total_x][y-1][z] - ey[total_x][y][z]);
                ey_x1[y][z] += ky_x1[y][z] * (ey[total_x - 1][y-1][z] - ey[total_x - 1][y][z]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_x1[y][z] += ky_x1[y][z] * (ey[total_x][y+1][z] - ey[total_x][y][z]);
                ey_x1[y][z] += ky_x1[y][z] * (ey[total_x - 1][y+1][z] - ey[total_x - 1][y][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ey_x1[y][z] += kz_x1[y][z] * (ey[total_x][y][z-1] - ey[total_x][y][z]);
                ey_x1[y][z] += kz_x1[y][z] * (ey[total_x - 1][y][z-1] - ey[total_x - 1][y][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ey_x1[y][z] += kz_x1[y][z] * (ey[total_x][y][z+1] - ey[total_x][y][z]);
                ey_x1[y][z] += kz_x1[y][z] * (ey[total_x - 1][y][z+1] - ey[total_x - 1][y][z]);
            }
        }
}

static void update_ez_x1 ()
{
    for (y = 0; y < total_y; y++)
        for (z = 0; z < total_z; z++)
        {
            ez[total_x][y][z] = ez_x1[y][z] - ez_x1_d[y][z] + k1_x1[y][z] * (ez[total_x - 1][y][z] + ez_x1_d[y][z]);

            ez_x1_dd[y][z] = ez_x1_d[y][z];
            ez_x1_d[y][z]  = ez[total_x][y][z];

            ez_x1[y][z] = k2_x1[y][z] * (ez[total_x][y][z] + ez[total_x - 1][y][z]);

            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ez_x1[y][z] += ky_x1[y][z] * (ez[total_x][y-1][z] - ez[total_x][y][z]);
                ez_x1[y][z] += ky_x1[y][z] * (ez[total_x - 1][y-1][z] - ez[total_x - 1][y][z]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ez_x1[y][z] += ky_x1[y][z] * (ez[total_x][y+1][z] - ez[total_x][y][z]);
                ez_x1[y][z] += ky_x1[y][z] * (ez[total_x - 1][y+1][z] - ez[total_x - 1][y][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ez_x1[y][z] += kz_x1[y][z] * (ez[total_x][y][z-1] - ez[total_x][y][z]);
                ez_x1[y][z] += kz_x1[y][z] * (ez[total_x - 1][y][z-1] - ez[total_x - 1][y][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ez_x1[y][z] += kz_x1[y][z] * (ez[total_x][y][z+1] - ez[total_x][y][z]);
                ez_x1[y][z] += kz_x1[y][z] * (ez[total_x - 1][y][z+1] - ez[total_x - 1][y][z]);
            }
        }
}

static void update_ex_y0 ()
{
    int x, z;

    for (x = 0; x < total_x; x++)
        for (z = 0; z < total_z; z++)
        {
            ex[x][0][z] = ex_y0[x][z] - ex_y0_d[x][z] + k1_y0[x][z] * (ex[x][1][z] + ex_y0_d[x][z]);

            ex_y0_dd[x][z] = ex_y0_d[x][z];
            ex_y0_d[x][z] = ex[x][0][z];

            ex_y0[x][z] = k2_y0[x][z] * (ex[x][0][z] + ex[x][1][z]);

            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ex_y0[x][z] += kx_y0[x][z] * (ex[x-1][0][z] - ex[x][0][z]);
                ex_y0[x][z] += kx_y0[x][z] * (ex[x-1][1][z] - ex[x][1][z]);
            }
            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ex_y0[x][z] += kx_y0[x][z] * (ex[x+1][0][z] - ex[x][0][z]);
                ex_y0[x][z] += kx_y0[x][z] * (ex[x+1][1][z] - ex[x][1][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ex_y0[x][z] += kz_y0[x][z] * (ex[x][0][z-1] - ex[x][0][z]);
                ex_y0[x][z] += kz_y0[x][z] * (ex[x][1][z-1] - ex[x][1][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ex_y0[x][z] += kz_y0[x][z] * (ex[x][0][z+1] - ex[x][0][z]);
                ex_y0[x][z] += kz_y0[x][z] * (ex[x][1][z+1] - ex[x][1][z]);
            }
        }
}

static void update_ez_y0 ()
{
    for (x = 0; x < total_x; x++)
        for (z = 0; z < total_z; z++)
        {
            ez[x][0][z] = ez_y0[x][z] - ez_y0_d[x][z] + k1_y0[x][z] * (ez[x][1][z] + ez_y0_d[x][z]);

            ez_y0_dd[x][z] = ez_y0_d[x][z];
            ez_y0_d[x][z]  = ez[x][0][z];

            ez_y0[x][z] = k2_y0[x][z] * (ez[x][0][z] + ez[x][1][z]);

            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ez_y0[x][z] += kx_y0[x][z] * (ez[x-1][0][z] - ez[x][0][z]);
                ez_y0[x][z] += kx_y0[x][z] * (ez[x-1][1][z] - ez[x][1][z]);
            }
            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ez_y0[x][z] += kx_y0[x][z] * (ez[x+1][0][z] - ez[x][0][z]);
                ez_y0[x][z] += kx_y0[x][z] * (ez[x+1][1][z] - ez[x][1][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ez_y0[x][z] += kz_y0[x][z] * (ez[x][0][z-1] - ez[x][0][z]);
                ez_y0[x][z] += kz_y0[x][z] * (ez[x][1][z-1] - ez[x][1][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ez_y0[x][z] += kz_y0[x][z] * (ez[x][0][z+1] - ez[x][0][z]);
                ez_y0[x][z] += kz_y0[x][z] * (ez[x][1][z+1] - ez[x][1][z]);
            }
        }
}

static void update_ex_y1 ()
{
    for (x = 0; x < total_x; x++)
        for (z = 0; z < total_z; z++)
        {
            ex[x][total_y][z] = ex_y1[x][z] - ex_y1_d[x][z] + k1_y1[x][z] * (ex[x][total_y - 1][z] + ex_y1_d[x][z]);

            ex_y1_dd[x][z] = ex_y1_d[x][z];
            ex_y1_d[x][z] = ex[x][total_y][z];

            ex_y1[x][z] = k2_y1[x][z] * (ex[x][total_y][z] + ex[x][total_y - 1][z]);

            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ex_y1[x][z] += kx_y1[x][z] * (ex[x-1][total_y][z] - ex[x][total_y][z]);
                ex_y1[x][z] += kx_y1[x][z] * (ex[x-1][total_y - 1][z] - ex[x][total_y - 1][z]);
            }
            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ex_y1[x][z] += kx_y1[x][z] * (ex[x+1][total_y][z] - ex[x][total_y][z]);
                ex_y1[x][z] += kx_y1[x][z] * (ex[x+1][total_y - 1][z] - ex[x][total_y - 1][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ex_y1[x][z] += kz_y1[x][z] * (ex[x][total_y][z-1] - ex[x][total_y][z]);
                ex_y1[x][z] += kz_y1[x][z] * (ex[x][total_y - 1][z-1] - ex[x][total_y - 1][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ex_y1[x][z] += kz_y1[x][z] * (ex[x][total_y][z+1] - ex[x][total_y][z]);
                ex_y1[x][z] += kz_y1[x][z] * (ex[x][total_y - 1][z+1] - ex[x][total_y - 1][z]);
            }
        }
}

static void update_ez_y1 ()
{
    int x, z;

    for (x = 0; x < total_x; x++)
        for (z = 0; z < total_z; z++)
        {
            ez[x][total_y][z] = ez_y1[x][z] - ez_y1_d[x][z] + k1_y1[x][z] * (ez[x][total_y - 1][z] + ez_y1_d[x][z]);

            ez_y1_dd[x][z] = ez_y1_d[x][z];
            ez_y1_d[x][z] = ez[x][total_y][z];

            ez_y1[x][z] = k2_y1[x][z] * (ez[x][total_y][z] + ez[x][total_y - 1][z]);

            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ez_y1[x][z] += kx_y1[x][z] * (ez[x-1][total_y][z] - ez[x][total_y][z]);
                ez_y1[x][z] += kx_y1[x][z] * (ez[x-1][total_y - 1][z] - ez[x][total_y - 1][z]);
            }
            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ez_y1[x][z] += kx_y1[x][z] * (ez[x+1][total_y][z] - ez[x][total_y][z]);
                ez_y1[x][z] += kx_y1[x][z] * (ez[x+1][total_y - 1][z] - ez[x][total_y - 1][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ez_y1[x][z] += kz_y1[x][z] * (ez[x][total_y][z-1] - ez[x][total_y][z]);
                ez_y1[x][z] += kz_y1[x][z] * (ez[x][total_y - 1][z-1] - ez[x][total_y - 1][z]);
            }
            if (z > 0 && z < total_z && mode != mode_tmz && mode != mode_tez)
            {
                ez_y1[x][z] += kz_y1[x][z] * (ez[x][total_y][z+1] - ez[x][total_y][z]);
                ez_y1[x][z] += kz_y1[x][z] * (ez[x][total_y - 1][z+1] - ez[x][total_y - 1][z]);
            }
        }
}

static void update_ex_z0 ()
{
    for (x = 0; x < total_x; x++)
        for (y = 0; y < total_y; y++)
        {
            ex[x][y][0] = ex_z0[x][y] - ex_z0_d[x][y] + k1_z0[x][y] * (ex[x][y][1] + ex_z0_d[x][y]);

            ex_z0_dd[x][y] = ex_z0_d[x][y];
            ex_z0_d[x][y]  = ex[x][y][0];

            ex_z0[x][y] = k2_z0[x][y] * (ex[x][y][0] + ex[x][y][1]);

            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ex_z0[x][y] += kx_z0[x][y] * (ex[x-1][y][0] - ex[x][y][0]);
                ex_z0[x][y] += kx_z0[x][y] * (ex[x-1][y][1] - ex[x][y][1]);
            }
            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ex_z0[x][y] += kx_z0[x][y] * (ex[x+1][y][0] - ex[x][y][0]);
                ex_z0[x][y] += kx_z0[x][y] * (ex[x+1][y][1] - ex[x][y][1]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y-1][0] - ex[x][y][0]);
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y-1][1] - ex[x][y][1]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y+1][0] - ex[x][y][0]);
                ex_z0[x][y] += ky_z0[x][y] * (ex[x][y+1][1] - ex[x][y][1]);
            }
        }
}

static void update_ey_z0 ()
{
    for (x = 0; x < total_x; x++)
        for (y = 0; y < total_y; y++)
        {
            ey[x][y][0] = ey_z0[x][y] - ey_z0_d[x][y] + k1_z0[x][y] * (ey[x][y][1] + ey_z0_d[x][y]);

            ey_z0_dd[x][y] = ey_z0_d[x][y];
            ey_z0_d[x][y]  = ey[x][y][0];

            ey_z0[x][y] = k2_z0[x][y] * (ey[x][y][0] + ey[x][y][1]);

            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ey_z0[x][y] += kx_z0[x][y] * (ey[x-1][y][0] - ey[x][y][0]);
                ey_z0[x][y] += kx_z0[x][y] * (ey[x-1][y][1] - ey[x][y][1]);
            }
            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ey_z0[x][y] += kx_z0[x][y] * (ey[x+1][y][0] - ey[x][y][0]);
                ey_z0[x][y] += kx_z0[x][y] * (ey[x+1][y][1] - ey[x][y][1]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y-1][0] - ey[x][y][0]);
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y-1][1] - ey[x][y][1]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y+1][0] - ey[x][y][0]);
                ey_z0[x][y] += ky_z0[x][y] * (ey[x][y+1][1] - ey[x][y][1]);
            }
        }
}

static void update_ex_z1 ()
{
    for (x = 0; x < total_x; x++)
        for (y = 0; y < total_y; y++)
        {
            ex[x][y][total_z] = ex_z1[x][y] - ex_z1_d[x][y] + k1_z1[x][y] * (ex[x][y][total_z - 1] + ex_z1_d[x][y]);

            ex_z1_dd[x][y] = ex_z1_d[x][y];
            ex_z1_d[x][y] = ex[x][y][total_z];

            ex_z1[x][y] = k2_z1[x][y] * (ex[x][y][total_z] + ex[x][y][total_z - 1]);

            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ex_z1[x][y] += kx_z1[x][y] * (ex[x-1][y][total_z] - ex[x][y][total_z]);
                ex_z1[x][y] += kx_z1[x][y] * (ex[x-1][y][total_z - 1] - ex[x][y][total_z - 1]);
            }
            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ex_z1[x][y] += kx_z1[x][y] * (ex[x+1][y][total_z] - ex[x][y][total_z]);
                ex_z1[x][y] += kx_z1[x][y] * (ex[x+1][y][total_z - 1] - ex[x][y][total_z - 1]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ex_z1[x][y] += ky_z1[x][y] * (ex[x][y-1][total_z] - ex[x][y][total_z]);
                ex_z1[x][y] += ky_z1[x][y] * (ex[x][y-1][total_z - 1] - ex[x][y][total_z - 1]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ex_z1[x][y] += ky_z1[x][y] * (ex[x][y+1][total_z] - ex[x][y][total_z]);
                ex_z1[x][y] += ky_z1[x][y] * (ex[x][y+1][total_z - 1] - ex[x][y][total_z - 1]);
            }
        }
}

static void update_ey_z1 ()
{
    for (x = 0; x < total_x; x++)
        for (y = 0; y < total_y; y++)
        {
            ey[x][y][total_z] = ey_z1[x][y] - ey_z1_d[x][y] + k1_z1[x][y] * (ey[x][y][total_z - 1] + ey_z1_d[x][y]);

            ey_z1_dd[x][y] = ey_z1_d[x][y];
            ey_z1_d[x][y] = ey[x][y][total_z];

            ey_z1[x][y] = k2_z1[x][y] * (ey[x][y][total_z] + ey[x][y][total_z - 1]);

            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ey_z1[x][y] += kx_z1[x][y] * (ey[x-1][y][total_z] - ey[x][y][total_z]);
                ey_z1[x][y] += kx_z1[x][y] * (ey[x-1][y][total_z - 1] - ey[x][y][total_z - 1]);
            }
            if (x > 0 && x < total_x && mode != mode_tmx && mode != mode_tex)
            {
                ey_z1[x][y] += kx_z1[x][y] * (ey[x+1][y][total_z] - ey[x][y][total_z]);
                ey_z1[x][y] += kx_z1[x][y] * (ey[x+1][y][total_z - 1] - ey[x][y][total_z - 1]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_z1[x][y] += ky_z1[x][y] * (ey[x][y-1][total_z] - ey[x][y][total_z]);
                ey_z1[x][y] += ky_z1[x][y] * (ey[x][y-1][total_z - 1] - ey[x][y][total_z - 1]);
            }
            if (y > 0 && y < total_y && mode != mode_tmy && mode != mode_tey)
            {
                ey_z1[x][y] += ky_z1[x][y] * (ey[x][y+1][total_z] - ey[x][y][total_z]);
                ey_z1[x][y] += ky_z1[x][y] * (ey[x][y+1][total_z - 1] - ey[x][y][total_z - 1]);
            }
        }
}

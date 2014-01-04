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
 *     File Name : pml.c
 * Last Modified : Thu 08 Nov 2012 12:37:22 PM EST
 */

#include "tools.h"
#include <math.h>
#include "mem.h"
#include "domain.h"
#include "h5io.h"
#include "pml.h"

double hxz;
double hyx;
double hzy;

double exz;
double eyx;
double ezy;

double *pml_c1x;
double *pml_c1y;
double *pml_c1z;
double *pml_d1x;
double *pml_d1y;
double *pml_d1z;

double *pml_c2x;
double *pml_c2y;
double *pml_c2z;

double *pml_d2x;
double *pml_d2y;
double *pml_d2z;

double ***hxy;
double ***hyz;
double ***hzx;

double ***exy;
double ***eyz;
double ***ezx;

double pml_c2x_0;
double pml_c2y_0;
double pml_c2z_0;
double pml_d2x_0;
double pml_d2y_0;
double pml_d2z_0;

int setup_pml (char* file_name)
{
    double pml_grading_order; 

    double pml_sigma_max;
    double pml_factor;
    double pml_sigma_e;
    double pml_sigma_h;

    int layer_index;
    int status;

    status = h5_get_attr(file_name, "boundaries", "pml_sigma_max", &pml_sigma_max);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "boundaries", "pml_grading_order", &pml_grading_order);
    inspect(status, "fail to get h5 attributes");

    hxy = bx;
    hyz = by;
    hzx = bz;
    exy = dx;
    eyz = dy;
    ezx = dz;

    pml_c1x = (double *)mem1(type_double, abc_size);
    inspect(pml_c1x, "fail to allocate mem1ory for field");
    pml_c1y = (double *)mem1(type_double, abc_size);
    inspect(pml_c1y, "fail to allocate mem1ory for field");
    pml_c1z = (double *)mem1(type_double, abc_size);
    inspect(pml_c1z, "fail to allocate mem1ory for field");
    pml_d1x = (double *)mem1(type_double, abc_size);
    inspect(pml_d1x, "fail to allocate mem1ory for field");
    pml_d1y = (double *)mem1(type_double, abc_size);
    inspect(pml_d1y, "fail to allocate mem1ory for field");
    pml_d1z = (double *)mem1(type_double, abc_size);
    inspect(pml_d1z, "fail to allocate mem1ory for field");

    pml_c2x = (double *)mem1(type_double, abc_size);
    inspect(pml_c2x, "fail to allocate mem1ory for field");
    pml_c2y = (double *)mem1(type_double, abc_size);
    inspect(pml_c2y, "fail to allocate mem1ory for field");
    pml_c2z = (double *)mem1(type_double, abc_size);
    inspect(pml_c2z, "fail to allocate mem1ory for field");

    pml_d2x = (double *)mem1(type_double, abc_size);
    inspect(pml_d2x, "fail to allocate mem1ory for field");
    pml_d2y = (double *)mem1(type_double, abc_size);
    inspect(pml_d2y, "fail to allocate mem1ory for field");
    pml_d2z = (double *)mem1(type_double, abc_size);
    inspect(pml_d2z, "fail to allocate mem1ory for field");


    pml_c2x_0 = d_t / (EPSILON0 * d_x);
    pml_c2y_0 = d_t / (EPSILON0 * d_y);
    pml_c2z_0 = d_t / (EPSILON0 * d_z);

    pml_d2x_0 = d_t / (MU0 * d_x);
    pml_d2y_0 = d_t / (MU0 * d_y);
    pml_d2z_0 = d_t / (MU0 * d_z);

    //boundaryThickness = (double )abc_size * dh;
    //pml_sigma_max = -log(reflectionCoefficient0) * (grading_order + 1.0) 
    //    * boundaryEpsilon * EPSILON0 * C0 / (2.0 * boundaryThickness);

    pml_factor = pml_sigma_max 
        / ((pml_grading_order + 1) * pow(2, pml_grading_order + 1) * pow(abc_size, pml_grading_order));

    for (layer_index = 0; layer_index < abc_size; layer_index++)
    {
        if (layer_index == 0)
            pml_sigma_e = pml_factor;
        else
            pml_sigma_e = pml_factor * (pow(2 * layer_index + 1, pml_grading_order + 1) 
                    - pow(2 * layer_index - 1, pml_grading_order + 1));

        pml_sigma_h = pml_sigma_e * MU0 / EPSILON0;

        pml_c1x[layer_index] = (1.0 - (pml_sigma_e * d_t) / (2.0 * EPSILON0)) / (1.0 + (pml_sigma_e * d_t) / (2.0 * EPSILON0));
        pml_c1y[layer_index] = (1.0 - (pml_sigma_e * d_t) / (2.0 * EPSILON0)) / (1.0 + (pml_sigma_e * d_t) / (2.0 * EPSILON0));
        pml_c1z[layer_index] = (1.0 - (pml_sigma_e * d_t) / (2.0 * EPSILON0)) / (1.0 + (pml_sigma_e * d_t) / (2.0 * EPSILON0));
        pml_d1x[layer_index] = (1.0 - (pml_sigma_h * d_t) / (2.0 * MU0)) / (1.0 + (pml_sigma_h * d_t) / (2.0 * MU0));
        pml_d1y[layer_index] = (1.0 - (pml_sigma_h * d_t) / (2.0 * MU0)) / (1.0 + (pml_sigma_h * d_t) / (2.0 * MU0));
        pml_d1z[layer_index] = (1.0 - (pml_sigma_h * d_t) / (2.0 * MU0)) / (1.0 + (pml_sigma_h * d_t) / (2.0 * MU0));
        pml_c2x[layer_index] = (d_t / d_x / EPSILON0) / (1.0 + pml_sigma_e * d_t / (2.0 * EPSILON0));
        pml_c2y[layer_index] = (d_t / d_y / EPSILON0) / (1.0 + pml_sigma_e * d_t / (2.0 * EPSILON0));
        pml_c2z[layer_index] = (d_t / d_z / EPSILON0) / (1.0 + pml_sigma_e * d_t / (2.0 * EPSILON0));
        pml_d2x[layer_index] = (d_t / d_x / MU0) / (1.0 + pml_sigma_h * d_t / (2.0 * MU0));
        pml_d2y[layer_index] = (d_t / d_y / MU0) / (1.0 + pml_sigma_h * d_t / (2.0 * MU0));
        pml_d2z[layer_index] = (d_t / d_z / MU0) / (1.0 + pml_sigma_h * d_t / (2.0 * MU0));
    }
    return 1;
}


void pml_get_h ()
{
    int x, y, z;
    double d1x, d1y, d1z;
    double d2x, d2y, d2z;

    for (x = 0; x < total_length_x; x++)
    {
        for (y = 0; y < total_length_y; y++)
        {
            for (z = 0; z < total_length_z; z++)
            {
                if (in_partition_main(x, y, z)) continue;

                d1x = 1.0;
                d1y = 1.0;
                d1z = 1.0;
                d2x = pml_d2x_0;
                d2y = pml_d2y_0;
                d2z = pml_d2z_0;

                if (in_partition_x0(x, y, z))
                {
                    d1x = pml_d1x[partition_data[partition_x0].x_stop - x - 1];
                    d2x = pml_d2x[partition_data[partition_x0].x_stop - x - 1];
                }
                if (in_partition_y0(x, y, z))
                {
                    d1y = pml_d1y[partition_data[partition_y0].y_stop - y - 1];
                    d2y = pml_d2y[partition_data[partition_y0].y_stop - y - 1];
                }
                if (in_partition_z0(x, y, z))
                {
                    d1z = pml_d1z[partition_data[partition_z0].z_stop - z - 1];
                    d2z = pml_d2z[partition_data[partition_z0].z_stop - z - 1];
                }
                if (in_partition_x1(x, y, z))
                {
                    d1x = pml_d1x[-partition_data[partition_x1].x_start + x];
                    d2x = pml_d2x[-partition_data[partition_x1].x_start + x];
                }
                if (in_partition_y1(x, y, z))
                {
                    d1y = pml_d1y[-partition_data[partition_y1].y_start + y];
                    d2y = pml_d2y[-partition_data[partition_y1].y_start + y];
                }
                if (in_partition_z1(x, y, z))
                {
                    d1z = pml_d1z[-partition_data[partition_z1].z_start + z];
                    d2z = pml_d2z[-partition_data[partition_z1].z_start + z];
                }

                switch (mode)
                {
                    case mode_full:
                        hxz = hx[x][y][z] - hxy[x][y][z];
                        hxy[x][y][z] = d1y * hxy[x][y][z] - d2y * (ez[x][y+1][z] - ez[x][y][z]);
                        hxz = d1z * hxz + d2z * (ey[x][y][z+1] - ey[x][y][z]);
                        hx[x][y][z] = hxy[x][y][z] + hxz;

                        hyx = hy[x][y][z] - hyz[x][y][z];
                        hyz[x][y][z] = d1z * hyz[x][y][z] - d2z * (ex[x][y][z+1] - ex[x][y][z]);
                        hyx = d1x * hyx + d2x * (ez[x+1][y][z] - ez[x][y][z]);
                        hy[x][y][z] = hyz[x][y][z] + hyx;

                        hzy = hz[x][y][z] - hzx[x][y][z];
                        hzx[x][y][z] = d1x * hzx[x][y][z] - d2x * (ey[x+1][y][z] - ey[x][y][z]);
                        hzy = d1y * hzy + d2y * (ex[x][y+1][z] - ex[x][y][z]);
                        hz[x][y][z] = hzx[x][y][z] + hzy;
                        break;
                    case mode_tmx:
                        hy[x][y][z] = d1z * hy[x][y][z] - d2z * (ex[x][y][z+1] - ex[x][y][z]);
                        hz[x][y][z] = d1y * hz[x][y][z] + d2y * (ex[x][y+1][z] - ex[x][y][z]);
                        break;
                    case mode_tmy:
                        hx[x][y][z] = d1z * hx[x][y][z] + d2z * (ey[x][y][z+1] - ey[x][y][z]);
                        hz[x][y][z] = d1x * hz[x][y][z] - d2x * (ey[x+1][y][z] - ey[x][y][z]);
                        break;
                    case mode_tmz:
                        hx[x][y][z] = d1y * hx[x][y][z] - d2y * (ez[x][y+1][z] - ez[x][y][z]);
                        hy[x][y][z] = d1x * hy[x][y][z] + d2x * (ez[x+1][y][z] - ez[x][y][z]);
                        break;
                    case mode_tex:
                        hxz = hx[x][y][z] - hxy[x][y][z];
                        hxy[x][y][z] = d1y * hxy[x][y][z] - d2y * (ez[x][y+1][z] - ez[x][y][z]);
                        hxz =  d1z * hxz + d2z * (ey[x][y][z+1] - ey[x][y][z]);
                        hx[x][y][z] = hxy[x][y][z] + hxz;
                        break;
                    case mode_tey:
                        hyx = hy[x][y][z] - hyz[x][y][z];
                        hyz[x][y][z] = d1z * hyz[x][y][z] - d2z * (ex[x][y][z+1] - ex[x][y][z]);
                        hyx = d1x * hyx + d2x * (ez[x+1][y][z] - ez[x][y][z]);
                        hy[x][y][z] = hyz[x][y][z] + hyx;
                        break;
                    case mode_tez:
                        hzy = hz[x][y][z] - hzx[x][y][z];
                        hzx[x][y][z] = d1x * hzx[x][y][z] - d2x * (ey[x+1][y][z] - ey[x][y][z]);
                        hzy = d1y * hzy + d2y * (ex[x][y+1][z] - ex[x][y][z]);
                        hz[x][y][z] = hzx[x][y][z] + hzy;
                        break;
                }
            }
        }
    }
}

void pml_get_e ()
{
    int x, y, z;
    double c1x, c1y, c1z;
    double c2x, c2y, c2z;

    for (x = 0; x < total_length_x; x++)
    {
        for (y = 0; y < total_length_y; y++)
        {
            for (z = 0; z < total_length_z; z++)
            {
                if (in_partition_main(x, y, z))  continue;

                c1x = 1.0;
                c1y = 1.0;
                c1z = 1.0;
                c2x = pml_c2x_0;
                c2y = pml_c2y_0;
                c2z = pml_c2z_0;

                if (in_partition_x0(x, y, z))
                {
                c1x = pml_c1x[partition_data[partition_x0].x_stop - x - 1];
                c2x = pml_c2x[partition_data[partition_x0].x_stop - x - 1];
                }
                if (in_partition_y0(x, y, z))
                {
                c1y = pml_c1y[partition_data[partition_y0].y_stop - y - 1];
                c2y = pml_c2y[partition_data[partition_y0].y_stop - y - 1];
                }
                if (in_partition_z0(x, y, z))
                {
                c1z = pml_c1z[partition_data[partition_z0].z_stop - z - 1];
                c2z = pml_c2z[partition_data[partition_z0].z_stop - z - 1];
                }
                if (in_partition_x1(x, y, z))
                {
                c1x = pml_c1x[-partition_data[partition_x1].x_start + x];
                c2x = pml_c2x[-partition_data[partition_x1].x_start + x];
                }
                if (in_partition_y1(x, y, z))
                {
                c1y = pml_c1y[-partition_data[partition_y1].y_start + y];
                c2y = pml_c2y[-partition_data[partition_y1].y_start + y];
                }
                if (in_partition_z1(x, y, z))
                {
                c1z = pml_c1z[-partition_data[partition_z1].z_start + z];
                c2z = pml_c2z[-partition_data[partition_z1].z_start + z];
                }

                switch (mode)
                {
                    case mode_full:
                        exz = ex[x][y][z] - exy[x][y][z];
                        if (y != 0)
                            exy[x][y][z] = c1y * exy[x][y][z] + c2y * (hz[x][y][z] - hz[x][y-1][z]);
                        else 
                            exy[x][y][z] = c1y * exy[x][y][z]; 
                        if (z != 0)
                            exz =  c1z * exz - c2z * (hy[x][y][z] - hy[x][y][z-1]);
                        else 
                            exz =  c1z * exz;
                        ex[x][y][z] = exy[x][y][z] + exz;

                        eyx = ey[x][y][z] - eyz[x][y][z];
                        if (z != 0)
                            eyz[x][y][z] = c1z * eyz[x][y][z] + c2z * (hx[x][y][z] - hx[x][y][z-1]);
                        else 
                            eyz[x][y][z] = c1z * eyz[x][y][z];
                        if (x != 0)
                            eyx = c1x * eyx - c2x * (hz[x][y][z] - hz[x-1][y][z]);
                        else 
                            eyx = c1x * eyx;
                        ey[x][y][z] = eyz[x][y][z] + eyx;

                        ezy = ez[x][y][z] - ezx[x][y][z];
                        if (x != 0)
                            ezx[x][y][z] = c1x * ezx[x][y][z] + c2x * (hy[x][y][z] - hy[x-1][y][z]);
                        else 
                            ezx[x][y][z] = c1x * ezx[x][y][z];
                        if (y != 0)
                            ezy = c1y * ezy - c2y * (hx[x][y][z] - hx[x][y-1][z]);
                        else 
                            ezy = c1y * ezy;
                        ez[x][y][z] = ezx[x][y][z] + ezy;
                        break;

                    case mode_tmx:
                        exz = ex[x][y][z] - exy[x][y][z];
                        if (y != 0)
                            exy[x][y][z] = c1y * exy[x][y][z] + c2y * (hz[x][y][z] - hz[x][y-1][z]);
                        else 
                            exy[x][y][z] = c1y * exy[x][y][z]; 
                        if (z != 0)
                            exz =  c1z * exz - c2z * (hy[x][y][z] - hy[x][y][z-1]);
                        else 
                            exz =  c1z * exz;
                        ex[x][y][z] = exy[x][y][z] + exz;
                        break;
                    case mode_tmy:
                        eyx = ey[x][y][z] - eyz[x][y][z];
                        if (z != 0)
                            eyz[x][y][z] = c1z * eyz[x][y][z] + c2z * (hx[x][y][z] - hx[x][y][z-1]);
                        else 
                            eyz[x][y][z] = c1z * eyz[x][y][z];
                        if (x != 0)
                            eyx = c1x * eyx - c2x * (hz[x][y][z] - hz[x-1][y][z]);
                        else 
                            eyx = c1x * eyx;
                        ey[x][y][z] = eyz[x][y][z] + eyx;
                        break;
                    case mode_tmz:
                        ezy = ez[x][y][z] - ezx[x][y][z];
                        if (x != 0)
                            ezx[x][y][z] = c1x * ezx[x][y][z] + c2x * (hy[x][y][z] - hy[x-1][y][z]);
                        else 
                            ezx[x][y][z] = c1x * ezx[x][y][z];
                        if (y != 0)
                            ezy = c1y * ezy - c2y * (hx[x][y][z] - hx[x][y-1][z]);
                        else 
                            ezy = c1y * ezy;
                        ez[x][y][z] = ezx[x][y][z] + ezy;
                        break;
                    case mode_tex:
                        if (z != 0)
                            ey[x][y][z] = c1z * ey[x][y][z] + c2z * (hx[x][y][z] - hx[x][y][z-1]);
                        if (y != 0)
                            ez[x][y][z] = c1y * ez[x][y][z] - c2y * (hx[x][y][z] - hx[x][y-1][z]);
                        break;
                    case mode_tey:
                        if (x != 0)
                            ez[x][y][z] = c1x * ez[x][y][z] + c2x * (hy[x][y][z] - hy[x-1][y][z]);
                        if (z != 0) 
                            ex[x][y][z] = c1z * ex[x][y][z] - c2z * (hy[x][y][z] - hy[x][y][z-1]);
                        break;
                    case mode_tez:
                        if (y != 0) 
                            ex[x][y][z] = c1y * ex[x][y][z] + c2y * (hz[x][y][z] - hz[x][y-1][z]);
                        if (z != 0)
                            ey[x][y][z] = c1z * ey[x][y][z] + c2z * (hx[x][y][z] - hx[x][y][z-1]);
                        break;
                }
            }
        }
    }
}





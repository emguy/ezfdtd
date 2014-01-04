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
 *     File Name : cpml.c
 * Last Modified : Thu 06 Dec 2012 10:39:02 AM EST
 */

#include "tools.h"
#include <math.h>
#include "mem.h"
#include "domain.h"
#include "h5io.h"
#include "cpml.h"

double *be;
double *bh;

double *cex;
double *cey;
double *cez;
double *chx;
double *chy;
double *chz;

double *pml_kappa_h;
double *pml_kappa_e;

CPMLFields psi[7];

int setup_cpml (char* file_name)
{
    double pml_grading_order; 
    double pml_sigma_max;
    double pml_alpha_max;
    double pml_kappa_max;
    double pml_sigma_e;
    double pml_sigma_h;
    double pml_alpha_e;
    double pml_alpha_h;

    int layer_index, partition_index;
    int x_start, x_stop, y_start, y_stop, z_start, z_stop;

    int status;

    status = h5_get_attr(file_name, "boundaries", "pml_sigma_max", &pml_sigma_max);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "boundaries", "pml_alpha_max", &pml_alpha_max);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "boundaries", "pml_kappa_max", &pml_kappa_max);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "boundaries", "pml_grading_order", &pml_grading_order);
    inspect(status, "fail to get h5 attributes");
    //printf("CPML: %e %e %e %e\n", pml_sigma_max, pml_alpha_max, pml_kappa_max, pml_grading_order);


    be  = (double *)mem1(type_double, abc_size);
    inspect(be, "fail to allocate memory for field");
    bh  = (double *)mem1(type_double, abc_size);
    inspect(bh, "fail to allocate memory for field");
    cex = (double *)mem1(type_double, abc_size);
    inspect(cex, "fail to allocate memory for field");
    cey = (double *)mem1(type_double, abc_size);
    inspect(cey, "fail to allocate memory for field");
    cez = (double *)mem1(type_double, abc_size);
    inspect(cez, "fail to allocate memory for field");
    chx = (double *)mem1(type_double, abc_size);
    inspect(chx, "fail to allocate memory for field");
    chy = (double *)mem1(type_double, abc_size);
    inspect(chy, "fail to allocate memory for field");
    chz = (double *)mem1(type_double, abc_size);
    inspect(chz, "fail to allocate memory for field");
    pml_kappa_e = (double *)mem1(type_double, abc_size);
    inspect(pml_kappa_e, "fail to allocate memory for field");
    pml_kappa_h = (double *)mem1(type_double, abc_size);
    inspect(pml_kappa_h, "fail to allocate memory for field");

    for (layer_index = 0; layer_index < abc_size; layer_index++)
    {
        pml_sigma_e = pml_sigma_max * pow((layer_index + 1.0) / (abc_size + 1.0), pml_grading_order);
        pml_sigma_h = pml_sigma_max * pow((layer_index + 1.0) / (abc_size + 1.0), pml_grading_order);
        pml_kappa_e[layer_index] = 1.0 + (pml_kappa_max - 1.0) * pow((layer_index + 1.0) / (abc_size + 1.0), pml_grading_order);
        pml_kappa_h[layer_index] = 1.0 + (pml_kappa_max - 1.0) * pow((layer_index + 1.0) / (abc_size + 1.0), pml_grading_order);
        pml_alpha_e = pml_alpha_max * (abc_size - layer_index - 1.0) / (abc_size - 1.0);
        pml_alpha_h = pml_alpha_max * (abc_size - layer_index - 1.0) / (abc_size - 1.0);

        //printf("layer %03d: %e, %e, %e, %e \n", layer_index, pml_sigma_e, pml_sigma_h, pml_alpha_e, pml_alpha_h);
        //printf("layer %03d: %e\n", layer_index, pml_kappa_e);

        be[layer_index]  = exp(-(pml_sigma_e / pml_kappa_e[layer_index] + pml_alpha_e) * d_t / EPSILON0);
        bh[layer_index]  = exp(-(pml_sigma_h / pml_kappa_h[layer_index] + pml_alpha_h) * d_t / EPSILON0);

        cex[layer_index] = pml_sigma_e * (be[layer_index] - 1) / (pml_sigma_e + pml_alpha_e * pml_kappa_e[layer_index]) / d_x / pml_kappa_e[layer_index];
        cey[layer_index] = pml_sigma_e * (be[layer_index] - 1) / (pml_sigma_e + pml_alpha_e * pml_kappa_e[layer_index]) / d_y / pml_kappa_e[layer_index];
        cez[layer_index] = pml_sigma_e * (be[layer_index] - 1) / (pml_sigma_e + pml_alpha_e * pml_kappa_e[layer_index]) / d_z / pml_kappa_e[layer_index];

        chx[layer_index] = pml_sigma_h * (be[layer_index] - 1) / (pml_sigma_h + pml_alpha_h * pml_kappa_h[layer_index]) / d_x / pml_kappa_h[layer_index];
        chy[layer_index] = pml_sigma_h * (be[layer_index] - 1) / (pml_sigma_h + pml_alpha_h * pml_kappa_h[layer_index]) / d_y / pml_kappa_h[layer_index];
        chz[layer_index] = pml_sigma_h * (be[layer_index] - 1) / (pml_sigma_h + pml_alpha_h * pml_kappa_h[layer_index]) / d_z / pml_kappa_h[layer_index];

        //printf("layer %03d: %e, %e, %e, %e, %e, %e, %e  \n", layer_index, pml_sigma_e, pml_kappa_e[layer_index], pml_alpha_e, be[layer_index], cex[layer_index], bh[layer_index], chx[layer_index]);
    }

    for (partition_index = 1; partition_index < 7; partition_index++)
    {
        x_start = partition_data[partition_index].x_start;
        x_stop = partition_data[partition_index].x_stop;
        y_start = partition_data[partition_index].y_start;
        y_stop = partition_data[partition_index].y_stop;
        z_start = partition_data[partition_index].z_start;
        z_stop = partition_data[partition_index].z_stop;

        if (partition_data[partition_index].boundary_type == boundary_pml)
        {
            switch (mode)
            {
                case mode_tmx:
                    psi[partition_index].hyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyz, "fail to allocate memory for field");
                    psi[partition_index].hzy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzy, "fail to allocate memory for field");
                    psi[partition_index].eyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyx, "fail to allocate memory for field");
                    psi[partition_index].ezx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezx, "fail to allocate memory for field");
                    break;
                case mode_tmy:
                    psi[partition_index].hxz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxz, "fail to allocate memory for field");
                    psi[partition_index].hzx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzx, "fail to allocate memory for field");
                    psi[partition_index].exy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exy, "fail to allocate memory for field");
                    psi[partition_index].ezy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezy, "fail to allocate memory for field");
                    break;
                case mode_tmz:
                    psi[partition_index].hyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyx, "fail to allocate memory for field");
                    psi[partition_index].hxy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxy, "fail to allocate memory for field");
                    psi[partition_index].eyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyz, "fail to allocate memory for field");
                    psi[partition_index].exz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exz, "fail to allocate memory for field");
                    break;
                case mode_tex:
                    psi[partition_index].hyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyx, "fail to allocate memory for field");
                    psi[partition_index].hzx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzx, "fail to allocate memory for field");
                    psi[partition_index].eyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyz, "fail to allocate memory for field");
                    psi[partition_index].ezy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezy, "fail to allocate memory for field");
                    break;
                case mode_tey:
                    psi[partition_index].hxy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxy, "fail to allocate memory for field");
                    psi[partition_index].hzy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzy, "fail to allocate memory for field");
                    psi[partition_index].exz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exz, "fail to allocate memory for field");
                    psi[partition_index].ezx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezx, "fail to allocate memory for field");
                    break;
                case mode_tez:
                    psi[partition_index].hxz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxz, "fail to allocate memory for field");
                    psi[partition_index].hyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyz, "fail to allocate memory for field");
                    psi[partition_index].exy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exy, "fail to allocate memory for field");
                    psi[partition_index].eyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyx, "fail to allocate memory for field");
                    break;
                case mode_full:
                    psi[partition_index].hxy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxy, "fail to allocate memory for field");
                    psi[partition_index].hxz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxz, "fail to allocate memory for field");
                    psi[partition_index].hyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyx, "fail to allocate memory for field");
                    psi[partition_index].hyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyz, "fail to allocate memory for field");
                    psi[partition_index].hzx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzx, "fail to allocate memory for field");
                    psi[partition_index].hzy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzy, "fail to allocate memory for field");
                    psi[partition_index].exy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exy, "fail to allocate memory for field");
                    psi[partition_index].exz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exz, "fail to allocate memory for field");
                    psi[partition_index].eyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyx, "fail to allocate memory for field");
                    psi[partition_index].eyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyz, "fail to allocate memory for field");
                    psi[partition_index].ezx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezx, "fail to allocate memory for field");
                    psi[partition_index].ezy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezy, "fail to allocate memory for field");
                    break;
            }
        }
    }
    return 1;
}

void cpml_get_d()
{
    int partition_index, layer_index;
    int x, y, z;
    int x_start, x_stop, y_start, y_stop, z_start, z_stop;

    for (partition_index = 1; partition_index < 7; partition_index++)
    {
        if (partition_data[partition_index].boundary_type == boundary_pml)
        {
            x_start = partition_data[partition_index].x_start;
            x_stop  = partition_data[partition_index].x_stop;
            y_start = partition_data[partition_index].y_start;
            y_stop  = partition_data[partition_index].y_stop;
            z_start = partition_data[partition_index].z_start;
            z_stop  = partition_data[partition_index].z_stop;

            for (x = x_start; x < x_stop; x++)
            {
                for (y = y_start; y < y_stop; y++)
                {
                    for (z = z_start; z < z_stop; z++)
                    {
                        switch (partition_index)
                        {
                            case partition_z0: layer_index = z_stop - z - 1; break;
                            case partition_z1: layer_index =-z_start + z;    break;
                            case partition_x0: layer_index = x_stop - x - 1; break;
                            case partition_x1: layer_index =-x_start + x;    break;
                            case partition_y0: layer_index = y_stop - y - 1; break;
                            case partition_y1: layer_index =-y_start + y;    break;
                            default: layer_index = 0;
                        }
                        if (psi[partition_index].hxy && x != 0)
                        {
                            psi[partition_index].hxy[x-x_start][y-y_start][z-z_start] 
                                = be[layer_index] * psi[partition_index].hxy[x-x_start][y-y_start][z-z_start]
                                + cex[layer_index] * (hy[x][y][z] - hy[x-1][y][z]);
                            /*
                            if (x == 20 && y == 49)
                            {
                                printf("x,y --- layer: %03d, %03d, %03d \n", x, y, layer_index);
                                printf("Delta_hy: %e\n", hy[x][y][z] - hy[x-1][y][z]);
                                printf("psi.hxy: %e\n", psi[partition_index].hxy[x-x_start][y-y_start][z-z_start]);
                                printf("be: %e\n", be[layer_index]);
                                printf("cex: %e\n", cex[layer_index]);
                            }
                            */
                        }
                        if (psi[partition_index].hxz && x != 0)
                            psi[partition_index].hxz[x-x_start][y-y_start][z-z_start] 
                                = be[layer_index] * psi[partition_index].hxz[x-x_start][y-y_start][z-z_start]
                                + cex[layer_index] * (hz[x][y][z] - hz[x-1][y][z]);
                        if (psi[partition_index].hyx && y != 0)
                            psi[partition_index].hyx[x-x_start][y-y_start][z-z_start] 
                                = be[layer_index] * psi[partition_index].hyx[x-x_start][y-y_start][z-z_start]
                                + cey[layer_index] * (hx[x][y][z] - hx[x][y-1][z]);
                        if (psi[partition_index].hyz && y != 0)
                            psi[partition_index].hyz[x-x_start][y-y_start][z-z_start] 
                                = be[layer_index] * psi[partition_index].hyz[x-x_start][y-y_start][z-z_start]
                                + cey[layer_index] * (hz[x][y][z] - hz[x][y-1][z]);
                        if (psi[partition_index].hzx && z != 0)
                            psi[partition_index].hzx[x-x_start][y-y_start][z-z_start] 
                                = be[layer_index] * psi[partition_index].hzx[x-x_start][y-y_start][z-z_start]
                                + cez[layer_index] * (hx[x][y][z] - hx[x][y][z-1]);
                        if (psi[partition_index].hzy && z != 0)
                            psi[partition_index].hzy[x-x_start][y-y_start][z-z_start] 
                                = be[layer_index] * psi[partition_index].hzy[x-x_start][y-y_start][z-z_start]
                                + cez[layer_index] * (hy[x][y][z] - hy[x][y][z-1]);


                        switch (mode)
                        {
                            case mode_full:
                                if (y != 0)
                                    dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]) / pml_kappa_e[layer_index];
                                if (z != 0)
                                    dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]) / pml_kappa_e[layer_index];
                                if (z != 0)
                                    dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]) / pml_kappa_e[layer_index];
                                if (x != 0)
                                    dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]) / pml_kappa_e[layer_index];
                                if (x != 0)
                                    dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]) / pml_kappa_e[layer_index];
                                if (y != 0)
                                    dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]) / pml_kappa_e[layer_index];
                                break;
                            case mode_tmx:
                                if (y != 0)
                                    dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]) / pml_kappa_e[layer_index];
                                if (z != 0)
                                    dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]) / pml_kappa_e[layer_index];
                                break;
                            case mode_tmy:
                                if (z != 0)
                                    dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]) / pml_kappa_e[layer_index];
                                if (x != 0)
                                    dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]) / pml_kappa_e[layer_index];
                                break;
                            case mode_tmz:
                                if (x != 0)
                                    dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]) / pml_kappa_e[layer_index];
                                if (y != 0)
                                    dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]) / pml_kappa_e[layer_index];
                                break;
                            case mode_tex:
                                if (z != 0)
                                    dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]) / pml_kappa_e[layer_index];
                                if (y != 0)
                                    dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]) / pml_kappa_e[layer_index];
                                break;
                            case mode_tey:
                                if (z != 0)
                                    dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]) / pml_kappa_e[layer_index];
                                if (x != 0)
                                    dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]) / pml_kappa_e[layer_index];
                                break;
                            case mode_tez:
                                if (y != 0)
                                    dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]) / pml_kappa_e[layer_index];
                                if (x != 0)
                                    dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]) / pml_kappa_e[layer_index];
                                break;
                        }

                        switch (mode)
                        {
                            case mode_full:
                                dx[x][y][z] = dx[x][y][z] + d_t * psi[partition_index].hyz[x-x_start][y-y_start][z-z_start];
                                dx[x][y][z] = dx[x][y][z] - d_t * psi[partition_index].hzy[x-x_start][y-y_start][z-z_start];
                                dy[x][y][z] = dy[x][y][z] + d_t * psi[partition_index].hzx[x-x_start][y-y_start][z-z_start];
                                dy[x][y][z] = dy[x][y][z] - d_t * psi[partition_index].hxz[x-x_start][y-y_start][z-z_start];
                                dz[x][y][z] = dz[x][y][z] + d_t * psi[partition_index].hxy[x-x_start][y-y_start][z-z_start];
                                dz[x][y][z] = dz[x][y][z] - d_t * psi[partition_index].hyx[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tmx:
                                dx[x][y][z] = dx[x][y][z] + d_t * psi[partition_index].hyz[x-x_start][y-y_start][z-z_start];
                                dx[x][y][z] = dx[x][y][z] - d_t * psi[partition_index].hzy[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tmy:
                                dy[x][y][z] = dy[x][y][z] + d_t * psi[partition_index].hzx[x-x_start][y-y_start][z-z_start];
                                dy[x][y][z] = dy[x][y][z] - d_t * psi[partition_index].hxz[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tmz:
                                dz[x][y][z] = dz[x][y][z] + d_t * psi[partition_index].hxy[x-x_start][y-y_start][z-z_start];
                                dz[x][y][z] = dz[x][y][z] - d_t * psi[partition_index].hyx[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tex:
                                dy[x][y][z] = dy[x][y][z] + d_t * psi[partition_index].hzx[x-x_start][y-y_start][z-z_start];
                                dz[x][y][z] = dz[x][y][z] - d_t * psi[partition_index].hyx[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tey:
                                dx[x][y][z] = dx[x][y][z] - d_t * psi[partition_index].hzy[x-x_start][y-y_start][z-z_start];
                                dz[x][y][z] = dz[x][y][z] + d_t * psi[partition_index].hxy[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tez:
                                dx[x][y][z] = dx[x][y][z] + d_t * psi[partition_index].hyz[x-x_start][y-y_start][z-z_start];
                                dy[x][y][z] = dy[x][y][z] - d_t * psi[partition_index].hxz[x-x_start][y-y_start][z-z_start];
                                break;
                        }
                    }
                }
            }
        }
    }
}

void cpml_get_b ()
{
    int partition_index, layer_index;
    int x, y, z;
    int x_start, x_stop, y_start, y_stop, z_start, z_stop;

    for (partition_index = 1; partition_index < 7; partition_index++)
    {
        if (partition_data[partition_index].boundary_type == boundary_pml)
        {
            x_start = partition_data[partition_index].x_start;
            x_stop  = partition_data[partition_index].x_stop;
            y_start = partition_data[partition_index].y_start;
            y_stop  = partition_data[partition_index].y_stop;
            z_start = partition_data[partition_index].z_start;
            z_stop  = partition_data[partition_index].z_stop;

            for (x = x_start; x < x_stop; x++)
            {
                for (y = y_start; y < y_stop; y++)
                {
                    for (z = z_start; z < z_stop; z++)
                    {
                        switch (partition_index)
                        {
                            case partition_z0: layer_index = z_stop - z - 1; break;
                            case partition_z1: layer_index =-z_start + z;    break;
                            case partition_x0: layer_index = x_stop - x - 1; break;
                            case partition_x1: layer_index =-x_start + x;    break;
                            case partition_y0: layer_index = y_stop - y - 1; break;
                            case partition_y1: layer_index =-y_start + y;    break;
                            default: layer_index = 0;
                        }

                        if (psi[partition_index].exy)
                            psi[partition_index].exy[x-x_start][y-y_start][z-z_start] 
                                = bh[layer_index] * psi[partition_index].exy[x-x_start][y-y_start][z-z_start]
                                + chx[layer_index] * (ey[x+1][y][z] - ey[x][y][z]);
                        if (psi[partition_index].exz)
                            psi[partition_index].exz[x-x_start][y-y_start][z-z_start] 
                                = bh[layer_index] * psi[partition_index].exz[x-x_start][y-y_start][z-z_start]
                                + chx[layer_index] * (ez[x+1][y][z] - ez[x][y][z]);
                        if (psi[partition_index].eyx)
                            psi[partition_index].eyx[x-x_start][y-y_start][z-z_start] 
                                = bh[layer_index] * psi[partition_index].eyx[x-x_start][y-y_start][z-z_start]
                                + chy[layer_index] * (ex[x][y+1][z] - ex[x][y][z]);
                        if (psi[partition_index].eyz)
                            psi[partition_index].eyz[x-x_start][y-y_start][z-z_start] 
                                = bh[layer_index] * psi[partition_index].eyz[x-x_start][y-y_start][z-z_start]
                                + chy[layer_index] * (ez[x][y+1][z] - ez[x][y][z]);
                        if (psi[partition_index].ezx)
                            psi[partition_index].ezx[x-x_start][y-y_start][z-z_start] 
                                = bh[layer_index] * psi[partition_index].ezx[x-x_start][y-y_start][z-z_start]
                                + chz[layer_index] * (ex[x][y][z+1] - ex[x][y][z]);
                        if (psi[partition_index].ezy)
                            psi[partition_index].ezy[x-x_start][y-y_start][z-z_start] 
                                = bh[layer_index] * psi[partition_index].ezy[x-x_start][y-y_start][z-z_start]
                                + chz[layer_index] * (ey[x][y][z+1] - ey[x][y][z]);

                        switch (mode)
                        {
                            case mode_full:
                                bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]) / pml_kappa_h[layer_index];
                                bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]) / pml_kappa_h[layer_index];
                                by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]) / pml_kappa_h[layer_index];
                                by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]) / pml_kappa_h[layer_index];
                                bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]) / pml_kappa_h[layer_index];
                                bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]) / pml_kappa_h[layer_index];
                                break;
                            case mode_tmx:
                                by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]) / pml_kappa_h[layer_index];
                                bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]) / pml_kappa_h[layer_index];
                                break;
                            case mode_tmy:
                                bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]) / pml_kappa_h[layer_index];
                                bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]) / pml_kappa_h[layer_index];
                                break;
                            case mode_tmz:
                                bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]) / pml_kappa_h[layer_index];
                                by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]) / pml_kappa_h[layer_index];
                                break;
                            case mode_tex:
                                bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]) / pml_kappa_h[layer_index];
                                bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]) / pml_kappa_h[layer_index];
                                break;
                            case mode_tey:
                                by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]) / pml_kappa_h[layer_index];
                                by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]) / pml_kappa_h[layer_index];
                                break;
                            case mode_tez:
                                bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]) / pml_kappa_h[layer_index];
                                bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]) / pml_kappa_h[layer_index];
                                break;
                        }

                        switch (mode)
                        {
                            case mode_full:
                                bx[x][y][z] = bx[x][y][z] - d_t * psi[partition_index].eyz[x-x_start][y-y_start][z-z_start];
                                bx[x][y][z] = bx[x][y][z] + d_t * psi[partition_index].ezy[x-x_start][y-y_start][z-z_start];
                                by[x][y][z] = by[x][y][z] - d_t * psi[partition_index].ezx[x-x_start][y-y_start][z-z_start];
                                by[x][y][z] = by[x][y][z] + d_t * psi[partition_index].exz[x-x_start][y-y_start][z-z_start];
                                bz[x][y][z] = bz[x][y][z] - d_t * psi[partition_index].exy[x-x_start][y-y_start][z-z_start];
                                bz[x][y][z] = bz[x][y][z] + d_t * psi[partition_index].eyx[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tmx:
                                by[x][y][z] = by[x][y][z] - d_t * psi[partition_index].ezx[x-x_start][y-y_start][z-z_start];
                                bz[x][y][z] = bz[x][y][z] + d_t * psi[partition_index].eyx[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tmy:
                                bx[x][y][z] = bx[x][y][z] + d_t * psi[partition_index].ezy[x-x_start][y-y_start][z-z_start];
                                bz[x][y][z] = bz[x][y][z] - d_t * psi[partition_index].exy[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tmz:
                                bx[x][y][z] = bx[x][y][z] - d_t * psi[partition_index].eyz[x-x_start][y-y_start][z-z_start];
                                by[x][y][z] = by[x][y][z] + d_t * psi[partition_index].exz[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tex:
                                bx[x][y][z] = bx[x][y][z] - d_t * psi[partition_index].eyz[x-x_start][y-y_start][z-z_start];
                                bx[x][y][z] = bx[x][y][z] + d_t * psi[partition_index].ezy[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tey:
                                by[x][y][z] = by[x][y][z] - d_t * psi[partition_index].ezx[x-x_start][y-y_start][z-z_start];
                                by[x][y][z] = by[x][y][z] + d_t * psi[partition_index].exz[x-x_start][y-y_start][z-z_start];
                                break;
                            case mode_tez:
                                bz[x][y][z] = bz[x][y][z] - d_t * psi[partition_index].exy[x-x_start][y-y_start][z-z_start];
                                bz[x][y][z] = bz[x][y][z] + d_t * psi[partition_index].eyx[x-x_start][y-y_start][z-z_start];
                                break;
                        }
                    }
                }
            }
        }
    }
}
